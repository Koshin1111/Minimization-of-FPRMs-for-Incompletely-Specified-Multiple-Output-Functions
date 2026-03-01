using Random


include("for_println.jl")   # 出力用 
include("read_benchmark.jl")

include("Alg_non_share.jl") # 共有なし
include("Alg_share.jl")     # 共有あり
#=
total = 10
k=5
idx = randperm(total)[1:k]

println("idx:$idx")
=#
#=
w_share = [
    [[2,4,5],[8,4],[3]],
    [[]],
    [[8,9],[4,5,5],[Inf]]
]

w_sum_n = sum( minimum.( length.(w_share) ) ) 
w_max_n = maximum( minimum.( length.(w_share) ) )

println("w_sum:$w_sum_n, w_max:$w_max_n")
=#
# ---------------------------------------------------------

# 割り当て(assign_d_sum)
function min_sum_terms(w_non_share, d_reduced, d_all) 
    w_sum_n = 0 # 重み(積項数)合計値
    assign_d_sum = [Vector{Int}() for _ in 1:length(d_all)]
    for i in 1:length(w_non_share)
        min_w = Int(minimum(w_non_share[i])) # 最小重み用
        min_num = findfirst(==(min_w), w_non_share[i])-1 # 最小重み要素番号用
        # min_num を2進数の割り当て用ベクトルに変換
        assign_d_i = Vector{Int}(undef, length(d_reduced[i]))
        for j in 0:length(d_reduced[i])-1
            assign_d_i[length(d_reduced[i])-j] = div(min_num, 2^j) % 2 # 2^j桁目の整数を抽出
        end
        if d_all[i] != d_reduced[i]
            # 削減された依存変数(0,1どちらを入れても積項数が変わらない依存変数)を assign_d_sum に追加
            # 積項数に影響しない依存変数には 0 を割り当てる
            assign_d_sum[i] = zeros(Int, length(d_all[i]))
            k = 1
            for j in 1:length(d_all[i])
                if d_all[i][j] in d_reduced[i]
                    assign_d_sum[i][j] = assign_d_i[k]
                    if (k < length(d_reduced[i])) k += 1 else break end
                end
            end
        else
            assign_d_sum[i] = copy(assign_d_i)
        end
        w_sum_n += min_w
    end
    return w_sum_n, assign_d_sum
end

# 割り当て(assign_d_sum_share)
# 再帰関数で全探索
function min_share_terms(w_share, d_all) 
    m = length(w_share) 
    best_weight = Ref(Inf) 
    best_choice = Ref(Vector{Int}()) 
    
    function dfs(f, used::Set{Int}, choice::Vector{Int}) # 全出力を選び終えた 
        if f > m 
            if length(used) < best_weight[] 
                best_weight[] = length(used) 
                best_choice[] = copy(choice) 
            end 
            return 
        end # f番目の候補を順番に試す 
        
        for k in 1:length(w_share[f])
             # [Inf] はスキップ 
            if w_share[f][k] == [Inf] continue end
            check = Int.(w_share[f][k]) 
            # 新しい積項集合 
            new_used = union(used, Set(check)) 
            # 枝刈り 
            if length(new_used) >= best_weight[] 
                continue 
            end # 選択履歴を追加して再帰 
            push!(choice, k) 
            dfs(f + 1, new_used, choice) 
            pop!(choice) 
        end 
    end 
    dfs(1, Set{Int}(), Int[])

    # best_coice の各要素を2進数の割り当て用ベクトルに変換
    #=
    assign_d = [Vector{Int}(undef, length(d_all[i])) for i in 1:m]
    for i in 1:m
        for j in 0:(length(d_all[i])-1)
            assign_d[i][length(d_all[i])-j] = div(best_choice[][i]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
    end
    =#
     
    return Int(best_weight[]), best_choice[]
end

function run_once(n, m, m_vec, threshold)

    # --- 共有なし ---
    w_non_share = Vector{Any}(undef, m)
    d_reduced   = Vector{Vector{Int}}(undef, m)
    T           = [ [Int[] for _ in 1:2^n] for _ in 1:m ]
    T_d         = [ [Int[] for _ in 1:2^n] for _ in 1:m ]

    # --- 共有あり ---
    w_share = Vector{Any}(undef, m)

    for i in 1:m
        # --- 共有なし ---
        w_non_share[i], d_reduced[i], T[i], T_d[i] = Alg_non_share(n, m_vec[i], threshold)
        # --- 共有あり ---
        w_share[i] = Alg_share(n, m_vec[i], threshold)
    end
    #for i in 1:m if w_non_share[i] == [Inf] error("No PPRM found within the threshold(", threshold, ") at f", i, ".") end end
    for i in 1:m if w_non_share[i] == [Inf] return -1,-1,-1,-1 end end # 閾値エラー


    # 全ての依存変数を記憶
    d_all = [Int[] for _ in 1:m]
    for i in 1:m
        for j in 1:2^n
            if m_vec[i][j] == 2
                push!(d_all[i], j-1)
            end
        end
    end

    # --- 共有なし ---
    # 合計
    #w_sum_n = sum(minimum.(w_non_share))
    w_sum_n, assign_d_sum = min_sum_terms(w_non_share, d_reduced, d_all) 
    # 最大
    w_max_n = maximum(minimum.(w_non_share))

    # --- 共有あり ---
    # 合計（best_choice も取得）
    w_sum_s, best_choice = min_share_terms(w_share, d_all)
    assign_d_share = [Vector{Int}(undef, length(d_all[i])) for i in 1:m]
    for i in 1:m
        for j in 0:(length(d_all[i])-1)
            assign_d_share[i][length(d_all[i])-j] = div(best_choice[i]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
    end

    # 最大（best_choice を使う）
    w_max_s = maximum([ length(w_share[i][best_choice[i]]) for i in 1:m ])


    println("--- PPRM of sum_n:", w_sum_n," ---")
    print_PPRM(n, m, T, T_d, assign_d_sum)
    println("")

    println("--- PPRM of sum_s:", w_sum_s," ---")
    print_PPRM(n, m, T, T_d, assign_d_share)

    return w_sum_n, w_sum_s, w_max_n, w_max_s
end

thtreshold = 5  # 閾値設定
n = 3 # 入力変数の個数
m = 3 # 出力ビット数
m_vec = [
    [1,0,2,1,1,2,1,0],   # 出力ビット0の真理値ベクトル
    [0,2,1,1,2,1,0,1],   # 出力ビット1の真理値ベクトル
    [2,2,1,0,2,0,1,0]    # 出力ビット2の真理値ベクトル
]
w_sum_n, w_sum_s, w_max_n, w_max_s = run_once(n, m, m_vec, thtreshold)
println("")
println("w_sum_n:$w_sum_n, w_sum_s:$w_sum_s, w_max_n:$w_max_n, w_max_s:$w_max_s")