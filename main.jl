# ------------------------------------------------------------------------
# 不完全定義多出力関数における FPRM の最小化を行うためのコード
# 積項を共有しない場合(non_shared)と共有する場合(shared)の両方の最小化結果を出力する。
# -------------------------------------------------------------------------

using Random

include("for_println.jl")   # 出力用 
include("read_benchmark.jl") # plaファイルから、真理値ベクトルを読み込むファイル

include("Alg_shared.jl") # 重みベクトル作成用ファイル

threshold = 8  # 閾値設定

#=
n = 3 # 入力変数の個数
m = 3 # 出力ビット数

m_vec = [
    [1,2,2,1,0,1,1,1],  # 出力ビット0の真理値ベクトル
    [2,1,1,0,0,1,2,2],  # 出力ビット1の真理値ベクトル
    [1,2,1,0,1,2,1,0]   # 出力ビット2の真理値ベクトル
]
=#

n, m, m_vec = pla_to_m("Sum.pla") # pla多出力関数読み込み

# ---------------------------------------------------------
# m_vec の中からランダムに k 個を 2 (don't-care) に書き換える
# ---------------------------------------------------------
function add_random_dontcare(m_vec::Vector{Vector{Int}}, k::Int)
    n_outputs = length(m_vec)
    n_inputs  = length(m_vec[1])
    total = n_outputs * n_inputs

    idx = randperm(total)[1:k]   # ランダムに k 個選ぶ
    new_m_vec = deepcopy(m_vec)

    for id in idx
        f = div(id-1, n_inputs) + 1
        t = mod(id-1, n_inputs) + 1
        new_m_vec[f][t] = 2
    end

    return new_m_vec
end

m_vec = add_random_dontcare(m_vec, 15) # m_vec の中からランダムに k 個を 2 (don't-care) に書き換える
# ---------------------------------------------------------

polarity_vector_ns = [Int[] for _ in 1:m] # 極性ベクトル(共有なし)
polarity_vector_s = [Int[] for _ in 1:m] # 極性ベクトル(共有あり)

assign_d_ns = [Int[] for _ in 1:m] # 割り当て(共有なし)
assign_d_s = [Int[] for _ in 1:m] # 割り当て(共有あり)

w = Vector{Any}(undef, m) # 重みベクトル

for i in 1:m
    # --- Alg_shared --- 
    w[i] = Alg_shared(n, m_vec[i], threshold)
end

# 閾値エラー確認
for i in 1:m for pola in 1:2^n if w[i][pola] == [Inf] error("Error: No FPRM found within the threshold(", threshold, ") for any polarity.") end end end # 閾値エラー

# 全ての依存変数を記憶
d_all = [Int[] for _ in 1:m]
for i in 1:m
    for j in 1:2^n
        if m_vec[i][j] == 2
            push!(d_all[i], j-1)
        end
    end
end

println("---------- truth_vector before minimization ----------")
for i  in 1:length(m_vec)
    println("f_$(i-1) : ", show_clean(m_vec[i]))
end
println("")

#=
println("---------- weight_vector before minimization ----------")
for i  in 1:length(w)
    println("--- f_$(i-1) ---")
    for j in 1:length(w[i])
        println("polarity_$(j-1): ", show_clean(w[i][j]))
    end
end
println("")
=#


# --- 共有しない場合の割り当て読み込み ---
# 割り当て(assign_d_sum)
function find_min_ns(n,w,d_all,threshold) 
    w_ns = 0 # 重み(積項数)合計値
    m = length(w) # 出力ビット数
    wei_vec = [ [Vector{Int}(undef, length(w[i][pola])) for pola in 1:2^n] for i in 1:m ] # 重みベクトル(積項数のみ)
    for i in 1:m
        for pola in 1:2^n
            for asg in 1:length(w[i][pola])
                if w[i][pola][asg] == [Inf]
                    wei_vec[i][pola][asg] = threshold + 1 # [Inf] は閾値+1(選択不可)の重みとする
                else
                    wei_vec[i][pola][asg] = length(w[i][pola][asg]) # 端点ベクトルを積項数(Int)に変換
                end
            end
        end
    end
    polarity_vector = [Vector{Int}(undef, n) for i in 1:m] # 極性ベクトル
    assign_d = [Vector{Int}(undef, length(d_all[i])) for i in 1:m] # 割り当て

    for i in 1:m
        min_w = threshold # 最小重み用
        min_num = 0
        for pola in 1:2^n
            if w[i][pola] == [Inf] continue end # [Inf] はスキップ
            min_w_pola = Int.(minimum(wei_vec[i][pola])) # 最小重み用
            if min_w_pola < min_w
                min_w = min_w_pola
                min_num = findfirst(==(min_w), wei_vec[i][pola])-1 # 最小重み要素番号用
                # pola を2進数の極性ベクトルに変換
                for j in 0:(n-1)   
                    polarity_vector[i][n-j] = div(pola-1, 2^j) % 2 # 極性ベクトル更新
                end
            end
        end
        # min_num を2進数の割り当て用ベクトルに変換
        assign_d_i = Vector{Int}(undef, length(d_all[i]))
        for j in 0:length(d_all[i])-1
            assign_d_i[length(d_all[i])-j] = div(min_num, 2^j) % 2 # 2^j桁目の整数を抽出
        end
        assign_d[i] = assign_d_i
        w_ns += min_w
    end
    return w_ns, polarity_vector, assign_d
end
 


# --- 共有する場合の割り当て読み込み ---

# 再帰関数で全探索
# 出力全体で極性を固定しない場合
function find_min_s_1(n, w, d_all) 
    m = length(w) # 出力ビット数
    best_weight = Ref(Inf) 
    max_weight = Ref(Inf) # max of best_weight
    best_choice = Ref(Vector{Int}()) 
    best_polarity = Ref(Vector{Int}())
    
    function dfs(f, used::Set{Int}, used_each::Vector{Vector{Int}}, polarity::Vector{Int}, choice::Vector{Int}) # 全出力を選び終えた 
        if f > m 
            if ( length(used) < best_weight[] ) || ( length(used) == best_weight[] && maximum(length.(used_each)) < max_weight[] )
                best_weight[] = length(used) 
                max_weight[] = maximum(length.(used_each))
                best_choice[] = copy(choice) 
                best_polarity[] = copy(polarity)
            end 
            return 
        end # f番目の候補を順番に試す 
        
        for k in 1:length(w[f][polarity[f]])
            # [Inf] はスキップ 
                if w[f][polarity[f]][k] == [Inf] continue end
            check = Int.(w[f][polarity[f]][k]) 
            # 新しい積項集合 
            new_used = union(used, Set(check)) 
            # 枝刈り 
            if length(new_used) >= best_weight[] 
                continue 
            end # 選択履歴を追加して再帰 

            # 極性の違いを確認 (極性が違う場合でも、影響がなければ OK とする)
            if f > 1
                shareable = true # 共有可能かどうかのフラグ

                for func in 1:f-1
                    diff_pola = 0 # 極性の違いを表すビット列
                    if polarity[func] != polarity[f] # 極性が違う場合
                        # xor 演算により、極性が異なる変数を求める
                        diff_pola = (polarity[func]-1) ⊻ (polarity[f]-1)
                    end
                    # 共有する積項を求める
                    shared_terms = intersect(used_each[func], check)
                    for term in shared_terms
                        # 共有したい積項の極性が違う場合 ( 極性が異なる変数を係数に含む場合 )
                        if (term-1) & diff_pola != 0 
                            shareable = false # 共有不可
                            break
                        end
                    end
                    if !shareable break end
                end
                if !shareable
                    continue # 共有不可な場合はスキップ
                end
            end
            new_used_each = [copy(v) for v in used_each]
            push!(new_used_each, check)
            
            push!(choice, k) 
            for pol in 1:2^n
                push!(polarity, pol)
                dfs(f + 1, new_used, new_used_each, polarity, choice) 
                pop!(polarity)
            end
            pop!(choice) 
        end 
    end 

    polarity = Int[]
    for pol in 1:2^n
        push!(polarity, pol)
        dfs(1, Set{Int}(), Vector{Vector{Int}}(), polarity, Int[])
        pop!(polarity)
    end

    # ---出力用配列作成---
    assign_d = [Vector{Int}(undef, length(d_all[i])) for i in 1:m]
    polarity_vector = [Vector{Int}(undef, n) for i in 1:m]
    for i in 1:m
        # best_coice の各要素を2進数の割り当て用ベクトルに変換
        for j in 0:(length(d_all[i])-1)
            assign_d[i][length(d_all[i])-j] = div(best_choice[][i]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
        # best_polarity の各要素を2進数の極性ベクトルに変換
        for j in 0:n-1
            polarity_vector[i][n-j] = div(best_polarity[][i]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
    end
     
    return Int(best_weight[]), polarity_vector, assign_d
end



# 出力全体で極性を固定する場合
function find_min_s_2(n, w, d_all) 
    m = length(w) # 出力ビット数
    best_weight = Ref(Inf) 
    max_weight = Ref(Inf) # max of best_weight
    best_choice = Ref(Vector{Int}()) 
    best_polarity = Ref(0)
    
    function dfs(f, used::Set{Int}, polarity::Int, choice::Vector{Int}) # 全出力を選び終えた 
        if f > m 
            if (length(used) < best_weight[]) || ( length(used) == best_weight[] && maximum([ length(w[i][polarity][choice[i]]) for i in 1:m ]) < max_weight[] )
                best_weight[] = length(used) 
                max_weight[] = maximum([ length(w[i][polarity][choice[i]]) for i in 1:m ])
                best_choice[] = copy(choice) 
                best_polarity[] = copy(polarity)
            end 
            return 
        end # f番目の候補を順番に試す 
        
        for k in 1:length(w[f][polarity])
            # [Inf] はスキップ 
                if w[f][polarity][k] == [Inf] continue end
            check = Int.(w[f][polarity][k]) 
            # 新しい積項集合 
            new_used = union(used, Set(check)) 
            # 枝刈り 
            if length(new_used) >= best_weight[] 
                continue 
            end # 選択履歴を追加して再帰 
            
            push!(choice, k) 
            dfs(f + 1, new_used, polarity, choice)
            pop!(choice) 
        end 
    end 

    polarity = Int[]
    for pol in 1:2^n
        dfs(1, Set{Int}(), pol, Int[])
    end

    # ---出力用配列作成---
    assign_d = [Vector{Int}(undef, length(d_all[i])) for i in 1:m]
    polarity_vector = [Vector{Int}(undef, n) for i in 1:m]
    for i in 1:m
        # best_coice の各要素を2進数の割り当て用ベクトルに変換
        for j in 0:(length(d_all[i])-1)
            assign_d[i][length(d_all[i])-j] = div(best_choice[][i]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
        # best_polarity の各要素を2進数の極性ベクトルに変換
        for j in 0:n-1
            polarity_vector[i][n-j] = div(best_polarity[]-1, 2^j) % 2 # 2^j桁目の整数を抽出
        end
    end
     
    return Int(best_weight[]), polarity_vector, assign_d
end

# --- 割り当て探索関数呼び出し ---
# find_min_ns : 積項を共有しない場合 (積項数のみで評価)
w_ns, polarity_vector_ns, assign_d_ns = find_min_ns(n, w, d_all, threshold)
# find_min_s_1 : 出力全体で極性を固定しない場合
# find_min_s_2 : 出力全体で極性を固定する場合
w_s, polarity_vector_s, assign_d_s = find_min_s_1(n, w, d_all)



println("---------- non_shared results ----------")
println("  weight : ", show_clean(w_ns))
println("polarity : ", show_clean(polarity_vector_ns))
println("assign_d : ", show_clean(assign_d_ns))
println("")
println("--- FPRM ---")
print_FPRM(n, m, w, assign_d_ns, polarity_vector_ns)
println("")

println("---------- shared results ----------")
println("  weight : ", show_clean(w_s))
println("polarity : ", show_clean(polarity_vector_s))
println("assign_d : ", show_clean(assign_d_s))
println("")
println("--- FPRM ---")
print_FPRM(n, m, w, assign_d_s, polarity_vector_s)
println("")