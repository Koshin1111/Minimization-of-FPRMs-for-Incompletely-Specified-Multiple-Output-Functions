# ------------------------------------------------------------------------
# 積項を共有しない場合と共有する場合の両方で、最小 FPRM の総積項数と最大積項数が 
# don`t care の個数に伴い、どのように推移するのか確認するための実験コード
# 
# 共有する場合の最小 FPRM 探索の計算コストが高すぎたため、後に
# 出力全体で極性を固定する手法を追加。
# それぞれの手法における最小 FPRM の探索時間についても計測する。
#
# ns : 積項を共有しない場合
# s1 : 積項を共有する場合 (出力全体で極性を固定しない)
# s2 : 積項を共有する場合 (出力全体で極性を固定する)
# -------------------------------------------------------------------------

using Random

include("for_println.jl")   # 出力用 
include("read_benchmark.jl")

include("Alg_shared.jl") 

threshold = 5  # 閾値設定


# --- 共有しない場合の割り当て読み込み ---
# 割り当て(assign_d_sum)
function find_min_ns(n,w,d_all,threshold) 
    w_ns = 0 # 重み(積項数)合計値
    w_max_n = 0 # 重み(積項数)最大値
    m = length(w) # 出力ビット数
    wei_vec = [ [Vector{Int}(undef, length(w[i][pola])) for pola in 1:2^n] for i in 1:m ] # 重みベクトル(積項数のみ)
    for i in 1:m
        for pola in 1:2^n
            for asg in 1:length(w[i][pola])

                if w[i][pola][asg] == [Inf] || w[i][pola][asg] == Inf
                    wei_vec[i][pola][asg] = threshold + 1 # [Inf] は閾値+1(選択不可)の重みとする
                else
                    wei_vec[i][pola][asg] = length(w[i][pola][asg]) # 端点ベクトルを積項数(Int)に変換
                end
            end
        end
    end
    #polarity_vector = [Vector{Int}(undef, n) for i in 1:m] # 極性ベクトル
    #assign_d = [Vector{Int}(undef, length(d_all[i])) for i in 1:m] # 割り当て

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
                #=
                for j in 0:(n-1)   
                    polarity_vector[i][n-j] = div(pola-1, 2^j) % 2 # 極性ベクトル更新
                end
                =#
            end
        end
        # min_num を2進数の割り当て用ベクトルに変換
        #=
        assign_d_i = Vector{Int}(undef, length(d_all[i]))
        for j in 0:length(d_all[i])-1
            assign_d_i[length(d_all[i])-j] = div(min_num, 2^j) % 2 # 2^j桁目の整数を抽出
        end
        assign_d[i] = assign_d_i
        =#
        w_ns += min_w
        w_max_n = max(min_w, w_max_n)
    end
    return w_ns, w_max_n
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
                if w[f][polarity[f]][k] == [Inf] || w[f][polarity[f]][k] == Inf continue end
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
                        # xor 演算により、極性が異なる部分を求める
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
            new_used_each = copy(used_each)
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
    #=
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
    =#
     
    return Int(best_weight[]), Int(max_weight[])
end

# --- 共有する場合の割り当て読み込み ---

# 再帰関数で全探索
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
                if w[f][polarity][k] == [Inf] || w[f][polarity][k] == Inf continue end
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
    #=
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
    =#
            
    return Int(best_weight[]), Int(max_weight[])
end



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


# ---------------------------------------------------------
# 1 回分の実験（共有なし / 共有あり）
# w_sum_n, w_sum_s, w_max_n, w_max_s を返す
# ---------------------------------------------------------
function run_once(n, m, m_vec, threshold)
    
    # --- 重みベクトル ---
    w = Vector{Any}(undef, m)

    for i in 1:m
        # --- Alg_shared --- 
        w[i] = Alg_shared(n, m_vec[i], threshold)
    end
    #for i in 1:m if w_non_share[i] == [Inf] error("No PPRM found within the threshold(", threshold, ") at f", i, ".") end end
    for i in 1:m for pola in 1:2^n if w[i][pola] == [Inf] return -1,-1,-1,-1,-1,-1,-1,-1,-1 end end end # 閾値エラー


    # 全ての依存変数を記憶
    d_all = [Int[] for _ in 1:m]
    for i in 1:m
        for j in 1:2^n
            if m_vec[i][j] == 2
                push!(d_all[i], j-1)
            end
        end
    end

    # --- 時間計測付き ---

    # 積項を共有しない場合
    t_n = @elapsed begin
        w_sum_n, w_max_n = find_min_ns(n, w, d_all, threshold)
    end
    # 積項を共有する場合
    # 出力全体で極性を固定しない場合
    t_s1 = @elapsed begin
        w_sum_s1, w_max_s1 = find_min_s_1(n, w, d_all)
    end
    # 出力全体で極性を固定する場合
    t_s2 = @elapsed begin
        w_sum_s2, w_max_s2 = find_min_s_2(n, w, d_all)
    end

    return w_sum_n, w_sum_s1, w_sum_s2, w_max_n, w_max_s1, w_max_s2, t_n, t_s1, t_s2
end


# ---------------------------------------------------------
# k 個の don't-care を入れて trials 回実験し平均を返す
# ---------------------------------------------------------
function experiment(n, m, m_vec, k, trials)
    sum_n = 0.0; sum_s1 = 0.0; sum_s2 = 0.0 
    max_n = 0.0; max_s1 = 0.0; max_s2 = 0.0

    time_n = 0.0; time_s1 = 0.0; time_s2 = 0.0

    adjustable_threshold = copy(threshold)

    trials = (k == 0 ? 1 : trials) # 完全定義のPPRMは一意に定まる

    count = 0 
    while count < trials
        m_vec_dc = add_random_dontcare(m_vec, k)
        w_n, w_s1, w_s2, wmn, wms1, wms2, t_n, t_s1, t_s2 = run_once(n, m, m_vec_dc, adjustable_threshold)
        if w_n == -1 # 閾値エラー
            adjustable_threshold += 1
            continue
        end

        sum_n += w_n; sum_s1 += w_s1; sum_s2 += w_s2
        max_n += wmn; max_s1 += wms1; max_s2 += wms2

        time_n += t_n; time_s1 += t_s1; time_s2 += t_s2

        count += 1
    end

    return sum_n / trials, sum_s1 / trials, sum_s2 / trials, 
           max_n / trials, max_s1 / trials, max_s2 / trials,
           time_n, time_s1, time_s2
end


# ---------------------------------------------------------
# メイン実行部
# ---------------------------------------------------------
n, m, m_vec = pla_to_m("Sum.pla")

ks = [0, 5, 10, 15]
results = Dict()

for k in ks
    avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2 = experiment(n, m, m_vec, k, 100)
    results[k] = (avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2)

    println("don't care = $k : avg(w_sum_n) = $avg_n,  avg(w_sum_s1) = $avg_s1,  avg(w_sum_s2) = $avg_s2")
    println("                  avg(w_max_n) = $avg_max_n, avg(w_max_s1) = $avg_max_s1, avg(w_max_s2) = $avg_max_s2")
    println("                  time_n = $time_n, time_s1 = $time_s1, time_s2 = $time_s2\n")
end

println("\n=== Sum 最終結果 ===")
for k in ks
    println("DC=$k :")
    println("  non_shared     : sum=$(results[k][1]), max=$(results[k][4]), time=$(results[k][7])")
    println("    shared_1     : sum=$(results[k][2]), max=$(results[k][5]), time=$(results[k][8])")
    println("    shared_2     : sum=$(results[k][3]), max=$(results[k][6]), time=$(results[k][9])\n")
end




n, m, m_vec = pla_to_m("Cmp.pla")

ks = [0, 5, 10, 15]
results = Dict()

for k in ks
    avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2 = experiment(n, m, m_vec, k, 100)
    results[k] = (avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2)

    println("don't care = $k : avg(w_sum_n) = $avg_n,  avg(w_sum_s1) = $avg_s1,  avg(w_sum_s2) = $avg_s2")
    println("                  avg(w_max_n) = $avg_max_n, avg(w_max_s1) = $avg_max_s1, avg(w_max_s2) = $avg_max_s2")
    println("                  time_n = $time_n, time_s1 = $time_s1, time_s2 = $time_s2\n")
end

println("\n=== Cmp 最終結果 ===")
for k in ks
    println("DC=$k :")
    println("  non_shared     : sum=$(results[k][1]), max=$(results[k][4]), time=$(results[k][7])")
    println("    shared_1     : sum=$(results[k][2]), max=$(results[k][5]), time=$(results[k][8])")
    println("    shared_2     : sum=$(results[k][3]), max=$(results[k][6]), time=$(results[k][9])\n")
end




n, m, m_vec = pla_to_m("Mul.pla")

ks = [0, 5, 10, 15]
results = Dict()

for k in ks
    avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2 = experiment(n, m, m_vec, k, 100)
    results[k] = (avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2)

    println("don't care = $k : avg(w_sum_n) = $avg_n,  avg(w_sum_s1) = $avg_s1,  avg(w_sum_s2) = $avg_s2")
    println("                  avg(w_max_n) = $avg_max_n, avg(w_max_s1) = $avg_max_s1, avg(w_max_s2) = $avg_max_s2")
    println("                  time_n = $time_n, time_s1 = $time_s1, time_s2 = $time_s2\n")
end

println("\n=== Mul 最終結果 ===")
for k in ks
    println("DC=$k :")
    println("  non_shared     : sum=$(results[k][1]), max=$(results[k][4]), time=$(results[k][7])")
    println("    shared_1     : sum=$(results[k][2]), max=$(results[k][5]), time=$(results[k][8])")
    println("    shared_2     : sum=$(results[k][3]), max=$(results[k][6]), time=$(results[k][9])\n")
end



n, m, m_vec = pla_to_m("ALU.pla")

ks = [0, 5, 10, 15]
results = Dict()

for k in ks
    avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2 = experiment(n, m, m_vec, k, 100)
    results[k] = (avg_n, avg_s1, avg_s2, avg_max_n, avg_max_s1, avg_max_s2, time_n, time_s1, time_s2)

    println("don't care = $k : avg(w_sum_n) = $avg_n,  avg(w_sum_s1) = $avg_s1,  avg(w_sum_s2) = $avg_s2")
    println("                  avg(w_max_n) = $avg_max_n, avg(w_max_s1) = $avg_max_s1, avg(w_max_s2) = $avg_max_s2")
    println("                  time_n = $time_n, time_s1 = $time_s1, time_s2 = $time_s2\n")
end

println("\n=== ALU 最終結果 ===")
for k in ks
    println("DC=$k :")
    println("  non_shared     : sum=$(results[k][1]), max=$(results[k][4]), time=$(results[k][7])")
    println("    shared_1     : sum=$(results[k][2]), max=$(results[k][5]), time=$(results[k][8])")
    println("    shared_2     : sum=$(results[k][3]), max=$(results[k][6]), time=$(results[k][9])\n")
end
