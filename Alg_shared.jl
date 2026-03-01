# ------------------------------------------------------------------------
# 不完全定義関数の入力から、積項の共有にも対応した重みベクトル w を出力するコード
#
# w : 不完全定義多出力関数における FPRM の最小化を行うための探索対象ベクトル
# w の構築には、拡張真理値ベクトル t と依存変数格納用 t_d を用いる。
# -------------------------------------------------------------------------

include("for_println.jl")   # 出力用

# ---------- Step 1: 拡張真理値ベクトルを作成 ----------
# m       :: Vector  → 通常の真理値表 (2^n 要素)
# n       :: Int          → 変数の個数
function build_extended_truth_vector(n::Int, m::Vector)
    N = 3^n
    t = [Int[] for _ in 1:N]  # Any型で Int またはネストされた配列を格納
    t_d = [Int[] for _ in 1:N]  # 依存変数格納用

    for i in 0:(N-1)
        digit = zeros(Int, n)

        for j in 0:(n-1)
            digit[j + 1] = div(i, 3^j) % 3
        end

        function extend_t(t::Vector, size::Int, target_index::Int)
            ex_t = zeros(Int, 2*length(t))
            loop_size = 2^(size - target_index)
            l = 1
            for m in 0:Int((length(ex_t) / (2*loop_size)))-1
                # 未定義：配列を拡張（左に 0 割当て、右に 1 割当て）
                # 左側
                for o in 1:2^(size - target_index)
                    ex_t[l] = t[o + m*loop_size]
                    l += 1
                end
                # 右側
                for o in 1:2^(size - target_index)
                    ex_t[l] = t[o + m*loop_size] ⊻ 1 # 分岐追加
                    l += 1
                end
            end

            return ex_t
        end
        
        function Construct_t(j::Int, digit::Vector{Int})
            if j == n && digit[n] !== 2
                # digit → m のインデックスに変換
                index = 0
                for k in 0:n-1
                    index += digit[k + 1] * 2^(k)
                end

                
                if m[index + 1] == 2
                    # 未定義：配列を拡張（左に 0 割当て、右に 1 割当て）
                    push!(t_d[i + 1], index)
                    sort!(t_d[i + 1])

                    for k in 1:length(t_d[i + 1]) # 追加した依存変数が何番目に分岐するか確認
                        if t_d[i + 1][k] == index
                            current_t = extend_t(current_t, length(t_d[i+1]), k) # 依存変数分岐を追加
                            break
                        end
                    end
                else
                    # 定義済み：すべての葉に XOR を適用
                    current_t = current_t .⊻ m[index + 1]
                end

            elseif digit[j] == 2
                digit[j] = 0
                Construct_t(j, digit)
                digit[j] = 1
                Construct_t(j, digit)
                digit[j] = 2
            else
                Construct_t(j + 1, digit)
            end
        end
       
        current_t = [0] # 初期値
        Construct_t(1, digit)   
        t[i + 1] = current_t
    end

    return t, t_d
end

# ---------- Step 2: 重みベクトルを作成 ----------
# t         :: Vector  → 拡張真理値ベクトル (3^n 要素)
# t_d       :: Vector  → 依存変数格納用  (3^n 要素)
# n         :: Int          → 変数の個数
# threshold :: Int          → 積項数の閾値
function build_weight_vector(t::Vector, t_d::Vector, n::Int, threshold::Int)
    w = [[] for _ in 1:(2^n)]   # 重みベクトル

    polarity_vector = zeros(Int, n) # 極性ベクトル
    digit = zeros(Int, n)

    for i in 0:2^n-1    # 極性ベクトル分ループ
        T = [Int[] for _ in 1:(2^n)]   # 加算対象の拡張真理値ベクトル
        T_d = [Int[] for _ in 1:(2^n)]  # 依存変数用

        # i を2進数の極性ベクトルに変換
        for j in 0:(n-1)   
            polarity_vector[n-j] = div(i, 2^j) % 2 # 2^j桁目の整数を抽出
        end 

        # 加算対象ベクトル T の生成
        for j in 0:2^n-1
            # j から対応する3進数を生成
            for k in 0:(n-1)   
                if div(j, 2^k)%2 == 0 # 2^j桁目の整数を抽出
                    digit[k+1] = polarity_vector[n-k]

                elseif div(j, 2^k)%2 == 1
                    digit[k+1] = 2
                end
            end   

            # digit → t の要素のインデックスに変換
            index = 0

            for k in 0:n-1
                index += digit[k+1] * 3^(k)
            end

            T[j+1] = copy(t[index+1]) # T[j] 格納
            T_d[j+1] = copy(t_d[index+1]) # T_d[j] 格納
        end

        function extend_w(w, size::Int, target_index::Int)
            ex_w = Vector{Any}(undef, 2*length(w))
            loop_size = 2^(size - target_index)
            for m in 0:Int((length(ex_w) / (2*loop_size)))-1
                for o in 0:1
                    for p in 1:2^(size - target_index)
                        ex_w[o*loop_size + p + m*2*loop_size] = copy(w[p + m*loop_size])
                    end
                end
            end
            return ex_w
        end

        function compute_w(current_w, t::Vector, checking::Int)
            for l in 1:length(current_w)
                if current_w[l] == [Inf]
                    continue 
                end
                if t[l] == 1 
                    push!(current_w[l], checking)
                end
                if length(current_w[l]) > threshold
                    current_w[l] = [Inf]
                end
            end
            return current_w
        end

        function Construct_w(T::Vector, T_d::Vector)
            current_w = [Any[]] # 初期値
            current_d = [] # 依存変数用
            unchecked = collect(1:length(T)) # 加算対象要素番号チェック用

            #定数関数削減
            for i in reverse(1:length(T))
                if length(T_d[i]) != 0 # 依存変数がある場合
                    continue # 定数関数でないためスキップ
                end
                deleteat!(unchecked, i) # 加算対象から削除
                if T[i] == [1]
                    push!(current_w[1], i) # 登場する積項を記憶
                    sort!(current_w[1])
                end
            end 

            if length(current_w[1]) > threshold
                return [Inf]
            end

            
            while !isempty(unchecked)
                min = count(!in(current_d), T_d[unchecked[1]])
                min_index = unchecked[1]
                for i in 2:length(unchecked)
                    if count(!in(current_d), T_d[unchecked[i]]) < min
                        min = count(!in(current_d), T_d[unchecked[i]])
                        min_index = unchecked[i]
                    end
                end

                if current_d == T_d[min_index]
                    current_w = compute_w(current_w,T[min_index],min_index) # MTBDDの算術加算
                    if all(==( [Inf] ), current_w)    # 探索中止
                        return [Inf]
                    end
                else     
                    len_T_d = length(T_d[min_index])

                    if len_T_d != 0 # T_d[i] が空でない場合

                        # current_d に含まれない T_d[i] の要素を追加
                        for j in 1:len_T_d
                            if !(T_d[min_index][j] in current_d)
                                push!(current_d, T_d[min_index][j])
                                sort!(current_d)
                                
                                for k in 1:length(current_d) # 追加した依存変数が何番目に分岐するか確認
                                    if current_d[k] == T_d[min_index][j]
                                        current_w = extend_w(current_w, length(current_d), k) # 依存変数分岐を追加
                                        break
                                    end
                                end
                        
                            end
                        end
                        
                    end

                    len_current_d = length(current_d)  
                    
                    if len_current_d != 0 # current_d が空でない場合
                        # T_d[min_index] に含まれない current_d の要素を追加
                        for j in 1:len_current_d
                            if !(current_d[j] in T_d[min_index])
                                push!(T_d[min_index], current_d[j])
                                sort!(T_d[min_index])

                                for k in 1:length(T_d[min_index]) # 追加した依存変数が何番目に分岐するか確認
                                    if T_d[min_index][k] == current_d[j]
                                        T[min_index] = extend_w(T[min_index], length(T_d[min_index]), k) # 依存変数分岐を追加
                                        break
                                    end
                                end
                            end
                        end

                    end
                    current_w = compute_w(current_w,T[min_index],min_index) # MTBDDの算術加算
                    if all(==( [Inf] ), current_w)    # 探索中止
                        return [Inf]
                    end
                end

                unchecked = deleteat!(unchecked, findfirst(==(min_index), unchecked)) # 加算対象から削除
            end

            return current_w
        end

        
        # 重みベクトル w_i の構築
        w[i+1] = Construct_w(T, T_d)

    end

    return w
end

function Alg_shared(n :: Int, m :: Vector, threshold :: Int)
    
    # ---------- Step 1: 拡張真理値ベクトルを作成 ----------
    t,t_d = build_extended_truth_vector(n, m)

    # ---------- Step 2: 重みベクトルを作成 ----------
    w = build_weight_vector(t, t_d, n, threshold)

    return w
end
#=
n = 3 # 変数の個数
m = [1,1,2,0,2,1,0,1]  # 真理値ベクトル
threshold = 4 # 閾値

# ---------- 重みベクトルを作成 ----------
w = Alg_share(n, m, threshold)

for i  in 1:length(w)
    println("w_$(i-1): ", show_clean(w[i]))
end
=#