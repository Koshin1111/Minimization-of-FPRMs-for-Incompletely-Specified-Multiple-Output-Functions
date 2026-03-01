# OffsetArray / ネスト配列を再帰的に普通の Vector と Int のみに変換（型注釈を除去）
function strip_offset(x)
    if x isa Int
        return x
    elseif x isa AbstractArray
        raw = x isa OffsetArray ? parent(x) : x
        return [strip_offset(xi) for xi in raw]
    else
        return x
    end
end

function show_clean(x)
    if x isa Int
        return string(x)
    elseif x isa AbstractArray
        return "[" * join(show_clean.(x), ", ") * "]"
    else
        return string(x)
    end
end

# ---------- PPRM可視化 ----------
function print_PPRM(n ::Int, m ::Int, T :: Vector{Vector{Vector{Int}}}, T_d :: Vector{Vector{Vector{Int}}}, assign_d :: Vector)
    PPRM = [zeros(Int, 2^n) for _ in 1:m]

    for k in 1:m
        all_d = copy(T_d[k][2^n])  # すべての依存変数
        for i in 1:2^n
            assign_d_i = []
            if !isempty(T_d[k][i])
                assign_d_i = zeros(Int, length(T_d[k][i]))
                for j in 1:length(T_d[k][i])
                    # 全ての依存変数を assign_d_i に変換
                    for l in 1:length(assign_d[k])
                        if T_d[k][i][j] == all_d[l]
                            assign_d_i[j] = assign_d[k][l]
                            break
                        end
                    end
                end
            end
            assign_index = 1
            for j in 1:length(assign_d_i)
                assign_index += assign_d_i[j] * 2^(j-1)
            end
            PPRM[k][i] = T[k][i][assign_index]
        end
    end

    X = fill("", n)
    for i in 1:n
        X[i] = "x$(i)"
    end

    for k in 1:m
        PPRM_term = ""
        for i in 1:2^n
            if PPRM[k][i] == 0
                continue
            else
                if !isempty(PPRM_term)
                    PPRM_term *= " ⊕ "
                end
                # i を2進数に変換
                if i == 1
                    PPRM_term *= "1"
                    continue
                end
                for j in reverse(0:n-1)
                    if div(i-1, 2^j) % 2 == 1# 2^j桁目の整数を抽出
                        PPRM_term *= X[n-j]
                    end
                end 
            end
        end

        println("f", k-1, ": ", PPRM_term)
    end
    
end


# ---------- FPRM可視化 ----------
function print_FPRM(n ::Int, m ::Int, w :: Vector{Any}, assign_d :: Vector{Vector{Int}}, polarity_vector :: Vector{Vector{Int}})
    FPRM = [zeros(Int, 2^n) for _ in 1:m]

    for f in 1:m
        # 極性ベクトルを配列用インデックスに変換
        pola_idx = 1
        for i in 1:n 
            pola_idx += polarity_vector[f][i] * 2^(n-i)
        end
        # 割り当てベクトルを配列用インデックスに変換
        asg_idx = 1
        for i in 1:length(assign_d[f])
            asg_idx += assign_d[f][i] * 2^(length(assign_d[f])-i)
        end
        for i in 1:length(w[f][pola_idx][asg_idx])
            # 登場する積項を 1 として FPRM に記憶
            FPRM[f][w[f][pola_idx][asg_idx][i]] = 1
        end
    end


    for f in 1:m
        FPRM_term = ""
        for i in 1:2^n
            if FPRM[f][i] == 0
                continue
            else
                if !isempty(FPRM_term)
                    FPRM_term *= " ⊕ "
                end
                # 変数を含まない積項
                if i == 1
                    FPRM_term *= "1"
                    continue
                end
                for j in reverse(0:n-1)
                    if div(i-1, 2^j) % 2 == 1# 2^j桁目の整数を抽出
                        if polarity_vector[f][n-j] == 0
                            FPRM_term *= "x$(n-j)"
                        elseif polarity_vector[f][n-j] == 1
                            FPRM_term *= "̄x$(n-j)"
                        end
                    end
                end 
            end
        end

        println("f", f-1, ": ", FPRM_term)
    end
    
end