function pla_to_m(filename::String)
    n_inputs  = 0
    n_outputs = 0
    m = Vector{Vector{Int}}()

    function expand_input(pattern)
        lists = [[]]
        for c in pattern
            if c == '-'
                lists = vcat(
                    [vcat(l, 0) for l in lists],
                    [vcat(l, 1) for l in lists]
                )
            elseif c == '0' || c == '1'
                lists = [vcat(l, parse(Int, string(c))) for l in lists]
            else
                error("Invalid character '$c' in input pattern")
            end
        end
        return lists
    end

    function binvec_to_index(v)
        idx = 0
        for b in v
            idx = 2 * idx + b
        end
        return idx + 1
    end

    for rawline in eachline(filename)
        line = strip(rawline)
        isempty(line) && continue

        # ---- 制御行 ----
        if startswith(line, ".i ")
            n_inputs = parse(Int, split(line)[2])
            continue
        elseif startswith(line, ".o ")
            n_outputs = parse(Int, split(line)[2])
            m = [fill(2, 2^n_inputs) for _ in 1:n_outputs]
            continue
        elseif startswith(line, ".e")
            break
        end

        # ---- データ行かを判定 ----
        parts = split(line)
        length(parts) == 2 || continue

        ins, outs = parts

        # 入力・出力が 0/1/- のみか確認
        all(c -> c in ('0','1','-'), ins) || continue
        all(c -> c in ('0','1','-'), outs) || continue

        # ---- 正当なデータ行のみここに到達 ----
        for v in expand_input(ins)
            idx = binvec_to_index(v)
            for j in 1:n_outputs
                c = outs[j]
                m[j][idx] = (c == '-' ? 2 : parse(Int, string(c)))
            end
        end
    end

    return n_inputs, n_outputs, m
end
