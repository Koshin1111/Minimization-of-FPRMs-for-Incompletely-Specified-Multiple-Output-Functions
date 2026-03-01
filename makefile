# 変数の定義
JULIA = julia
MAIN = main.jl

# 実行されるデフォルトのターゲット
all:
	run

# 実行コマンド
run:
	$(JULIA) $(mAIN)

# パッケージのインストール
# setup:
# 	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'

# 一義ファイルの削除
# clean:
# 	rm -rf *= *.tmp