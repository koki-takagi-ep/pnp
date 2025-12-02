# Claude Code プロジェクトルール

このファイルは Claude Code がプロジェクトで作業する際のルールを定義します。

## 必須ルール

### README の更新

**図やプロットを追加・更新した場合は、必ず README.md も更新すること。**

- 新しい図を `results/` に追加したら、README の該当セクションにも追加
- 図には通し番号（図1、図2、...）をつける
- 図のキャプションには内容の説明を含める

### コミット前のチェックリスト

1. コード変更をコミット
2. **README.md が最新の状態か確認**
3. 図やデータファイルも忘れずにコミット
4. ブランチにプッシュ

## プロジェクト構成

- `src/`: C++ ソースコード
- `include/`: ヘッダファイル
- `scripts/`: Python スクリプト（プロット生成等）
- `results/`: 出力データ、図

## ビルド

```bash
make        # ビルド
make run    # 実行
```

## Python 環境

macOS では以下の Python を使用：
```bash
/opt/homebrew/opt/python@3.11/bin/python3.11
```
