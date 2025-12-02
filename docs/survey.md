# 非定常PNP方程式の数値解法サーベイ（安定性重視）

本メモでは手元にある論文PDFと`.bib`情報を基に、過去約10年の非定常（時間依存）Poisson–Nernst–Planck (PNP) 系に対する数値解法を安定性の観点で整理する。主眼は解の非負性・エネルギー減衰・質量保存などの構造保持と、時間積分の取り扱い（陽／陰／半陰・ETD など）。

正値保存（positivity-preserving）とは、本来負にならない量（濃度・電荷密度など）を離散化・時間発展の各ステップでも非負に保つ性質を指す。PNPでは非物理的な負濃度を防ぎ、離散最大値原理やエネルギー減衰と並んで安定性を特徴付ける重要な指標となる。

## 有限差分系（FD）
- **ding_2019_jcp (JCP 397, 108864, 2019, DOI: 10.1016/j.jcp.2019.108864)**  
  Steric 効果を含む修正PNP。調和平均フラックスによる有限差分で濃度の非負性を保持。スキームはポテンシャルと濃度を適合させつつ、数値的安定性を重視。
- **liu_2023_jsc (J. Sci. Comput. 97:23, 2023, DOI: 10.1007/s10915-023-02345-9)**  
  2次精度の陽・陰を組み合わせた構造保持型スキーム。濃度の非負性とエネルギー安定性（収束解析付き）を理論保証。
- **shen_2021_nummath (Numer. Math. 148:671–697, 2021, DOI: 10.1007/s00211-021-01203-w)**  
  非負性・エネルギー減衰を両立する“unconditionally positivity preserving and energy dissipative”差分スキーム。Slotboom 変数変換を用いた半陰的離散化が特徴。

## 有限要素・時空間高次法（FEM / space–time）
- **fu_2022_cma (Comput. Methods Appl. Mech. Eng. 395:115031, 2022, DOI: 10.1016/j.cma.2022.115031)**  
  時空間高次FEM (space–time FEM) による任意次数の構造保持スキーム。質量保存・非負性・無条件エネルギー安定を証明。時間方向はDG型の高次積分で、最低次は後退オイラーに一致。

## 指数時間積分（ETD）・構造保存型
- **guo_2024_arxiv (arXiv:2410.00306, DOI: 10.48550/arXiv.2410.00306)**  
  Maxwell–Ampère–PNP を対象とした構造保存型 ETD スキーム。Slotboom 変換とETDを組み合わせ、濃度非負性とエネルギー減衰をねらった半陰的構成。

## まとめと今後の補完
- 手元PDFに基づき、FD系3本、space–time FEM 1本、ETD系1本を整理。いずれも非負性・エネルギー安定性を強調する近年の代表例。
- references.bib に挙がるが未取得の論文（例: qian_2024_sjna, hu_2024_cnsns など）は、アクセス可能になり次第追記予定。

## 対応するファイル
- ding_2019_jcp.pdf / ding_2019_jcp.bib
- liu_2023_jsc.pdf / liu_2023_jsc.bib
- shen_2021_nummath.pdf / shen_2021_nummath.bib
- fu_2022_cma.pdf / fu_2022_cma.bib
- guo_2024_arxiv.pdf / guo_2024_arxiv.bib
