# 1次元 Poisson-Nernst-Planck ソルバー

イオン液体中の電気二重層（EDL: Electric Double Layer）をシミュレーションするための C++ 数値解析ライブラリ。

## 概要

本プロジェクトは、Poisson-Nernst-Planck (PNP) 方程式を用いて、帯電した界面近傍のイオン分布と電位分布を計算する1次元ソルバーを実装しています。

### 主な機能

- **Newton-Raphson ソルバー**: 非線形 Poisson-Boltzmann 方程式の高速求解
- **非一様グリッド**: 界面付近にグリッドを集中配置して解像度を向上
- **Gouy-Chapman 解析解**: 検証用の理論解との比較
- **Bikerman モデル**: 有限イオンサイズによる立体効果（steric effects）
- **過渡解析ソルバー**: 陰的 Scharfetter-Gummel スキームによる時間発展計算
- **誤差解析**: L2/L∞ ノルムによる精度評価

## 支配方程式

### Poisson-Nernst-Planck 方程式

電解質中のイオン輸送は以下の連立方程式で記述される：

**Poisson 方程式**（静電ポテンシャル）:

$$\nabla^2 \phi = -\frac{\rho}{\varepsilon} = -\frac{e}{\varepsilon}(z_+ c_+ + z_- c_-)$$

ここで：
- $\phi$ : 電位 [V]
- $\rho$ : 電荷密度 [C/m³]
- $\varepsilon = \varepsilon_r \varepsilon_0$ : 誘電率 [F/m]
- $e = 1.602 \times 10^{-19}$ C : 素電荷
- $z_\pm$ : イオン価数
- $c_\pm$ : イオン濃度 [mol/m³]

**Nernst-Planck 方程式**（イオンフラックス）:

$$\mathbf{J}_i = -D_i \left( \nabla c_i + \frac{z_i e c_i}{k_B T} \nabla \phi \right)$$

ここで：
- $\mathbf{J}_i$ : イオン種 $i$ のフラックス [mol/(m²·s)]
- $D_i$ : 拡散係数 [m²/s]
- $k_B = 1.381 \times 10^{-23}$ J/K : Boltzmann 定数
- $T$ : 温度 [K]

**連続の式**:

$$\frac{\partial c_i}{\partial t} = -\nabla \cdot \mathbf{J}_i$$

### 1次元定式化

$x$ 方向のみを考慮した1次元系では：

$$\frac{d^2 \phi}{d x^2} = -\frac{e}{\varepsilon}(z_+ c_+ + z_- c_-)$$

$$\frac{\partial c_\pm}{\partial t} = D_\pm \frac{\partial}{\partial x} \left( \frac{\partial c_\pm}{\partial x} \pm \frac{e c_\pm}{k_B T} \frac{\partial \phi}{\partial x} \right)$$

### 定常状態解（Boltzmann 分布）

定常状態かつフラックスゼロ（$\mathbf{J} = 0$）の条件下では、イオン濃度は Boltzmann 分布に従う：

$$c_\pm = c_0 \exp\left( \mp \frac{e \phi}{k_B T} \right)$$

これを Poisson 方程式に代入すると **Poisson-Boltzmann 方程式** が得られる：

$$\frac{d^2 \phi}{d x^2} = \frac{2 e N_A c_0}{\varepsilon} \sinh\left( \frac{e \phi}{k_B T} \right)$$

### 無次元化

**特性スケール**:

| 量 | 定義 | 物理的意味 |
|---|---|---|
| Debye 長 | $\lambda_D = \sqrt{\dfrac{\varepsilon k_B T}{2 e^2 c_0 N_A}}$ | 遮蔽長さ |
| 熱電圧 | $\phi_T = \dfrac{k_B T}{e} \approx 25.7$ mV (298 K) | 熱エネルギー |
| バルク濃度 | $c_0$ | 基準濃度 |

**無次元変数**:

$$\xi = \frac{x}{\lambda_D}, \quad \psi = \frac{\phi}{\phi_T} = \frac{e \phi}{k_B T}$$

**無次元 Poisson-Boltzmann 方程式**:

$$\frac{d^2 \psi}{d \xi^2} = \sinh(\psi)$$

### Gouy-Chapman 解析解

1:1 電解質において、表面電位 $\psi_0$ が与えられたときの解析解：

$$\tanh\left( \frac{\psi}{4} \right) = \tanh\left( \frac{\psi_0}{4} \right) \exp(-\xi)$$

次元量に戻すと：

$$\phi(x) = \frac{4 k_B T}{e} \tanh^{-1} \left[ \tanh\left( \frac{e \phi_0}{4 k_B T} \right) \exp\left( -\frac{x}{\lambda_D} \right) \right]$$

**導出**:

Poisson-Boltzmann 方程式を一度積分すると：

$$\left( \frac{d\psi}{d\xi} \right)^2 = 2(\cosh\psi - 1) = 4\sinh^2\left( \frac{\psi}{2} \right)$$

したがって：

$$\frac{d\psi}{d\xi} = -2\sinh\left( \frac{\psi}{2} \right)$$

（負号は電位が界面から離れるにつれ減少することを反映）

変数分離して積分：

$$\int_{\psi_0}^{\psi} \frac{d\psi'}{2\sinh(\psi'/2)} = -\int_0^{\xi} d\xi'$$

$$\ln\left| \tanh\left( \frac{\psi}{4} \right) \right| - \ln\left| \tanh\left( \frac{\psi_0}{4} \right) \right| = -\xi$$

これより Gouy-Chapman 解が得られる。

### Bikerman モデル（修正 Poisson-Boltzmann）

標準 Poisson-Boltzmann 方程式は、高電位領域でイオン濃度が非物理的に高くなる問題がある。Bikerman モデルは有限イオンサイズを考慮することでこの問題を解決する。

**修正 Boltzmann 分布**:

$$c_\pm = \frac{c_0 \exp(\mp \psi)}{g(\psi)}$$

**混雑関数（crowding function）**:

$$g(\psi) = 1 - \nu + \nu \cosh(\psi)$$

**充填率（packing fraction）**:

$$\nu = 2 a^3 c_0 N_A$$

ここで $a$ はイオン直径（イオン液体では典型的に 0.5〜1.0 nm）。

**修正 Poisson-Boltzmann 方程式**:

$$\frac{d^2 \psi}{d \xi^2} = \frac{\sinh(\psi)}{g(\psi)}$$

**物理的解釈**:

- 低電位（$|\psi| \ll 1$）: $g(\psi) \approx 1$ となり標準 PB に帰着
- 高電位: $g(\psi)$ により最大濃度が $\sim c_0/\nu$ に制限
- 非物理的な濃度 $c > 1/(a^3 N_A)$ を防止

**参考文献**: Kilic, Bazant & Ajdari, *Phys. Rev. E* 75, 021502 (2007)

## 数値解法

### 空間離散化

**非一様グリッド**:

界面付近（$x = 0$）にグリッド点を集中配置：

$$x_i = L \left[ 1 - (1 - \xi_i)^\beta \right]$$

ここで $\xi_i = i/(N-1) \in [0, 1]$、$\beta > 1$ はストレッチング係数。

**2階微分（非一様グリッド）**:

内点における離散化：

$$\frac{d^2 \phi}{d x^2} \bigg|_i \approx \frac{(\phi_{i+1} - \phi_i)/\Delta x^+ - (\phi_i - \phi_{i-1})/\Delta x^-}{\Delta x_{\text{avg}}}$$

ここで：
- $\Delta x^+ = x_{i+1} - x_i$
- $\Delta x^- = x_i - x_{i-1}$
- $\Delta x_{\text{avg}} = (\Delta x^+ + \Delta x^-) / 2$

### Newton-Raphson 法

非線形 Poisson-Boltzmann 方程式を Newton-Raphson 反復で解く。

**残差**:

$$F(\phi) = \frac{d^2 \phi}{d x^2} - \kappa^2 \phi_T \sinh\left( \frac{\phi}{\phi_T} \right) = 0$$

ここで $\kappa = 1/\lambda_D$。

**Jacobian**:

$$J = \frac{dF}{d\phi} = \frac{d^2}{d x^2} - \kappa^2 \cosh\left( \frac{\phi}{\phi_T} \right)$$

**Newton 更新**:

$$J \cdot \delta\phi = -F$$

$$\phi^{n+1} = \phi^n + \alpha \cdot \delta\phi$$

ここで $\alpha \in (0, 1]$ は安定性のための適応的ダンピング係数。

**三重対角行列の解法（Thomas アルゴリズム）**:

Newton 法で現れる連立一次方程式 $J \cdot \delta\phi = -F$ は三重対角形式であり、$O(N)$ の計算量で効率的に解ける。

### 過渡解析ソルバー

過渡解析では陰的スキームを使用：

**時間離散化**:

$$\frac{c_i^{n+1} - c_i^n}{\Delta t} = D \frac{\partial}{\partial x} \left[ \frac{\partial c}{\partial x} + \frac{z e c}{k_B T} \frac{\partial \phi}{\partial x} \right]$$

**Scharfetter-Gummel スキーム**:

ドリフト-拡散問題に対して数値的に安定な離散化を提供：

$$J_{i+1/2} = \frac{D}{\Delta x} \left[ B(v \Delta x) c_{i+1} - B(-v \Delta x) c_i \right]$$

ここで $B(x) = x/(\exp(x) - 1)$ は Bernoulli 関数、$v = zeE/(k_B T)$。

**Bernoulli 関数の性質**:

- $B(0) = 1$
- $B(x) + B(-x) = x$
- $x \to 0$ で $B(x) \approx 1 - x/2 + x^2/12$

**各時間ステップのアルゴリズム**:

1. 準静的 Poisson 方程式を解いて $\phi$ を更新
2. 陰的 Nernst-Planck で $c_+$ を更新（Scharfetter-Gummel）
3. 陰的 Nernst-Planck で $c_-$ を更新（Scharfetter-Gummel）
4. 定常状態判定（$\Delta c / c_0 < $ 許容誤差）

**特性時間スケール**:

$$\tau_D = \frac{\lambda_D^2}{D} \approx 0.1 - 1 \text{ ns}$$（典型的なイオン液体）

## 検証結果

### 格子収束性解析

数値解の精度を検証するため、Gouy-Chapman 解析解との比較による格子収束性解析を実施した。

#### 誤差評価指標

**L2 誤差ノルム**（二乗平均平方根誤差）:

$$L_2 = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (\phi_i^{\text{num}} - \phi_i^{\text{exact}})^2}$$

**収束次数** $p$ の算出:

$$p = \frac{\log(L_2^{(1)} / L_2^{(2)})}{\log(N^{(2)} / N^{(1)})}$$

ここで上付き添字 $(1), (2)$ は粗いグリッドと細かいグリッドをそれぞれ示す。

#### 収束性テスト条件

- バルク濃度: $c_0 = 0.01$ mol/L（Debye 長 ≈ 3.76 nm）
- 表面電位: $\phi_0 = 100$ mV
- 計算領域: $L = 100$ nm
- グリッド: 一様グリッド（`--stretch 1.0`）

#### 結果

<div align="center">

![Grid Convergence](results/grid_convergence.png)

*図: 格子収束性解析。横軸はセル数 N、縦軸は L2 誤差 [mV]。破線は1次および2次精度の参照勾配。*

</div>

| グリッド点数 N | グリッド幅 Δx [nm] | L2 誤差 [mV] | 収束次数 p |
|:--------------:|:------------------:|:------------:|:----------:|
| 51 | 2.0000 | 0.8638 | — |
| 101 | 1.0000 | 0.5046 | 0.78 |
| 201 | 0.5000 | 0.2120 | 1.25 |
| 401 | 0.2500 | 0.0694 | 1.61 |
| 801 | 0.1250 | 0.0193 | 1.85 |
| 1601 | 0.0625 | 0.0050 | 1.95 |

**平均収束次数: 1.49**（漸近的に2次精度に収束）

#### 考察

- グリッドを2倍に細かくすると、誤差は約4分の1に減少（2次精度の特徴）
- 粗いグリッドでは漸近領域に達していないため収束次数が低い
- $N \geq 400$ で理論的な2次精度（$p \approx 2$）が達成される
- 最終的な L2 誤差は 0.005 mV（熱電圧の約0.02%）と極めて高精度

### Gouy-Chapman 理論との比較

**テスト条件**: $c_0 = 0.1$ M, $\phi_0 = 100$ mV, $\varepsilon_r = 12$

| パラメータ | 値 |
|:-----------|------:|
| Debye 長 | 0.376 nm |
| 熱電圧 | 25.7 mV |
| 無次元電位 $\psi_0$ | 3.89 |
| L2 誤差 | ~1.5 mV |
| 相対 L2 誤差 | ~3.5% |
| L∞ 誤差 | ~4 mV |
| 収束反復数 | 4 回 |

## ビルドと実行

### 必要環境

- C++17 対応コンパイラ（g++ 推奨）
- Python 3 + NumPy + Matplotlib（可視化用）

### ビルド

```bash
make
```

### 実行

```bash
# デフォルトパラメータ（1 M, 100 mV, 定常状態）
make run

# 標準 Poisson-Boltzmann
./build/pnp_solver --c0 0.1 --phi0 100 --output results/standard.dat

# Bikerman モデル（立体効果あり）
./build/pnp_solver --c0 0.1 --phi0 100 --model bikerman --ion-size 0.7 --output results/bikerman.dat

# 一様グリッド
./build/pnp_solver --stretch 1.0 --output results/uniform.dat

# 過渡解析
./build/pnp_solver --transient --dt 0.1 --t-final 1.0 --output results/transient.dat
```

### コマンドラインオプション

| オプション | 説明 | デフォルト値 |
|-----------|------|-------------|
| `--phi0 <value>` | 表面電位 [mV] | 100 |
| `--c0 <value>` | バルク濃度 [mol/L] | 1.0 |
| `--eps <value>` | 比誘電率 | 12 |
| `--L <value>` | 計算領域長 [nm] | 100 |
| `--N <value>` | グリッド点数 | 1001 |
| `--stretch <value>` | グリッドストレッチング係数 | 3.0 |
| `--model <type>` | モデル種類（standard/bikerman） | standard |
| `--ion-size <value>` | イオン直径 [nm]（Bikerman用） | 0.7 |
| `--transient` | 過渡解析モード | -- |
| `--dt <value>` | 時間刻み [ns] | 0.1 |
| `--t-final <value>` | 終了時間 [µs] | 1.0 |
| `--output <file>` | 出力ファイル名 | results/pnp_results.dat |

### 可視化

```bash
# 結果のプロット
python3 scripts/plot_results.py

# 収束性解析
bash scripts/run_convergence.sh
python3 scripts/plot_convergence.py
```

## ファイル構成

```
pnp/
├── include/
│   └── pnp_solver.hpp      # ヘッダファイル（クラス定義）
├── src/
│   ├── pnp_solver.cpp      # ソルバー実装
│   └── main.cpp            # メインプログラム
├── scripts/
│   ├── plot_results.py     # 結果可視化スクリプト
│   ├── plot_convergence.py # 収束性プロット
│   ├── plot_parametric.py  # パラメトリックスタディ用
│   └── run_convergence.sh  # 収束性テスト実行スクリプト
├── results/                # 出力データ・図
├── Makefile
└── README.md
```

## 今後の課題

1. **Carnahan-Starling モデル**: 硬球の状態方程式をより正確に記述
2. **2D/3D 拡張**: より複雑な形状への対応
3. **印加電圧**: 電流が流れる非平衡定常状態
4. **非対称イオン**: カチオンとアニオンの異なるサイズ
5. **多成分系**: 複数のイオン種を含む電解質

## 参考文献

1. Newman, J., & Thomas-Alyea, K. E. (2004). *Electrochemical Systems* (3rd ed.). Wiley.

2. Bazant, M. Z., Kilic, M. S., Storey, B. D., & Ajdari, A. (2009). Towards an understanding of induced-charge electrokinetics at large applied voltages in concentrated solutions. *Advances in Colloid and Interface Science*, 152(1-2), 48-88.

3. Kilic, M. S., Bazant, M. Z., & Ajdari, A. (2007). Steric effects in the dynamics of electrolytes at large applied voltages. *Physical Review E*, 75(2), 021502.

4. Bazant, M. Z., Storey, B. D., & Kornyshev, A. A. (2011). Double layer in ionic liquids: Overscreening versus crowding. *Physical Review Letters*, 106(4), 046102.

## ライセンス

MIT License
