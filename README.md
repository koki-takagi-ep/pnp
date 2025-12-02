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

## モデルの仮定と適用範囲

本ソルバーは以下の仮定に基づいている。これらの仮定の妥当性を理解した上で使用することが重要である。

### 物理的仮定

| 仮定 | 説明 | 影響・制限 |
|------|------|-----------|
| **1. 連続体仮定** | イオンおよび溶媒を連続体として扱う | Debye長が分子サイズ（〜0.5 nm）より十分大きい系に適用可能。極めて高濃度の電解質（> 5 M）では分子論的効果が顕著になる |
| **2. 平均場近似** | イオン間の相互作用を平均電場を通じてのみ考慮 | イオン間の直接的な相関（静電相関、硬球排除）を無視。高価数イオンや高濃度系では誤差が大きくなる |
| **3. 希薄溶液仮定**（標準PB） | イオンを点電荷として扱う | 高電位領域で非物理的な高濃度を予測。Bikermanモデルで部分的に改善 |
| **4. 誘電率一定** | 電解質全体で誘電率が均一かつ一定 | 電気二重層内での誘電飽和（高電場による誘電率低下）を無視 |
| **5. 温度一定** | 等温条件 | Joule加熱や温度勾配を無視。急速な充放電過程では不正確になる可能性 |
| **6. 1:1対称電解質** | カチオンとアニオンが同じ価数（±1）と拡散係数を持つ | 非対称電解質（異なる価数、サイズ、拡散係数）には直接適用不可 |
| **7. 準静的電場** | 過渡解析でもPoisson方程式は各時刻で定常状態 | 電荷緩和が電場の変化より十分速い場合に妥当（RC時定数 << 拡散時定数） |

### 境界条件の仮定

| 境界 | 条件 | 物理的意味 |
|------|------|-----------|
| **左端（x = 0）** | Dirichlet（電位固定）、ゼロフラックス（濃度） | 理想的な阻止電極、イオンは電極に吸着・脱着しない |
| **右端（x = L）** | Dirichlet（φ = 0, c = c₀） | 十分遠方でバルク条件が成立、L >> λ_D が必要 |

### 数値的仮定

| 仮定 | 説明 |
|------|------|
| **1次元系** | 電極は無限平面、界面に垂直な方向のみ変化 |
| **単一界面** | 一方の電極のみを考慮、対向電極からの影響なし |

### 各モデルの適用範囲

**標準 Poisson-Boltzmann**:
- 低〜中程度の電位（$|\phi| \lesssim 100$ mV）
- 低〜中程度の濃度（$c_0 \lesssim 0.1$ M）
- Gouy-Chapman理論と同等

**Bikerman モデル**:
- より高電位の系に適用可能
- 有限イオンサイズ効果を考慮
- ただし、静電相関は考慮されていない

### 本モデルでは考慮されていない効果

1. **静電相関**: 高価数イオンや高濃度系で重要
2. **誘電飽和**: 高電場での誘電率低下
3. **比吸着**: イオンの電極表面への特異吸着
4. **溶媒構造**: 溶媒分子の配向や構造化
5. **画像電荷効果**: 電極界面でのイオンの鏡像電荷
6. **非電気化学的相互作用**: van der Waals力など

### 妥当性の目安

| パラメータ | 推奨範囲 | 備考 |
|-----------|---------|------|
| 濃度 $c_0$ | 0.001〜1 M | 高濃度では平均場近似が破綻 |
| 表面電位 $\phi_0$ | $\lesssim 4\phi_T \approx 100$ mV | 高電位ではBikermanモデル推奨 |
| イオン価数 | ±1 | 多価イオンでは静電相関が重要 |
| 計算領域 $L$ | $\geq 10 \lambda_D$ | バルク条件の成立に必要 |

## 数値解法

### 離散化手法

本ソルバーでは、**2次精度中心差分法**（2nd-order central finite difference method）を用いて空間離散化を行う。

#### 解く方程式

定常状態の Poisson-Boltzmann 方程式を離散化する：

$$\frac{d^2 \phi}{dx^2} = \kappa^2 \phi_T \sinh\left(\frac{\phi}{\phi_T}\right)$$

ここで $\kappa = 1/\lambda_D$ は逆 Debye 長、$\phi_T = k_B T/e$ は熱電圧。

#### 空間離散化

**非一様グリッド**:

界面付近（$x = 0$）にグリッド点を集中配置：

$$x_i = L \left[ 1 - (1 - \xi_i)^\beta \right]$$

ここで $\xi_i = i/(N-1) \in [0, 1]$、$\beta > 1$ はストレッチング係数。$\beta = 1$ で一様グリッドとなる。

**2階微分の離散化（2次精度）**:

非一様グリッド上での2階微分は、Taylor展開に基づく以下の式で離散化される：

$$\frac{d^2 \phi}{d x^2} \bigg|_i = \frac{2}{\Delta x^- + \Delta x^+} \left[ \frac{\phi_{i+1} - \phi_i}{\Delta x^+} - \frac{\phi_i - \phi_{i-1}}{\Delta x^-} \right] + O(h^2)$$

ここで：
- $\Delta x^+ = x_{i+1} - x_i$（前方グリッド幅）
- $\Delta x^- = x_i - x_{i-1}$（後方グリッド幅）

**2次精度の証明**:

Taylor展開により：
$$\phi_{i+1} = \phi_i + \Delta x^+ \phi'_i + \frac{(\Delta x^+)^2}{2}\phi''_i + \frac{(\Delta x^+)^3}{6}\phi'''_i + O(h^4)$$
$$\phi_{i-1} = \phi_i - \Delta x^- \phi'_i + \frac{(\Delta x^-)^2}{2}\phi''_i - \frac{(\Delta x^-)^3}{6}\phi'''_i + O(h^4)$$

これらを組み合わせると、上記の離散化公式の打ち切り誤差が $O(h^2)$ となることが示される。

**離散化された方程式系**:

内点（$i = 1, 2, \ldots, N-2$）において：

$$\frac{2}{\Delta x^- + \Delta x^+} \left[ \frac{\phi_{i+1} - \phi_i}{\Delta x^+} - \frac{\phi_i - \phi_{i-1}}{\Delta x^-} \right] = \kappa^2 \phi_T \sinh\left(\frac{\phi_i}{\phi_T}\right)$$

境界条件（Dirichlet）：
- $\phi_0 = \phi_{\text{surface}}$（表面電位）
- $\phi_{N-1} = 0$（バルク、電気的中性）

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

- バルク濃度: $c_0 = 0.001$ mol/L（Debye 長 ≈ 11.9 nm）
- 表面電位: $\phi_0 = 50$ mV
- 計算領域: $L = 100$ nm
- グリッド: 一様グリッド（`--stretch 1.0`）

#### 結果

<div align="center">

![Grid Convergence](results/grid_convergence.png)

*図: 格子収束性解析。横軸はグリッド幅 Δx [nm]、縦軸は L2 誤差 [mV]。破線は1次精度、点線は2次精度の参照勾配を示す。*

</div>

| グリッド点数 N | グリッド幅 Δx [nm] | L2 誤差 [mV] | 収束次数 p |
|:--------------:|:------------------:|:------------:|:----------:|
| 51 | 2.0000 | 1.064e-01 | — |
| 101 | 1.0000 | 2.994e-02 | 1.86 |
| 201 | 0.5000 | 7.758e-03 | 1.96 |
| 401 | 0.2500 | 1.959e-03 | 1.99 |
| 801 | 0.1250 | 4.911e-04 | 2.00 |
| 1601 | 0.0625 | 1.229e-04 | 2.00 |

**収束次数（最後の3点平均）: 2.00**（理論値と完全に一致）

#### 考察

- グリッドを2倍に細かくすると、誤差は約4分の1に減少（2次精度の特徴）
- $N \geq 201$ で理論的な2次精度（$p \approx 2$）を達成
- 最終的な L2 誤差は 0.123 µV（熱電圧の約0.0005%）と極めて高精度
- 2次精度中心差分法の理論通りの収束性を確認

### 計算結果

#### 電気二重層構造

<div align="center">

![Combined Results](results/combined_results.png)

*図: 電気二重層の数値解析結果。(a) 電位分布（PNP数値解とGouy-Chapman解析解の比較）、(b) 無次元電位、(c) イオン濃度分布（対数スケール）、(d) 空間電荷密度。*

</div>

#### 電位・濃度プロファイル

| プロット | 説明 |
|----------|------|
| ![Potential](results/potential_profile.png) | 電位プロファイル：数値解と解析解の比較 |
| ![Concentration](results/concentration_profiles.png) | イオン濃度：カチオン減少、アニオン増加 |

### 過渡解析（電気二重層の形成過程）

陰的 Scharfetter-Gummel スキームによる真の時間発展計算で、**100 mV ステップ応答**の EDL 形成過程をシミュレーションする。

<div align="center">

![EDL Evolution](results/edl_evolution.gif)

*図: EDL時間発展のアニメーション（0〜200 ns）。左: 電位分布（青: 数値解、赤破線: Gouy-Chapman解析解）、右: イオン濃度分布（赤: カチオン c+、青: アニオン c-）。t=0 で 100 mV をステップ印加し、初期状態（一様濃度）から定常状態への緩和過程を表示。*

</div>

**シミュレーション条件**:
- 表面電位: 100 mV ステップ印加（t=0 で瞬時に印加）
- バルク濃度: $c_0 = 1.0$ mol/L
- 時間刻み: $\Delta t = 10$ ps
- 計算時間: 200 ns（20,000ステップ）
- 空間離散化: Scharfetter-Gummel スキーム
- 時間積分: 陰的 Euler 法 + Gummel 反復

**物理パラメータ**:
- Debye長: $\lambda_D = 0.12$ nm
- 熱電圧: $\phi_T = 25.7$ mV
- 無次元電位: $\psi_0 = \phi_0/\phi_T \approx 3.9$
- Debye時間: $\tau_D = \lambda_D^2/D \approx 140$ ps

**定常状態の濃度**（約70 ns で収束）:
- 表面での c+/c₀ ≈ 0.022（カチオン排除）
- 表面での c-/c₀ ≈ 37（アニオン蓄積）

電位印加直後から急速にイオン再分布が始まり、約70 ns で定常状態に収束する。高濃度イオン液体（1 M）では Debye 長が極めて短く（0.12 nm）、EDL 内の電場が非常に強い（~1 GV/m）ため、ps オーダーのタイムステップが必要となる。

### Gouy-Chapman 理論との比較

**テスト条件**: $c_0 = 0.001$ M, $\phi_0 = 50$ mV, $\varepsilon_r = 12$

| パラメータ | 値 |
|:-----------|------:|
| Debye 長 | 11.9 nm |
| 熱電圧 | 25.7 mV |
| 無次元電位 $\psi_0$ | 1.95 |
| L2 誤差 | 0.0066 mV |
| 相対 L2 誤差 | 0.16 % |
| L∞ 誤差 | 0.050 mV |
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
| `--animation` | アニメーション用スナップショット保存 | -- |
| `--snapshot-dir <dir>` | スナップショット出力先 | results/snapshots |
| `--snapshot-interval <n>` | スナップショット間隔（ステップ数） | 10 |
| `--output <file>` | 出力ファイル名 | results/pnp_results.dat |

### 可視化

```bash
# 結果のプロット
python3 scripts/plot_results.py

# 収束性解析
bash scripts/run_convergence.sh
python3 scripts/plot_convergence.py

# 過渡解析アニメーション生成
./build/pnp_solver --animation --phi0 100 --c0 0.1 --dt 1.0 --t-final 5 --snapshot-interval 50
python3 scripts/create_animation.py
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
│   ├── create_animation.py # GIFアニメーション生成
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
