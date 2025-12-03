# 理論と数値解法

本ドキュメントでは、1D PNP ソルバーで使用している理論的背景と数値解法について詳細に解説する。

## 目次

1. [支配方程式](#支配方程式)
2. [Poisson-Boltzmann方程式の導出](#poisson-boltzmann方程式の導出)
3. [Gouy-Chapman解析解](#gouy-chapman解析解)
4. [Bikermanモデル](#bikermanモデル)
5. [数値解法](#数値解法)
6. [境界条件](#境界条件)

---

## 支配方程式

### Poisson-Nernst-Planck 方程式

電解質中のイオン輸送は以下の連立方程式で記述される：

**Poisson 方程式**（静電ポテンシャル）:

$$\nabla^2 \phi = -\frac{\rho}{\varepsilon} = -\frac{e}{\varepsilon}(z_+ c_+ + z_- c_-)$$

ここで：
- φ : 電位 [V]
- ρ : 電荷密度 [C/m³]
- ε = εᵣε₀ : 誘電率 [F/m]
- e = 1.602 × 10⁻¹⁹ C : 素電荷
- z± : イオン価数
- c± : イオン濃度 [mol/m³]

**Nernst-Planck 方程式**（イオンフラックス）:

$$\mathbf{J}_i = -D_i \left( \nabla c_i + \frac{z_i e c_i}{k_B T} \nabla \phi \right)$$

ここで：
- **J**ᵢ : イオン種 i のフラックス [mol/(m²·s)]
- Dᵢ : 拡散係数 [m²/s]
- kB = 1.381 × 10⁻²³ J/K : Boltzmann 定数
- T : 温度 [K]

**連続の式**:

$$\frac{\partial c_i}{\partial t} = -\nabla \cdot \mathbf{J}_i$$

### 1次元定式化

x 方向のみを考慮した1次元系では：

$$\frac{d^2 \phi}{d x^2} = -\frac{e}{\varepsilon}(z_+ c_+ + z_- c_-)$$

$$\frac{\partial c_\pm}{\partial t} = D_\pm \frac{\partial}{\partial x} \left( \frac{\partial c_\pm}{\partial x} \pm \frac{e c_\pm}{k_B T} \frac{\partial \phi}{\partial x} \right)$$

---

## Poisson-Boltzmann方程式の導出

**Poisson-Boltzmann (PB) 方程式は、PNP 方程式系の定常状態解として導出される。**

### 導出の流れ

1. **定常状態条件**: ∂c/∂t = 0 より、連続の式から ∇·**J** = 0
2. **ゼロフラックス条件**: 阻止電極（イオンが透過しない境界）では **J** = 0
3. **Nernst-Planck 方程式にゼロフラックス条件を適用**:

$$\mathbf{J} = -D \left( \nabla c + \frac{ze c}{k_B T} \nabla \phi \right) = 0$$

これを整理すると：

$$\frac{\nabla c}{c} = -\frac{ze}{k_B T} \nabla \phi$$

4. **積分して Boltzmann 分布を得る**:

$$\ln c - \ln c_0 = -\frac{ze}{k_B T}(\phi - \phi_{\text{bulk}})$$

バルク条件（φ_bulk = 0, c = c₀）を基準として：

$$c_\pm = c_0 \exp\left( \mp \frac{e \phi}{k_B T} \right)$$

5. **Poisson 方程式に代入して PB 方程式を得る**:

$$\frac{d^2 \phi}{d x^2} = -\frac{e NA}{\varepsilon}(c_+ - c_-) = \frac{2 e NA c_0}{\varepsilon} \sinh\left( \frac{e \phi}{k_B T} \right)$$

### 無次元化

**特性スケール**:

| 量 | 定義 | 物理的意味 |
|---|---|---|
| Debye 長 | λD = √(εkBT / 2e²c₀NA) | 遮蔽長さ |
| 熱電圧 | φT = kBT/e ≈ 25.7 mV (298 K) | 熱エネルギー |
| バルク濃度 | c₀ | 基準濃度 |

**無次元変数**:

$$\xi = \frac{x}{\lambda_D}, \quad \psi = \frac{\phi}{\phi_T} = \frac{e \phi}{k_B T}$$

**無次元 Poisson-Boltzmann 方程式**:

$$\frac{d^2 \psi}{d \xi^2} = \sinh(\psi)$$

---

## Gouy-Chapman解析解

1:1 電解質において、表面電位 ψ₀ が与えられたときの解析解：

$$\tanh\left( \frac{\psi}{4} \right) = \tanh\left( \frac{\psi_0}{4} \right) \exp(-\xi)$$

次元量に戻すと：

$$\phi(x) = \frac{4 k_B T}{e} \tanh^{-1} \left[ \tanh\left( \frac{e \phi_0}{4 k_B T} \right) \exp\left( -\frac{x}{\lambda_D} \right) \right]$$

### 導出

Poisson-Boltzmann 方程式を一度積分すると：

$$\left( \frac{d\psi}{d\xi} \right)^2 = 2(\cosh\psi - 1) = 4\sinh^2\left( \frac{\psi}{2} \right)$$

したがって：

$$\frac{d\psi}{d\xi} = -2\sinh\left( \frac{\psi}{2} \right)$$

（負号は電位が界面から離れるにつれ減少することを反映）

変数分離して積分：

$$\int_{\psi_0}^{\psi} \frac{d\psi'}{2\sinh(\psi'/2)} = -\int_0^{\xi} d\xi'$$

$$\ln\left| \tanh\left( \frac{\psi}{4} \right) \right| - \ln\left| \tanh\left( \frac{\psi_0}{4} \right) \right| = -\xi$$

### 表面電荷密度

Gouy-Chapman 理論による表面電荷密度：

$$\sigma = \sqrt{8\varepsilon\varepsilon_0 k_B T c_0 N_A} \sinh\left(\frac{e\phi_0}{2k_B T}\right)$$

- 低電位（eφ << kT）: σ ∝ φ（線形）
- 高電位（eφ >> kT）: σ ∝ exp(eφ/2kT)（指数関数的増加）

---

## Bikermanモデル（Poisson-Fermi方程式）

標準 PB 方程式は、高電位領域でイオン濃度が非物理的に高くなる問題がある。Bikerman モデルは有限イオンサイズを考慮することでこの問題を解決する。

### Bikerman モデルと Poisson-Fermi 方程式の関係

**Bikerman モデル**（1942）と **Poisson-Fermi 方程式**は本質的に同一のものを指す。

| 名称 | 由来・強調点 |
|------|-------------|
| Bikerman モデル | 提唱者 J.J. Bikerman（1942）に由来する歴史的名称 |
| Poisson-Fermi 方程式 | 統計力学的解釈を強調した現代的名称 |

**「Fermi」の由来**：
- Fermi-Dirac 統計では、フェルミ粒子は同じ量子状態を占有できない（Pauli の排他原理）
- 同様に、有限サイズのイオンは同じ空間を占有できない（体積排除効果）
- この類似性から、Bikerman の修正 Boltzmann 分布を「Fermi 型」と呼ぶ

**数学的同値性**：
両者とも同一の修正 Boltzmann 分布を使用：

$$c_\pm = \frac{c_0 \exp(\mp \psi)}{1 - \nu + \nu \cosh(\psi)}$$

文献によっては「modified Poisson-Boltzmann (mPB)」「lattice-gas model」とも呼ばれる。

### 修正 Boltzmann 分布

$$c_\pm = \frac{c_0 \exp(\mp \psi)}{g(\psi)}$$

### 混雑関数（crowding function）

$$g(\psi) = 1 - \nu + \nu \cosh(\psi)$$

### 充填率（packing fraction）

$$\nu = 2 a^3 c_0 NA$$

ここで *a* はイオン直径（イオン液体では典型的に 0.5〜1.0 nm）。

### 修正 Poisson-Boltzmann 方程式

$$\frac{d^2 \psi}{d \xi^2} = \frac{\sinh(\psi)}{g(\psi)}$$

### 物理的解釈

- 低電位（|ψ| ≪ 1）: g(ψ) ≈ 1 となり標準 PB に帰着
- 高電位: g(ψ) により最大濃度が c₀/ν 程度に制限
- 非物理的な濃度 c > 1/(a³NA) を防止

**参考文献**: Kilic, Bazant & Ajdari, *Phys. Rev. E* 75, 021502 (2007)

---

## 数値解法

### 離散化手法

**2次精度中心差分法**を用いて空間離散化を行う。

#### 非一様グリッド

界面付近（x = 0）にグリッド点を集中配置：

$$x_i = L \left[ 1 - (1 - \xi_i)^\beta \right]$$

ここで ξ_i = i/(N-1) ∈ [0, 1]、β > 1 はストレッチング係数。

#### 2階微分の離散化

非一様グリッド上での2階微分：

$$\frac{d^2 \phi}{d x^2} \bigg|_i = \frac{2}{\Delta x^- + \Delta x^+} \left[ \frac{\phi_{i+1} - \phi_i}{\Delta x^+} - \frac{\phi_i - \phi_{i-1}}{\Delta x^-} \right] + O(h^2)$$

ここで：
- Δx⁺ = x_{i+1} - x_i（前方グリッド幅）
- Δx⁻ = x_i - x_{i-1}（後方グリッド幅）

### Newton-Raphson 法

非線形 PB 方程式を Newton-Raphson 反復で解く。

**残差**:

$$F(\phi) = \frac{d^2 \phi}{d x^2} - \kappa^2 \phi_T \sinh\left( \frac{\phi}{\phi_T} \right) = 0$$

**Jacobian**:

$$J = \frac{dF}{d\phi} = \frac{d^2}{d x^2} - \kappa^2 \cosh\left( \frac{\phi}{\phi_T} \right)$$

**Newton 更新**:

$$J \cdot \delta\phi = -F$$
$$\phi^{n+1} = \phi^n + \alpha \cdot \delta\phi$$

ここで α ∈ (0, 1] は適応的ダンピング係数。

### 過渡解析：Scharfetter-Gummel スキーム

過渡解析では**陰的 Scharfetter-Gummel スキーム**を使用。

#### フラックスの離散化

$$J_{i+1/2} = \frac{D}{\Delta x_{i+1/2}} \left[ B(\eta_{i+1/2}) \, n_{i+1} - B(-\eta_{i+1/2}) \, n_i \right]$$

ここで：
- η_{i+1/2} = -ze(φ_{i+1} - φ_i)/(k_BT)（無次元電位差）
- B(η) = η/(e^η - 1) は **Bernoulli 関数**

#### Bernoulli 関数の性質

| 性質 | 数式 |
|------|------|
| 正値性 | B(η) > 0 (∀η) |
| 連続性 | B(0) = 1 |
| 相反関係 | B(-η) = B(η) + η |
| Taylor 展開 | B(η) ≈ 1 - η/2 + η²/12 |

---

## 境界条件

### 開放系（デフォルト）

| 境界 | 電位 φ | 濃度 c |
|------|--------|--------|
| 左端 (x = 0) | Dirichlet: φ = φ₀ | ゼロフラックス: J = 0 |
| 右端 (x = L) | Dirichlet: φ = 0 | Dirichlet: c = c₀ |

### 閉鎖系（キャパシタ）

| 境界 | 電位 φ | 濃度 c |
|------|--------|--------|
| 左端 (x = 0) | Dirichlet: φ = φ_L | ゼロフラックス: J = 0 |
| 右端 (x = L) | Dirichlet: φ = φ_R | ゼロフラックス: J = 0 |

閉鎖系では**電荷中性条件**によりバルク電位が自己無撞着に決定される：

$$\int_0^L \rho \, dx = 0$$

対称な1:1電解質の場合：

$$\phi_{\text{bulk}} = \frac{\phi_L + \phi_R}{2}$$

---

## 参考文献

1. Gouy, L. G. (1910). *J. Phys. Theor. Appl.* 9, 457-468.
2. Chapman, D. L. (1913). *Philos. Mag.* 25, 475-481.
3. Bikerman, J. J. (1942). *Philos. Mag.* 33, 384-397.
4. Scharfetter, D. L. & Gummel, H. K. (1969). *IEEE Trans. Electron Devices* 16, 64-77.
5. Bazant, M. Z. et al. (2009). *Adv. Colloid Interface Sci.* 152, 48-88.
