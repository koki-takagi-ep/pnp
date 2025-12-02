# 完全陰的Newton法による過渡PNPソルバー設計書

## 1. 支配方程式

### Poisson方程式
```
∇²φ = -ρ/ε = -(e·NA/ε)(z₊c₊ + z₋c₋)
```

### Nernst-Planck方程式（カチオン）
```
∂c₊/∂t = ∇·J₊ = ∇·(D₊∇c₊ + (D₊z₊e/kT)c₊∇φ)
```

### Nernst-Planck方程式（アニオン）
```
∂c₋/∂t = ∇·J₋ = ∇·(D₋∇c₋ + (D₋z₋e/kT)c₋∇φ)
```

## 2. 離散化スキーム

### 2.1 時間離散化（後退Euler法）

時刻 n+1 での値を陰的に評価：

```
(c₊^(n+1) - c₊^n) / Δt = ∇·J₊^(n+1)
(c₋^(n+1) - c₋^n) / Δt = ∇·J₋^(n+1)
```

### 2.2 空間離散化（Scharfetter-Gummelスキーム）

セル面 i+1/2 でのフラックス：

```
J₊,i+1/2 = (D₊/Δx) · [B(v₊)·c₊,i+1 - B(-v₊)·c₊,i]
```

ここで：
- `v₊ = z₊·e·(φ_{i+1} - φᵢ) / (kT)` （正規化電場）
- `B(x) = x / (exp(x) - 1)` （Bernoulli関数）

Bernoulli関数の性質：
- `B(0) = 1`
- `B(x) + B(-x) = x`
- `x → ∞`: `B(x) → 0`, `B(-x) → -x`

### 2.3 離散化された方程式系

各格子点 i で3つの方程式：

**Poisson (F_φ,i = 0)**:
```
F_φ,i = (φ_{i+1} - 2φᵢ + φ_{i-1})/Δx² + (e·NA/ε)(z₊c₊,i + z₋c₋,i) = 0
```

**NP+ (F_+,i = 0)**:
```
F_+,i = (c₊,i - c₊,i^n)/Δt - (J₊,i+1/2 - J₊,i-1/2)/Δx = 0
```

**NP- (F_-,i = 0)**:
```
F_-,i = (c₋,i - c₋,i^n)/Δt - (J₋,i+1/2 - J₋,i-1/2)/Δx = 0
```

## 3. Newton法

### 3.1 残差ベクトルとヤコビアン

未知数ベクトル: `U = [φ₀, c₊,₀, c₋,₀, φ₁, c₊,₁, c₋,₁, ..., φ_{N-1}, c₊,N-1, c₋,N-1]`

残差ベクトル: `F(U) = [F_φ,0, F_+,0, F_-,0, ..., F_φ,N-1, F_+,N-1, F_-,N-1]`

Newton更新: `J·δU = -F(U)`, `U^(k+1) = U^(k) + δU`

### 3.2 ヤコビアン行列の構造

3N × 3N のブロック三重対角行列：

```
J = | A₀  B₀   0   0  ... |
    | C₁  A₁  B₁   0  ... |
    |  0  C₂  A₂  B₂  ... |
    | ...                  |
```

各ブロック Aᵢ, Bᵢ, Cᵢ は 3×3 行列。

### 3.3 ヤコビアン要素の導出

**Poisson方程式の微分**:
```
∂F_φ,i/∂φᵢ = -2/Δx²
∂F_φ,i/∂φ_{i±1} = 1/Δx²
∂F_φ,i/∂c₊,i = e·NA·z₊/ε
∂F_φ,i/∂c₋,i = e·NA·z₋/ε
```

**NP方程式の微分**:

フラックスの微分（Scharfetter-Gummel）:
```
∂J₊,i+1/2/∂c₊,i = -(D₊/Δx)·B(-v₊)
∂J₊,i+1/2/∂c₊,i+1 = (D₊/Δx)·B(v₊)
∂J₊,i+1/2/∂φᵢ = (D₊/Δx)·[B'(v₊)·c₊,i+1 + B'(-v₊)·c₊,i]·(z₊e/kT)/Δx
∂J₊,i+1/2/∂φ_{i+1} = -(D₊/Δx)·[B'(v₊)·c₊,i+1 + B'(-v₊)·c₊,i]·(z₊e/kT)/Δx
```

ここで `B'(x) = dB/dx = (B(x) - 1 + B(x)·exp(x)) / x`

## 4. 境界条件

### 左境界 (x = 0): Dirichlet電位、Zero-flux濃度
```
φ₀ = φ_left  (固定電位)
J₊,1/2 = 0   (ブロッキング電極)
J₋,1/2 = 0
```

### 右境界 (x = L): Dirichletすべて
```
φ_{N-1} = φ_right = 0  (接地)
c₊,N-1 = c₀           (バルク濃度)
c₋,N-1 = c₀
```

## 5. アルゴリズム

```
初期化:
  c₊,i = c₀, c₋,i = c₀ for all i
  φᵢ = φ_left·(1 - xᵢ/L)  (線形プロファイル)

時間ループ: for n = 0, 1, 2, ..., N_time

  Newton反復: for k = 0, 1, 2, ...

    1. 残差ベクトル F(U^k) を計算

    2. ヤコビアン行列 J(U^k) を構築

    3. 線形系を解く: J·δU = -F
       (ブロック三重対角なのでThomas法の拡張を使用)

    4. 更新: U^(k+1) = U^k + α·δU
       (α: ダンピング係数、必要に応じてline searchで調整)

    5. 収束判定: ||δU|| < tol または ||F|| < tol
       収束したら break

  end Newton

  U^(n+1) = U^k
  スナップショット保存（必要に応じて）

end 時間ループ
```

## 6. 実装上の注意点

### 6.1 Bernoulli関数の数値的取り扱い

`x ≈ 0` のとき Taylor展開を使用：
```
B(x) ≈ 1 - x/2 + x²/12 - x⁴/720 + ...  (|x| < 0.01)
```

### 6.2 濃度の正値性保証

Newton更新後に濃度が負になる場合：
- ダンピング係数 α を減少
- または下限 `c_min = 10⁻¹² × c₀` でクランプ

### 6.3 収束加速

- 良い初期推定（前時刻の解を使用）
- 適応的ダンピング（α を反復ごとに調整）
- 悪条件時はGummel反復にフォールバック

## 7. 参考文献

1. Scharfetter & Gummel (1969) - IEEE Trans. Electron Devices
2. Bank et al. - SIAM J. Numer. Anal. (Finite Element for drift-diffusion)
3. Liu & Wang (2020) - Numer. Math. (Positivity-preserving schemes)
4. Bousquet et al. (2018) - SIAM J. Sci. Comput. (Newton solvers)
