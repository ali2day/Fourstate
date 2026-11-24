# Spin Hamiltonian and conventions

## 1. Hamiltonian (arXiv:2407.20489, Eqs. 1–2) — used by `paper_convention.py`

$$
H=\sum_i \mathbf S_i\cdot\mathbf A_{ii}\cdot\mathbf S_i
+\sum_{i<j}\mathbf S_i\cdot\mathcal J_{ij}\cdot\mathbf S_j
+\sum_{i<j}K_{ij}\,(\mathbf S_i\cdot\mathbf S_j)^2
$$

$$
\mathbf S_i\cdot\mathcal J_{ij}\cdot\mathbf S_j
=J_{ij}\,\mathbf S_i\cdot\mathbf S_j+\mathbf D_{ij}\cdot(\mathbf S_i\times\mathbf S_j)+\mathbf S_i\cdot\Gamma_{ij}\cdot\mathbf S_j
$$

| quantity | definition | sign meaning |
|---|---|---|
| spins | **unit vectors**, every pair counted **once** (i < j) | – |
| J | Tr 𝒥 / 3 | J < 0 ferromagnetic, J > 0 antiferromagnetic |
| D | D_x = (𝒥_yz − 𝒥_zy)/2, D_y = (𝒥_zx − 𝒥_xz)/2, D_z = (𝒥_xy − 𝒥_yx)/2 | zero with an inversion centre on the bond |
| Γ | (𝒥 + 𝒥ᵀ)/2 − J·1 (symmetric, traceless) | anisotropic (Kitaev-like) exchange |
| K | coefficient of +(S_i·S_j)² | K < 0 favours collinear spins |
| A | single-ion tensor, +S·A·S | A^zz > 0 easy plane, A^zz < 0 easy axis |

Units: meV. Frame: x, y, z of the VASP cell (the spin quantisation axes in M_CONSTR, SAXIS = 0 0 1).

### Relation to the other convention in this code
`analyze.py`, `analyze_sia.py`, `analyze_bq.py` use the TB2J sign, H = −Σ S_i·J·S_j − K (S_i·S_j)²:

| | paper (`paper_convention.py`) | TB2J (`analyze*.py`) |
|---|---|---|
| exchange tensor | 𝒥 | −J |
| biquadratic | K | −K |
| single-ion A | A | A (same) |
| ferromagnet | J < 0 | J > 0 |

## 2. Four-state formulas (paper sign)

Bilinear term, component (a, b): S_i = ±a, S_j = ±b, all other spins along c ⊥ a, b
(e.g. jxx: c = z; dx = (y, z), c = x):

$$
\mathcal J_{ab}=\tfrac14\left[E(+a,+b)+E(-a,-b)-E(+a,-b)-E(-a,+b)\right]
$$

Everything linear in S_i or S_j (coupling to the other spins, Zeeman-like terms), the single-ion
terms and K cancel, because K depends on (S_i·S_j)² = 0 for a ≠ b and = 1 in all four states for a = b.

Script labels → components (1, 2, 3, 4 = ++, +−, −+, −−):

| label | S_i, S_j | others | gives |
|---|---|---|---|
| jxx, jyy, jzz | (x,x), (y,y), (z,z) | z, z, x | 𝒥_xx, 𝒥_yy, 𝒥_zz |
| dx / dxb | (y,z) / (z,y) | x | 𝒥_yz / 𝒥_zy |
| dy / dyb | (x,z) / (z,x) | y | 𝒥_xz / 𝒥_zx |
| dz / dzb | (x,y) / (y,x) | z | 𝒥_xy / 𝒥_yx |

**Biquadratic** (k1–k4): S_i = x, others along z, S_j at 0°, 180°, 45°, 225° in the xy plane, so
c = S_i·S_j = 1, −1, 1/√2, −1/√2. With E = E₀ + 𝒥_xx-type bilinear terms + K c² + single-ion of S_j:

$$
E_1+E_2-(E_3+E_4)=2K-K+2(a_x-a_{45})
\;\Rightarrow\;
K=E_1+E_2-E_3-E_4-2(a_x-a_{45})
$$

with a_x = A_xx and a_45 = (A_xx + A_yy)/2 + A_xy (single-ion energy of S_j along x and at 45°).
For a site with only A^zz (pristine trigonal), a_x = a_45 and the correction is zero.
Check from the same runs: E₁ − E₂ = 2𝒥_xx.

**Single-ion** (only S_i rotated, others ⊥): a1–4 = +z, −z, +x, −x (others y);
ayyxx1–4 = +y, −y, +x, −x (others z); axy/axz/ayz1–4 = (±a ±b)/√2.

$$
A^{zz}-A^{xx}=\tfrac12[E(+z)+E(-z)-E(+x)-E(-x)],\qquad
A^{ab}=\tfrac14[E(++)+E(--)-E(+-)-E(-+)]
$$

## 3. Extraction

```bash
python fourstate/extract.py all_sub_9_6  > eng_1NN.dat      # one line per run: energy(sigma->0), E_p, P, mag
python fourstate/extract.py all_sub_9_12 > eng_bq.dat
python fourstate/paper_convention.py --table eng_1NN.dat eng_bq.dat
```
* Energy = `energy(sigma->0)` of the converged ionic step (`energy  without entropy=` line).
* `--minus-penalty` subtracts the constraint energy E_p (OSZICAR); `--first` uses the first ionic step
  when the ions relax in every state.
* Check that every run converged (`conv` column) and that the total moments (`Mag_x..z`) are what the
  state requires: e.g. k1 (4, 0, 32), k2 (0, 0, 32), k3 (3.41, 1.41, 32), k4 (0.59, −1.41, 32) for
  2 μB spins with 16 background spins along z. With gen.py folders use `--runs`, which checks every atom.
* Two states that should be equal by symmetry (e.g. jxx1 and jxx4) are a good test of convergence.

## 4. Worked example: NiI2 monolayer, nearest neighbours (6×3 cell, Ni 9 – Ni 6)

```
𝒥 (meV)    x        y        z
  x      -3.916    0.000    0.000
  y       0.000   -2.977   -0.588
  z       0.000   -0.588   -3.585
```
J = −3.49 meV (FM), D = 0, Γ^xx = −0.42, Γ^yy = +0.52, Γ^zz = −0.09, Γ^yz = −0.59 meV, Γ^yz/J = 0.17.
Biquadratic (Ni 9 – Ni 12, k1–k4): K = E₁ + E₂ − E₃ − E₄ = −0.70 meV, K/J = 0.20; 𝒥_xx from k1/k2 = −3.33 meV.
These reproduce Table 1 of arXiv:2407.20489 (J −3.49, Γ −0.42/0.52/−0.09/−0.59, K −0.67 meV).
