# Spin Hamiltonian and four-state formulas

## 1. Spin Hamiltonian

We use the Hamiltonian of Ghojavand *et al.*, arXiv:2407.20489 (Eqs. 1–2):

```math
H \;=\; \sum_{i} \mathbf{S}_i \cdot \mathbf{A}_{ii} \cdot \mathbf{S}_i
\;+\; \sum_{i<j} \mathbf{S}_i \cdot \mathcal{J}_{ij} \cdot \mathbf{S}_j
\;+\; \sum_{i<j} K_{ij} \left( \mathbf{S}_i \cdot \mathbf{S}_j \right)^2 ,
```

where the exchange tensor $`\mathcal{J}_{ij}`$ is decomposed as

```math
\mathbf{S}_i \cdot \mathcal{J}_{ij} \cdot \mathbf{S}_j
= J_{ij}\, \mathbf{S}_i \cdot \mathbf{S}_j
+ \mathbf{D}_{ij} \cdot \left( \mathbf{S}_i \times \mathbf{S}_j \right)
+ \mathbf{S}_i \cdot \Gamma_{ij} \cdot \mathbf{S}_j ,
```

with

```math
J_{ij} = \tfrac{1}{3}\,\mathrm{Tr}\,\mathcal{J}_{ij}, \qquad
\Gamma_{ij} = \tfrac{1}{2}\left( \mathcal{J}_{ij} + \mathcal{J}_{ij}^{T} \right) - J_{ij}\,\mathbb{1}, \qquad
D^{x}_{ij} = \tfrac{1}{2}\left( \mathcal{J}^{yz}_{ij} - \mathcal{J}^{zy}_{ij} \right),\;
D^{y}_{ij} = \tfrac{1}{2}\left( \mathcal{J}^{zx}_{ij} - \mathcal{J}^{xz}_{ij} \right),\;
D^{z}_{ij} = \tfrac{1}{2}\left( \mathcal{J}^{xy}_{ij} - \mathcal{J}^{yx}_{ij} \right).
```

**Conventions**

* $`\mathbf{S}_i`$ is a **unit vector**; every pair is counted **once** ($`i<j`$); energies in meV.
* $`J_{ij} < 0`$: ferromagnetic, $`J_{ij} > 0`$: antiferromagnetic.
* $`K_{ij} < 0`$ favours collinear spins.
* $`A^{zz} > 0`$: easy plane, $`A^{zz} < 0`$: easy axis.
* $`x, y, z`$ are the Cartesian axes of the VASP cell (the directions written in `MAGMOM` / `M_CONSTR`, `SAXIS = 0 0 1`).

**Relation to the TB2J sign** (used by `analyze.py`, `analyze_sia.py`, `analyze_bq.py`,
$`H = -\sum_{i<j} \mathbf{S}_i \cdot \mathbf{J}_{ij} \cdot \mathbf{S}_j - \sum_{i<j} K^{\mathrm{TB2J}}_{ij} (\mathbf{S}_i \cdot \mathbf{S}_j)^2`$):

```math
\mathcal{J}_{ij} = -\,\mathbf{J}^{\mathrm{TB2J}}_{ij}, \qquad K_{ij} = -\,K^{\mathrm{TB2J}}_{ij}, \qquad \mathbf{A}_{ii} \text{ unchanged.}
```

`paper_convention.py` reports everything in the convention of this section.

---

## 2. Exchange tensor $`\mathcal{J}_{ij}`$

For one component $`\mathcal{J}^{\alpha\beta}_{ij}`$ ($`\alpha, \beta = x, y, z`$) four DFT calculations are made.
Spin $`i`$ points along $`\pm\hat{\boldsymbol{\alpha}}`$, spin $`j`$ along $`\pm\hat{\boldsymbol{\beta}}`$, and **all other
spins** point along the third axis $`\hat{\boldsymbol{\gamma}}`$, perpendicular to both:

| state | $`\mathbf{S}_i`$ | $`\mathbf{S}_j`$ | other spins |
|:---:|:---:|:---:|:---:|
| 1 | $`+\hat{\boldsymbol{\alpha}}`$ | $`+\hat{\boldsymbol{\beta}}`$ | $`\hat{\boldsymbol{\gamma}}`$ |
| 2 | $`+\hat{\boldsymbol{\alpha}}`$ | $`-\hat{\boldsymbol{\beta}}`$ | $`\hat{\boldsymbol{\gamma}}`$ |
| 3 | $`-\hat{\boldsymbol{\alpha}}`$ | $`+\hat{\boldsymbol{\beta}}`$ | $`\hat{\boldsymbol{\gamma}}`$ |
| 4 | $`-\hat{\boldsymbol{\alpha}}`$ | $`-\hat{\boldsymbol{\beta}}`$ | $`\hat{\boldsymbol{\gamma}}`$ |

```math
\mathcal{J}^{\alpha\beta}_{ij} = \frac{E_1 + E_4 - E_2 - E_3}{4}
```

All terms that are linear in $`\mathbf{S}_i`$ or in $`\mathbf{S}_j`$ (coupling to the other spins), the single-ion
terms (equal for $`\pm\mathbf{S}`$) and the biquadratic term ($`(\mathbf{S}_i\cdot\mathbf{S}_j)^2`$ is the same in the
four states) cancel in this combination.

The nine components, with the folder names of the shell scripts (`1, 2, 3, 4` = states above):

| folders | $`\alpha\beta`$ | other spins | component |
|:---:|:---:|:---:|:---:|
| `jxx1–4` | $`xx`$ | $`\hat{z}`$ | $`\mathcal{J}^{xx}`$ |
| `jyy1–4` | $`yy`$ | $`\hat{z}`$ | $`\mathcal{J}^{yy}`$ |
| `jzz1–4` | $`zz`$ | $`\hat{x}`$ | $`\mathcal{J}^{zz}`$ |
| `dx1–4` / `dxb1–4` | $`yz`$ / $`zy`$ | $`\hat{x}`$ | $`\mathcal{J}^{yz}`$ / $`\mathcal{J}^{zy}`$ |
| `dy1–4` / `dyb1–4` | $`xz`$ / $`zx`$ | $`\hat{y}`$ | $`\mathcal{J}^{xz}`$ / $`\mathcal{J}^{zx}`$ |
| `dz1–4` / `dzb1–4` | $`xy`$ / $`yx`$ | $`\hat{z}`$ | $`\mathcal{J}^{xy}`$ / $`\mathcal{J}^{yx}`$ |

---

## 3. Biquadratic exchange $`K_{ij}`$

$`\mathbf{S}_i`$ is fixed along $`\hat{x}`$, all other spins along $`\hat{z}`$, and $`\mathbf{S}_j`$ is rotated in the $`xy`$ plane
by the angle $`\theta`$ to $`\mathbf{S}_i`$:

| state | $`\mathbf{S}_i`$ | $`\mathbf{S}_j`$ | $`\theta`$ | $`\mathbf{S}_i\cdot\mathbf{S}_j`$ | $`(\mathbf{S}_i\cdot\mathbf{S}_j)^2`$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| `k1` | $`\hat{x}`$ | $`\hat{x}`$ | $`0^\circ`$ | $`1`$ | $`1`$ |
| `k2` | $`\hat{x}`$ | $`-\hat{x}`$ | $`180^\circ`$ | $`-1`$ | $`1`$ |
| `k3` | $`\hat{x}`$ | $`\tfrac{1}{\sqrt2}(\hat{x}+\hat{y})`$ | $`45^\circ`$ | $`\tfrac{1}{\sqrt2}`$ | $`\tfrac12`$ |
| `k4` | $`\hat{x}`$ | $`-\tfrac{1}{\sqrt2}(\hat{x}+\hat{y})`$ | $`225^\circ`$ | $`-\tfrac{1}{\sqrt2}`$ | $`\tfrac12`$ |

The energy of each state is

```math
E(\mathbf{S}_j) = E_0 + \mathbf{S}_i\cdot\mathcal{J}_{ij}\cdot\mathbf{S}_j + K_{ij}\,(\mathbf{S}_i\cdot\mathbf{S}_j)^2
+ \mathbf{S}_j\cdot\mathbf{A}_{jj}\cdot\mathbf{S}_j + \mathbf{S}_j\cdot\mathbf{h} ,
```

where $`\mathbf{h}`$ collects the coupling of $`\mathbf{S}_j`$ to all other (fixed) spins. In $`E_{k1}+E_{k2}`$ and in
$`E_{k3}+E_{k4}`$ the two spin directions of $`j`$ are opposite, so the bilinear term and $`\mathbf{S}_j\cdot\mathbf{h}`$
cancel and

```math
\left(E_{k1}+E_{k2}\right) - \left(E_{k3}+E_{k4}\right) = 2K_{ij} - K_{ij} + 2\left(a_{0} - a_{45}\right),
```

```math
K_{ij} = E_{k1} + E_{k2} - E_{k3} - E_{k4} - 2\left(a_{0} - a_{45}\right)
```

with the single-ion energies of $`\mathbf{S}_j`$ along $`\hat{x}`$ and at $`45^\circ`$:

```math
a_{0} = A^{xx}_{jj}, \qquad a_{45} = \tfrac12\left( A^{xx}_{jj} + A^{yy}_{jj} \right) + A^{xy}_{jj} .
```

For a site with only $`A^{zz}`$ (trigonal symmetry, e.g. pristine NiI$`_2`$) $`a_0 = a_{45}`$ and the correction vanishes.
The same four runs also give $`E_{k1} - E_{k2} = 2\,\mathcal{J}^{xx}_{ij}`$ as a check.

---

## 4. Single-ion anisotropy $`\mathbf{A}_{ii}`$

Only $`\mathbf{S}_i`$ is rotated; all other spins point along a fixed axis perpendicular to the rotation plane.

| folders | $`\mathbf{S}_i`$ (states 1–4) | other spins | gives |
|:---:|:---:|:---:|:---:|
| `a1–4` | $`+\hat{z},\,-\hat{z},\,+\hat{x},\,-\hat{x}`$ | $`\hat{y}`$ | $`A^{zz}-A^{xx}`$ |
| `ayyxx1–4` | $`+\hat{y},\,-\hat{y},\,+\hat{x},\,-\hat{x}`$ | $`\hat{z}`$ | $`A^{yy}-A^{xx}`$ |
| `axy1–4` | $`\tfrac{1}{\sqrt2}(\pm\hat{x}\pm\hat{y})`$ | $`\hat{z}`$ | $`A^{xy}`$ |
| `axz1–4` | $`\tfrac{1}{\sqrt2}(\pm\hat{x}\pm\hat{z})`$ | $`\hat{y}`$ | $`A^{xz}`$ |
| `ayz1–4` | $`\tfrac{1}{\sqrt2}(\pm\hat{y}\pm\hat{z})`$ | $`\hat{x}`$ | $`A^{yz}`$ |

```math
A^{zz}-A^{xx} = \frac{E_1 + E_2 - E_3 - E_4}{2}, \qquad
A^{\alpha\beta} = \frac{E_{++} + E_{--} - E_{+-} - E_{-+}}{4}\quad(\alpha\neq\beta).
```

$`A^{zz}`$ is reported with $`A^{xx} = A^{yy} = 0`$, as in Table 1 of the paper.

---

## 5. From run folders to parameters

```bash
python fourstate/extract.py all_sub_9_6  > eng_1NN.dat     # exchange runs (jxx1 ... dzb4)
python fourstate/extract.py all_sub_9_12 > eng_bq.dat      # biquadratic runs (k1 ... k4)
python fourstate/paper_convention.py --table eng_1NN.dat eng_bq.dat
```

* $`E`$ = `energy(sigma->0)` of the converged ionic step (the `energy  without entropy=` line of the OUTCAR).
* `--minus-penalty` subtracts the constraint energy $`E_p`$ (OSZICAR); `--first` takes the first ionic step
  when the ions relax in every state.
* Check that every run converged (`conv`) and that the total moment is the one of the state. For
  $`2\,\mu_B`$ spins and 16 other spins along $`\hat{z}`$: `k1` $`(4,0,32)`$, `k2` $`(0,0,32)`$, `k3` $`(3.41,1.41,32)`$,
  `k4` $`(0.59,-1.41,32)\,\mu_B`$.
* States that are equal by symmetry (e.g. `jxx1` and `jxx4`) must have the same energy; a difference
  signals an unconverged run or a different electronic state.

---

## 6. Example: NiI$`_2`$ monolayer, nearest neighbours

6×3 cell, Ni 9 – Ni 6 (exchange) and Ni 9 – Ni 12 (biquadratic), PBE+U, SOC:

```math
\mathcal{J}_{\mathrm{NN}} =
\begin{pmatrix}
-3.916 & 0.000 & 0.000 \\
0.000 & -2.977 & -0.588 \\
0.000 & -0.588 & -3.585
\end{pmatrix}\ \mathrm{meV}
```

| $`J`$ | $`\mathbf{D}`$ | $`\Gamma^{xx}`$ | $`\Gamma^{yy}`$ | $`\Gamma^{zz}`$ | $`\Gamma^{yz}`$ | $`\Gamma^{yz}/J`$ | $`K`$ | $`K/J`$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $`-3.49`$ | $`0`$ | $`-0.42`$ | $`0.52`$ | $`-0.09`$ | $`-0.59`$ | $`0.17`$ | $`-0.70`$ | $`0.20`$ |

(meV). These values reproduce Table 1 of arXiv:2407.20489.
