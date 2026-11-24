# fourstate

Four-state method with VASP for any magnetic material:

* **energy mapping** → exchange tensor (isotropic $`J`$, Dzyaloshinskii–Moriya $`\mathbf D`$, anisotropic $`\Gamma`$),
  single-ion anisotropy $`\mathbf A`$ and biquadratic exchange $`K`$;
* **polarization mapping** (Berry phase, same spin states) → spin-induced polarization (magnetoelectric)
  tensor of the generalized spin-current (GSC) model.

References: Xiang *et al.*, PRB **84**, 224429 (2011); Dalton Trans. **42**, 823 (2013); PRL **107**, 157202 (2011).
Hamiltonian and sign convention: Ghojavand *et al.*, arXiv:2407.20489.

Requires Python 3 with `numpy`, `ase`, `spglib`.

## Spin Hamiltonian

```math
H=\sum_i \mathbf S_i\cdot\mathbf A_{ii}\cdot\mathbf S_i
+\sum_{i<j}\mathbf S_i\cdot\mathcal J_{ij}\cdot\mathbf S_j
+\sum_{i<j}K_{ij}\,(\mathbf S_i\cdot\mathbf S_j)^2
```

```math
\mathbf S_i\cdot\mathcal J_{ij}\cdot\mathbf S_j
=J_{ij}\,\mathbf S_i\cdot\mathbf S_j+\mathbf D_{ij}\cdot(\mathbf S_i\times\mathbf S_j)+\mathbf S_i\cdot\Gamma_{ij}\cdot\mathbf S_j
```

Unit spins, each pair counted once, meV. $`J<0`$ ferromagnetic, $`K<0`$ favours collinear spins,
$`A^{zz}>0`$ easy plane. Derivations of every formula: **[CONVENTIONS.md](CONVENTIONS.md)**.

> `analyze.py`, `analyze_sia.py`, `analyze_bq.py` print the TB2J sign
> ($`H=-\sum_{i<j}\mathbf S_i\cdot\mathbf J_{ij}\cdot\mathbf S_j`$, $`J>0`$ ferromagnetic);
> `paper_convention.py` prints the convention above ($`\mathcal J=-\mathbf J^{\mathrm{TB2J}}`$, $`K=-K^{\mathrm{TB2J}}`$).

## Four-state formulas

Exchange component $`\alpha\beta`$: $`\mathbf S_i=\pm\hat{\boldsymbol\alpha}`$, $`\mathbf S_j=\pm\hat{\boldsymbol\beta}`$,
all other spins along $`\hat{\boldsymbol\gamma}\perp\hat{\boldsymbol\alpha},\hat{\boldsymbol\beta}`$:

```math
\mathcal J^{\alpha\beta}_{ij}=\frac{E_{++}+E_{--}-E_{+-}-E_{-+}}{4}
```

Biquadratic: $`\mathbf S_i=\hat x`$, other spins $`\hat z`$, $`\mathbf S_j`$ at $`0^\circ,180^\circ,45^\circ,225^\circ`$ in the $`xy`$ plane (k1–k4):

```math
K_{ij}=E_{k1}+E_{k2}-E_{k3}-E_{k4}-2\,(a_0-a_{45})
```

with $`a_0-a_{45}`$ the single-ion energy difference of $`\mathbf S_j`$ (zero for a trigonal site).

Magnetoelectric matrix (mode `me`, other spins along $`+\hat z`$; k1–4: $`(\pm\hat x,\pm\hat y)`$, k5–8: $`(\pm\hat x,\pm\hat z)`$,
k9–12: $`(\pm\hat y,\pm\hat z)`$), with $`\mathbf P=\mathbf M\,(\mathbf S_i\times\mathbf S_j)`$:

```math
\mathbf M_{\cdot z}=\tfrac14(\mathbf P_1+\mathbf P_4-\mathbf P_2-\mathbf P_3),\qquad
\mathbf M_{\cdot y}=\tfrac14(\mathbf P_6+\mathbf P_7-\mathbf P_5-\mathbf P_8),\qquad
\mathbf M_{\cdot x}=\tfrac14(\mathbf P_9+\mathbf P_{12}-\mathbf P_{10}-\mathbf P_{11})
```

## Modes

| mode | runs | result | analysis |
|---|---|---|---|
| `iso` | 4 / bond, collinear | $`J`$ | `analyze.py` |
| `tensor` | 36 / bond, SOC | $`\mathcal J_{ij}`$ (3×3), $`J`$, $`\mathbf D`$, $`\Gamma`$ | `analyze.py` |
| `gsc` | 36 / bond, SOC + Berry phase | full ME tensor, $`\mathbf P=\mathbf S_i\cdot\mathbf M\cdot\mathbf S_j`$ | `analyze_gsc.py` |
| `me` | 12 / bond, SOC + Berry phase | spin-current ME matrix, $`\mathbf P=\mathbf M\,(\mathbf S_i\times\mathbf S_j)`$ | `analyze_gsc.py` |
| `sia` | 20 / site, SOC | $`\mathbf A`$ | `analyze_sia.py` |
| `bq` | 4 / bond, SOC | $`K`$ | `analyze_bq.py` |

The pair $`(i,j)`$ is chosen either automatically (`gen.py list`, then `build --bonds`: symmetry-inequivalent
bonds and a supercell with no periodic image of $`j`$ close to $`i`$) or by hand (`build --atoms 9 6`: atoms of
the cell you give, POSCAR numbering from 1). Every run is analysed with its **actual** final spins (the
constraint fixes only the axis, so spins can flip); runs with other spins off target are rejected.

## Use

```bash
# generate
python fourstate/gen.py list  --poscar POSCAR --mag Ni --rmax 8
python fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 6 --mode tensor \
       --incar templates/INCAR.tensor --moment 2 --potcar-dir POTCARS --out runs_exchange
bash   fourstate/submit.sh runs_exchange --exe vasp_ncl --dry      # then without --dry

# analyse
python fourstate/paper_convention.py --runs runs_exchange runs_sia runs_bq
python fourstate/analyze_gsc.py runs_me

# runs made with other scripts (folders jxx1 ... dzb4, a1 ..., k1 ... k4)
python fourstate/extract.py all_sub_9_6 > eng.dat
python fourstate/paper_convention.py --table eng.dat
```

* INCAR templates: set ENCUT, $`+U`$, RWIGS for your system; keep the `{MAGMOM}` line. Species order = POSCAR order.
* Energy = `energy(sigma->0)` of the converged ionic step; the constraint energy $`E_p`$ (OSZICAR) is
  subtracted by `analyze*.py` (`--ignore-penalty` to keep it).
* `--potcar-dir DIR` expects `DIR/<Element>/POTCAR`; `--moment Cr=3 Mn=4` for several species.
* Ions relaxing in every state (`NSW > 0`): `--step first` gives the input geometry, `--step last` the relaxed one.
  Polarization modes need one fixed geometry.
* `LORBIT = 11` is required.

## Example

`examples/NiI2_exchange_and_ME`: monolayer NiI$`_2`$ (6×3 cell), Ni 9 – Ni 6 → exchange tensor and ME matrix,
Ni 9 → single-ion anisotropy, Ni 9 – Ni 12 → biquadratic $`K`$.

## Tests

```bash
python tests/test_energy.py    # known J tensor + spin flip
python tests/test_gsc.py       # known M tensor + polarization quantum
python tests/test_me.py        # 12-run ME scheme, states k1-k12
python tests/test_sia.py       # known A tensor + reversed spin
python tests/test_bq.py        # known K + single-ion correction
python tests/test_paper.py     # arXiv:2407.20489 convention from a table
python tests/test_extract.py   # run folders -> table -> K
```

## License

MIT
