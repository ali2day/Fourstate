# fourstate

Four-state method with VASP for any magnetic material:
* **energy mapping** → isotropic exchange, the full exchange tensor (DMI, anisotropic exchange),
  the single-ion anisotropy and the biquadratic exchange;
* **polarization mapping** (Berry phase, same four spin states) → the spin-induced polarization
  (magnetoelectric) tensor of the generalized spin-current (GSC) model.

* Xiang et al., PRB 84, 224429 (2011); Dalton Trans. 42, 823 (2013)
* GSC model: Xiang et al., PRL 107, 157202 (2011)

Requires Python 3 with `numpy`, `ase`, `spglib`.

## What it does
1. Choose the pair (i, j):
   * automatically: `gen.py list` finds the symmetry-inequivalent bonds, and `build --bonds` makes a
     supercell where no periodic image of j is close to i;
   * by hand: `build --atoms 9 6` uses atoms 9 and 6 (POSCAR numbering, from 1) of the cell you give,
     as it is (e.g. a supercell with a defect).
2. Write one VASP folder per spin state: MAGMOM and M_CONSTR for every magnetic atom.
3. Read the results, put every run in the slot of its **actual** final spins (the constraint fixes
   only the axis, so spins can flip), reject runs with background spins off target, give the tensors.

| mode | runs | result | analysis |
|---|---|---|---|
| `iso` | 4 / bond (collinear) | J_iso | `analyze.py` |
| `tensor` | 36 / bond (SOC) | 3×3 J, J_iso, DMI, anisotropic exchange | `analyze.py` |
| `gsc` | 36 / bond (SOC + Berry phase) | full ME tensor M, P = S_i · M · S_j | `analyze_gsc.py` |
| `me` | 12 / bond (SOC + Berry phase) | spin-current ME matrix, P = M (S_i × S_j) | `analyze_gsc.py` |
| `sia` | 20 / site (SOC) | single-ion anisotropy A, H = S · A · S | `analyze_sia.py` |
| `bq` | 4 / bond (SOC) | biquadratic K, H = −J S_i·S_j − K (S_i·S_j)² | `analyze_bq.py` |

For component (a, b): S_i = ±a, S_j = ±b, all other spins along the third axis c.
`me`: background +z; k1–4 (S_i, S_j) = (±x, ±y), k5–8 (±x, ±z), k9–12 (±y, ±z);
M[:,z] = (P1+P4−P2−P3)/4, M[:,y] = (P6+P7−P5−P8)/4, M[:,x] = (P9+P12−P10−P11)/4. k5, k7, k9, k11 have
S_j = +z = background, so they are the same for every j of the same i.
`bq`: S_i = x, S_j at 0°, 180°, 45°, 225° in the xy plane, others z:
K = −[(E1 + E2) − (E3 + E4)] + 2(a_x − a_45), where the single-ion term of j comes from a `sia` run
(`analyze_bq.py runs_bq --sia runs_sia/results_sia.json`).
Convention of `analyze*.py`: H = −Σ S_i·J·S_j (unit spins, J > 0 ferromagnetic), as in TB2J.

### Convention of arXiv:2407.20489 (`paper_convention.py`)
H = Σ_i S_i·A·S_i + Σ_{i<j} S_i·𝒥_ij·S_j + Σ_{i<j} K_ij (S_i·S_j)², unit spins, meV, J < 0 ferromagnetic,
J = Tr 𝒥/3, D from (𝒥 − 𝒥ᵀ)/2, Γ = (𝒥 + 𝒥ᵀ)/2 − J·1, K < 0 favours collinear spins, A^zz > 0 easy plane.
Full explanation, derivations and a NiI2 worked example: **[CONVENTIONS.md](CONVENTIONS.md)**.
```bash
python fourstate/extract.py all_sub_9_6 > eng.dat                   # run folders -> one line per run
python fourstate/paper_convention.py --table eng.dat                # labels jxx1.., dx1.., a1.., k1..
python fourstate/paper_convention.py --runs runs_exchange runs_sia runs_bq
```
Formulas (S_i = ±a, S_j = ±b, others ⟂): 𝒥_ab = [E(++) + E(−−) − E(+−) − E(−+)]/4;
K = [E(k1) + E(k2)] − [E(k3) + E(k4)] − 2(a_x − a_45), with a_x = A_xx, a_45 = (A_xx + A_yy)/2 + A_xy.

## Use
```bash
python fourstate/gen.py list  --poscar POSCAR --mag Ni --rmax 8
python fourstate/gen.py build --poscar POSCAR --mag Ni --bonds 0 --mode tensor \
       --incar templates/INCAR.tensor --moment 2 --potcar-dir POTCARS --out runs
python fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 6 --mode gsc \
       --incar templates/INCAR.gsc --moment 2 --out runs_me
bash   fourstate/submit.sh runs --exe vasp_ncl --dry      # then without --dry
python fourstate/analyze.py runs                         # J
python fourstate/analyze_gsc.py runs_me                  # ME tensor
python fourstate/table.py runs                           # one line per run: E, E_p, spins, P
```
* INCAR templates: set ENCUT, +U, RWIGS for your system; keep the `{MAGMOM}` line. Species order in
  the runs = order in your POSCAR.
* Energy = `energy(sigma->0)` of the converged ionic step; the constraint penalty `E_p` (OSZICAR) is
  subtracted (`--ignore-penalty` to keep it).
* `--potcar-dir DIR` expects `DIR/<Element>/POTCAR`. `--moment Cr=3 Mn=4` for several species.
* Ions relaxing in every state (`NSW > 0` in the INCAR or `--relax NSW`): `--step first` gives the
  result at the input geometry, `--step last` after relaxation. `gsc` needs one fixed geometry.
* `LORBIT = 11` is required (the analysis checks every moment).

## Examples
`examples/NiI2_exchange_and_ME`: monolayer NiI2, Ni 9 – Ni 6 chosen by hand → exchange tensor,
ME tensor, single-ion anisotropy and biquadratic exchange (Ni 9 – Ni 12). `examples/exchange_tensor_VBr3`: bulk VBr3, automatic bonds.

## Tests
```bash
python tests/test_energy.py   # known J tensor + spin flip
python tests/test_gsc.py      # known M tensor + polarization quantum
python tests/test_me.py       # 12-run ME scheme, states k1-k12
python tests/test_sia.py      # known A tensor + reversed spin, --atoms
python tests/test_bq.py       # known K + single-ion correction + reversed spin
python tests/test_paper.py    # arXiv:2407.20489 convention from a table
python tests/test_extract.py  # run folders -> table -> K
```

## License
MIT
