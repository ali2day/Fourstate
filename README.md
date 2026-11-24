# fourstate

Four-state method with VASP for any magnetic material:
* **energy mapping** → isotropic exchange and the full exchange tensor (DMI, anisotropic exchange);
* **polarization mapping** (Berry phase, same four spin states) → the spin-induced polarization
  tensor of the generalized spin-current (GSC) model.

* Xiang et al., PRB 84, 224429 (2011); Dalton Trans. 42, 823 (2013)
* GSC model: Xiang et al., PRL 107, 157202 (2011)

Requires Python 3 with `numpy`, `ase`, `spglib`.

## What it does
1. Finds all symmetry-inequivalent bonds between the magnetic atoms (`gen.py list`).
2. Builds a supercell where no periodic image of j is close to i, and writes one VASP folder per
   spin state (`gen.py build`).
3. Reads the results, puts every run in the slot of its **actual** final spins (the moment
   constraint can flip them), rejects runs with background spins off target, and gives the tensors.

| mode | runs / bond | result |
|---|---|---|
| `iso` | 4 (collinear) | J_iso |
| `tensor` | 36 (SOC) | 3×3 J, J_iso, DMI, anisotropy |
| `gsc` | 36 (SOC + Berry phase) | polarization mapping: M tensor of P = S_i · M · S_j |

Convention: H = −Σ S_i·J·S_j (unit spins, J > 0 ferromagnetic), as in TB2J.

## Use
```bash
python fourstate/gen.py list  --poscar POSCAR --mag Cr --rmax 8
python fourstate/gen.py build --poscar POSCAR --mag Cr --bonds 0 1 --mode tensor \
       --incar templates/INCAR.tensor --moment 3 --potcar-dir POTCARS --out runs
bash   fourstate/submit.sh runs --exe vasp_ncl --dry      # then without --dry
python fourstate/analyze.py runs                         # J  (mode iso / tensor)
python fourstate/analyze_gsc.py runs                     # M  (mode gsc)
```
* Edit the INCAR templates (ENCUT, +U, RWIGS) for your system; keep the `{MAGMOM}` line. Species
  order in the runs = order in your input POSCAR.
* `--potcar-dir DIR` expects `DIR/<Element>/POTCAR`.
* `--moment Cr=3 Mn=4` for several magnetic species; `--supercell`, `--kpts` to override the defaults.
* `--relax NSW` relaxes the ions in every state; `analyze.py --step first` then gives J at the input
  geometry. The standard method uses one fixed geometry.
* `LORBIT = 11` is required (the analysis checks every moment).

## Tests
```bash
python tests/test_energy.py   # known J tensor + spin flip
python tests/test_gsc.py      # known M tensor + polarization quantum
```

## License
MIT
