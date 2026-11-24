# Examples

| folder | system | what | runs |
|---|---|---|---|
| `NiI2_exchange_and_ME` | NiI2 monolayer, 6×3 cell (18 Ni) | 1NN bond Ni 9 – Ni 6 chosen by hand (`--atoms 9 6`): exchange tensor (`INCAR.exchange`), ME matrix by Berry-phase polarization mapping (12 runs k1–k12, `INCAR.me`, fixed geometry), single-ion anisotropy of Ni 9, biquadratic K of Ni 9 – Ni 12 | 36 + 12 + 20 + 4 |
| `exchange_tensor_VBr3` | bulk VBr3, 5 K experimental structure | full exchange tensor (J_iso, DMI, anisotropy) of the three J1 bonds | 108 |

Run from inside a folder:
```bash
export POTCARS=/path/to/potcars      # contains <Element>/POTCAR
bash run.sh                          # generates runs/, dry-run submit
```
Adapt `templates/job.sh` (modules, partition) to your cluster first.
