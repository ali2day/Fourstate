# Examples

| folder | system | what | runs |
|---|---|---|---|
| `NiI2_exchange_and_ME` | NiI2 monolayer, 5×5 cell | exchange tensor **and** ME (GSC) tensor of the 1NN Ni–Ni bond from the same 36 SOC runs (`--mode gsc`): energies → `analyze.py`, Berry-phase polarizations → `analyze_gsc.py` | 36 |
| `exchange_tensor_VBr3` | bulk VBr3, 5 K experimental structure | full exchange tensor (J_iso, DMI, anisotropy) of the three J1 bonds | 108 |

Run from inside a folder:
```bash
export POTCARS=/path/to/potcars      # contains <Element>/POTCAR
bash run.sh                          # generates runs/, dry-run submit
```
Adapt `templates/job.sh` (modules, partition) to your cluster first.
