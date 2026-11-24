# Example

`NiI2_exchange_and_ME`: monolayer NiI2, 6×3 cell (18 Ni), pairs chosen by hand.

| set | atoms | mode | runs |
|---|---|---|---|
| exchange tensor | Ni 9 – Ni 6 | `tensor` (`INCAR.exchange`) | 36 |
| magnetoelectric matrix | Ni 9 – Ni 6 | `me` (`INCAR.me`, Berry phase, fixed geometry) | 12 |
| single-ion anisotropy | Ni 9 | `sia` | 20 |
| biquadratic exchange | Ni 9 – Ni 12 | `bq` | 4 |

```bash
cd NiI2_exchange_and_ME
export POTCARS=/path/to/potcars      # contains Ni/POTCAR, I/POTCAR
bash run.sh                          # generates the four sets, dry-run submit
```
Adapt `templates/job.sh` (modules, partition) to your cluster first.
