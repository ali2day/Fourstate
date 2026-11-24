#!/bin/bash
# Monolayer NiI2 (6x3 cell, 18 Ni), nearest-neighbour bond between Ni 9 and Ni 6 (POSCAR numbering),
# chosen by hand as in the original 4states scripts (atom1 = 9, atom2 = 6, s = 2 muB).
#   exchange tensor  : energy mapping        (36 SOC runs, INCAR.exchange)
#   ME matrix        : polarization mapping  (12 SOC runs k1-k12, INCAR.me, Berry phase, fixed geometry)
#                      (--mode gsc instead of me gives the full 36-run tensor)
#   single-ion aniso : Ni 9 alone            (20 SOC runs)
#   biquadratic K    : Ni 9 - Ni 12          (4 SOC runs, k1-k4)
set -e
R=../..
P=${POTCARS:?set POTCARS=folder with Ni/POTCAR I/POTCAR}
python $R/fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 6 --mode tensor \
       --incar INCAR.exchange --moment 2 --kpts 2 2 1 --potcar-dir $P --out runs_exchange
python $R/fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 6 --mode me \
       --incar INCAR.me --moment 2 --kpts 2 2 1 --potcar-dir $P --out runs_me
python $R/fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 --mode sia \
       --incar INCAR.exchange --moment 2 --kpts 2 2 1 --potcar-dir $P --out runs_sia
python $R/fourstate/gen.py build --poscar POSCAR --mag Ni --atoms 9 12 --mode bq \
       --incar INCAR.exchange --moment 2 --kpts 2 2 1 --potcar-dir $P --out runs_bq
for r in runs_exchange runs_me runs_sia runs_bq; do bash $R/fourstate/submit.sh $r --exe vasp_ncl --dry; done
# when finished:
# python $R/fourstate/analyze.py     runs_exchange   # J tensor (J_iso, DMI, anisotropy)
# python $R/fourstate/analyze_gsc.py runs_me         # ME matrix, P = M (S_i x S_j)
# python $R/fourstate/analyze_sia.py runs_sia        # single-ion anisotropy A
# python $R/fourstate/analyze_bq.py  runs_bq --sia runs_sia/results_sia.json   # biquadratic K
# python $R/fourstate/table.py       runs_exchange   # one line per run (E, E_p, spins, P)
