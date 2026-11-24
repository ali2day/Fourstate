#!/bin/bash
# Monolayer NiI2: exchange tensor (energy mapping) AND magnetoelectric GSC tensor (Berry-phase
# polarization mapping) of the nearest-neighbour Ni-Ni bond, from the SAME 36 SOC runs.
# 5x5 supercell (75 atoms), Gamma 2x2x1, initial Ni moment 1.25 muB (avoids a high-moment SCF trap).
set -e
R=../..
python $R/fourstate/gen.py list  --poscar POSCAR --mag Ni --rmax 8
python $R/fourstate/gen.py build --poscar POSCAR --mag Ni --bonds 0 --mode gsc        --incar INCAR --moment 1.25 --supercell 5 5 1 --kpts 2 2 1        --potcar-dir ${POTCARS:?set POTCARS=folder with Ni/POTCAR I/POTCAR} --out runs
bash $R/fourstate/submit.sh runs --exe vasp_ncl --dry      # remove --dry to submit
# when finished:
# python $R/fourstate/analyze.py runs       # J tensor: J_iso, DMI, anisotropy  (energies)
# python $R/fourstate/analyze_gsc.py runs   # ME / GSC tensor M, P = S_i M S_j   (polarizations)
