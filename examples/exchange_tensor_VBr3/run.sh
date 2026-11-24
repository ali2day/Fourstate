#!/bin/bash
# Full exchange tensor (J_iso, DMI, anisotropy) of the three nearest-neighbour V-V bonds of VBr3.
# 3 bonds x 36 SOC runs = 108 runs in a 2x2x2 supercell (64 atoms).
set -e
R=../..
python $R/fourstate/gen.py list  --poscar POSCAR --mag V --rmax 7.5
python $R/fourstate/gen.py build --poscar POSCAR --mag V --bonds 0 1 2 --mode tensor \
       --incar INCAR.tensor --moment 2 --supercell 2 2 2 --kpts 2 2 2 \
       --potcar-dir ${POTCARS:?set POTCARS=folder with V/POTCAR Br/POTCAR} --out runs
bash $R/fourstate/submit.sh runs --exe vasp_ncl --dry      # remove --dry to submit
# when finished:
# python $R/fourstate/analyze.py runs
