#!/bin/bash
# Example Slurm job for one run. Edit partition / modules for your cluster.
#SBATCH --job-name=4s
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=2-00:00:00

module load VASP          # <- your VASP module

EXE=${VASP_EXE:-vasp_ncl} # vasp_std for --mode iso
srun $EXE > vasp.out 2>&1
if grep -q "^ *NSW *= *[1-9]" INCAR; then OK="reached required accuracy"; else OK="EDIFF is reached"; fi
grep -q "$OK" OUTCAR && echo DONE > STATUS || echo NOT_CONVERGED > STATUS
rm -f WAVECAR CHG CHGCAR PROCAR vaspout.h5
