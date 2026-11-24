#!/usr/bin/env python
"""Self-test: energies from a spin model with a known J tensor (incl. DMI) -> analyze.py must return it.
Extra couplings to all background spins, single-ion anisotropy and one sign-flipped run are included."""
import glob, json, os, subprocess, sys, tempfile
import numpy as np
from ase.build import bulk
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J0 = np.array([[2.0, 0.3, -0.1], [-0.5, 1.5, 0.2], [0.4, 0.0, 1.0]])     # meV


def main():
    out = tempfile.mkdtemp()
    bulk('Fe', 'bcc', a=2.87, cubic=True).write(f'{out}/POSCAR', format='vasp')
    subprocess.check_call([sys.executable, f'{ROOT}/fourstate/gen.py', 'build', '--poscar', f'{out}/POSCAR',
                           '--mag', 'Fe', '--bonds', '0', '--mode', 'tensor', '--incar',
                           f'{ROOT}/templates/INCAR.tensor', '--supercell', '2', '2', '2', '--out', f'{out}/runs'])
    rng = np.random.default_rng(1)
    K, A = rng.normal(0, 0.5, (3, 3)), np.diag([0.0, 0.1, 0.4])
    for d in glob.glob(f'{out}/runs/bond*/*/'):
        st = json.load(open(f'{d}/state.json'))
        t = np.array(st['target'], float)
        I, J = st['I'], st['J']
        if os.path.basename(os.path.normpath(d)) == 'xy_pp':      # spin flip of S_J
            t[J] *= -1
        E = t[I] @ J0 @ t[J]
        for k in range(len(t)):
            if np.linalg.norm(t[k]) > 0 and k not in (I, J):
                E -= t[I] @ K @ t[k] + t[J] @ K @ t[k]
        E -= t[I] @ A @ t[I] + t[J] @ A @ t[J]
        outcar(d, E / 1000 - 100.0, 2.0 * t)
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze.py', f'{out}/runs'], check=True,
                   capture_output=True)
    J = np.array(json.load(open(f'{out}/runs/results.json'))[0]['J'])
    assert np.isnan(J[0, 1]), 'flipped run not detected'
    J[0, 1] = J0[0, 1]
    assert np.allclose(J, J0, atol=1e-6), J
    print('test_energy PASSED (J tensor recovered, spin flip detected)')


if __name__ == '__main__':
    main()
