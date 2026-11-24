#!/usr/bin/env python
"""Self-test: polarizations from a known GSC tensor M (P = S_i M S_j) plus background and single-site
terms and one jump by a polarization quantum -> analyze_gsc.py must return M."""
import glob, json, os, subprocess, sys, tempfile
import numpy as np
from ase.build import bulk
from ase.io import read
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    out = tempfile.mkdtemp()
    bulk('Fe', 'bcc', a=2.87, cubic=True).write(f'{out}/POSCAR', format='vasp')
    subprocess.check_call([sys.executable, f'{ROOT}/fourstate/gen.py', 'build', '--poscar', f'{out}/POSCAR',
                           '--mag', 'Fe', '--bonds', '0', '--mode', 'gsc', '--incar',
                           f'{ROOT}/templates/INCAR.gsc', '--supercell', '2', '2', '2', '--out', f'{out}/runs'])
    rng = np.random.default_rng(2)
    M0 = rng.normal(0, 1e-3, (3, 3, 3))          # e*A
    K = rng.normal(0, 1e-3, (3, 3, 3))           # bond spins - background
    Q = rng.normal(0, 1e-3, (3, 3, 3))           # single site
    for n, d in enumerate(sorted(glob.glob(f'{out}/runs/bond*/*/'))):
        st = json.load(open(f'{d}/state.json'))
        t = np.array(st['target'], float)
        I, J = st['I'], st['J']
        P = np.einsum('a,abk,b->k', t[I], M0, t[J])
        for k in range(len(t)):
            if np.linalg.norm(t[k]) > 0 and k not in (I, J):
                P += np.einsum('a,abk,b->k', t[I], K, t[k]) + np.einsum('a,abk,b->k', t[J], K, t[k])
        P += np.einsum('a,abk,b->k', t[I], Q, t[I]) + np.einsum('a,abk,b->k', t[J], Q, t[J])
        if n == 5:                                 # Berry-phase branch jump: + one lattice vector
            P = P + read(f'{d}/POSCAR').cell.array[1]
        outcar(d, -100.0, 2.0 * t, dip=P + np.array([0.3, -0.2, 0.1]) - np.array([1.0, 2.0, 3.0]))
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_gsc.py', f'{out}/runs'], check=True,
                   capture_output=True)
    M = np.array(json.load(open(f'{out}/runs/results_gsc.json'))[0]['M'])
    assert np.allclose(M, M0, atol=1e-9), np.abs(M - M0).max()
    print('test_gsc PASSED (GSC tensor recovered, polarization quantum removed)')


if __name__ == '__main__':
    main()
