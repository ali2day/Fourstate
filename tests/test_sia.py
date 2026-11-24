#!/usr/bin/env python
"""Self-test: single-ion anisotropy A (H = S.A.S) with exchange to the background and one reversed run;
atom i chosen by hand with --atoms in the NiI2 example cell (Ali's 6x3 cell, atom 9)."""
import glob, json, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSCAR = f'{ROOT}/examples/NiI2_exchange_and_ME/POSCAR'


def main():
    out = tempfile.mkdtemp()
    subprocess.check_call([sys.executable, f'{ROOT}/fourstate/gen.py', 'build', '--poscar', POSCAR, '--mag', 'Ni',
                           '--atoms', '9', '--mode', 'sia', '--incar', f'{ROOT}/templates/INCAR.tensor',
                           '--moment', '2', '--out', f'{out}/runs'])
    rng = np.random.default_rng(3)
    A0 = rng.normal(0, 0.3, (3, 3)); A0 = (A0 + A0.T) / 2; A0 -= np.eye(3) * np.trace(A0) / 3      # meV
    K = rng.normal(0, 1.0, (3, 3))
    for d in glob.glob(f'{out}/runs/*/*/'):
        st = json.load(open(f'{d}/state.json'))
        t = np.array(st['target'], float)
        I = st['I']
        if os.path.basename(os.path.normpath(d)) == 'zx_xp':      # spin i reversed -> slot xm
            t[I] *= -1
        E = t[I] @ A0 @ t[I]
        for k in range(len(t)):
            if np.linalg.norm(t[k]) > 0 and k != I:
                E -= t[I] @ K @ t[k]
        outcar(d, E / 1000 - 50.0, 2.0 * t)
    r = subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_sia.py', f'{out}/runs'],
                       capture_output=True, text=True, check=True)
    A = np.array(json.load(open(f'{out}/runs/results_sia.json'))[0]['A'])
    # the reversed xp run fills xm twice -> xp missing -> A_zz-A_xx unknown
    assert np.isnan(A[0, 0]), 'reversed run not detected'
    assert np.allclose([A[0, 1], A[0, 2], A[1, 2]], [A0[0, 1], A0[0, 2], A0[1, 2]], atol=1e-6)
    # now the same run without the reversal: full tensor
    d = glob.glob(f'{out}/runs/*/zx_xp/')[0]
    st = json.load(open(f'{d}/state.json')); t = np.array(st['target'], float); I = st['I']
    E = t[I] @ A0 @ t[I] - sum(t[I] @ K @ t[k] for k in range(len(t)) if np.linalg.norm(t[k]) > 0 and k != I)
    outcar(d, E / 1000 - 50.0, 2.0 * t)
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_sia.py', f'{out}/runs'], capture_output=True, check=True)
    A = np.array(json.load(open(f'{out}/runs/results_sia.json'))[0]['A'])
    assert np.allclose(A, A0, atol=1e-6), (A, A0)
    print('test_sia PASSED (off-diagonal A recovered, reversed spin detected; --atoms works)')


if __name__ == '__main__':
    main()
