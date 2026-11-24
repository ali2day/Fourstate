#!/usr/bin/env python
"""Self-test of the 12-run ME scheme (--mode me, background +z): P = M_sc (S_i x S_j) + terms linear in
S_i or S_j (incl. coupling to the background) + one jump by a polarization quantum + one spin flip.
Also checks that the states are those of the NiI2 V_I scripts (k1: S_i=+x S_j=+y ... k12: -y, -z)."""
import glob, json, os, subprocess, sys, tempfile
import numpy as np
from ase.io import read
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSCAR = f'{ROOT}/examples/NiI2_exchange_and_ME/POSCAR'
EXPECT = {1: ('+x', '+y'), 2: ('+x', '-y'), 3: ('-x', '+y'), 4: ('-x', '-y'), 5: ('+x', '+z'), 6: ('+x', '-z'),
          7: ('-x', '+z'), 8: ('-x', '-z'), 9: ('+y', '+z'), 10: ('+y', '-z'), 11: ('-y', '+z'), 12: ('-y', '-z')}


def vec(s):
    return (1 if s[0] == '+' else -1) * np.eye(3)['xyz'.index(s[1])]


def main():
    out = tempfile.mkdtemp()
    subprocess.check_call([sys.executable, f'{ROOT}/fourstate/gen.py', 'build', '--poscar', POSCAR, '--mag', 'Ni',
                           '--atoms', '9', '6', '--mode', 'me', '--incar', f'{ROOT}/examples/NiI2_exchange_and_ME/INCAR.me',
                           '--moment', '2', '--out', f'{out}/runs'])
    rng = np.random.default_rng(5)
    Msc = rng.normal(0, 1e-3, (3, 3))            # P = Msc @ (S_i x S_j)
    Li, Lj = rng.normal(0, 1e-3, (3, 3)), rng.normal(0, 1e-3, (3, 3))
    for d in glob.glob(f'{out}/runs/*/k*/'):
        st = json.load(open(f'{d}/state.json'))
        t = np.array(st['target'], float)
        I, J = st['I'], st['J']
        n = int(os.path.basename(os.path.normpath(d))[1:])
        assert np.allclose(t[I], vec(EXPECT[n][0])) and np.allclose(t[J], vec(EXPECT[n][1])), (n, t[I], t[J])
        assert np.allclose(t[[k for k in range(18) if k not in (I, J)]], [0, 0, 1]), 'background must be +z'
        if n == 2:
            t[J] *= -1                              # flip: this run is state k1
        P = Msc @ np.cross(t[I], t[J]) + Li @ t[I] + Lj @ t[J]
        if n == 10:
            P = P + read(f'{d}/POSCAR').cell.array[0]
        outcar(d, -100.0, 2.0 * t, dip=P)
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_gsc.py', f'{out}/runs'], capture_output=True, check=True)
    r = json.load(open(f'{out}/runs/results_gsc.json'))[0]
    M = np.array(r['M_spin_current'])
    assert np.isnan(M[:, 2]).all(), 'flipped run not detected'
    assert np.allclose(M[:, :2], Msc[:, :2], atol=1e-9), (M, Msc)
    print('test_me PASSED (NiI2 k1-k12 states, ME matrix recovered, flip + polarization quantum handled)')


if __name__ == '__main__':
    main()
