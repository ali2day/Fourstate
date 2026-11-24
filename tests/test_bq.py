#!/usr/bin/env python
"""Self-test: biquadratic K with bilinear J, single-ion anisotropy of j (corrected with --sia),
coupling to the background and one run where S_j ended reversed."""
import glob, json, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSCAR = f'{ROOT}/examples/NiI2_exchange_and_ME/POSCAR'
J0, K0 = 1.7, 0.35                                        # meV


def main():
    out = tempfile.mkdtemp()
    subprocess.check_call([sys.executable, f'{ROOT}/fourstate/gen.py', 'build', '--poscar', POSCAR, '--mag', 'Ni',
                           '--atoms', '9', '12', '--mode', 'bq', '--incar', f'{ROOT}/templates/INCAR.tensor',
                           '--moment', '2', '--out', f'{out}/runs'])
    rng = np.random.default_rng(4)
    A = rng.normal(0, 0.3, (3, 3)); A = (A + A.T) / 2; A -= np.eye(3) * np.trace(A) / 3
    Kb = rng.normal(0, 1.0, (3, 3))
    json.dump([dict(A=A.tolist())], open(f'{out}/sia.json', 'w'))
    for d in glob.glob(f'{out}/runs/*/k?/'):
        st = json.load(open(f'{d}/state.json'))
        t = np.array(st['target'], float)
        I, J = st['I'], st['J']
        if os.path.basename(os.path.normpath(d)) == 'k3':     # S_j reversed -> this run is state k4
            t[J] *= -1
        c = t[I] @ t[J]
        E = -J0 * c - K0 * c ** 2 + t[J] @ A @ t[J] + t[I] @ A @ t[I]
        for k in range(len(t)):
            if np.linalg.norm(t[k]) > 0 and k not in (I, J):
                E -= t[I] @ Kb @ t[k] + t[J] @ Kb @ t[k]
        outcar(d, E / 1000 - 70.0, 2.0 * t)
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_bq.py', f'{out}/runs', '--sia', f'{out}/sia.json'],
                   capture_output=True, check=True)
    r = json.load(open(f'{out}/runs/results_bq.json'))[0]
    assert 'K' not in r, 'reversed run not detected (k3 missing, k4 duplicated)'
    d = glob.glob(f'{out}/runs/*/k3/')[0]                    # redo k3 properly -> K must be exact
    st = json.load(open(f'{d}/state.json')); t = np.array(st['target'], float); I, J = st['I'], st['J']
    c = t[I] @ t[J]
    E = -J0 * c - K0 * c ** 2 + t[J] @ A @ t[J] + t[I] @ A @ t[I] - sum(
        t[I] @ Kb @ t[k] + t[J] @ Kb @ t[k] for k in range(len(t)) if np.linalg.norm(t[k]) > 0 and k not in (I, J))
    outcar(d, E / 1000 - 70.0, 2.0 * t)
    subprocess.run([sys.executable, f'{ROOT}/fourstate/analyze_bq.py', f'{out}/runs', '--sia', f'{out}/sia.json'],
                   capture_output=True, check=True)
    r = json.load(open(f'{out}/runs/results_bq.json'))[0]
    assert abs(r['K'] - K0) < 1e-6, (r['K'], K0)
    print('test_bq PASSED (K recovered with the single-ion correction, reversed spin detected)')


if __name__ == '__main__':
    main()
