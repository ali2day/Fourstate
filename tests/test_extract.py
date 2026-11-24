#!/usr/bin/env python
"""Self-test: extract.py on run folders (nested like all_sub_9_12/6-k/k1) -> paper_convention.py.
Uses the NiI2 1NN biquadratic energies (k1-k4); K must be -0.699 meV and J_xx -3.330 meV."""
import json, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fake import outcar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
E = {'k1': -127.29046333, 'k2': -127.28380273, 'k3': -127.28920773, 'k4': -127.28435978}


def main():
    top = tempfile.mkdtemp()
    for k, e in E.items():
        d = os.path.join(top, 'all_sub_9_12', '6-k', k)
        os.makedirs(d)
        outcar(d, e, np.zeros((2, 3)), ep=0.0)
    tab = subprocess.run([sys.executable, f'{ROOT}/fourstate/extract.py', os.path.join(top, 'all_sub_9_12')],
                         capture_output=True, text=True, check=True).stdout
    open(f'{top}/eng_bq.dat', 'w').write(tab)
    assert len(tab.strip().split('\n')) == 5, tab
    subprocess.run([sys.executable, f'{ROOT}/fourstate/paper_convention.py', '--table', f'{top}/eng_bq.dat',
                    '--out', f'{top}/r.json'], check=True, capture_output=True)
    r = json.load(open(f'{top}/r.json'))
    assert abs(r['K'] - (-0.69855)) < 1e-4, r
    assert abs(r['J_xx_from_bq'] - (-3.3303)) < 1e-3, r
    print('test_extract PASSED (folders -> table -> K = %.3f meV, J_xx = %.3f meV)' % (r['K'], r['J_xx_from_bq']))


if __name__ == '__main__':
    main()
