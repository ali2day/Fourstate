#!/usr/bin/env python
"""Self-test of paper_convention.py (arXiv:2407.20489 convention) with a table in the 4states-script format.
Energies from H = S_i.A.S_i + S_j.A.S_j + S_i.J.S_j + K (S_i.S_j)^2 + coupling to the background (spins of 2 muB
in the scripts, unit vectors in H).  J, Gamma, A, K must come back exactly."""
import json, os, subprocess, sys, tempfile
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X, Y, Z = np.eye(3)
r2 = 1 / np.sqrt(2)


def main():
    rng = np.random.default_rng(7)
    J0 = np.array([[-3.9, 0.0, 0.0], [0.0, -3.0, -0.6], [0.0, -0.6, -3.6]])
    A0 = np.array([[0.0, 0.03, -0.02], [0.03, 0.05, 0.01], [-0.02, 0.01, 0.07]])
    K0 = -0.67
    Kb = rng.normal(0, 1, (3, 3))                       # i, j - background (same for all 16 others)

    def E(si, sj, bg):
        e = si @ A0 @ si + si @ J0 @ sj + K0 * (si @ sj) ** 2 + 16 * (si @ Kb @ bg + bg @ Kb @ sj)
        if sj is not None:
            e += sj @ A0 @ sj
        return e

    lines = []
    pairs = {'jxx': (X, X, Z), 'jyy': (Y, Y, Z), 'jzz': (Z, Z, X), 'dx': (Y, Z, X), 'dxb': (Z, Y, X),
             'dy': (X, Z, Y), 'dyb': (Z, X, Y), 'dz': (X, Y, Z), 'dzb': (Y, X, Z)}
    for lab, (a, b, c) in pairs.items():
        for n, (s, t) in zip('1234', ((1, 1), (1, -1), (-1, 1), (-1, -1))):
            lines.append((f'{lab}{n}', E(s * a, t * b, c)))
    ej = lambda si, bg: si @ A0 @ si + 16 * (si @ Kb @ bg)          # single ion: j = background
    for n, v in zip('1234', (Z, -Z, X, -X)):
        lines.append((f'a{n}', ej(v, Y)))
    for n, v in zip('1234', (Y, -Y, X, -X)):
        lines.append((f'ayyxx{n}', ej(v, Z)))
    for lab, (a, b, c) in {'axy': (X, Y, Z), 'axz': (X, Z, Y), 'ayz': (Y, Z, X)}.items():
        for n, (s, t) in zip('1234', ((1, 1), (1, -1), (-1, 1), (-1, -1))):
            lines.append((f'{lab}{n}', ej(r2 * (s * a + t * b), c)))
    t45 = r2 * (X + Y)
    for n, sj in zip('1234', (X, -X, t45, -t45)):
        lines.append((f'k{n}', E(X, sj, Z)))
    out = tempfile.mkdtemp()
    with open(f'{out}/eng.dat', 'w') as f:
        f.write('i Eng-2\n')
        for lab, e in lines:
            f.write(f'{lab} {e / 1000 - 264.0:.12f}\n')
    subprocess.run([sys.executable, f'{ROOT}/fourstate/paper_convention.py', '--table', f'{out}/eng.dat',
                    '--out', f'{out}/r.json'], check=True, capture_output=True)
    r = json.load(open(f'{out}/r.json'))
    assert np.allclose(r['J_tensor'], J0, atol=1e-6), r['J_tensor']
    assert abs(r['J'] - np.trace(J0) / 3) < 1e-6
    assert abs(r['A_zz'] - (A0[2, 2] - A0[0, 0])) < 1e-6, r['A_zz']
    assert abs(r['K'] - K0) < 1e-6, (r['K'], K0)
    print('test_paper PASSED (J tensor, A^zz and K in the arXiv:2407.20489 convention)')


if __name__ == '__main__':
    main()
