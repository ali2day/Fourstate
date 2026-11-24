#!/usr/bin/env python
"""Exchange from 4-state energies.

  python analyze.py RUNS [--step first|last] [--ignore-penalty] [--cosmin 0.9]

H = - sum_<ij> S_i^T J_ij S_j  (unit spins, J > 0 = FM, same as TB2J)
J_ab = [E(+a,-b) + E(-a,+b) - E(+a,+b) - E(-a,-b)] / 4
Output per bond: 3x3 J (meV), J_iso = Tr/3, DMI D = (J_yz-J_zy, J_zx-J_xz, J_xy-J_yx)/2,
symmetric traceless anisotropy.  Written to RUNS/results.json.
Every run is placed by the ACTUAL signs of its bond spins (constraints can flip them);
runs with background spins off target are rejected.  --step first uses the first ionic step
of relaxation runs (= the input geometry).
"""
import argparse, json
import numpy as np
from common import collect, AX


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--step', choices=['first', 'last'], default='last')
    ap.add_argument('--cosmin', type=float, default=0.9)
    ap.add_argument('--ignore-penalty', action='store_true', help='keep the constraint energy E_p')
    a = ap.parse_args()
    res = []
    for b, meta, slots, notes, _ in collect(a.runs, a.cosmin, 0 if a.step == 'first' else -1):
        E = {k: o['E'] - (0 if a.ignore_penalty else o['Ep']) for k, (o, _) in slots.items()}
        J = np.full((3, 3), np.nan)
        for x in 'xyz':
            for y in 'xyz':
                k = [(x, y, s, t) for s, t in ((1, -1), (-1, 1), (1, 1), (-1, -1))]
                if all(q in E for q in k):
                    J[AX[x], AX[y]] = 1000 * (E[k[0]] + E[k[1]] - E[k[2]] - E[k[3]]) / 4
        r = dict(id=b['id'], d=b['d'], J=J.tolist(), notes=notes)
        print(f'\nbond {b["id"]:2d}  d = {b["d"]:.3f} A')
        if meta['mode'] == 'iso':
            r['J_iso'] = J[2, 2]
            print('   J_iso = %+.3f meV' % J[2, 2])
        else:
            if not np.isnan(J).any():
                r['J_iso'] = float(np.trace(J) / 3)
                r['DMI'] = [float((J[1, 2] - J[2, 1]) / 2), float((J[2, 0] - J[0, 2]) / 2),
                            float((J[0, 1] - J[1, 0]) / 2)]
                S = (J + J.T) / 2
                r['J_ani'] = (S - np.eye(3) * np.trace(J) / 3).tolist()
                print('   J_iso = %+.3f meV' % r['J_iso'])
                print('   D     = (%+.3f, %+.3f, %+.3f) meV  |D| = %.3f' % (*r['DMI'], np.linalg.norm(r['DMI'])))
            print('   J (meV):'); print(np.array2string(J, precision=3, suppress_small=True, prefix='   '))
        for n in notes:
            print('   !', n)
        res.append(r)
    json.dump(res, open(f'{a.runs}/results.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
