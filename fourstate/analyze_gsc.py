#!/usr/bin/env python
"""GSC (generalized spin-current) tensor by four-state POLARIZATION mapping (Berry phase, not energies).

  python analyze_gsc.py RUNS [--cosmin 0.9]

Model (Xiang et al., PRL 107, 157202 (2011)):  P_ij = S_i^T M S_j, i.e. P_ij = sum_ab M_ab S_ia S_jb
with M_ab a 3-vector (3x3x3 tensor).  Runs made with  gen.py build --mode gsc  (36 per bond,
one fixed geometry, LCALCPOL in the INCAR).
    M_ab = [P(+a,+b) + P(-a,-b) - P(+a,-b) - P(-a,+b)] / 4
P = p[elc] + p[ion] from the OUTCAR (e*A).  Each P is first brought to the branch of a reference
run by removing polarization quanta e*R (R = supercell lattice vectors).
Also printed: the spin-current part M_anti = (M_ab - M_ba)/2 written as P = M_sc (S_i x S_j)
(3x3: column k = component (S_i x S_j)_k).  Written to RUNS/results_gsc.json (units e*A).
"""
import argparse, json
import numpy as np
from common import collect, AX


def fold(p, ref, cell):
    """p moved by lattice quanta to the branch closest to ref."""
    f = (p - ref) @ np.linalg.inv(cell)
    return ref + (f - np.round(f)) @ cell


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--cosmin', type=float, default=0.9)
    a = ap.parse_args()
    res = []
    for b, meta, slots, notes, cell in collect(a.runs, a.cosmin, -1):
        P = {k: o['dip'] for k, (o, _) in slots.items() if o['dip'] is not None}
        if len(P) < len(slots):
            notes.append('some runs have no p[elc] (LCALCPOL missing?)')
        M = np.full((3, 3, 3), np.nan)                  # M[a, b, :] = vector
        if P:
            ref = next(iter(P.values()))
            P = {k: fold(v, ref, cell) for k, v in P.items()}
        for x in 'xyz':
            for y in 'xyz':
                k = [(x, y, s, t) for s, t in ((1, 1), (-1, -1), (1, -1), (-1, 1))]
                if all(q in P for q in k):
                    M[AX[x], AX[y]] = (P[k[0]] + P[k[1]] - P[k[2]] - P[k[3]]) / 4
        r = dict(id=b['id'], d=b['d'], M=M.tolist(), notes=notes)
        print(f'\nbond {b["id"]:2d}  d = {b["d"]:.3f} A   (units 1e-5 e*A)')
        for x in 'xyz':
            for y in 'xyz':
                print('   M_%s%s = (%+9.1f %+9.1f %+9.1f)' % (x, y, *(1e5 * M[AX[x], AX[y]])))
        if not np.isnan(M).any():
            # P = M_sc (S_i x S_j):  (S_i x S_j)_k = eps_kab S_ia S_jb  ->  M_sc[:, k] = (M_ab - M_ba)/2
            sc = np.array([(M[1, 2] - M[2, 1]) / 2, (M[2, 0] - M[0, 2]) / 2, (M[0, 1] - M[1, 0]) / 2]).T
            r['M_spin_current'] = sc.tolist()
            print('   spin-current part (P = M_sc . S_i x S_j), columns = (S_i x S_j)_x,y,z:')
            print(np.array2string(1e5 * sc, precision=1, suppress_small=True, prefix='   '))
        for n in notes:
            print('   !', n)
        res.append(r)
    json.dump(res, open(f'{a.runs}/results_gsc.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
