#!/usr/bin/env python
"""GSC (generalized spin-current) tensor by four-state POLARIZATION mapping (Berry phase, not energies).

  python analyze_gsc.py RUNS [--cosmin 0.9]

Model (Xiang et al., PRL 107, 157202 (2011)):  P_ij = S_i^T M S_j, i.e. P_ij = sum_ab M_ab S_ia S_jb
with M_ab a 3-vector (3x3x3 tensor).  Runs made with  gen.py build --mode gsc  (36 per bond,
one fixed geometry, LCALCPOL in the INCAR).
    M_ab = [P(+a,+b) + P(-a,-b) - P(+a,-b) - P(-a,+b)] / 4
P = p[elc] + p[ion] from the OUTCAR (e*A).  Each P is first brought to the branch of a reference
run by removing polarization quanta e*R (R = supercell lattice vectors).
--mode me (12 runs, background +z): only M_xy, M_xz, M_yz; printed as the 3x3 matrix
P = M (S_i x S_j) with M[:,z] = (P1+P4-P2-P3)/4, M[:,y] = (P6+P7-P5-P8)/4, M[:,x] = (P9+P12-P10-P11)/4.
Also printed (gsc): the spin-current part M_anti = (M_ab - M_ba)/2 written as P = M_sc (S_i x S_j)
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
    ap.add_argument('--elc-only', action='store_true',
                    help='use p[elc] only (as output_vasp_polar.sh); the ionic part is the same in every '
                         'state at fixed geometry, so M does not change')
    a = ap.parse_args()
    res = []
    for b, meta, slots, notes, cell in collect(a.runs, a.cosmin, -1):
        key = 'dip_elc' if a.elc_only else 'dip'
        P = {k: o[key] for k, (o, _) in slots.items() if o[key] is not None}
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
                if not np.isnan(M[AX[x], AX[y]]).any():
                    print('   M_%s%s = (%+9.1f %+9.1f %+9.1f)' % (x, y, *(1e5 * M[AX[x], AX[y]])))
        if meta['mode'] == 'me':
            # 12-run scheme: only (x,y), (x,z), (y,z) are computed; with P = M_sc (S_i x S_j):
            #   column z <- M_xy,  column y <- -M_xz (x cross z = -y),  column x <- M_yz
            sc = np.array([M[1, 2], -M[0, 2], M[0, 1]]).T
            r['M_spin_current'] = sc.tolist()
            print('   ME matrix (P = M . S_i x S_j), rows P_x,y,z, columns (S_i x S_j)_x,y,z:')
            print(np.array2string(1e5 * sc, precision=1, suppress_small=True, prefix='   '))
        elif not np.isnan(M).any():
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
