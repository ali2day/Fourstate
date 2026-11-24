#!/usr/bin/env python
"""Single-ion anisotropy tensor from the 20 runs of  gen.py build --mode sia.

  python analyze_sia.py RUNS [--step first|last] [--ignore-penalty] [--cosmin 0.9]

H_SIA = S_i^T A S_i (unit spin; only atom i is rotated, all other spins along the third axis c).
  A_zz - A_xx = [E(+z) + E(-z) - E(+x) - E(-x)] / 2      (background y)
  A_yy - A_xx = [E(+y) + E(-y) - E(+x) - E(-x)] / 2      (background z)
  A_ab        = [E(++) + E(--) - E(+-) - E(-+)] / 4      (S_i = (+-a +-b)/sqrt 2, background c)
Exchange with the background is linear in S_i and cancels.  A is reported traceless (meV);
the eigenvector with the lowest eigenvalue is the easy axis.  Written to RUNS/results_sia.json.
A run whose spin i ended antiparallel to its target is put in the opposite slot.
"""
import argparse, glob, json, os
import numpy as np
from common import read_outcar

NEG = {'zp': 'zm', 'zm': 'zp', 'yp': 'ym', 'ym': 'yp', 'xp': 'xm', 'xm': 'xp',
       'pp': 'mm', 'mm': 'pp', 'pm': 'mp', 'mp': 'pm'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--step', choices=['first', 'last'], default='last')
    ap.add_argument('--cosmin', type=float, default=0.9)
    ap.add_argument('--ignore-penalty', action='store_true')
    a = ap.parse_args()
    meta = json.load(open(f'{a.runs}/meta.json'))
    res = []
    for b in meta['bonds']:
        bd = os.path.join(a.runs, os.path.basename(b['dir']))
        E, notes = {}, []
        for d in sorted(glob.glob(f'{bd}/*')):
            if not os.path.isfile(f'{d}/state.json'):
                continue
            st = json.load(open(f'{d}/state.json'))
            tag = os.path.basename(d)
            t = np.array(st['target'])
            o = read_outcar(f'{d}/OUTCAR', len(t), 0 if a.step == 'first' else -1) \
                if os.path.exists(f'{d}/OUTCAR') else None
            if o is None:
                notes.append(f'{tag}: not finished'); continue
            if not o['conv']:
                notes.append(f'{tag}: not converged')
            m = o['m']
            bad = [k for k in np.where(np.linalg.norm(t, axis=1) > 0)[0] if k != st['I'] and
                   (np.linalg.norm(m[k]) < 0.5 or m[k] @ t[k] / np.linalg.norm(m[k]) < a.cosmin)]
            if bad:
                notes.append(f'{tag}: {len(bad)} background spins off target'); continue
            n = np.linalg.norm(m[st['I']])
            c = m[st['I']] @ t[st['I']] / n if n > 0 else 0
            group, lab = tag.split('_')
            if c <= -a.cosmin:
                lab = NEG[lab]
                notes.append(f'{tag}: spin i reversed (relabelled {group}_{lab})')
            elif c < a.cosmin:
                notes.append(f'{tag}: spin i off target (cos {c:.2f})'); continue
            E[(group, lab)] = 1000 * (o['E'] - (0 if a.ignore_penalty else o['Ep']))   # meV
        g = lambda *k: all(q in E for q in k)
        A = np.full((3, 3), np.nan)
        d1 = (E[('zx', 'zp')] + E[('zx', 'zm')] - E[('zx', 'xp')] - E[('zx', 'xm')]) / 2 \
            if g(('zx', 'zp'), ('zx', 'zm'), ('zx', 'xp'), ('zx', 'xm')) else np.nan
        d2 = (E[('yx', 'yp')] + E[('yx', 'ym')] - E[('yx', 'xp')] - E[('yx', 'xm')]) / 2 \
            if g(('yx', 'yp'), ('yx', 'ym'), ('yx', 'xp'), ('yx', 'xm')) else np.nan
        axx = -(d1 + d2) / 3
        A[0, 0], A[1, 1], A[2, 2] = axx, axx + d2, axx + d1
        for (p, q), grp in (((0, 1), 'txy'), ((0, 2), 'txz'), ((1, 2), 'tyz')):
            k = [(grp, s) for s in ('pp', 'mm', 'pm', 'mp')]
            if g(*k):
                A[p, q] = A[q, p] = (E[k[0]] + E[k[1]] - E[k[2]] - E[k[3]]) / 4
        r = dict(site=b.get('i'), A=A.tolist(), Azz_minus_Axx=d1, Ayy_minus_Axx=d2, notes=notes)
        print(f'\n{os.path.basename(bd)}   H = S.A.S  (meV, traceless)')
        print(np.array2string(A, precision=4, suppress_small=True, prefix='   '))
        if not np.isnan(A).any():
            w, v = np.linalg.eigh(A)
            r['easy_axis'] = v[:, 0].tolist()
            r['eigenvalues'] = w.tolist()
            print('   eigenvalues', np.round(w, 4), ' easy axis', np.round(v[:, 0], 3))
        for n_ in notes:
            print('   !', n_)
        res.append(r)
    json.dump(res, open(f'{a.runs}/results_sia.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
