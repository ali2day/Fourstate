#!/usr/bin/env python
"""Biquadratic exchange from the 4 runs of  gen.py build --mode bq.

  python analyze_bq.py RUNS [--sia RUNS_SIA/results_sia.json] [--step first|last] [--ignore-penalty]

Model: H = -J S_i.S_j - K (S_i.S_j)^2 + S_j.A.S_j + ...   (unit spins; K > 0 favours collinear spins)
States: S_i = x, others along z, S_j at 0 (k1), 180 (k2), 45 (k3), 225 deg (k4) in the xy plane,
so c = S_i.S_j = 1, -1, 1/sqrt2, -1/sqrt2.  Coupling to the background is linear in S_j and cancels:
    (E1 + E2) - (E3 + E4) = -K + 2 [a_x - a_45]
    a_x = A_xx,  a_45 = (A_xx + A_yy)/2 + A_xy       (single-ion anisotropy of j)
    K   = -[(E1 + E2) - (E3 + E4)] + 2 (a_x - a_45)
The anisotropy term is taken from --sia (results of  --mode sia  for atom j, or an equivalent site);
without it, it is assumed zero and a warning is printed.  Also  J ~ (E2 - E1)/2  (= J_xx).
A run whose spin j ended reversed is relabelled (k1<->k2, k3<->k4); reversed spin i -> rejected.
"""
import argparse, glob, json, os
import numpy as np
from common import read_outcar

SWAP = {'k1': 'k2', 'k2': 'k1', 'k3': 'k4', 'k4': 'k3'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--sia', help='results_sia.json of the j site (single-ion anisotropy correction)')
    ap.add_argument('--step', choices=['first', 'last'], default='last')
    ap.add_argument('--cosmin', type=float, default=0.9)
    ap.add_argument('--ignore-penalty', action='store_true')
    a = ap.parse_args()
    corr, cnote = 0.0, 'no single-ion correction (--sia not given): K contains 2(a_x - a_45) of atom j'
    if a.sia:
        A = np.array(json.load(open(a.sia))[0]['A'])
        corr = 2 * (A[0, 0] - ((A[0, 0] + A[1, 1]) / 2 + A[0, 1]))
        cnote = f'single-ion correction 2(a_x - a_45) = {corr:+.4f} meV from {a.sia}'
    meta = json.load(open(f'{a.runs}/meta.json'))
    res = []
    for b in meta['bonds']:
        bd = os.path.join(a.runs, os.path.basename(b['dir']))
        E, notes = {}, []
        for d in sorted(glob.glob(f'{bd}/k?')):
            st = json.load(open(f'{d}/state.json'))
            t = np.array(st['target'])
            tag = os.path.basename(d)
            o = read_outcar(f'{d}/OUTCAR', len(t), 0 if a.step == 'first' else -1) \
                if os.path.exists(f'{d}/OUTCAR') else None
            if o is None:
                notes.append(f'{tag}: not finished'); continue
            if not o['conv']:
                notes.append(f'{tag}: not converged')
            m = o['m']
            I, J = st['I'], st['J']
            bad = [k for k in np.where(np.linalg.norm(t, axis=1) > 0)[0] if k not in (I, J) and
                   (np.linalg.norm(m[k]) < 0.5 or m[k] @ t[k] / np.linalg.norm(m[k]) < a.cosmin)]
            if bad:
                notes.append(f'{tag}: {len(bad)} background spins off target'); continue
            sg = []
            for k in (I, J):
                c = m[k] @ t[k] / max(np.linalg.norm(m[k]), 1e-9)
                sg.append(1 if c >= a.cosmin else (-1 if c <= -a.cosmin else 0))
            if 0 in sg:
                notes.append(f'{tag}: spin i or j off its direction'); continue
            if sg[0] == -1:            # S_i must stay along +x (its coupling to the background must cancel)
                notes.append(f'{tag}: spin i reversed (run rejected)'); continue
            lab = tag if sg[1] == 1 else SWAP[tag]
            if lab != tag:
                notes.append(f'{tag}: one spin reversed (relabelled {lab})')
            E[lab] = 1000 * (o['E'] - (0 if a.ignore_penalty else o['Ep']))
        r = dict(id=b['id'], d=b['d'], notes=notes)
        print(f'\n{os.path.basename(bd)}   d = {b["d"]:.3f} A')
        if all(k in E for k in ('k1', 'k2', 'k3', 'k4')):
            K = -((E['k1'] + E['k2']) - (E['k3'] + E['k4'])) + corr
            Jxx = (E['k2'] - E['k1']) / 2
            r.update(K=K, J_xx=Jxx, sia_correction=corr)
            print('   K (biquadratic) = %+.4f meV     J_xx = %+.3f meV' % (K, Jxx))
        print('   ' + cnote)
        for n in notes:
            print('   !', n)
        res.append(r)
    json.dump(res, open(f'{a.runs}/results_bq.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
