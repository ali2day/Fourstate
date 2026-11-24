#!/usr/bin/env python
"""One line per run (like an eng_press.dat summary): energy, penalty, spins of i and j, polarization.

  python table.py RUNS [--step first|last] > RUNS/table.dat
"""
import argparse, glob, json, os
import numpy as np
from common import read_outcar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs')
    ap.add_argument('--step', choices=['first', 'last'], default='last')
    a = ap.parse_args()
    print('%-28s %16s %10s %4s %6s  %-22s %-22s %s' % ('run', 'E(sigma->0) eV', 'E_p eV', 'conv', 'steps',
                                                    'S_i (muB)', 'S_j (muB)', 'P = p_elc + p_ion (|e| A)'))
    for st_file in sorted(glob.glob(f'{a.runs}/*/*/state.json')):
        d = os.path.dirname(st_file)
        st = json.load(open(st_file))
        tag = '/'.join(d.replace('\\', '/').split('/')[-2:])
        o = read_outcar(f'{d}/OUTCAR', len(st['target']), 0 if a.step == 'first' else -1) \
            if os.path.exists(f'{d}/OUTCAR') else None
        if o is None:
            print('%-28s  not finished' % tag); continue
        f = lambda v: '(%+.2f %+.2f %+.2f)' % tuple(v)
        mj = f(o['m'][st['J']]) if st.get('J') is not None else '-'
        p = f(o['dip']) if o['dip'] is not None else '-'
        print('%-28s %16.6f %10.6f %4s %6d  %-22s %-22s %s' % (tag, o['E'], o['Ep'], 'y' if o['conv'] else 'n',
                                                             o['nsteps'], f(o['m'][st['I']]), mj, p))


if __name__ == '__main__':
    main()
