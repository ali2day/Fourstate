#!/usr/bin/env python
"""Exchange tensor, single-ion anisotropy and biquadratic exchange in the convention of
A. Ghojavand et al., "Strain-tunable magnetic and magnonic states in Ni-dihalide monolayers",
arXiv:2407.20489, Eqs. (1)-(2):

    H = sum_i S_i.A.S_i  +  sum_{i<j} S_i.J_ij.S_j  +  sum_{i<j} K_ij (S_i.S_j)^2
    S_i.J_ij.S_j = J S_i.S_j + D.(S_i x S_j) + S_i.Gamma.S_j
    J = Tr(J_ij)/3,   D from (J_ij - J_ij^T)/2,   Gamma = (J_ij + J_ij^T)/2 - J*1

* S_i are UNIT vectors, every pair counted once (i<j), energies in meV.
* J < 0 ferromagnetic, J > 0 antiferromagnetic (opposite sign to TB2J and to analyze.py).
* K < 0 favours collinear spins.  A^zz > 0 easy plane, A^zz < 0 easy axis.

Four-state formulas in this convention (S_i = +-a, S_j = +-b, all other spins along c, perpendicular):
    J_ab = [E(+a,+b) + E(-a,-b) - E(+a,-b) - E(-a,+b)] / 4
Biquadratic (S_i = x, S_j at 0, 180, 45, 225 deg in the xy plane, others along z; S_i.S_j = 1, -1, 1/sqrt2,
-1/sqrt2): the bilinear terms and the coupling to the background cancel in
    K = [E(k1) + E(k2)] - [E(k3) + E(k4)] - 2 (a_x - a_45),   a_x = A_xx,  a_45 = (A_xx + A_yy)/2 + A_xy
(a_x - a_45 is the single-ion anisotropy of S_j, which points along x in k1/k2 but at 45 deg in k3/k4).
Single-ion (only S_i rotated, background perpendicular):
    A^zz - A^xx = [E(+z) + E(-z) - E(+x) - E(-x)] / 2,  A^yy - A^xx likewise,
    A^ab = [E(++) + E(--) - E(+-) - E(-+)] / 4 with S_i = (+-a +-b)/sqrt2.
A^zz is reported with A^xx = A^yy = 0 (uniaxial form, as in the paper's Table 1).

Input, either
  --runs DIR [DIR ...]   folders made by gen.py (modes tensor / sia / bq): OUTCARs are read and every run is
                         placed by its ACTUAL final spins (constraints can flip them);
  --table FILE [...]     a summary with one line per run: label and energy in eV in the first two columns, labels
                         as in the 4states shell scripts: jxx1-4, jyy1-4, jzz1-4, dx/dxb, dy/dyb, dz/dzb 1-4
                         (1,2,3,4 = ++, +-, -+, --), a1-4, ayyxx1-4, axy/axz/ayz1-4 (single ion), k1-4 (biquadratic).
                         NOTE: a table has no spin moments, so sign flips cannot be checked.

  python paper_convention.py --table eng_press.dat
  python paper_convention.py --runs runs_exchange runs_sia runs_bq
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
AX = 'xyz'
# shell-script labels -> (a, b) of S_i = a, S_j = b
PAIR = {'jxx': 'xx', 'jyy': 'yy', 'jzz': 'zz', 'dx': 'yz', 'dxb': 'zy', 'dy': 'xz', 'dyb': 'zx',
        'dz': 'xy', 'dzb': 'yx'}
SG = {'1': (1, 1), '2': (1, -1), '3': (-1, 1), '4': (-1, -1)}


def from_table(files):
    E = {}
    for f in files:
        for ln in open(f):
            p = ln.split()
            if len(p) < 2:
                continue
            try:
                E[p[0]] = 1000 * float(p[1])        # meV
            except ValueError:
                continue                            # header
    pair, sia, bq = {}, {}, {}
    for lab, e in E.items():
        base, n = lab.rstrip('1234'), lab[len(lab.rstrip('1234')):]
        if base in PAIR and n in SG:
            pair[(PAIR[base][0], PAIR[base][1]) + SG[n]] = e
        elif base == 'k' and n in SG:
            bq['k' + n] = e
        elif base == 'a' and n in SG:               # a1 +z, a2 -z, a3 +x, a4 -x
            sia[('zx', {'1': 'zp', '2': 'zm', '3': 'xp', '4': 'xm'}[n])] = e
        elif base == 'ayyxx' and n in SG:           # +y, -y, +x, -x
            sia[('yx', {'1': 'yp', '2': 'ym', '3': 'xp', '4': 'xm'}[n])] = e
        elif base in ('axy', 'axz', 'ayz') and n in SG:
            sia[('t' + base[1:], {'1': 'pp', '2': 'pm', '3': 'mp', '4': 'mm'}[n])] = e
    return pair, sia, bq, ['table input: spin signs not checked']


def from_runs(dirs, step):
    """Use the analysis scripts on gen.py folders (actual spin signs)."""
    from common import collect, read_outcar
    pair, sia, bq, notes = {}, {}, {}, []
    for d in dirs:
        mode = json.load(open(f'{d}/meta.json'))['mode']
        if mode in ('tensor', 'iso', 'gsc'):
            for b, meta, slots, nts, _ in collect(d, 0.9, step):
                for k, (o, _) in slots.items():
                    pair[k] = 1000 * (o['E'] - o['Ep'])
                notes += nts
                break                               # first bond of the folder
        elif mode in ('sia', 'bq'):
            import subprocess
            tool = 'analyze_sia.py' if mode == 'sia' else 'analyze_bq.py'
            subprocess.run([sys.executable, f'{HERE}/{tool}', d, '--step', 'first' if step == 0 else 'last'],
                           check=True, capture_output=True)
            r = json.load(open(f'{d}/results_{mode}.json'))[0]
            notes += r['notes']
            if mode == 'sia':
                sia['_A'] = np.array(r['A'])        # H = S.A.S already, traceless
            else:
                bq['_K_tb'] = r.get('K')            # analyze_bq.py convention (H = -K c^2), no SIA correction
    return pair, sia, bq, notes


def sia_tensor(sia):
    if '_A' in sia:
        return sia['_A']
    g = lambda *k: all(q in sia for q in k)
    A = np.full((3, 3), np.nan)
    d1 = (sia[('zx', 'zp')] + sia[('zx', 'zm')] - sia[('zx', 'xp')] - sia[('zx', 'xm')]) / 2 \
        if g(('zx', 'zp'), ('zx', 'zm'), ('zx', 'xp'), ('zx', 'xm')) else np.nan
    d2 = (sia[('yx', 'yp')] + sia[('yx', 'ym')] - sia[('yx', 'xp')] - sia[('yx', 'xm')]) / 2 \
        if g(('yx', 'yp'), ('yx', 'ym'), ('yx', 'xp'), ('yx', 'xm')) else np.nan
    A[0, 0], A[1, 1], A[2, 2] = 0.0, d2, d1         # A_xx = 0 reference
    for (p, q), grp in (((0, 1), 'txy'), ((0, 2), 'txz'), ((1, 2), 'tyz')):
        k = [(grp, s) for s in ('pp', 'mm', 'pm', 'mp')]
        if g(*k):
            A[p, q] = A[q, p] = (sia[k[0]] + sia[k[1]] - sia[k[2]] - sia[k[3]]) / 4
    return A


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--table', nargs='+')
    src.add_argument('--runs', nargs='+')
    ap.add_argument('--step', choices=['first', 'last'], default='last')
    ap.add_argument('--out', default='paper_convention.json')
    a = ap.parse_args()
    pair, sia, bq, notes = from_table(a.table) if a.table else from_runs(a.runs, 0 if a.step == 'first' else -1)
    out = {'convention': 'arXiv:2407.20489: H = S.A.S + sum_{i<j} S_i.J.S_j + K (S_i.S_j)^2, unit spins, meV, '
                         'J<0 FM'}

    # ---- exchange tensor
    J = np.full((3, 3), np.nan)
    for x in AX:
        for y in AX:
            k = [(x, y, s, t) for s, t in ((1, 1), (-1, -1), (1, -1), (-1, 1))]
            if all(q in pair for q in k):
                J[AX.index(x), AX.index(y)] = (pair[k[0]] + pair[k[1]] - pair[k[2]] - pair[k[3]]) / 4
    if not np.isnan(J).all():
        np.set_printoptions(precision=3, suppress=True)
        print('Exchange tensor J_ij (meV):'); print(J)
        out['J_tensor'] = J.tolist()
        if not np.isnan(J).any():
            Jiso = np.trace(J) / 3
            D = [(J[1, 2] - J[2, 1]) / 2, (J[2, 0] - J[0, 2]) / 2, (J[0, 1] - J[1, 0]) / 2]
            G = (J + J.T) / 2 - Jiso * np.eye(3)
            out.update(J=Jiso, D=D, Gamma=G.tolist())
            print(f'J     = {Jiso:+.3f} meV   ({"FM" if Jiso < 0 else "AFM"})')
            print('D     = (%+.3f, %+.3f, %+.3f) meV' % tuple(D))
            print('Gamma =  xx %+.3f  yy %+.3f  zz %+.3f  xy %+.3f  xz %+.3f  yz %+.3f  meV'
                  % (G[0, 0], G[1, 1], G[2, 2], G[0, 1], G[0, 2], G[1, 2]))
            print(f'Gamma_yz / J = {G[1, 2] / Jiso:+.3f}')

    # ---- single-ion anisotropy
    A = sia_tensor(sia) if sia else None
    if A is not None:
        Azz = A[2, 2] - A[0, 0]
        out['A'] = A.tolist(); out['A_zz'] = Azz
        print(f'\nSingle-ion: A^zz = {Azz:+.3f} meV ({"easy plane" if Azz > 0 else "easy axis"}),  full A (A_xx = 0):')
        print(A)

    # ---- biquadratic
    if all(f'k{n}' in bq for n in '1234') or bq.get('_K_tb') is not None:
        raw = (bq['k1'] + bq['k2'] - bq['k3'] - bq['k4']) if 'k1' in bq else -bq['_K_tb']
        corr = 0.0
        if A is not None and not np.isnan(A[[0, 1, 0], [0, 1, 1]]).any():
            corr = 2 * (A[0, 0] - ((A[0, 0] + A[1, 1]) / 2 + A[0, 1]))
        else:
            notes.append('no single-ion data: K not corrected for 2(a_x - a_45)')
        K = raw - corr
        out.update(K=K, K_uncorrected=raw)
        print(f'\nBiquadratic: K = {K:+.3f} meV  (uncorrected {raw:+.3f}, single-ion correction {-corr:+.3f})')
        if 'J' in out:
            print(f'K / J = {K / out["J"]:+.3f}')
    for n in notes:
        print('!', n)
    json.dump(out, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
