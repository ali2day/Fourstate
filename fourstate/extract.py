#!/usr/bin/env python
"""Collect the VASP results of a 4-state set into one table (one line per run folder).

  python extract.py all_sub_9_6 [more folders ...] > eng_press.dat
  python paper_convention.py --table eng_press.dat

Works for any folder whose sub-folders are VASP runs (jxx1, dx1, a1, k1, ... as made by the 4states shell
scripts, or the folders made by gen.py).  Columns, as in the eng_press scripts:
  i        run folder name (the label used by paper_convention.py)
  Eng-2    energy(sigma->0) of the last converged ionic step, eV  (the 'energy  without entropy' line)
  E_p      constraint penalty of that step (OSZICAR), eV
  P_x..z   Berry-phase p[elc] (|e| A), if LCALCPOL was used
  Mag_x..z total moment of the last step (OSZICAR 'mag=')
  conv     y if 'EDIFF is reached' (and 'reached required accuracy' when ions relax)
--first takes the first ionic step instead (= the input geometry when ions relax in every state).
--minus-penalty writes Eng-2 - E_p.
"""
import argparse, os, re, sys


def read(d, first):
    o = os.path.join(d, 'OUTCAR')
    if not os.path.exists(o):
        return None
    txt = open(o, errors='replace').read()
    e = re.findall(r'energy  without entropy=\s*\S+\s+energy\(sigma->0\)\s*=\s*(\S+)', txt)
    if not e:
        return None
    k = 0 if first else -1
    pe = re.findall(r'p\[elc\]=\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)', txt)
    ep, mag = '0', ['-', '-', '-']
    osz = os.path.join(d, 'OSZICAR')
    if os.path.exists(osz):
        cur, eps, mags = '0', [], []
        for ln in open(osz, errors='replace'):
            m = re.search(r'E_p\s*=\s*(\S+)', ln)
            if m:
                cur = m.group(1)
            if ' F= ' in ln:
                eps.append(cur)
                mm = ln.split('mag=')
                mags.append(mm[1].split() if len(mm) > 1 else ['-'])
        if eps:
            ep, mag = eps[k], (mags[k] + ['-', '-', '-'])[:3]
    conv = 'aborting loop because EDIFF is reached' in txt
    if re.search(r'^\s*NSW\s*=\s*[1-9]', txt, re.M):
        conv = conv and 'reached required accuracy' in txt
    return float(e[k]), float(ep), (pe[k] if pe else ('-', '-', '-')), mag, conv


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--first', action='store_true')
    ap.add_argument('--minus-penalty', action='store_true')
    a = ap.parse_args()
    print('i Eng-2 E_p P_x P_y P_z Mag_x Mag_y Mag_z conv')
    for top in a.dirs:
        runs = sorted(r for r, _, f in os.walk(top) if 'OUTCAR' in f or 'INCAR' in f)   # also 1-j/jxx1, ...
        for path in runs:
            if os.path.normpath(path) == os.path.normpath(top) and len(runs) > 1:
                continue                                    # the template INCAR of the parent folder
            s = os.path.basename(os.path.normpath(path))
            r = read(path, a.first)
            if r is None:
                print(f'# {s}: not finished', file=sys.stderr); continue
            E, ep, p, mag, conv = r
            if a.minus_penalty:
                E -= ep
            print(s, '%.8f' % E, '%.6f' % ep, *p, *mag, 'y' if conv else 'n')
            if not conv:
                print(f'# {s}: NOT converged', file=sys.stderr)


if __name__ == '__main__':
    main()
