"""Shared helpers: read VASP OUTCARs of 4-state runs and put every run in the slot of its ACTUAL spins."""
import glob, json, os, re
import numpy as np

AX = {'x': 0, 'y': 1, 'z': 2}
SLOTS = ((1, 1), (1, -1), (-1, 1), (-1, -1))       # pp, pm, mp, mm


def read_outcar(path, nat, step=-1):
    """Energy, constraint penalty, moments, dipole of ionic step `step` (0 = first, -1 = last)."""
    txt = open(path, errors='replace').read()
    # converged energy of every ionic step: 'energy  without entropy=  ...  energy(sigma->0) = E'
    # (the same line the eng_press scripts read; per-SCF lines are written differently)
    fin = list(re.finditer(r'energy  without entropy=\s*\S+\s+energy\(sigma->0\)\s*=\s*(\S+)', txt))
    if not fin:
        fin = list(re.finditer(r'free  energy   TOTEN  =\s+(\S+) eV', txt))
    if not fin:
        return None
    e = [m.group(1) for m in fin]
    ep = ['0']                         # constraint penalty: only in OSZICAR, last E_p before ionic step `step`
    osz = os.path.join(os.path.dirname(path), 'OSZICAR')
    if os.path.exists(osz):
        cur, per_step = None, []
        for ln in open(osz, errors='replace'):
            m = re.search(r'E_p\s*=\s*(\S+)', ln)
            if m:
                cur = m.group(1)
            if ' F= ' in ln:
                per_step.append(cur or '0')
        if per_step and len(per_step) >= (step + 1 if step >= 0 else -step):
            ep = [per_step[step]]
    tables = {}
    for comp in 'xyz':
        blocks = txt.split(f' magnetization ({comp})')[1:]
        rows_all = []
        for blk in blocks:
            rows = [float(l.split()[-1]) for l in blk.split('\n')[4:4 + nat] if l.split() and l.split()[0].isdigit()]
            if len(rows) == nat:
                rows_all.append(rows)
        if rows_all:
            tables[comp] = rows_all
    m = np.zeros((nat, 3))
    if 'y' in tables:                                   # noncollinear
        for comp in 'xyz':
            m[:, AX[comp]] = tables[comp][step]
    elif 'x' in tables:                                 # collinear: one table = m_z
        m[:, 2] = tables['x'][step]
    nsteps = len(re.findall(r'F=', open(os.path.join(os.path.dirname(path), 'OSZICAR')).read())) \
        if os.path.exists(os.path.join(os.path.dirname(path), 'OSZICAR')) else len(e)
    num = r'\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)'
    pe = re.findall(r'p\[elc\]=' + num, txt)
    pi = re.findall(r'p\[ion\]=' + num, txt)
    dip = dip_elc = None
    if pe:
        dip_elc = np.array(pe[step], float)
        dip = dip_elc + (np.array(pi[step], float) if pi else 0)
    conv = 'aborting loop because EDIFF is reached' in txt
    if re.search(r'^\s*NSW\s*=\s*[1-9]', txt, re.M):
        conv = conv and 'reached required accuracy' in txt
    return dict(E=float(e[step]), Ep=float(ep[0]), m=m, dip=dip, dip_elc=dip_elc, conv=conv,
                nsteps=nsteps)


def spin_signs(st, m, cosmin):
    """Actual signs of the two bond spins along their axes, and background atoms off target."""
    t = np.array(st['target'])
    bad = []
    for k in np.where(np.linalg.norm(t, axis=1) > 0)[0]:
        if k in (st['I'], st['J']):
            continue
        n = np.linalg.norm(m[k])
        if n < 0.5 or np.dot(m[k], t[k]) / n < cosmin:
            bad.append(int(k))
    sgn = []
    for k, ax in ((st['I'], st['a']), (st['J'], st['b'])):
        n = np.linalg.norm(m[k])
        c = m[k][AX[ax]] / n if n > 0 else 0
        sgn.append(int(np.sign(c)) if abs(c) >= cosmin and n >= 0.5 else 0)
    return sgn, bad


def collect(runs, cosmin=0.9, step=-1):
    """Yield (bond, {(a, b, si, sj): (out, state)}, notes, cell) for every bond of a runs folder."""
    meta = json.load(open(f'{runs}/meta.json'))
    for b in meta['bonds']:
        bd = os.path.join(runs, os.path.basename(b['dir']))
        slots, notes, cell = {}, [], None
        for d in sorted(glob.glob(f'{bd}/*')):
            if not os.path.isfile(f'{d}/state.json'):
                continue
            st = json.load(open(f'{d}/state.json'))
            tag = os.path.basename(d)
            o = read_outcar(f'{d}/OUTCAR', len(st['target']), step) if os.path.exists(f'{d}/OUTCAR') else None
            if o is None:
                notes.append(f'{tag}: not finished'); continue
            if cell is None and os.path.exists(f'{d}/POSCAR'):
                from ase.io import read
                cell = read(f'{d}/POSCAR').cell.array
            if not o['conv']:
                notes.append(f'{tag}: not converged')
            sgn, bad = spin_signs(st, o['m'], cosmin)
            if 0 in sgn:
                notes.append(f'{tag}: bond spin off axis or collapsed'); continue
            if bad:
                notes.append(f'{tag}: {len(bad)} background spins off target'); continue
            if tuple(sgn) != (st['si'], st['sj']):
                notes.append(f'{tag}: spin flip, target {st["si"]:+d}{st["sj"]:+d} actual {sgn[0]:+d}{sgn[1]:+d} (relabelled)')
            key = (st['a'], st['b'], sgn[0], sgn[1])
            if key in slots:
                notes.append(f'{tag}: duplicate of state {key} (kept the first)'); continue
            slots[key] = (o, st)
        yield b, meta, slots, notes, cell
