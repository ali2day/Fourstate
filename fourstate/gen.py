#!/usr/bin/env python
"""Four-state method (energy mapping for exchange, polarization mapping for the GSC model):
VASP input generator for any magnetic material.

Refs: Xiang et al., PRB 84, 224429 (2011); Dalton Trans. 42, 823 (2013);
      GSC model: Xiang et al., PRL 107, 157202 (2011).

Steps
  1. find all pairs of magnetic atoms up to --rmax in the input cell
  2. group them into symmetry-inequivalent bonds (spglib)
  3. choose the smallest diagonal supercell with no periodic image of j within --rclean of i
  4. write one VASP folder per spin state (+ state.json with the target spins)

Modes
  iso    : collinear, 4 runs per bond             -> isotropic J
  tensor : noncollinear + SOC, 36 runs per bond   -> full 3x3 J (J_iso, DMI, anisotropy)
  gsc    : same 36 states as tensor, with Berry-phase polarization -> GSC tensor M (P = S_i M S_j)
  me     : 12 runs, background +z: (S_i,S_j) along (x,y), (x,z), (y,z) -> spin-current ME matrix
           P = M_sc (S_i x S_j); k5, k7, k9, k11 (S_j = +z = background) are the same for every j
  sia    : only atom i rotated (20 runs)              -> single-ion anisotropy A (H = S A S)
  bq     : S_i = x, S_j at 0/180/45/225 deg, others z (4 runs) -> biquadratic K (H = S_i.J.S_j + K (S_i.S_j)^2)
Pair: --bonds (symmetry + automatic supercell) or --atoms i j (by hand, in the given cell)
For component (a, b): S_i = +-a, S_j = +-b, all other magnetic atoms along c, perpendicular to a and b.

Usage
  python gen.py list  --poscar POSCAR --mag Cr --rmax 8
  python gen.py build --poscar POSCAR --mag Cr --bonds 0 1 --mode tensor \
                      --incar templates/INCAR.tensor --moment 3 --out runs
"""
import argparse, json, os, itertools
import numpy as np
import spglib
from ase.io import read, write
from ase.build import make_supercell

AX = {'x': np.array([1., 0, 0]), 'y': np.array([0, 1., 0]), 'z': np.array([0, 0, 1.])}


# ---------------------------------------------------------------- bonds
def pairs(atoms, mag, rmax):
    """All (i, j, R, d) with i, j magnetic and 0 < |r_j + R - r_i| <= rmax."""
    cell, pos = atoms.cell.array, atoms.positions
    idx = [k for k, s in enumerate(atoms.get_chemical_symbols()) if s in mag]
    n = [int(np.ceil(rmax / (atoms.cell.volume / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3]))))) + 1
         for k in range(3)]
    out = []
    for i in idx:
        for j in idx:
            for R in itertools.product(*[range(-m, m + 1) for m in n]):
                d = np.linalg.norm(pos[j] + np.array(R) @ cell - pos[i])
                if 1e-3 < d <= rmax:
                    out.append((i, j, tuple(R), d))
    return out


def classes(atoms, mag, rmax, symprec=1e-3):
    """Symmetry-inequivalent bonds, sorted by distance."""
    sym = spglib.get_symmetry((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=symprec)
    frac, nat = atoms.get_scaled_positions(), len(atoms)

    def image(k, f):
        for m in range(nat):
            dR = f - frac[m]
            if np.allclose(dR, np.round(dR), atol=1e-3) and atoms.numbers[m] == atoms.numbers[k]:
                return m, tuple(np.round(dR).astype(int))
        raise RuntimeError('symmetry image not found (try another --symprec)')

    P = pairs(atoms, mag, rmax)
    keys = {(p[0], p[1], p[2]) for p in P}
    seen, out = set(), []
    for p in sorted(P, key=lambda x: (round(x[3], 4), x[0], x[1], x[2])):
        if (p[0], p[1], p[2]) in seen:
            continue
        orbit = set()
        for Rm, t in zip(sym['rotations'], sym['translations']):
            mi, Ri = image(p[0], Rm @ frac[p[0]] + t)
            mj, Rj = image(p[1], Rm @ (frac[p[1]] + np.array(p[2])) + t)
            R = tuple(np.array(Rj) - np.array(Ri))
            orbit.add((mi, mj, R))
            orbit.add((mj, mi, tuple(-np.array(R))))          # (i,j) ~ (j,i)
        seen |= orbit
        out.append(dict(i=p[0], j=p[1], R=[int(x) for x in p[2]], d=round(p[3], 4),
                        mult=len(orbit & keys)))
    for n, c in enumerate(out):
        c['id'] = n
    return out


# ---------------------------------------------------------------- supercell
def supercell_ok(atoms, P, bond, rclean):
    """No image of j (other than the bond partner) within rclean of i."""
    L = P @ atoms.cell.array
    pi = atoms.positions[bond['i']]
    pj = atoms.positions[bond['j']] + np.array(bond['R']) @ atoms.cell.array
    return all(np.linalg.norm(pj + np.array(n) @ L - pi) >= rclean
               for n in itertools.product(range(-2, 3), repeat=3) if n != (0, 0, 0))


def choose_supercell(atoms, bonds, rclean, nmax=4):
    ok = [n for n in itertools.product(range(1, nmax + 1), repeat=3)
          if all(supercell_ok(atoms, np.diag(n), b, rclean) for b in bonds)]
    if not ok:
        raise SystemExit('no diagonal supercell up to %d found; lower --rclean or use --supercell' % nmax)
    return np.diag(min(ok, key=np.prod))


def map_bond(atoms, P, bond):
    """Supercell and the indices (I, J) of the bond atoms in it."""
    sc = make_supercell(atoms, P)
    pi = atoms.positions[bond['i']]
    pj = atoms.positions[bond['j']] + np.array(bond['R']) @ atoms.cell.array
    sf = sc.get_scaled_positions()
    IJ = []
    for f in sc.cell.scaled_positions(np.array([pi, pj])):
        d = sf - f
        d -= np.round(d)
        IJ.append(int(np.argmin(np.linalg.norm(d @ sc.cell.array, axis=1))))
    I, J = IJ
    assert abs(sc.get_distance(I, J, mic=True) - bond['d']) < 1e-3, 'bond mapping failed'
    return sc, I, J


# ---------------------------------------------------------------- states
R2 = 1 / np.sqrt(2)


def states(mode):
    """(label, a, b, c, si, sj): S_i = si*a, S_j = sj*b, background along c.
    sia: (label, vi, None, c, 0, 0) with vi the unit vector of S_i (only atom i is rotated)."""
    sg = [(1, 1, 'pp'), (1, -1, 'pm'), (-1, 1, 'mp'), (-1, -1, 'mm')]
    if mode == 'me':       # spin-current ME scheme (12 runs): background always +z, as in the NiI2 V_I study
        n = 0              # k1-4: (x, y)  k5-8: (x, z)  k9-12: (y, z);  each ++ +- -+ --
        for a, b in (('x', 'y'), ('x', 'z'), ('y', 'z')):
            for si, sj, t in sg:
                n += 1
                yield f'k{n}', a, b, 'z', si, sj
        return
    if mode == 'bq':       # biquadratic: S_i = x, S_j at angle 0, 180, 45, 225 deg in the xy plane, others z
        t = R2 * (AX['x'] + AX['y'])
        for lab, vj in (('k1', AX['x']), ('k2', -AX['x']), ('k3', t), ('k4', -t)):
            yield lab, AX['x'], vj, 'z', 0, 0
        return
    if mode == 'sia':
        for a, b, c in (('z', 'x', 'y'), ('y', 'x', 'z')):          # A_aa - A_bb, background c
            for lab, v in ((f'{a}p', AX[a]), (f'{a}m', -AX[a]), (f'{b}p', AX[b]), (f'{b}m', -AX[b])):
                yield f'{a}{b}_{lab}', v, None, c, 0, 0
        for a, b, c in (('x', 'y', 'z'), ('x', 'z', 'y'), ('y', 'z', 'x')):   # A_ab: tilted 45 deg
            for si, sj, t in sg:
                yield f't{a}{b}_{t}', R2 * (si * AX[a] + sj * AX[b]), None, c, si, sj
        return
    if mode == 'iso':
        for si, sj, t in sg:
            yield t, 'z', 'z', 'z', si, sj
        return
    for a in 'xyz':
        for b in 'xyz':
            c = [x for x in 'xyz' if x not in (a, b)][0]
            for si, sj, t in sg:
                yield f'{a}{b}_{t}', a, b, c, si, sj


def write_run(d, sc, mag, mom, I, J, a, b, c, si, sj, mode, incar, potcar_dir, kpts, relax):
    os.makedirs(d, exist_ok=True)
    spec = list(dict.fromkeys(sc.get_chemical_symbols()))          # species order of the input POSCAR
    order = np.argsort([spec.index(s) for s in sc.get_chemical_symbols()], kind='stable')
    s2 = sc[order]
    inv = {int(o): k for k, o in enumerate(order)}
    I2, J2 = inv[I], (inv[J] if J is not None else None)
    write(f'{d}/POSCAR', s2, format='vasp', direct=True, sort=False)
    syms = s2.get_chemical_symbols()
    vec = []
    for k, s in enumerate(syms):
        v = np.zeros(3)
        if s in mag:
            v = AX[c] * mom[s]
            if mode == 'sia':
                if k == I2: v = a * mom[s]
            elif mode == 'bq':
                if k == I2: v = a * mom[s]
                if k == J2: v = b * mom[s]
            else:
                if k == I2: v = si * AX[a] * mom[s]
                if k == J2: v = sj * AX[b] * mom[s]
        vec.append(v)
    if mode == 'iso':
        line = 'MAGMOM = ' + ' '.join('%.3f' % v[2] for v in vec)
    else:
        f = ' '.join('%.3f %.3f %.3f' % tuple(v) for v in vec)
        line = f'MAGMOM = {f}\nM_CONSTR = {f}'
    txt = open(incar).read()
    if '{MAGMOM}' not in txt:
        raise SystemExit('INCAR template needs a {MAGMOM} line')
    txt = txt.replace('{MAGMOM}', line)
    if relax:
        txt += ('\n# ionic relaxation in this spin state, cell fixed\nIBRION = 2\nISIF = 2\n'
                f'NSW = {relax}\nEDIFFG = -0.01\n')
    open(f'{d}/INCAR', 'w', newline='\n').write(txt)
    open(f'{d}/KPOINTS', 'w', newline='\n').write('Gamma-centred\n0\nGamma\n%d %d %d\n0 0 0\n' % tuple(kpts))
    if potcar_dir:
        with open(f'{d}/POTCAR', 'w', newline='\n') as fh:
            for s in spec:
                fh.write(open(os.path.join(potcar_dir, s, 'POTCAR')).read())
    target = [list(v / np.linalg.norm(v)) if np.linalg.norm(v) > 0 else [0, 0, 0] for v in vec]
    if mode in ('sia', 'bq'):
        a = [float(x) for x in a]
        b = [float(x) for x in b] if b is not None else None
    json.dump(dict(I=I2, J=J2, a=a, b=b, c=c, si=si, sj=sj, mode=mode, target=target),
              open(f'{d}/state.json', 'w'), indent=1)


def moments(mag, spec):
    """--moment 3  or  --moment Cr=3 Mn=4"""
    out = {}
    for t in spec:
        if '=' in t:
            k, v = t.split('=')
            out[k] = float(v)
        else:
            out.update({m: float(t) for m in mag})
    missing = [m for m in mag if m not in out]
    if missing:
        raise SystemExit(f'no --moment for {missing}')
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['list', 'build'])
    ap.add_argument('--poscar', required=True, help='structure (any format ASE reads)')
    ap.add_argument('--mag', nargs='+', required=True, help='magnetic element(s), e.g. Cr  or  Fe Co')
    ap.add_argument('--rmax', type=float, default=8.0, help='longest bond to list (A)')
    ap.add_argument('--symprec', type=float, default=1e-3)
    ap.add_argument('--bonds', type=int, nargs='*', help='bond ids from "list"')
    ap.add_argument('--atoms', type=int, nargs='+', metavar='N',
                    help='choose i (and j) yourself: atom numbers (1-based, POSCAR order) in the GIVEN cell, '
                         'which is used as it is (e.g. a defect supercell); one atom for --mode sia')
    ap.add_argument('--rclean', type=float, default=None,
                    help='no image of j closer than this to i (A); default = longest selected bond')
    ap.add_argument('--supercell', type=int, nargs=3, help='force a diagonal supercell')
    ap.add_argument('--mode', choices=['iso', 'tensor', 'gsc', 'me', 'sia', 'bq'], default='tensor',
                    help='sia = single-ion anisotropy of atom i (20 runs); bq = biquadratic K of a pair (4 runs)')
    ap.add_argument('--incar', help='INCAR template with a {MAGMOM} line')
    ap.add_argument('--potcar-dir', help='folder with <Element>/POTCAR')
    ap.add_argument('--moment', nargs='+', default=['3.0'], help='initial moment: 3  or  Cr=3 Mn=4')
    ap.add_argument('--kspacing', type=float, default=0.25, help='Gamma grid density (1/A)')
    ap.add_argument('--kpts', type=int, nargs=3, help='force the Gamma k-grid')
    ap.add_argument('--relax', type=int, default=0, metavar='NSW',
                    help='relax ions in every state (not standard; not allowed for gsc)')
    ap.add_argument('--out', default='runs')
    a = ap.parse_args()

    atoms = read(a.poscar)
    cl = classes(atoms, a.mag, a.rmax, a.symprec)
    if a.cmd == 'list':
        sg = spglib.get_spacegroup((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), a.symprec)
        print(f'space group {sg};  {len(cl)} inequivalent bonds between {"/".join(a.mag)} up to {a.rmax} A')
        print(' id   d(A)   i   j   R             mult')
        for c in cl:
            print('%3d  %6.3f %3d %3d   %-13s %3d' % (c['id'], c['d'], c['i'], c['j'], c['R'], c['mult']))
        return

    if not ((a.bonds or a.atoms) and a.incar):
        raise SystemExit('build needs --incar and --bonds or --atoms')
    if a.mode in ('gsc', 'me') and a.relax:
        raise SystemExit('gsc needs one fixed geometry for all states (no --relax)')
    mom = moments(a.mag, a.moment)

    if a.atoms:                                        # user-chosen atoms in the given cell
        P = np.diag(a.supercell or [1, 1, 1])
        sc = make_supercell(atoms, P)
        I = a.atoms[0] - 1
        J = a.atoms[1] - 1 if len(a.atoms) > 1 else None
        for k in [I] + ([J] if J is not None else []):
            if sc.get_chemical_symbols()[k] not in a.mag:
                raise SystemExit(f'atom {k + 1} is {sc.get_chemical_symbols()[k]}, not in --mag')
        if a.mode != 'sia' and J is None:
            raise SystemExit('--atoms needs two atoms (i j) for this mode')
        d = float(sc.get_distance(I, J, mic=True)) if J is not None else 0.0
        if J is not None:
            near = sorted(np.linalg.norm(sc.positions[J] + np.array(n) @ sc.cell.array - sc.positions[I])
                          for n in itertools.product(range(-2, 3), repeat=3))
            if len(near) > 1 and near[1] < (a.rclean or d) + 1e-6:
                print('WARNING: another periodic image of j is only %.2f A from i (bond %.2f A)' % (near[1], d))
        jobs = [(dict(id=-1, i=I + 1, j=(J + 1 if J is not None else None), d=round(d, 4)), sc, I, J,
                 f'atoms{I + 1}_{J + 1}_d{d:.3f}' if J is not None else f'site{I + 1}')]
    else:
        sel = [c for c in cl if c['id'] in a.bonds]
        rclean = a.rclean or max(b['d'] for b in sel)
        P = np.diag(a.supercell) if a.supercell else choose_supercell(atoms, sel, rclean)
        for b in sel:
            if not supercell_ok(atoms, P, b, rclean):
                print('WARNING bond %d: an image of j is closer than %.2f A in this supercell' % (b['id'], rclean))
        jobs = []
        for b in sel:
            sc, I, J = map_bond(atoms, P, b)
            if a.mode == 'sia':
                J = None
            jobs.append((b, sc, I, J, f'bond{b["id"]:02d}_d{b["d"]:.3f}' if a.mode != 'sia' else f'site{I + 1}'))
    sc0 = jobs[0][1]
    rec = np.linalg.norm(np.linalg.inv(sc0.cell.array).T, axis=1)            # 1/A (no 2 pi)
    kpts = a.kpts or [max(1, int(np.ceil(x / a.kspacing))) for x in rec]
    print('cell', np.diag(P).tolist(), '| atoms', len(sc0), '| k-grid', kpts)
    os.makedirs(a.out, exist_ok=True)
    meta = dict(poscar=os.path.abspath(a.poscar), mag=a.mag, supercell=np.diag(P).tolist(), mode=a.mode,
                relax=a.relax, kpts=kpts, bonds=[])
    for b, sc, I, J, bd in jobs:
        n = 0
        for lab, x, y, z, si, sj in states(a.mode):
            write_run(f'{a.out}/{bd}/{lab}', sc, a.mag, mom, I, J, x, y, z, si, sj, a.mode, a.incar,
                      a.potcar_dir, kpts, a.relax)
            n += 1
        meta['bonds'].append(dict(b, dir=bd, runs=n))
        print(f'  {bd}: {n} runs -> {a.out}/{bd}')
    json.dump(meta, open(f'{a.out}/meta.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
