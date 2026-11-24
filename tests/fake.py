"""Write fake OUTCAR/OSZICAR files for the self-tests, in the format VASP prints them.

The energy written is E + E_p (VASP includes the constraint penalty); E_p goes to OSZICAR only.
Per-SCF energy lines with a wrong value are added: the reader must take the converged line.
"""


def outcar(d, E, m, dip=None, ep=0.0123):
    with open(f'{d}/OUTCAR', 'w') as f:
        f.write('  free energy    TOTEN  =       999.0 eV\n'
                '  energy without entropy =      999.0  energy(sigma->0) =      999.0\n')
        for c, comp in enumerate('xyz'):
            f.write(f' magnetization ({comp})\n \n# of ion       s       p       d       tot\n'
                    '------------------------------------------\n')
            for k in range(len(m)):
                f.write('    %d   0.0 0.0 %.3f %.3f\n' % (k + 1, m[k, c], m[k, c]))
        if dip is not None:
            f.write(' Total electronic dipole moment: p[elc]=( %.10f %.10f %.10f ) |e| Angst\n' % tuple(dip))
            f.write('            Ionic dipole moment: p[ion]=(    1.00000    2.00000    3.00000 ) |e| Angst\n')
        f.write('  free  energy   TOTEN  =  %.10f eV\n' % (E + ep))
        f.write('  energy  without entropy=  %.10f  energy(sigma->0) =  %.10f\n' % (E + ep, E + ep))
        f.write(' aborting loop because EDIFF is reached\n')
    open(f'{d}/OSZICAR', 'w').write('RMM:  10  -0.1E+03\n E_p =  %.6E  lambda =  0.100E+02\n'
                                    '   1 F= %.6E E0= %.6E  d E =0.0  mag= 0 0 0\n' % (ep, E + ep, E + ep))
