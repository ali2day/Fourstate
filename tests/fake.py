"""Write fake OUTCAR/OSZICAR files for the self-tests (format as VASP prints them)."""


def outcar(d, E, m, dip=None):
    with open(f'{d}/OUTCAR', 'w') as f:
        for c, comp in enumerate('xyz'):
            f.write(f' magnetization ({comp})\n \n# of ion       s       p       d       tot\n'
                    '------------------------------------------\n')
            for k in range(len(m)):
                f.write('    %d   0.0 0.0 %.3f %.3f\n' % (k + 1, m[k, c], m[k, c]))
        if dip is not None:
            f.write(' Total electronic dipole moment: p[elc]=( %.10f %.10f %.10f ) electrons Angst\n' % tuple(dip))
            f.write(' Ionic dipole moment: p[ion]=( 1.0 2.0 3.0 ) electrons Angst\n')
        f.write('  free  energy   TOTEN  =  %.10f eV\n aborting loop because EDIFF is reached\n' % E)
    open(f'{d}/OSZICAR', 'w').write('   1 F= %.6E\n' % E)
