from __future__ import print_function
import sys

import numpy as np
from ase import Atoms
from ase.io import write, read
from ase.calculators.singlepoint import SinglePointCalculator

a = 5.0
d = 1.9
c = a / 2
atoms = Atoms('AuH',
              positions=[(c, c, 0), (c, c, d)],
              cell=(a, a, 2 * d),
              pbc=(0, 0, 1))
extra = np.array([2.3, 4.2])
atoms.set_array('extra', extra)
atoms *= (1, 1, 2)
images = [atoms.copy(), atoms.copy()]
r = ['xyz', 'traj', 'cube', 'pdb', 'cfg', 'struct',
     'cif', 'gen', 'extxyz', 'res']

# attach some results to the Atoms. These are serialised by the extxyz writer.
spc = SinglePointCalculator(atoms,
                            energy=-1.0,
                            stress=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                            forces=-1.0*atoms.get_positions())
atoms.set_calculator(spc)

try:
    import json
except ImportError:
    pass
else:
    r += ['json', 'db']

try:
    import Scientific
    version = Scientific.__version__.split('.')
    print('Found ScientificPython version: ', Scientific.__version__)
    if list(map(int, version)) < [2, 8]:
        print('ScientificPython 2.8 or greater required for numpy support')
        raise ImportError
except ImportError:
    print('No Scientific python found. Check your PYTHONPATH')
else:
    r += ['etsf']

w = r + ['xsf',
         'findsym']
try:
    import matplotlib
except ImportError:
    pass
else:
    w += ['png', 'eps']

if sys.version_info[0] == 3:
    r.remove('cif')
    
only_one_image = ['cube', 'png', 'eps', 'cfg', 'struct', 'etsf', 'gen',
                  'json', 'db', 'res']

for format in w:
    print(format, 'O', end=' ')
    fname1 = 'io-test.1.' + format
    fname2 = 'io-test.2.' + format
    if format == 'xsf':
        # The xsf format supports only pbc=000, 100, 110, or 111.
        # No crazy ideas like 001.  So let's try all the combinations
        # writing all the files on top of each other.
        atoms1 = atoms.copy()
        for pbc in [1, 0, 0], [1, 1, 0], [1, 1, 1]:
            atoms1.pbc = pbc
            images1 = [atoms1] * 2
            write(fname2, images1, format='xsf')
            atoms2 = read(fname2)
    else:
        write(fname1, atoms, format=format)
        if format not in only_one_image:
            write(fname2, images, format=format)

    if format in r:
        print('I')
        a1 = read(fname1)
        assert np.all(np.abs(a1.get_positions() -
                             atoms.get_positions()) < 1e-6)
        if format in ['traj', 'cube', 'cfg', 'struct', 'gen', 'extxyz']:
            assert np.all(np.abs(a1.get_cell() - atoms.get_cell()) < 1e-6)
        if format in ['cfg', 'extxyz']:
            assert np.all(np.abs(a1.get_array('extra') -
                                 atoms.get_array('extra')) < 1e-6)
        if format in ['extxyz']:
            assert np.all(a1.get_pbc() == atoms.get_pbc())
            assert np.all(a1.get_potential_energy() == atoms.get_potential_energy())
            assert np.all(a1.get_stress() == atoms.get_stress())
            assert np.all(abs(a1.get_forces() - atoms.get_forces()) < 1e-6)

        if format not in only_one_image:
            a2 = read(fname2)
            a3 = read(fname2, index=0)
            a4 = read(fname2, index=slice(None))
            if format in ['cif'] and sys.platform in ['win32']:
                # Fails on Windows:
                # https://trac.fysik.dtu.dk/projects/ase/ticket/62
                pass
            else:
                assert len(a4) == 2
    else:
        print()
