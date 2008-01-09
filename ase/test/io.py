from ase import *

a = 5.0
d = 1.9
c = a / 2
atoms = Atoms(positions=[(c, c, 0), (c, c, d)],
              symbols='AuH',
              cell=(a, a, 2 * d),
              pbc=(0, 0, 1))
atoms *= (1, 1, 2)
images = [atoms, atoms]

rw = ['traj']
try:
    import Scientific.IO.NetCDF
except ImportError:
    pass
else:
    rw += ['nc']

w = ['xyz', 'cube', 'png', 'eps']
for format in rw + w:
    write('io-test.x', atoms, format=format)
