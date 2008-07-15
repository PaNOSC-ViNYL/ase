from ase import *

a = 5.0
d = 1.9
c = a / 2
atoms = Atoms('AuH',
              positions=[(c, c, 0), (c, c, d)],
              cell=(a, a, 2 * d),
              pbc=(0, 0, 1))
atoms *= (1, 1, 2)
images = [atoms.copy(), atoms.copy()]

r = ['xyz', 'traj', 'cube']
w = r + ['pdb', 'png', 'eps', 'xsf']
for format in w:
    print format
    write('io-test.1', atoms, format=format)
    if format not in ['cube', 'png', 'eps']:
        write('io-test.2', images, format=format)

    if format in r:
        a1 = read('io-test.1')
        if format != 'cube':
            a2 = read('io-test.2')
            a3 = read('io-test.2', index=0)
            a4 = read('io-test.2', index=slice(None))
