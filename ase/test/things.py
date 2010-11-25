from ase.dft import monkhorst_pack

assert [0, 0, 0] in  monkhorst_pack((1, 3, 5)).tolist()
assert [0, 0, 0] not in  monkhorst_pack((1, 3, 6)).tolist()
assert len(monkhorst_pack((3, 4, 6))) == 3 * 4 * 6

from ase.units import Hartree, Bohr, kJ, mol, kcal, kB, fs
print Hartree, Bohr, kJ/mol, kcal/mol, kB*300, fs, 1/fs

from ase.structure import bulk
ru = bulk('Ru', 'hcp', a=2.7) * (2, 2, 1)
assert abs(ru.get_distance(0, 7, mic=True) - ru.get_distance(1, 6)) < 1e-14
assert abs(ru.get_distance(0, 5, mic=True) - 2.7) < 1e-14
