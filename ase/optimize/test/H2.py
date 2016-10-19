from ase import Atoms
cell = (5, 5, 5)
atoms = Atoms('H2', [(0, 0, 0), (0, 0, 1.4)], cell=cell)
atoms.center()
