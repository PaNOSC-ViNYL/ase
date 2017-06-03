from ase.build import molecule
from ase.io.mpl import plot_atoms
system = molecule('H2O')

ax = plot_atoms(system)
print(ax)
