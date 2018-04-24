import os

if 'DISPLAY' not in os.environ:
    from ase.test.testsuite import NotAvailable
    raise NotAvailable('No $DISPLAY on which to plot')

import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.lattice.cubic import FaceCenteredCubic

slab = FaceCenteredCubic('Au', size=(2, 2, 2))

fig, ax = plt.subplots()
plot_atoms(slab, ax, radii=0.5, rotation=('10x,10y,10z'))

assert len(ax.patches) == len(slab)
print(ax)
