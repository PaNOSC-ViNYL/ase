import matplotlib.pyplot as plt
from ase.io.matplotlib import write_matplotlib
from ase.lattice.cubic import FaceCenteredCubic

slab = FaceCenteredCubic('Au', size=(2, 2, 2))

fig, ax = plt.subplots()
write_matplotlib(slab, ax, radii=0.5, rotation=('10x,10y,10z'))

assert len(ax.patches) == len(slab)
