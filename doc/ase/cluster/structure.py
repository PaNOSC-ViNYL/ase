# creates: culayer.png truncated.png
from ase.io import write
from ase.cluster.cubic import FaceCenteredCubic
import numpy as np

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [6, 9, 5]
lc = 3.61000
culayer = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)
culayer.rotate('x', 0.1, rotate_cell=True)
culayer.rotate('y', 0.04, rotate_cell=True)
write('culayer.pov', culayer, show_unit_cell=2, display=False, run_povray=True)

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [6, 5, -1]
trunc = FaceCenteredCubic('Cu', surfaces, layers)
trunc.rotate('x', 0.1, rotate_cell=True)
trunc.rotate('y', 0.04, rotate_cell=True)
write('truncated.pov', trunc, show_unit_cell=2, display=False, run_povray=True)

