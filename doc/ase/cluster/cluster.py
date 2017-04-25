# creates: culayer.png truncated.png

from ase.io import write
from ase.cluster.cubic import FaceCenteredCubic
#from ase.cluster.hexagonal import HexagonalClosedPacked
#import numpy as np

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [6, 9, 5]
lc = 3.61000
culayer = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)
culayer.rotate(6, 'x', rotate_cell=True)
culayer.rotate(2, 'y', rotate_cell=True)
write('culayer.pov', culayer, show_unit_cell=0, display=False, run_povray=True)

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [6, 5, -1]
trunc = FaceCenteredCubic('Cu', surfaces, layers)
trunc.rotate(6, 'x', rotate_cell=True)
trunc.rotate(2, 'y', rotate_cell=True)
write('truncated.pov', trunc, show_unit_cell=0, display=False, run_povray=True)

# This does not work!
#surfaces = [(0, 0, 0, 1), (1, 1, -2, 0), (1, 0, -1, 1)]
#layers = [6, 6, 6]
#graphite = Graphite('C', surfaces, layers, latticeconstant=(2.461, 6.708))
#write('graphite.pov', graphite, show_unit_cell=2, display=False, run_povray=True)

# surfaces = [(0, 0, 0, 1), (1, 1, -2, 0), (1, 0, -1, 1)]
# layers = [6, 6, 6]
# magn = HexagonalClosedPacked('Mg', surfaces, layers)
# magn.rotate('x', np.pi/2 - 0.1, rotate_cell=True)
# magn.rotate('y', 0.04, rotate_cell=True)
# write('magnesium.pov', magn, show_unit_cell=2, display=False, run_povray=True)
