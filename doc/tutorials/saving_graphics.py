from ase import *

atoms = Atoms('Ag', cell=(2.7, 2.7, 2.7), pbc=True) * (18, 8, 8)

# view with ag
#view(atoms)
rotation = '-70x, -20y, -2z' # found using 'view -> rotate'

#Make colors
from ase.utils import hsv
colors = hsv(atoms.positions[:, 0])

# Textures
tex = ['jmol',] * 432 + ['vmd',] * 432+ ['ase3',] * 288

# Make graphics files
write('flat.png', atoms, rotation=rotation, colors=colors, show_unit_cell=2,
      radii=None)
write('nice.pov', atoms, rotation=rotation, colors=colors, show_unit_cell=2,
      radii=None, textures=tex, display=False)

import os; os.system('povray nice.ini')

extra_keywords = { # For povray files only
'showunitcell' : 2,
'display'      : True,  # Display while rendering
'pause'        : True,  # Pause when done rendering (only if display)
'transparent'  : True,  # Transparent background
'canvas_width' : None,  # Width of canvas in pixels
'canvas_height': None,  # Height of canvas in pixels 
'camera_dist'  : 10.,   # Distance from camera to image plane
'camera_type'  : 'orthographic', # perspective, ultra_wide_angle
'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
'area_light'   : [(2., 3., 40.) ,# location
                  'White',       # color
                  .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
'background'   : 'White',        # color
'textures'     : None, # Length of atoms list of texture names
}
