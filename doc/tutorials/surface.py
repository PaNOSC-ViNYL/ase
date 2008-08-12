# -*- coding: utf-8 -*-
# creates:  surface.png
import os
from ase import *
execfile('N2Cu.py')
image = read('N2Cu.traj@-1')
write('surface.pov', image, pause=False)
os.system('povray surface.ini')
