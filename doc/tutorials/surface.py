# -*- coding: utf-8 -*-
# creates:  surface.png
import os
from ase import *
from ase.io.pov import povpng
execfile('N2Cu.py')
image = read('N2Cu.traj@-1')
povpng('surface.pov', image)
