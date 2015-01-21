# -*- coding: utf-8 -*-
# creates: surface.png
from ase.io import read, write
exec(compile(open('N2Cu.py').read(), 'N2Cu.py', 'exec'))
image = read('N2Cu.traj@-1')
write('surface.pov', image, transparent=False, display=False, run_povray=True)
