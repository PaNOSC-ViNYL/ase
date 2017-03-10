# -*- coding: utf-8 -*-
# creates: diffusion-I.png, diffusion-T.png, diffusion-F.png
# creates: diffusion-barrier.png

from ase.io import read, write
from ase.neb import NEBTools

if 1:
    exec(compile(open('diffusion1.py').read(), 'diffusion1.py', 'exec'))
    exec(compile(open('diffusion2.py').read(), 'diffusion2.py', 'exec'))
    exec(compile(open('diffusion4.py').read(), 'diffusion4.py', 'exec'))
    exec(compile(open('diffusion5.py').read(), 'diffusion5.py', 'exec'))

images = read('neb.traj@-5:')
for name, a in zip('ITF', images[::2]):
    cell = a.get_cell()
    del a.constraints
    a = a * (2, 2, 1)
    a.set_cell(cell)
    write('diffusion-%s.pov' % name, a, show_unit_cell=True,
          transparent=False, display=False, run_povray=True)

nebtools = NEBTools(images)
assert abs(nebtools.get_barrier()[0] - 0.374) < 1e-3
