# creates: spacegroup-al.png spacegroup-fe.png spacegroup-rutile.png spacegroup-cosb3.png spacegroup-mg.png spacegroup-skutterudite.png spacegroup-diamond.png spacegroup-nacl.png

import ase.io

for name in ['al', 'mg', 'fe', 'diamond', 'nacl', 'rutile', 'skutterudite']:
    py = 'spacegroup-{0}.py'.format(name)
    exec(compile(open(py).read(), py, 'exec'))
    atoms = globals()[name]
    ase.io.write('spacegroup-%s.pov' % name,
                 atoms,
                 transparent=False,
                 display=False,
                 run_povray=True,
                 # canvas_width=128,
                 show_unit_cell=2,
                 rotation='10x,-10y',
                 # celllinewidth=0.02,
                 celllinewidth=0.05)
    
exec(compile(open('spacegroup-cosb3.py').read(),
             'spacegroup-cosb3.py', 'exec'))
