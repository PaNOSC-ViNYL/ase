# creates: spacegroup-al.png spacegroup-fe.png spacegroup-rutile.png spacegroup-cosb3.png spacegroup-mg.png spacegroup-skutterudite.png spacegroup-diamond.png spacegroup-nacl.png

import ase.io

exec(compile(open('spacegroup-al.py').read(), 'spacegroup-al.py', 'exec'))
exec(compile(open('spacegroup-mg.py').read(), 'spacegroup-mg.py', 'exec'))
exec(compile(open('spacegroup-fe.py').read(), 'spacegroup-fe.py', 'exec'))
exec(compile(open('spacegroup-diamond.py').read(), 'spacegroup-diamond.py', 'exec'))
exec(compile(open('spacegroup-nacl.py').read(), 'spacegroup-nacl.py', 'exec'))
exec(compile(open('spacegroup-rutile.py').read(), 'spacegroup-rutile.py', 'exec'))
exec(compile(open('spacegroup-skutterudite.py').read(), 'spacegroup-skutterudite.py', 'exec'))

for name in ['al', 'mg', 'fe', 'diamond', 'nacl', 'rutile', 'skutterudite']:
    atoms = globals()[name]
    ase.io.write('spacegroup-%s.pov'%name, 
                 atoms, 
                 transparent=False, 
                 display=False, 
                 run_povray=True,
                 #canvas_width=128,
                 show_unit_cell=2,
                 rotation='10x,-10y', 
                 #celllinewidth=0.02,
                 celllinewidth=0.05,
                 )

exec(compile(open('spacegroup-cosb3.py').read(), 'spacegroup-cosb3.py', 'exec'))
    
