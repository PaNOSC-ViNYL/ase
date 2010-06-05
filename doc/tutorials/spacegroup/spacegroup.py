# creates: spacegroup-*.png 

import ase

execfile('spacegroup-al.py')
execfile('spacegroup-mg.py')
execfile('spacegroup-fe.py')
execfile('spacegroup-diamond.py')
execfile('spacegroup-nacl.py')
execfile('spacegroup-rutile.py')
execfile('spacegroup-skutterudite.py')

for name in ['al', 'mg', 'fe', 'diamond', 'nacl', 'rutile', 'skutterudite']:
    atoms = globals()[name]
    ase.write('spacegroup-%s.pov'%name, 
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

execfile('spacegroup-cosb3.py')
    
