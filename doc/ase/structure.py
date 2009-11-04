# creates: a1.png a2.png a3.png
from ase.io import write
from ase.structure import bulk
for i, a in enumerate([
    bulk('Cu', 'fcc', a=3.6),
    bulk('Cu', 'fcc', a=3.6, orthorhombic=True),
    bulk('Cu', 'fcc', a=3.6, cubic=True)]):
    write('a%d.pov' % (i + 1), a,
          show_unit_cell=2, display=False, run_povray=True)
