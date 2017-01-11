import sys

from ase import Atoms
from ase.io import write


write('x.json', Atoms('X'))

# Make sure ase-gui can run in terminal mode without $DISPLAY and tkinter:
sys.argv = ['ase-gui', '--verbose', '--terminal', 'x.json@id=1']
from ase.gui.ag import main
main()
assert 'tkinter' not in sys.modules
assert 'Tkinter' not in sys.modules  # legacy Python
