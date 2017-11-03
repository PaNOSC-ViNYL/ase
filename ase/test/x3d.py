from ase import Atoms
import os

try:
    from IPython.display import HTML
    from ase.visualize import x3d
except:
    pass
else:
    a = 3.6
    b = a / 2
    atoms = Atoms('Cu4',
                    positions=[(0, 0, 0),
                               (0, b, b),
                               (b, 0, b),
                               (b, b, 0)],
                    cell=(a, a, a),
                    pbc=True)
    my_obj, my_tempfile = x3d.view_x3d(atoms, return_path=True)
    assert isinstance(my_obj, HTML)
    assert not os.path.isfile(my_tempfile)