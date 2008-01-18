import os
import tempfile

from ase.io.xyz import write_xyz
from ase.io.cube import write_cube
from ase.io.plt import write_plt
import ase.parallel as parallel


def view(atoms, data=None, viewer=None):
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    viewers = ['ase.gui', 'gopenmol', 'vmd', 'rasmol', 'nanolab']
    if viewer is not None:
        viewer = viewer.lower()
        viewers.remove(viewer)
        viewers.insert(0, viewer)
    for viewer in viewers:
        try:
            if viewer == 'ase.gui':
                from ase.gui import gui
                gui(atoms)
                break
            if viewer == 'gopenmol':
                fd, filename = tempfile.mkstemp('.xyz', 'ag-')
                os.close(fd)
                write_xyz(filename, atoms)
                if data is not None:
                    write_plt('data.plt', atoms, data)
                os.system('(rungOpenMol %s &); (sleep 5; rm %s) &' %
                          (filename, filename))
                break
            if viewer == 'vmd':
                fd, filename = tempfile.mkstemp('.cube', 'ag-')
                os.close(fd)
                write_cube(filename, atoms, data)
                os.system('(vmd %s &); (sleep 5; rm %s) &' %
                          (filename, filename))
                break
        except:
            pass


#def g2():
#    pass

#def vmd, rasmol, xmakemol
