import os
import pickle
import tempfile

from ase.io.xyz import write_xyz
from ase.io.cube import write_cube
from ase.io.plt import write_plt
import ase.parallel as parallel

try:
    from ase.old import OldASEListOfAtomsWrapper
except:
    oldase = False
else:
    oldase = True

def view(atoms, data=None, viewer=None, repeat=None):
    print repeat
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    if hasattr(atoms, 'GetUnitCell'):
        # Convert old ASE ListOfAtoms to new style.
        if oldase:
            atoms = OldASEListOfAtomsWrapper(atoms).copy()
        else:
            raise RuntimeError('conversion to old ASE not available')

    viewers = ['ase.gui', 'gopenmol', 'vmd', 'rasmol']
    if viewer is not None:
        viewer = viewer.lower()
        viewers.remove(viewer)
        viewers.insert(0, viewer)
    for viewer in viewers:
        try:
            if viewer == 'ase.gui':
                from ase.io.trajectory import write_trajectory
                filename = tempfile.mktemp('.traj', 'ag-')
                calc = atoms.get_calculator()
                atoms.set_calculator(None)
                write_trajectory(filename, atoms)
                atoms.set_calculator(calc)
                if repeat is None:
                    option = ''
                else:
                    option = '--repeat=%d,%d,%d ' % tuple(repeat)
                print option
                os.system('(ag %s%s &); (sleep 15; rm %s) &' %
                          (option, filename, filename))
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
