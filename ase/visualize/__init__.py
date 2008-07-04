import os
import tempfile

from ase.io import write
import ase.parallel as parallel

try:
    from ase.old import OldASEListOfAtomsWrapper
except:
    oldase = False
else:
    oldase = True

def view(atoms, data=None, viewer='ag', repeat=None, block=False):
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    if hasattr(atoms, 'GetUnitCell'):
        # Convert old ASE ListOfAtoms to new style.
        if oldase:
            atoms = OldASEListOfAtomsWrapper(atoms).copy()
        else:
            raise RuntimeError('conversion to old ASE not available')

    if viewer == 'ag':
        format = 'traj'
        if repeat is None:
            command = 'ag'
        else:
            command = 'ag --repeat=%d,%d,%d' % tuple(repeat)
            repeat = None
    elif viewer == 'vmd':
        format = 'cube'
        command = 'vmd'
    elif viewer == 'rasmol':
        format = 'pdb'
        command = 'rasmol -pdb'
    elif viewer == 'xmakemol':
        format = 'xyz'
        command = 'xmakemol -f'
    elif viewer == 'gopenmol':
        format = 'xyz'
        command = 'rungOpenMol'
    else:
        raise RuntimeError('Unknown viewer: ' + viewer)

    fd, filename = tempfile.mkstemp('.' + format, 'ase-')
    fd = os.fdopen(fd, 'w')
    if repeat is not None:
        atoms = atoms.repeat()
    write(fd, atoms, format=format, data=data)
    fd.close()
    if block:
        os.system('%s %s' % (command, filename))
    else:
        os.system('%s %s &' % (command, filename))
    os.system('(sleep 60; rm %s) &' % filename)
