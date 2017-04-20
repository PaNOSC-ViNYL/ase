import os
import subprocess
import sys
import tempfile

from ase.io import write
import ase.parallel as parallel


def view(atoms, data=None, viewer='ase', repeat=None, block=False):
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    vwr = viewer.lower()

    if vwr == 'ase':
        format = 'traj'
        command = sys.executable + ' -m ase gui'
        if repeat is not None:
            command += ' --repeat={},{},{}'.format(*repeat)
            repeat = None
    elif vwr == 'vmd':
        format = 'cube'
        command = 'vmd'
    elif vwr == 'rasmol':
        format = 'pdb'
        command = 'rasmol -pdb'
    elif vwr == 'xmakemol':
        format = 'xyz'
        command = 'xmakemol -f'
    elif vwr == 'gopenmol':
        format = 'xyz'
        command = 'rungOpenMol'
    elif vwr == 'avogadro':
        format = 'cube'
        command = 'avogadro'
    elif vwr == 'sage':
        from ase.visualize.sage import view_sage_jmol
        view_sage_jmol(atoms)
        return
    else:
        raise RuntimeError('Unknown viewer: ' + viewer)

    fd, filename = tempfile.mkstemp('.' + format, 'ase-')
    if repeat is not None:
        atoms = atoms.repeat()
    if data is None:
        write(filename, atoms, format=format)
    else:
        write(filename, atoms, format=format, data=data)
    if block:
        subprocess.call(command.split() + [filename])
        os.remove(filename)
    else:
        subprocess.Popen(command.split() + [filename])
        subprocess.Popen(['sleep 60; rm {0}'.format(filename)], shell=True)
