import os
import tempfile
import subprocess

from ase.io import write
import ase.parallel as parallel


def view(atoms, data=None, viewer='ase', repeat=None, block=False):
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    vwr = viewer.lower()

    if vwr == 'ase':
        format = 'traj'
        if repeat is None:
            command = 'ase gui'
        else:
            command = 'ase gui --repeat=%d,%d,%d' % tuple(repeat)
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
    elif vwr == 'paraview':
        # macro for showing atoms in paraview
        macro = """\
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
source = GetActiveSource()
renderView1 = GetActiveViewOrCreate('RenderView')
atoms = Glyph(Input=source,
              GlyphType='Sphere',
              GlyphMode='All Points',
              Scalars='radii',
              ScaleMode='scalar',
              ScaleFactor=0.8)
RenameSource('Atoms', atoms)
atomsDisplay = Show(atoms, renderView1)
ColorBy(atomsDisplay, 'atomic numbers')
atomsDisplay.SetScalarBarVisibility(renderView1, True)
renderView1.ResetCamera()
        """
        script_name = os.path.join(tempfile.gettempdir(), 'draw_atoms.py')
        with open(script_name, 'w') as f:
            f.write(macro)
        format = 'vtu'
        command = 'paraview --script=' + script_name
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
