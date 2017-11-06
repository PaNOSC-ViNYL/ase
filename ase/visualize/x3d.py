"""Inline viewer for jupyter notebook using X3D."""

from io import BytesIO
from IPython.display import HTML

def view_x3d(atoms):
    """View atoms inline in a jupyter notbook. This command
    should only be used within a jupyter/ipython notebook.
    
    Args:
        atoms - ase.Atoms, atoms to be rendered"""
    
    with BytesIO() as output:
        atoms.write(output, format='html')
        return HTML(output.getvalue())
