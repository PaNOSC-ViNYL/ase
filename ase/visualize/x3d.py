"""Inline viewer for jupyter notebook using X3D."""

from tempfile import mkstemp
from IPython.display import HTML
import os

def view_x3d(atoms, return_path=False):
    """View atoms inline in a jupyter notbook. This command
    should only be used within a jupyter/ipython notebook.
    
    Args:
        return_path - bool, whether to return the path to the tempfile 
            created during this operation [debug purposes]"""
    
    ntf_fp, ntf = mkstemp(suffix=".html")
    try:
        atoms.write(ntf, format='html')
        with open(ntf, 'r') as fp:
            html_atoms = fp.read()
    finally:
        os.close(ntf_fp)
        os.remove(ntf)
        
    if return_path:
        # Debug purposes
        return HTML(html_atoms), ntf
    return HTML(html_atoms)
