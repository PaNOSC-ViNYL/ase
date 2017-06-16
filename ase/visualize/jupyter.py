"""Inline viewer for jupyter notebook."""

from tempfile import NamedTemporaryFile
from IPython.display import HTML


def jupview(atoms):
    """View atoms inline in a jupyter notbook. This command
    should only be used within a jupyter/ipython notebook."""
    with NamedTemporaryFile('r+', suffix='.html') as ntf:
        atoms.write(ntf.name, format='html')
        ntf.seek(0)
        html_atoms = ntf.read()
    return HTML(html_atoms)
