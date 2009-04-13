try:
    import vtk
    hasvtk = True
except ImportError:
    hasvtk = False

def requirevtk(code=0):
    from ase.test import NotAvailable
    if not hasvtk:
        # VTK required but not installed, force termination
        # with exit status determined by the code argument.
        raise NotAvailable('VTK is not installed.', code)

