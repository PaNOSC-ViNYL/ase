try:
    import vtk
    hasvtk = True
except ImportError:
    hasvtk = False

def testvtk(code=0):
    from ase.test import NotAvailable
    if not hasvtk:
        # VTK required but not installed, force exit
        # exit status is determined by code argument
        raise NotAvailable('VTK is not installed.', code)

