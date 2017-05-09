# coding: utf-8

from ase import Atoms
import nglview 

def view_ngl( atoms ):
    if isinstance(atoms[0],Atoms):
        # Assume this is a trajectory or struct list
        view = nglview.show_asetraj( atoms )
    else :
        view = nglview.show_ase( atoms )
    view.add_unitcell()
    view.add_spacefill(scale=0.6)
    view.center()
    return view
    

