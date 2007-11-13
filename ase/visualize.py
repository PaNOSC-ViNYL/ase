
def view(atoms, repeat=(1, 1, 1), viewer=None, filename=None):
    viewers = ['ase.gui', 'vmd', 'rasmol', 'nanolab']
    if viewer is not None:
        viewers = [viewer]
    for viewer in viewers:
        if viewer == 'ase.gui':
            from ase.gui import gui
            gui(atoms)



def g2():
    pass

#def vmd, rasmol, xmakemol
