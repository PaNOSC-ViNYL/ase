"Base class for simulation windows"

import gtk
from ase.gui.widgets import oops, pack
from ase import Atoms

class Simulation(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui

    def packtext(self, vbox, text, label=None):
        "Pack an text frame into the window."
        pack(vbox, gtk.Label(""))
        txtframe = gtk.Frame(label)
        txtlbl = gtk.Label(text)
        txtframe.add(txtlbl)
        txtlbl.show()
        pack(vbox, txtframe)
        pack(vbox, gtk.Label(""))

    def makebutbox(self, vbox):
        self.buttons = gtk.HButtonBox()
        runbut = gtk.Button("Run")
        runbut.connect('clicked', self.run)
        closebut = gtk.Button(stock=gtk.STOCK_CLOSE)
        closebut.connect('clicked', lambda x: self.destroy())
        for w in (runbut, closebut):
            self.buttons.pack_start(w, 0, 0)
            w.show()
        pack(vbox, [self.buttons], end=True, bottom=True)

    def setup_atoms(self):
        self.atoms = self.get_atoms()
        if self.atoms is None:
            return False
        try:
            self.calculator = self.gui.simulation['calc']
        except KeyError:
            oops("No calculator: Use Calculate/Set Calculator on the menu.")
            return False
        self.atoms.set_calculator(self.calculator())
        return True
    
    def get_atoms(self):
        "Make an atoms object from the active image"
        images = self.gui.images
        if images.natoms < 1:
            oops("No atoms present")
            return None
        return Atoms(positions=images.P[0], symbols=images.Z,
                           cell=images.A[0])


