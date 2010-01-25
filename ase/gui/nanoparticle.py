# encoding: utf-8
"""nanoparticle.py - Window for setting up crystalline nanoparticles.
"""

import gtk
from ase.gui.widgets import pack, cancel_apply_ok, oops
from ase.gui.setupwindow import SetupWindow
from ase.gui.pybutton import PyButton
import ase
import numpy as np
# Delayed imports:
# ase.cluster.data

introtext = """\
Use this dialog to create crystalline nanoparticles.

1. Select the element, the desired crystal structure and the lattice
   constant.  The normal crystal structure and lattice constant can be
   looked up by pressing the button.

2. Specify the size of the particle by specifying the number of atomic
   layers in the different low-index crystal directions.  The number
   of layers is normally specified for a family of directions, but
   individual directions can be selected and specified for example for
   illustrating supported clusters.


*** NOT YET FUNCTIONAL! ***
"""

class SetupNanoparticle(SetupWindow):
    "Window for setting up a nanoparticle."
    families = {'fcc': [(0,0,1), (0,1,1), (1,1,1)]}
    
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title("Nanoparticle")
        self.atoms = None
        import ase.cluster.data
        self.data_module = ase.cluster.data

        vbox = gtk.VBox()

        # Intoductory text
        self.packtext(vbox, introtext)
           
        # Choose the element
        label = gtk.Label("Element: ")
        label.set_alignment(0.0, 0.2)
        element = gtk.Entry(max=3)
        self.element = element
        lattice_button = gtk.Button("Get structure")
        lattice_button.connect('clicked', self.get_structure)
        self.elementinfo = gtk.Label(" ")
        vbox2 = gtk.VBox()
        pack(vbox2, [element, gtk.Label("     "), lattice_button], end=True)
        pack(vbox2, [self.elementinfo])
        pack(vbox, [label, vbox2])
        self.element.connect('activate', self.update)
        self.legal_element = False

        # The structure
        label = gtk.Label("Structure: ")
        self.structure = gtk.combo_box_new_text()
        self.allowed_structures = ('fcc',)
        for struct in self.allowed_structures:
            self.structure.append_text(struct)
        self.structure.set_active(0)
        pack(vbox, [label, self.structure])
        
        # Lattice constant
        label = gtk.Label("Lattice constant: ")
        self.lattice_const = gtk.Adjustment(3.0, 0.0, 1000.0, 0.01)
        lattice_box = gtk.SpinButton(self.lattice_const, 10.0, 3)
        lattice_box.numeric = True
        pack(vbox, [label, lattice_box])
        self.lattice_const.connect('value-changed', self.update)
        pack(vbox, gtk.Label(""))

        # The number of layers
        pack(vbox, [gtk.Label("Number of layers:")])
        self.layerbox = gtk.VBox()
        pack(vbox, self.layerbox)
        self.make_layer_gui()
        
        # Finalize setup
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def update(self, *args):
        self.update_element()

    def get_structure(self, *args):
        if not self.update_element():
            oops("Invalid element.")
            return
        z = ase.atomic_numbers[self.legal_element]
        ref = ase.data.reference_states[z]
        if ref is None:
            structure = None
        else:
            structure = ref['symmetry'].lower()
                
        if ref is None or not structure in self.allowed_structures:
            oops("Unsupported or unknown structure",
                 "Element = %s,  structure = %s" % (self.legal_element,
                                                    structure))
            return
        for i, s in enumerate(self.allowed_structures):
            if structure == s:
                self.structure.set_active(i)
        a = ref['a']
        self.lattice_const.set_value(a)

    def make_layer_gui(self):
        "Make the part of the gui specifying the layers of the particle"
        # Clear the box
        children = self.layerbox.get_children()
        for c in children:
            self.layerbox.remove(c)
        del children

        # Get the crystal structure
        struct = self.structure.get_active_text()
        # Get the surfaces in the order the ase.cluster module expects
        surfaces = self.data_module.lattice[struct]['surface_names']
        # Get the surface families
        families = self.families[struct]

        # Empty array for the gtk.Adjustments for the layer numbers
        self.layers = [None] * len(surfaces)
        self.layer_lbl = [None] * len(surfaces)
        self.layer_spin = [None] * len(surfaces)
        self.famlayers = [None] * len(families)
        
        # Now, make a box for each family of surfaces
        frames = []
        for i, family in enumerate(families):
            frames.append(self.make_layer_family(i, family, surfaces))
        for a in self.layers:
            assert a is not None

        pack(self.layerbox, frames)
        self.layerbox.show_all()

    def make_layer_family(self, n, family, surfaces):
        """Make a frame box for a single family of surfaces.

        The layout is a frame containing a table.  For example

        {0,0,1}, SPIN, EMPTY, EMPTY
        -- empty line --
        (0,0,1), SPIN, Label(actual), Checkbox
        ...
        """
        tbl = gtk.Table(2, 4)
        lbl = gtk.Label("{%i,%i,%i}: " % family)
        lbl.set_alignment(1, 0.5)
        tbl.attach(lbl, 0, 1, 0, 1)
        famlayers = gtk.Adjustment(1, 1, 100, 1)
        tbl.attach(gtk.SpinButton(famlayers, 0, 0),
                   1, 2, 0, 1)
        tbl.attach(gtk.Label(" "), 0, 1, 1, 2)
        assert self.famlayers[n] is None
        self.famlayers[n] = famlayers
        row = 2
        myspin = []
        for i, s in enumerate(surfaces):
            s2 = [abs(x) for x in s]
            s2.sort()
            if tuple(s2) == family:
                tbl.resize(row+1, 4)
                lbl = gtk.Label("(%i,%i,%i): " % s)
                lbl.set_alignment(1, 0.5)
                tbl.attach(lbl, 0, 1, row, row+1)
                lay = gtk.Adjustment(1, 1, 100, 1)
                spin = gtk.SpinButton(lay, 0, 0)
                spin.set_sensitive(False)
                tbl.attach(spin, 1, 2, row, row+1)
                assert self.layers[i] is None
                self.layers[i] = lay
                self.layer_lbl[i] = lbl
                self.layer_spin[i] = spin
                myspin.append(spin)
                # Add actual number later
                # XXX
                chkbut = gtk.CheckButton()
                tbl.attach(chkbut, 3, 4, row, row+1)
                chkbut.connect("toggled", self.toggle_surface, i)
                row += 1
        famlayers.connect('value-changed', self.changed_family_layers, myspin)
        vbox = gtk.VBox()
        vbox.pack_start(tbl, False, False, 0)
        fr = gtk.Frame()
        fr.add(vbox)
        fr.show_all()
        return fr

    def toggle_surface(self, widget, number):
        "Toggle whether a layer in a family can be specified."
        self.layer_spin[number].set_sensitive(widget.get_active())
        
    def changed_family_layers(self, widget, myspin):
        "Change the number of layers in inactive members of a family."
        x = widget.value
        for s in myspin:
            if s.state == gtk.STATE_INSENSITIVE:
                adj = s.get_adjustment()
                if adj.value != x:
                    adj.value = x
                    
                
    
