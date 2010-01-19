# encoding: utf-8
"""surfaceslab.py - Window for setting up surfaces

This is part of the Visual ASE for teaching (vase) GUI.
"""

import gtk
from ase.gui.widgets import pack, cancel_apply_ok

introtext = """\
  Use this dialog to create surface slabs.  Select the element by
writing the chemical symbol or the atomic number in the box.  Then
select the desired surface structure.  Note that some structures can
be created with an othogonal or a non-orthogonal unit cell, in these
cases the non-orthogonal unit cell will contain fewer atoms.

  If the structure matches the experimental crystal structure, you can
look up the lattice constant, otherwise you have to specify it
yourself."""

import ase.lattice.surface as _surf
import ase
import numpy as np

# Name, structure, orthogonal, support-nonorthogonal, function
surfaces = [('FCC(100)', 'fcc', True, False, _surf.fcc100),
            ('FCC(110)', 'fcc', True, False, _surf.fcc110),
            ('FCC(111) non-orthogonal', 'fcc', False, True, _surf.fcc111),
            ('FCC(111) orthogonal', 'fcc', True, True, _surf.fcc111),
            ('BCC(100)', 'bcc', True, False, _surf.bcc100),
            ('BCC(110) non-orthogonal', 'bcc', False, True, _surf.bcc110),
            ('BCC(110) orthogonal', 'bcc', True, True, _surf.bcc110),
            ('BCC(111) non-orthogonal', 'bcc', False, True, _surf.bcc111),
            ('BCC(111) orthogonal', 'bcc', True, True, _surf.bcc111),
            ]

class SetupSurfaceSlab(gtk.Window):
    """Window for setting up a surface."""
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title("Surface")
        self.atoms = None

        vbox = gtk.VBox()

        # Intoductory text
        pack(vbox, gtk.Label(""))
        txtframe = gtk.Frame()
        txtlbl = gtk.Label(introtext)
        txtframe.add(txtlbl)
        txtlbl.show()
        pack(vbox, txtframe)
        pack(vbox, gtk.Label(""))
             
        # Choose the element
        label = gtk.Label("Element: ")
        element = gtk.Entry(max=3)
        self.element = element
        self.elementinfo = gtk.Label("")
        pack(vbox, [label, element, self.elementinfo])
        self.element.connect('activate', self.update)
        self.legal_element = False
        
        # Choose the surface structure
        label = gtk.Label("Structure: ")
        self.structchoice = gtk.combo_box_new_text()
        self.surfinfo = {}
        for s in surfaces:
            assert len(s) == 5
            self.structchoice.append_text(s[0])
            self.surfinfo[s[0]] = s
        pack(vbox, [label, self.structchoice])
        self.structchoice.connect('changed', self.update)

        # Choose the lattice constant
        label = gtk.Label("Lattice constant: ")
        self.lattice_const = gtk.Adjustment(3.0, 0.0, 1000.0, 0.01)
        lattice_box = gtk.SpinButton(self.lattice_const, 10.0, 3)
        lattice_box.numeric = True
        lattice_button = gtk.Button("Get from database")
        pack(vbox, [label, lattice_box, lattice_button], end=True)
        self.lattice_const.connect('value-changed', self.update)
        lattice_button.connect('clicked', self.get_lattice_const)
        pack(vbox, gtk.Label(""))

        # System size
        self.size = [gtk.Adjustment(1, 1, 100, 1) for i in range(3)]
        buttons = [gtk.SpinButton(s, 0, 0) for s in self.size]
        self.vacuum = gtk.Adjustment(10.0, 0, 100.0, 0.1)
        vacuum_box = gtk.SpinButton(self.vacuum, 0.0, 1)
        pack(vbox, [gtk.Label("Size: \tx: "), buttons[0],
                    gtk.Label(" unit cells")])
        pack(vbox, [gtk.Label("\t\ty: "), buttons[1],
                    gtk.Label(" unit cells")])
        pack(vbox, [gtk.Label("      \t\tz: "), buttons[2],
                    gtk.Label(" layers,  "),
                    vacuum_box, gtk.Label(" Å vacuum")])
        self.nosize = "\t\tNo size information yet."
        self.sizelabel = gtk.Label(self.nosize)
        pack(vbox, [self.sizelabel])
        for s in self.size:
            s.connect('value-changed', self.update)
        self.vacuum.connect('value-changed', self.update)
        pack(vbox, gtk.Label(""))

        # Buttons
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True)
        
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def update_element(self, *args):
        "Called when a new element may have been entered."
        elem = self.element.get_text()
        if not elem:
            self.invalid_element("  No element specified!")
            return False
        try:
            z = int(elem)
        except ValueError:
            # Probably a symbol
            try:
                z = ase.data.atomic_numbers[elem]
            except KeyError:
                self.invalid_element()
                return False
        try:
            symb = ase.data.chemical_symbols[z]
        except KeyError:
            self.invalid_element()
            return False
        name = ase.data.atomic_names[z]
        ref = ase.data.reference_states[z]
        if ref is None:
            struct = "No crystal structure data"
        else:
            struct = ref['symmetry'].lower()
            if struct == 'fcc' or struct == 'bcc':
                struct = "%s (a=%.3f Å)" % (struct, ref['a'])
        
        txt = "  %s: %s, Z=%i, %s" % (name, symb, z, struct)
        self.elementinfo.set_text(txt)
        self.legal_element = symb
        return True

    def update(self, *args):
        "Called when something has changed."
        struct = self.structchoice.get_active_text()
        if not (self.update_element() and struct):
            self.sizelabel.set_text(self.nosize)
            self.atoms = None
            return False
        # Make the atoms
        assert self.legal_element
        structinfo = self.surfinfo[struct]
        xtra_keywords = {}
        if structinfo[3]:  # Support othogonal keyword?
            xtra_keywords['orthogonal'] = structinfo[2]
        self.atoms = structinfo[4](symbol=self.legal_element,
                                   size=[int(s.value) for s in self.size],
                                   a=self.lattice_const.value,
                                   vacuum=self.vacuum.value,
                                   **xtra_keywords)
        # Find the heights of the unit cell
        h = np.zeros(3)
        uc = self.atoms.get_cell()
        for i in range(3):
            norm = np.cross(uc[i-1], uc[i-2])
            norm /= np.sqrt(np.dot(norm, norm))
            h[i] = np.abs(np.dot(norm, uc[i]))
        natoms = len(self.atoms)
        txt = ("\t\t%.2f Å x %.2f Å x %.2f Å,  %i atoms."
               % (h[0], h[1], h[2], natoms))
        self.sizelabel.set_text(txt)
        return True
    
    def invalid_element(self, txt="  ERROR: Invalid element!"):
        self.legal_element = False
        self.elementinfo.set_text(txt)

    def get_lattice_const(self, *args):
        if not self.update_element():
            self.oops("Invalid element.")
            return
        z = ase.atomic_numbers[self.legal_element]
        ref = ase.data.reference_states[z]
        surface = self.structchoice.get_active_text()
        if not surface:
            self.oops("No structure specified!")
            return
        struct = self.surfinfo[surface][1]
        if ref is None or ref['symmetry'].lower() != struct:
            self.oops(struct.upper() + " lattice constant unknown for "
                      + self.legal_element + ".")
            return
        a = ref['a']
        self.lattice_const.set_value(a)

    def apply(self, *args):
        if self.atoms is not None:
            self.gui.new_atoms(self.atoms)
            return True
        else:
            self.oops("No valid atoms.")
            return False

    def ok(self, *args):
        if self.apply():
            self.destroy()
            
    def oops(self, message):
        dialog = gtk.MessageDialog(flags=gtk.DIALOG_MODAL,
                                   type=gtk.MESSAGE_WARNING,
                                   buttons=gtk.BUTTONS_CLOSE,
                                   message_format=message)
        dialog.connect('response', lambda x, y: dialog.destroy())
        dialog.show()
        
