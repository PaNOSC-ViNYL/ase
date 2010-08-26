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
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster.data.fcc import surface_names as fcc_surface_names
from ase.cluster import wulff_construction

introtext = """\
You specify the size of the particle by specifying the number of atomic layers
in the different low-index crystal directions.  Often, the number of layers is
specified for a family of directions, but they can be given individually.

When the particle is created, the actual numbers of layers are printed, they
may be less than specified if a surface is cut of by other surfaces."""

py_template = """
from ase.cluster import Cluster
import ase

layers = %(layers)s
atoms = Cluster(symbol='%(element)s', layers=layers, latticeconstant=%(a).5f,
                symmetry='%(structure)s')

# OPTIONAL: Cast to ase.Atoms object, discarding extra information:
# atoms = ase.Atoms(atoms)
"""

class attribute_collection:
    pass

class SetupNanoparticle(SetupWindow):
    "Window for setting up a nanoparticle."
    families = {'fcc': [(0,0,1), (0,1,1), (1,1,1)]}
    defaults = {'fcc': [6, 9, 5]}
    
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title("Nanoparticle")
        self.atoms = None
        #import ase.cluster.data
        #self.data_module = ase.cluster.data
        #import ase.cluster
        #self.Cluster = ase.cluster.Cluster
        #self.wulffconstruction = ase.cluster.wulff_construction
        self.no_update = True
        
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
        pack(vbox, [label, element, self.elementinfo, lattice_button], end=True)
        self.element.connect('activate', self.update)
        self.legal_element = False

        # The structure and lattice constant
        label = gtk.Label("Structure: ")
        self.structure = gtk.combo_box_new_text()
        self.allowed_structures = ('fcc',)
        for struct in self.allowed_structures:
            self.structure.append_text(struct)
        self.structure.set_active(0)
        self.structure.connect('changed', self.update)
        
        label2 = gtk.Label("   Lattice constant: ")
        self.lattice_const = gtk.Adjustment(3.0, 0.0, 1000.0, 0.01)
        lattice_box = gtk.SpinButton(self.lattice_const, 10.0, 3)
        lattice_box.numeric = True
        pack(vbox, [label, self.structure, label2, lattice_box])
        self.lattice_const.connect('value-changed', self.update)

        # Choose specification method
        label = gtk.Label("Method: ")
        self.method = gtk.combo_box_new_text()
        for meth in ("Layer specification", "Wulff construction"):
            self.method.append_text(meth)
        self.method.set_active(0)
        self.method.connect('changed', self.update_gui)
        pack(vbox, [label, self.method])
        pack(vbox, gtk.Label(""))
        
        # The number of layers
        self.layerbox = gtk.VBox()
        self.layerdata = attribute_collection()
        pack(vbox, self.layerbox)
        self.make_layer_gui(self.layerbox, self.layerdata, 0)

        # The Wulff construction
        self.wulffbox = gtk.VBox()
        self.wulffdata = attribute_collection()
        pack(vbox, self.wulffbox)
        self.make_layer_gui(self.wulffbox, self.wulffdata, 1)
        label = gtk.Label("Particle size: ")
        self.size_n_radio = gtk.RadioButton(None, "Number of atoms: ")
        self.size_n_radio.set_active(True)
        self.size_n_adj = gtk.Adjustment(100, 1, 100000, 1)
        self.size_n_spin = gtk.SpinButton(self.size_n_adj, 0, 0)
        self.size_dia_radio = gtk.RadioButton(self.size_n_radio,
                                              "Volume: ")
        self.size_dia_adj = gtk.Adjustment(1.0, 0, 100.0, 0.1)
        self.size_dia_spin = gtk.SpinButton(self.size_dia_adj, 10.0, 2)
        pack(self.wulffbox, [label, self.size_n_radio, self.size_n_spin,
                    gtk.Label("   "), self.size_dia_radio, self.size_dia_spin,
                    gtk.Label("Å³")])
        self.size_n_radio.connect("toggled", self.update_gui)
        self.size_dia_radio.connect("toggled", self.update_gui)
        self.size_n_adj.connect("value-changed", self.update_size_n)
        self.size_dia_adj.connect("value-changed", self.update_size_dia)
        self.update_size_n()
        label = gtk.Label(
            "Rounding: If exact size is not possible, choose the size")
        pack(self.wulffbox, [label])
        self.round_above = gtk.RadioButton(None, "above  ")
        self.round_below = gtk.RadioButton(self.round_above, "below  ")
        self.round_closest = gtk.RadioButton(self.round_above, "closest  ")
        self.round_closest.set_active(True)
        buts = [self.round_above, self.round_below, self.round_closest]
        pack(self.wulffbox, buts)
        for b in buts:
            b.connect("toggled", self.update)
        
        # Information
        pack(vbox, gtk.Label(""))
        infobox = gtk.VBox()
        label1 = gtk.Label("Number of atoms: ")
        self.natoms_label = gtk.Label("-")
        label2 = gtk.Label("   Approx. diameter: ")
        self.dia1_label = gtk.Label("-")
        pack(infobox, [label1, self.natoms_label, label2, self.dia1_label])
        pack(infobox, gtk.Label(""))
        infoframe = gtk.Frame("Information about the created cluster:")
        infoframe.add(infobox)
        infobox.show()
        pack(vbox, infoframe)
        
        # Buttons
        self.pybut = PyButton("Creating a nanoparticle.")
        self.pybut.connect('clicked', self.makeatoms)
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [self.pybut, buts], end=True, bottom=True)
        self.auto = gtk.CheckButton("Automatic Apply")
        fr = gtk.Frame()
        fr.add(self.auto)
        fr.show_all()
        pack(vbox, [fr], end=True, bottom=True)
        
        # Finalize setup
        self.update_gui()
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui
        self.no_update = False

    def update_gui(self, widget=None):
        method = self.method.get_active()
        if method == 0:
            self.wulffbox.hide()
            self.layerbox.show()
        elif method == 1:
            self.layerbox.hide()
            self.wulffbox.show()
            self.size_n_spin.set_sensitive(self.size_n_radio.get_active())
            self.size_dia_spin.set_sensitive(self.size_dia_radio.get_active())

    def update_size_n(self, widget=None):
        if not self.size_n_radio.get_active():
            return
        at_vol = self.get_atomic_volume()
        dia = 2.0 * (3 * self.size_n_adj.value * at_vol / (4 * np.pi))**(1.0/3)
        self.size_dia_adj.value = dia
        self.update()

    def update_size_dia(self, widget=None):
        if not self.size_dia_radio.get_active():
            return
        at_vol = self.get_atomic_volume()
        n = round(np.pi / 6 * self.size_dia_adj.value**3 / at_vol)
        self.size_n_adj.value = n
        self.update()
                
    def update(self, *args):
        if self.no_update:
            return
        self.update_gui()
        self.update_element()
        if self.auto.get_active():
            self.makeatoms()
            if self.atoms is not None:
                self.gui.new_atoms(self.atoms)
        else:
            self.clearatoms()
        self.makeinfo()

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

    def make_layer_gui(self, box, data, method):
        "Make the part of the gui specifying the layers of the particle"
        assert method in (0,1)
        
        # Clear the box
        children = box.get_children()
        for c in children:
            box.remove(c)
        del children

        # Make the label
        if method == 0:
            txt = "Number of layers:"
        else:
            txt = "Surface energies (unit: energy/area, i.e. J/m<sup>2</sup> or eV/nm<sup>2</sup>, <i>not</i> eV/atom):"
        label = gtk.Label()
        label.set_markup(txt)
        pack(box, [label])

        # Get the crystal structure
        struct = self.structure.get_active_text()
        # Get the surfaces in the order the ase.cluster module expects
        surfaces = fcc_surface_names
        # Get the surface families
        families = self.families[struct]
        if method == 0:
            defaults = self.defaults[struct]
        else:
            defaults = [1.0] * len(self.defaults[struct])
        
        # Empty array for the gtk.Adjustments for the layer numbers
        data.layers = [None] * len(surfaces)
        data.layer_lbl = [None] * len(surfaces)
        data.layer_spin = [None] * len(surfaces)
        data.layer_owner = [None] * len(surfaces)
        data.layer_label = [None] * len(surfaces)
        data.famlayers = [None] * len(families)
        data.infamily = [None] * len(families)
        data.family_label = [None] * len(families)
        
        # Now, make a box for each family of surfaces
        frames = []
        for i in range(len(families)):
            family = families[i]
            default = defaults[i]
            frames.append(self.make_layer_family(data, i, family, surfaces,
                                                 default, method))
        for a in data.layers:
            assert a is not None

        pack(box, frames)
        box.show_all()

    def make_layer_family(self, data, n, family, surfaces, default, method):
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
        if method == 0:
            famlayers = gtk.Adjustment(default, 1, 100, 1)
            sbut = gtk.SpinButton(famlayers, 0, 0)
        else:
            flimit = 1000.0
            famlayers = gtk.Adjustment(default, 0.0, flimit, 0.1)
            sbut = gtk.SpinButton(famlayers, 10.0, 3)
        tbl.attach(sbut, 2, 3, 0, 1)
        tbl.attach(gtk.Label(" "), 0, 1, 1, 2)
        assert data.famlayers[n] is None
        data.famlayers[n] = famlayers
        data.infamily[n] = []
        data.family_label[n] = gtk.Label("")
        tbl.attach(data.family_label[n], 1, 2, 0, 1)
        row = 2
        myspin = []
        for i, s in enumerate(surfaces):
            s2 = [abs(x) for x in s]
            s2.sort()
            if tuple(s2) == family:
                data.infamily[n].append(i)
                tbl.resize(row+1, 4)
                lbl = gtk.Label("(%i,%i,%i): " % s)
                lbl.set_alignment(1, 0.5)
                tbl.attach(lbl, 0, 1, row, row+1)
                label = gtk.Label("    ")
                tbl.attach(label, 1, 2, row, row+1)
                data.layer_label[i] = label
                if method == 0:
                    lay = gtk.Adjustment(default, -100, 100, 1)
                    spin = gtk.SpinButton(lay, 0, 0)
                else:
                    lay = gtk.Adjustment(default, -flimit, flimit, 0.1)
                    spin = gtk.SpinButton(lay, 10.0, 3)
                lay.connect('value-changed', self.update)
                spin.set_sensitive(False)
                tbl.attach(spin, 2, 3, row, row+1)
                assert data.layers[i] is None
                data.layers[i] = lay
                data.layer_lbl[i] = lbl
                data.layer_spin[i] = spin
                data.layer_owner[i] = n
                myspin.append(spin)
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
        active = widget.get_active()
        data = self.get_data()
        data.layer_spin[number].set_sensitive(active)
        if not active:
            data.layers[number].value = \
                data.famlayers[data.layer_owner[number]].value
        
    def changed_family_layers(self, widget, myspin):
        "Change the number of layers in inactive members of a family."
        self.no_update = True
        x = widget.value
        for s in myspin:
            if s.state == gtk.STATE_INSENSITIVE:
                adj = s.get_adjustment()
                if adj.value != x:
                    adj.value = x
        self.no_update = False
        self.update()

    def get_data(self):
        return self.layerdata
    
    def makeatoms(self, *args):
        "Make the atoms according to the current specification."
        if not self.update_element():
            self.clearatoms()
            self.makeinfo()
            return False
        assert self.legal_element is not None
        struct = self.structure.get_active_text()
        lc = self.lattice_const.value
        if self.method.get_active() == 0:
            # Layer-by-layer specification
            layers = [int(x.value) for x in self.layerdata.layers]
            self.atoms = FaceCenteredCubic(self.legal_element, fcc_surface_names,
                                           layers=layers, latticeconstant=lc)
            self.pybut.python = py_template % {'element': self.legal_element,
                                               'layers': str(layers),
                                               'structure': struct,
                                               'a': lc}
        else:
            # Wulff construction
            assert self.method.get_active() == 1
            surfaceenergies = [x.value for x in self.wulffdata.layers]
            self.update_size_dia()
            if self.round_above.get_active():
                rounding = "above"
            elif self.round_below.get_active():
                rounding = "below"
            elif self.round_closest.get_active():
                rounding = "closest"
            else:
                raise RuntimeError("No rounding!")
            self.atoms = wulffconstruction(self.legal_element,
                                           surfaceenergies,
                                           self.size_n_adj.value,
                                           rounding, struct, lc)
                                   
        self.makeinfo()

    def clearatoms(self):
        self.atoms = None
        self.pybut.python = None

    def get_atomic_volume(self):
        s = self.structure.get_active_text()
        a = self.lattice_const.value
        if s == 'fcc':
            return a**3 / 4
        else:
            raise RuntimeError("Unknown structure: "+s)

    def makeinfo(self):
        "Fill in information field about the atoms."
        #data = self.get_data()
        if self.atoms is None:
            self.natoms_label.set_label("-")
            self.dia1_label.set_label("-")
            for d in (self.layerdata, self.wulffdata):
                for label in d.layer_label+d.family_label:
                    label.set_text("    ")
        else:
            self.natoms_label.set_label(str(len(self.atoms)))
            at_vol = self.get_atomic_volume()
            dia = 2 * (3 * len(self.atoms) * at_vol / (4 * np.pi))**(1.0/3.0)
            self.dia1_label.set_label("%.1f Å" % (dia,))
            actual = self.atoms.get_layers()
            for i, a in enumerate(actual):
                self.layerdata.layer_label[i].set_text("%2i " % (a,))
                self.wulffdata.layer_label[i].set_text("%2i " % (a,))
            for d in (self.layerdata, self.wulffdata):
                for i, label in enumerate(d.family_label):
                    relevant = actual[d.infamily[i]]
                    if relevant.min() == relevant.max():
                        label.set_text("%2i " % (relevant[0]))
                    else:
                        label.set_text("-- ")
            
    def apply(self, *args):
        self.makeatoms()
        if self.atoms is not None:
            self.gui.new_atoms(self.atoms)
            return True
        else:
            oops("No valid atoms.",
                 "You have not (yet) specified a consistent set of parameters.")
            return False

    def ok(self, *args):
        if self.apply():
            self.destroy()
            
