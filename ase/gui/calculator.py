# encoding: utf-8
"calculator.py - module for choosing a calculator."

import gtk
import numpy as np
from copy import copy
from ase.gui.setupwindow import SetupWindow
from ase.gui.widgets import pack, oops, cancel_apply_ok
from ase import Atoms
from ase.data import chemical_symbols

# Asap and GPAW may be imported if selected.

introtext = """\
To make most calculations on the atoms, a Calculator object must first
be associated with it.  ASE supports a number of calculators, supporting
different elements, and implementing different physical models for the
interatomic interactions.\
"""

# Informational text about the calculators
lj_info_txt = """\
The Lennard-Jones pair potential is one of the simplest
possible models for interatomic interactions, mostly
suitable for noble gasses and model systems.

Interactions are described by an interaction length and an
interaction strength.\
"""

emt_info_txt = """\

The EMT potential is a many-body potential, giving a
good description of the late transition metals crystalling
in the FCC crystal structure.  The elements described by the
main set of EMT parameters are Al, Ni, Cu, Pd, Ag, Pt, and
Au, the Al potential is however not suitable for materials
science application, as the stacking fault energy is wrong.

A number of parameter sets are provided.

<b>Default parameters:</b>

The default EMT parameters, as published in K. W. Jacobsen,
P. Stoltze and J. K. Nørskov, <i>Surf. Sci.</i> <b>366</b>, 394 (1996).

<b>Alternative Cu, Ag and Au:</b>

An alternative set of parameters for Cu, Ag and Au,
reoptimized to experimental data including the stacking
fault energies by Torben Rasmussen (partly unpublished).

<b>Ruthenium:</b>

Parameters for Ruthenium, as published in J. Gavnholt and
J. Schiøtz, <i>Phys. Rev. B</i> <b>77</b>, 035404 (2008).

<b>Metallic glasses:</b>

Parameters for MgCu and CuZr metallic glasses. MgCu
parameters are in N. P. Bailey, J. Schiøtz and
K. W. Jacobsen, <i>Phys. Rev. B</i> <b>69</b>, 144205 (2004).
CuZr in A. Paduraru, A. Kenoufi, N. P. Bailey and
J. Schiøtz, <i>Adv. Eng. Mater.</i> <b>9</b>, 505 (2007).
"""

emt_parameters = (
    ("Default (Al, Ni, Cu, Pd, Ag, Pt, Au)", None),
    ("Alternative Cu, Ag and Au", "EMTRasmussenParameters"),
    ("Ruthenium", "EMThcpParameters"),
    ("CuMg and CuZr metallic glass", "EMTMetalGlassParameters")
    )

class SetCalculator(SetupWindow):
    "Window for selecting a calculator."
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title("Select calculator")
        vbox = gtk.VBox()
        
        # Intoductory text
        self.packtext(vbox, introtext)
        
        pack(vbox, [gtk.Label("Calculator:")])
        
        # No calculator (the default)
        self.none_radio = gtk.RadioButton(None, "None")
        pack(vbox, [self.none_radio])

        # Lennard-Jones
        self.lj_radio = gtk.RadioButton(self.none_radio, "Lennard-Jones (ASAP)")
        self.lj_setup = gtk.Button("Setup")
        self.lj_info = InfoButton(lj_info_txt)
        self.lj_setup.connect("clicked", self.lj_setup_window)
        self.pack_line(vbox, self.lj_radio, self.lj_setup, self.lj_info)

        # EMT
        self.emt_radio = gtk.RadioButton(
            self.none_radio, "EMT - Effective Medium Theory (ASAP)")
        self.emt_setup = gtk.combo_box_new_text()
        self.emt_param_info = {}
        for p in emt_parameters:
            self.emt_setup.append_text(p[0])
            self.emt_param_info[p[0]] = p[1]
        self.emt_setup.set_active(0)
        self.emt_info = InfoButton(emt_info_txt)
        self.pack_line(vbox, self.emt_radio, self.emt_setup, self.emt_info)

        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True, bottom=True)
        self.check = gtk.CheckButton("Check that the calculator is reasonable.")
        self.check.set_active(True)
        fr = gtk.Frame()
        fr.add(self.check)
        fr.show_all()
        pack(vbox, [fr], end=True, bottom=True)
        
        # Finalize setup
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui
        
        
    def pack_line(self, box, radio, setup, info):
        hbox = gtk.HBox()
        hbox.pack_start(radio, 0, 0)
        hbox.pack_start(gtk.Label("  "), 0, 0)
        hbox.pack_end(info, 0, 0)
        if setup is not None:
            radio.connect("toggled", self.radio_toggled, setup)
            setup.set_sensitive(False)
            hbox.pack_end(setup, 0, 0)
        hbox.show_all()
        box.pack_start(hbox, 0, 0)

    def radio_toggled(self, radio, button):
        button.set_sensitive(radio.get_active())

    def lj_setup_window(self, widget):
        if not self.get_atoms():
            return
        lj_param = getattr(self, "lj_parameters", None)
        LJ_Window(self, lj_param, "lj_parameters")
        # When control is retuned, self.lj_parameters has been set.
        
    def get_atoms(self):
        "Make an atoms object from the active image"
        images = self.gui.images
        if images.natoms < 1:
            oops("No atoms present")
            return False
        self.atoms = Atoms(positions=images.P[0], symbols=images.Z,
                           cell=images.A[0], pbc=images.pbc)
        return True

    def apply(self, *widget):
        nochk = not self.check.get_active()
        if self.none_radio.get_active():
            self.gui.simulation['calc'] = None
            return True
        elif self.lj_radio.get_active():
            if nochk or self.lj_check():
                self.choose_lj()
                return True
        elif self.emt_radio.get_active():
            if nochk or self.emt_check():
                self.choose_emt()
                return True
        return False

    def ok(self, *widget):
        if self.apply():
            self.destroy()
            
    def lj_check(self):
        try:
            import asap3
        except ImportError:
            oops("ASAP is not installed. (Failed to import asap3)")
            return False
        if not hasattr(self, "lj_parameters"):
            oops("You must set up the Lennard-Jones parameters")
            return False
        try:
            self.atoms.set_calculator(asap3.LennardJones(**self.lj_parameters))
        except (asap3.AsapError, TypeError, ValueError), e:
            oops("Could not create useful Lennard-Jones calculator.",
                 str(e))
            return False
        return True

    def choose_lj(self):
        # Define a function on the fly!
        import asap3
        def lj_factory(p=self.lj_parameters, lj=asap3.LennardJones):
            return lj(**p)
        self.gui.simulation["calc"] = lj_factory

    def emt_get(self):
        import asap3
        provider_name = self.emt_setup.get_active_text()
        provider =  self.emt_param_info[provider_name]
        if provider is not None:
            provider = getattr(asap3, provider)
        return (asap3.EMT, provider, asap3)
                                      
    def emt_check(self):
        if not self.get_atoms():
            return False
        try:
            emt, provider, asap3 = self.emt_get()
        except ImportError:
            oops("ASAP is not installed. (Failed to import asap3)")
            return False
        try:
            if provider is not None:
                self.atoms.set_calculator(emt(provider()))
            else:
                self.atoms.set_calculator(emt())
        except (asap3.AsapError, TypeError, ValueError), e:
            oops("Could not attach EMT calculator to the atoms.",
                 str(e))
            return False
        return True

    def choose_emt(self):
        emt, provider, asap3 = self.emt_get()
        if provider is None:
            emt_factory = emt
        else:
            def emt_factory(emt=emt, prov=provider):
                return emt(prov())
        self.gui.simulation["calc"] = emt_factory

    
class InfoButton(gtk.Button):
    def __init__(self, txt):
        gtk.Button.__init__(self, "Info")
        self.txt = txt
        self.connect('clicked', self.run)

    def run(self, widget):
        dialog = gtk.MessageDialog(flags=gtk.DIALOG_MODAL,
                                   type=gtk.MESSAGE_INFO,
                                   buttons=gtk.BUTTONS_CLOSE)
        dialog.set_markup(self.txt)
        dialog.connect('response', lambda x, y: dialog.destroy())
        dialog.show()
        

class LJ_Window(gtk.Window):
    def __init__(self, owner, param, attrname):
        gtk.Window.__init__(self)
        self.set_title("Lennard-Jones parameters")
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        atnos = atoms.get_atomic_numbers()
        found = {}
        for z in atnos:
            found[z] = True
        self.present = found.keys()
        self.present.sort()  # Sorted list of atomic numbers
        nelem = len(self.present)
        vbox = gtk.VBox()
        label = gtk.Label("Specify the Lennard-Jones parameters here")
        pack(vbox, [label])
        pack(vbox, gtk.Label(""))
        pack(vbox, [gtk.Label("Epsilon (eV):")])
        tbl, self.epsilon_adj = self.makematrix(self.present)
        pack(vbox, [tbl])
        pack(vbox, gtk.Label(""))
        pack(vbox, [gtk.Label("Sigma (Å):")])
        tbl, self.sigma_adj = self.makematrix(self.present)
        pack(vbox, [tbl])
        self.modif = gtk.CheckButton("Shift to make smooth at cutoff")
        self.modif.set_active(True)
        pack(vbox, gtk.Label(""))
        pack(vbox, self.modif)
        pack(vbox, gtk.Label(""))
        butbox = gtk.HButtonBox()
        cancel_but = gtk.Button(stock=gtk.STOCK_CANCEL)
        cancel_but.connect('clicked', lambda widget: self.destroy())
        ok_but = gtk.Button(stock=gtk.STOCK_OK)
        ok_but.connect('clicked', self.ok)
        butbox.pack_start(cancel_but, 0, 0)
        butbox.pack_start(ok_but, 0, 0)
        butbox.show_all()
        pack(vbox, [butbox], end=True, bottom=True)
        vbox.show()
        self.add(vbox)

        # Now, set the parameters
        if param and param['elements'] == self.present:
            self.set_param(self.epsilon_adj, param["epsilon"], nelem)
            self.set_param(self.sigma_adj, param["sigma"], nelem)
            self.modif.set_active(param["modified"])

        self.show()
        self.grab_add()  # Lock all other windows
        
    def makematrix(self, present):
        nelem = len(present)
        adjdict = {}
        tbl = gtk.Table(2+nelem, 2+nelem)
        for i in range(nelem):
            s = chemical_symbols[present[i]]
            tbl.attach(gtk.Label(" " + str(present[i])), 0, 1, i, i+1)
            tbl.attach(gtk.Label("  "+s+" "), 1, 2, i, i+1)
            tbl.attach(gtk.Label(str(present[i])), i+2, i+3, 1+nelem, 2+nelem)
            tbl.attach(gtk.Label(s), i+2, i+3, nelem, 1+nelem)
            for j in range(i+1):
                adj = gtk.Adjustment(1.0, 0.0, 100.0, 0.1)
                spin = gtk.SpinButton(adj, 0.1, 3)
                tbl.attach(spin, 2+j, 3+j, i, i+1)
                adjdict[(i,j)] = adj
        tbl.show_all()
        return tbl, adjdict
    
    def set_param(self, adj, params, n):
        for i in range(n):
            for j in range(n):
                if j <= i:
                    adj[(i,j)].value = params[i,j]

    def get_param(self, adj, params, n):
        for i in range(n):
            for j in range(n):
                if j <= i:
                    params[i,j] = params[j,i] = adj[(i,j)].value


    def destroy(self):
        self.grab_remove()
        gtk.Window.destroy(self)

    def ok(self, *args):
        params = {}
        params["elements"] = copy(self.present)
        n = len(self.present)
        eps = np.zeros((n,n))
        self.get_param(self.epsilon_adj, eps, n)
        sigma = np.zeros((n,n))
        self.get_param(self.sigma_adj, sigma, n)
        params["epsilon"] = eps
        params["sigma"] = sigma
        params["modified"] = self.modif.get_active()
        setattr(self.owner, self.attrname, params)
        print params
        self.destroy()
        
