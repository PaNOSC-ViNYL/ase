# encoding: utf-8
"calculator.py - module for choosing a calculator."

import gtk
import numpy as np
from copy import copy
from ase.gui.setupwindow import SetupWindow
from ase.gui.progress import DefaultProgressIndicator, GpawProgressIndicator
from ase.gui.widgets import pack, oops, cancel_apply_ok
from ase import Atoms
from ase.data import chemical_symbols
import ase

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

aseemt_info_txt = """\
The EMT potential is a many-body potential, giving a
good description of the late transition metals crystalling
in the FCC crystal structure.  The elements described by the
main set of EMT parameters are Al, Ni, Cu, Pd, Ag, Pt, and
Au.  In addition, this implementation allows for the use of
H, N, O and C adatoms, although the description of these is
most likely not very good.

<b>This is the ASE implementation of EMT.</b> For large
simulations the ASAP implementation is more suitable; this
implementation is mainly to make EMT available when ASAP is
not installed.
"""

gpaw_info_txt = """\
GPAW implements Density Functional Theory using a
<b>G</b>rid-based real-space representation of the wave
functions, and the <b>P</b>rojector <b>A</b>ugmented <b>W</b>ave
method for handling the core regions.  
"""

emt_parameters = (
    ("Default (Al, Ni, Cu, Pd, Ag, Pt, Au)", None),
    ("Alternative Cu, Ag and Au", "EMTRasmussenParameters"),
    ("Ruthenium", "EMThcpParameters"),
    ("CuMg and CuZr metallic glass", "EMTMetalGlassParameters")
    )

class SetCalculator(SetupWindow):
    "Window for selecting a calculator."

    # List the names of the radio button attributes
    radios = ("none", "lj", "emt", "aseemt", "gpaw")
    # List the names of the parameter dictionaries
    paramdicts = ("lj_parameters",)
    # The name used to store parameters on the gui object
    classname = "SetCalculator"
    
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

        # EMT (ASE implementation)
        self.aseemt_radio = gtk.RadioButton(
            self.none_radio, "EMT - Effective Medium Theory (ASE)")
        self.aseemt_info = InfoButton(aseemt_info_txt)
        self.pack_line(vbox, self.aseemt_radio, None, self.aseemt_info)

        # GPAW
        self.gpaw_radio = gtk.RadioButton(self.none_radio,
                                          "Density Functional Theory (GPAW)")
        self.gpaw_setup = gtk.Button("Setup")
        self.gpaw_info = InfoButton(gpaw_info_txt)
        self.gpaw_setup.connect("clicked", self.gpaw_setup_window)
        self.pack_line(vbox, self.gpaw_radio, self.gpaw_setup, self.gpaw_info)
        
        # Buttons etc.
        pack(vbox, gtk.Label(""))
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
        self.load_state()
        
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
        
    def gpaw_setup_window(self, widget):
        if not self.get_atoms():
            return
        gpaw_param = getattr(self, "gpaw_parameters", None)
        GPAW_Window(self, gpaw_param, "gpaw_parameters")
        # When control is retuned, self.gpaw_parameters has been set.
        
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
        if self.do_apply():
            self.save_state()
            return True
        else:
            return False
        
    def do_apply(self):
        nochk = not self.check.get_active()
        self.gui.simulation["progress"] = DefaultProgressIndicator()
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
        elif self.aseemt_radio.get_active():
            if nochk or self.aseemt_check():
                self.choose_aseemt()
                return True
        elif self.gpaw_radio.get_active():
            if nochk or self.gpaw_check():
                self.choose_gpaw()
                return True
        return False

    def ok(self, *widget):
        if self.apply():
            self.destroy()

    def save_state(self):
        state = {}
        for r in self.radios:
            radiobutton = getattr(self, r+"_radio")
            if radiobutton.get_active():
                state["radio"] = r
        state["emtsetup"] = self.emt_setup.get_active()
        state["check"] = self.check.get_active()
        for p in self.paramdicts:
            if hasattr(self, p):
                state[p] = getattr(self, p)
        self.gui.module_state[self.classname] = state

    def load_state(self):
        try:
            state = self.gui.module_state[self.classname]
        except KeyError:
            return
        r = state["radio"]
        radiobutton = getattr(self, r+"_radio")
        radiobutton.set_active(True)
        self.emt_setup.set_active(state["emtsetup"])
        self.check.set_active(state["check"])
        for p in self.paramdicts:
            if state.has_key(p):
                setattr(self, p, state[p])
            
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

    def aseemt_check(self):
        return self.element_check("ASE EMT", ['H', 'Al', 'Cu', 'Ag', 'Au',
                                              'Ni', 'Pd', 'Pt', 'C', 'N', 'O'])

    def choose_aseemt(self):
        self.gui.simulation["calc"] = ase.EMT
        ase.EMT.disabled = False  # In case Asap has been imported.

    def gpaw_check(self):
        try:
            import gpaw
        except ImportError:
            oops("GPAW is not installed. (Failed to import gpaw)")
            return False
        if not hasattr(self, "gpaw_parameters"):
            oops("You must set up the GPAW parameters")
            return False
        return True

    def choose_gpaw(self):
        # This reuses the same GPAW object.
        try:
            import gpaw
        except ImportError:
            oops("GPAW is not installed. (Failed to import gpaw)")
            return False
        p = self.gpaw_parameters
        use = ["xc", "kpts", "mode"]
        if p["use_h"]:
            use.append("h")
        else:
            use.append("gpts")
        if p["mode"] == "lcao":
            use.append("basis")
        gpaw_param = {}
        for s in use:
            gpaw_param[s] = p[s]
        if p["use mixer"]:
            mx = getattr(gpaw, p["mixer"])
            mx_args = {}
            mx_arg_n = ["beta", "nmaxold", "weight"]
            if p["mixer"] == "MixerDiff":
                mx_arg_n.extend(["beta_m", "nmaxold_m", "weight_m"])
            for s in mx_arg_n:
                mx_args[s] = p[s]
            gpaw_param["mixer"] = mx(**mx_args)
        progress = GpawProgressIndicator()
        self.gui.simulation["progress"] = progress
        gpaw_param["txt"] = progress.get_gpaw_stream()
        gpaw_calc = gpaw.GPAW(**gpaw_param)
        def gpaw_factory(calc = gpaw_calc):
            return calc
        self.gui.simulation["calc"] = gpaw_factory
                
    def element_check(self, name, elements):
        "Check that all atoms are allowed"
        elements = [ase.data.atomic_numbers[s] for s in elements]
        elements_dict = {}
        for e in elements:
            elements_dict[e] = True
        if not self.get_atoms():
            return False
        try:
            for e in self.atoms.get_atomic_numbers():
                elements_dict[e]
        except KeyError:
            oops("Element %s is not allowed by the '%s' calculator"
                 % (ase.data.chemical_symbols[e], name))
            return False
        return True
    
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
        self.destroy()


class GPAW_Window(gtk.Window):
    gpaw_xc_list = ['LDA', 'PBE', 'RPBE', 'revPBE']
    gpaw_xc_default = 'PBE'
    def __init__(self, owner, param, attrname):
        gtk.Window.__init__(self)
        self.set_title("GPAW parameters")
        self.owner = owner
        self.attrname = attrname
        atoms = owner.atoms
        self.ucell = atoms.get_cell()
        self.size = tuple([self.ucell[i,i] for i in range(3)])
        self.pbc = atoms.get_pbc()
        self.orthogonal = self.isorthogonal(self.ucell)
        self.natoms = len(atoms)
        
        vbox = gtk.VBox()
        #label = gtk.Label("Specify the GPAW parameters here")
        #pack(vbox, [label])

        # Print some info
        txt = "%i atoms.\n" % (self.natoms,)
        if self.orthogonal:
            txt += "Orthogonal unit cell: %.2f x %.2f x %.2f Å." % self.size
        else:
            txt += "Non-orthogonal unit cell:\n"
            txt += str(self.ucell)
        pack(vbox, [gtk.Label(txt)])
        
        # XC potential
        self.xc = gtk.combo_box_new_text()
        for i, x in enumerate(self.gpaw_xc_list):
            self.xc.append_text(x)
            if x == self.gpaw_xc_default:
                self.xc.set_active(i)
        pack(vbox, [gtk.Label("Exchange-correlation functional: "),
                    self.xc])
        
        # Grid spacing
        self.radio_h = gtk.RadioButton(None, "Grid spacing")
        self.h = gtk.Adjustment(0.18, 0.0, 1.0, 0.01)
        self.h_spin = gtk.SpinButton(self.h, 0, 2)
        pack(vbox, [self.radio_h, gtk.Label(" h = "), self.h_spin,
                    gtk.Label("Å")])
        self.radio_gpts = gtk.RadioButton(self.radio_h, "Grid points")
        self.gpts = []
        self.gpts_spin = []
        for i in range(3):
            g = gtk.Adjustment(4, 4, 1000, 4)
            s = gtk.SpinButton(g, 0, 0)
            self.gpts.append(g)
            self.gpts_spin.append(s)
        self.gpts_hlabel = gtk.Label("")
        self.gpts_hlabel_format = "h<sub>eff</sub> = (%.3f, %.3f, %.3f) Å"
        pack(vbox, [self.radio_gpts, gtk.Label(" gpts = ("), self.gpts_spin[0],
                    gtk.Label(", "), self.gpts_spin[1], gtk.Label(", "),
                    self.gpts_spin[2], gtk.Label(")  "), self.gpts_hlabel])
        self.radio_h.connect("toggled", self.radio_grid_toggled)
        self.radio_gpts.connect("toggled", self.radio_grid_toggled)
        self.radio_grid_toggled(None)
        for g in self.gpts:
            g.connect("value-changed", self.gpts_changed)
        self.h.connect("value-changed", self.h_changed)
        
        # K-points
        self.kpts = []
        self.kpts_spin = []
        for i in range(3):
            if self.pbc[i] and self.orthogonal:
                default = np.ceil(20.0 / self.size[i])
            else:
                default = 1
            g = gtk.Adjustment(default, 1, 100, 1)
            s = gtk.SpinButton(g, 0, 0)
            self.kpts.append(g)
            self.kpts_spin.append(s)
            if not self.pbc[i]:
                s.set_sensitive(False)
            g.connect("value-changed", self.k_changed)
        pack(vbox, [gtk.Label("k-points  k = ("), self.kpts_spin[0],
                    gtk.Label(", "), self.kpts_spin[1], gtk.Label(", "),
                    self.kpts_spin[2], gtk.Label(")")])
        self.kpts_label = gtk.Label("")
        self.kpts_label_format = "k-points x size:  (%.1f, %.1f, %.1f) Å"
        pack(vbox, [self.kpts_label])
        self.k_changed()
        
        # Spin polarized
        self.spinpol = gtk.CheckButton("Spin polarized")
        pack(vbox, [self.spinpol])
        pack(vbox, gtk.Label(""))

        # Mode and basis functions
        self.mode = gtk.combo_box_new_text()
        self.mode.append_text("FD - Finite Difference (grid) mode")
        self.mode.append_text("LCAO - Linear Combination of Atomic Orbitals")
        self.mode.set_active(0)
        pack(vbox, [gtk.Label("Mode: "), self.mode])
        self.basis = gtk.combo_box_new_text()
        self.basis.append_text("sz - Single Zeta")
        self.basis.append_text("szp - Single Zeta polarized")
        self.basis.append_text("dzp - Double Zeta polarized")
        self.basis.set_active(2) # dzp
        pack(vbox, [gtk.Label("Basis functions: "), self.basis])
        pack(vbox, gtk.Label(""))
        self.mode.connect("changed", self.mode_changed)
        self.mode_changed()
        
        # Mixer
        self.use_mixer = gtk.CheckButton("Non-standard mixer parameters")
        pack(vbox, [self.use_mixer])
        self.radio_mixer = gtk.RadioButton(None, "Mixer   ")
        self.radio_mixersum = gtk.RadioButton(self.radio_mixer, "MixerSum   ")
        self.radio_mixerdiff = gtk.RadioButton(self.radio_mixer, "MixerDiff")
        pack(vbox, [self.radio_mixer, self.radio_mixersum,
                    self.radio_mixerdiff])
        self.beta_adj = gtk.Adjustment(0.25, 0.0, 1.0, 0.05)
        self.beta_spin = gtk.SpinButton(self.beta_adj, 0, 2)
        self.nmaxold_adj = gtk.Adjustment(3, 1, 10, 1)
        self.nmaxold_spin = gtk.SpinButton(self.nmaxold_adj, 0, 0)
        self.weight_adj = gtk.Adjustment(50, 1, 500, 1)
        self.weight_spin = gtk.SpinButton(self.weight_adj, 0, 0)
        pack(vbox, [gtk.Label("beta = "), self.beta_spin,
                    gtk.Label("  nmaxold = "), self.nmaxold_spin,
                    gtk.Label("  weight = "), self.weight_spin])
        self.beta_m_adj = gtk.Adjustment(0.70, 0.0, 1.0, 0.05)
        self.beta_m_spin = gtk.SpinButton(self.beta_m_adj, 0, 2)
        self.nmaxold_m_adj = gtk.Adjustment(2, 1, 10, 1)
        self.nmaxold_m_spin = gtk.SpinButton(self.nmaxold_m_adj, 0, 0)
        self.weight_m_adj = gtk.Adjustment(10, 1, 500, 1)
        self.weight_m_spin = gtk.SpinButton(self.weight_m_adj, 0, 0)
        pack(vbox, [gtk.Label("beta_m = "), self.beta_m_spin,
                    gtk.Label("  nmaxold_m = "), self.nmaxold_m_spin,
                    gtk.Label("  weight_m = "), self.weight_m_spin])
        for but in (self.spinpol, self.use_mixer, self.radio_mixer,
                    self.radio_mixersum, self.radio_mixerdiff):
            but.connect("clicked", self.mixer_changed)
        self.mixer_changed()
        
        # Eigensolver
        # Poisson-solver
        
        vbox.show()
        self.add(vbox)

        # Buttons at the bottom
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

        # Set stored parameters
        if param:
            self.xc.set_active(param["xc#"])
            if param["use_h"]:
                self.radio_h.set_active(True)
            else:
                self.radio_gpts.set_active(True)
            for i in range(3):
                self.gpts[i].value = param["gpts"][i]
                self.kpts[i].value = param["kpts"][i]
            self.spinpol.set_active(param["spinpol"])
            self.mode.set_active(param["mode#"])
            self.basis.set_active(param["basis#"])
            self.use_mixer.set_active(param["use mixer"])
            getattr(self, "radio_"+param["mixer"].lower()).set_active(True)
            for t in ("beta", "nmaxold", "weight", "beta_m", "nmaxold_m",
                      "weight_m"):                    
                getattr(self, t+"_adj").value = param[t]

        self.show()
        self.grab_add()  # Lock all other windows

    def radio_grid_toggled(self, widget):
        hmode = self.radio_h.get_active()
        self.h_spin.set_sensitive(hmode)
        for s in self.gpts_spin:
            s.set_sensitive(not hmode)
        self.gpts_changed()

    def gpts_changed(self, *args):
        if self.radio_gpts.get_active():
            g = np.array([int(g.value) for g in self.gpts])
            size = np.array([self.ucell[i,i] for i in range(3)])
            txt = self.gpts_hlabel_format % tuple(size / g)
            self.gpts_hlabel.set_markup(txt)
        else:
            self.gpts_hlabel.set_markup("")

    def h_changed(self, *args):
        h = self.h.value
        for i in range(3):
            g = 4 * round(self.ucell[i,i] / (4*h))
            self.gpts[i].value = g

    def k_changed(self, *args):
        if self.orthogonal:
            size = [self.kpts[i].value * self.size[i] for i in range(3)]
        self.kpts_label.set_text(self.kpts_label_format % tuple(size))

    def mode_changed(self, *args):
        self.basis.set_sensitive(self.mode.get_active() == 1)

    def mixer_changed(self, *args):
        radios = (self.radio_mixer, self.radio_mixersum, self.radio_mixerdiff)
        spin1 = (self.beta_spin, self.nmaxold_spin, self.weight_spin)
        spin2 = (self.beta_m_spin, self.nmaxold_m_spin, self.weight_m_spin)
        if self.use_mixer.get_active():
            # Mixer parameters can be specified.
            if self.spinpol.get_active():
                self.radio_mixer.set_sensitive(False)
                self.radio_mixersum.set_sensitive(True)
                self.radio_mixerdiff.set_sensitive(True)
                if self.radio_mixer.get_active():
                    self.radio_mixersum.set_active(True)
            else:
                self.radio_mixer.set_sensitive(True)
                self.radio_mixersum.set_sensitive(False)
                self.radio_mixerdiff.set_sensitive(False)
                self.radio_mixer.set_active(True)
            if self.radio_mixerdiff.get_active():
                active = spin1 + spin2
                passive = ()
            else:
                active = spin1
                passive = spin2
            for widget in active:
                widget.set_sensitive(True)
            for widget in passive:
                widget.set_sensitive(False)
        else:
            # No mixer parameters
            for widget in radios + spin1 + spin2:
                widget.set_sensitive(False)
                
    def isorthogonal(self, matrix):
        ortho = True
        for i in range(3):
            for j in range(3):
                if i != j and matrix[i][j] != 0.0:
                    ortho = False
        return ortho

    def ok(self, *args):
        param = {}
        param["xc"] = self.xc.get_active_text()
        param["xc#"] = self.xc.get_active()
        param["use_h"] = self.radio_h.get_active()
        param["h"] = self.h.value
        param["gpts"] = [int(g.value) for g in self.gpts]
        param["kpts"] = [int(k.value) for k in self.kpts]
        param["spinpol"] = self.spinpol.get_active()
        param["mode"] = self.mode.get_active_text().split()[0].lower()
        param["mode#"] = self.mode.get_active()
        param["basis"] = self.basis.get_active_text().split()[0].lower()
        param["basis#"] = self.basis.get_active()
        param["use mixer"] = self.use_mixer.get_active()
        if self.radio_mixer.get_active():
            m = "Mixer"
        elif self.radio_mixersum.get_active():
            m = "MixerSum"
        else:
            assert self.radio_mixerdiff.get_active()
            m = "MixerDiff"
        param["mixer"] = m
        for t in ("beta", "nmaxold", "weight", "beta_m", "nmaxold_m",
                  "weight_m"):
            param[t] = getattr(self, t+"_adj").value
        setattr(self.owner, self.attrname, param)
        self.destroy()
