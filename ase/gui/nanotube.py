# encoding: utf-8
"""nanotube.py - Window for setting up Carbon nanotubes and similar tubes.
"""
from gettext import gettext as _

import numpy as np

import ase.gui.ui as ui
from ase.gui.widgets import pack, cancel_apply_ok, oops
from ase.gui.setupwindow import SetupWindow
from ase.gui.pybutton import pybutton
from ase.gui.status import formula
from ase.build import nanotube
import ase


introtext = _("""\
Set up a Carbon nanotube by specifying the (n,m) roll-up vector.
Please note that m <= n.

Nanotubes of other elements can be made by specifying the element
and bond length.\
""")

py_template = """\
from ase.build import nanotube
atoms = nanotube({n}, {m}, length={length}, bond={bl:.3f}, symbol='{symb}')
"""

label_template = _('{natoms} atoms: {symbols}, diameter: {diameter:.3f} Å, '
                   'cell volume: {volume:.3f} Å<sup>3</sup>')


class SetupNanotube:
    "Window for setting up a (Carbon) nanotube."
    def __init__(self, gui):
        win = ui.Window(_("Nanotube"))
        win.add(ui.Text(introtext))

        self.element = ui.Entry('C', 3, self.make)
        self.bondlength = ui.SpinBox(1.42, 0.0, 10.0, 0.01, self.make)

        win.add([_("Element: "),
                 self.element,
                 _("  Bond length: "),
                 self.bondlength,
                 _(u"Å")])

        self.elementinfo = ui.Label("")
        self.elementinfo.modify_fg(ui.STATE_NORMAL,
                                   '#FF0000')
        pack(vbox, [self.elementinfo])
        pack(vbox, ui.Label(""))

        # Choose the structure.
        pack(vbox, [ui.Label(_("Select roll-up vector (n,m) "
                                "and tube length:"))])
        label1 = ui.Label("n: ")
        label2 = ui.Label("  m: ")
        self.n = ui.Adjustment(5, 1, 100, 1)
        self.m = ui.Adjustment(5, 0, 100, 1)
        spinn = ui.SpinButton(self.n, 0, 0)
        spinm = ui.SpinButton(self.m, 0, 0)
        label3 = ui.Label(_("  Length: "))
        self.length = ui.Adjustment(1, 1, 100, 1)
        spinl = ui.SpinButton(self.length, 0, 0)
        pack(vbox, [label1, spinn, label2, spinm, label3, spinl])
        self.err = ui.Label("")
        self.err.modify_fg(ui.STATE_NORMAL, '#FF0000')
        pack(vbox, [self.err])
        pack(vbox, ui.Label(""))

        self.status = ui.Label("")
        pack(vbox,[self.status])
        pack(vbox,[ui.Label("")])

        # Buttons
        self.pybut = PyButton(_("Creating a nanoparticle."))
        self.pybut.connect('clicked', self.makeatoms)
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [self.pybut, buts], end=True, bottom=True)

        # Finalize setup
        self.makeatoms()
        self.bondlength.connect('value-changed', self.makeatoms)
        self.m.connect('value-changed', self.makeatoms)
        self.n.connect('value-changed', self.makeatoms)
        self.length.connect('value-changed', self.makeatoms)
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def update_element(self, *args):
        "Called when a new element may have been entered."
        # Assumes the element widget is self.element and that a label
        # for errors is self.elementinfo.  The chemical symbol is
        # placed in self.legalelement - or None if the element is
        # invalid.
        elem = self.element.get_text()
        if not elem:
            self.invalid_element(_("  No element specified!"))
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
        self.elementinfo.set_text("")
        self.legal_element = symb
        return True

    def makeatoms(self, *args):
        self.update_element()
        if self.legal_element is None:
            self.atoms = None
            self.pybut.python = None
        else:
            n = int(self.n.value)
            m = int(self.m.value)
            symb = self.legal_element
            length = int(self.length.value)
            bl = self.bondlength.value
            self.atoms = nanotube(n, m, length=length, bond=bl, symbol=symb)
            # XXX can this be translated?
            self.pybut.python = py_template % {'n': n, 'm':m, 'length':length,
                                               'symb':symb, 'bl':bl}
            h = np.zeros(3)
            uc = self.atoms.get_cell()
            for i in range(3):
                norm = np.cross(uc[i-1], uc[i-2])
                norm /= np.sqrt(np.dot(norm, norm))
                h[i] = np.abs(np.dot(norm, uc[i]))
            label = label_template % {'natoms'   : len(self.atoms),
                                      'symbols'  : formula(self.atoms.get_atomic_numbers()),
                                      'volume'   : self.atoms.get_volume(),
                                      'diameter' : self.atoms.get_cell()[0][0]/2.0}
            self.status.set_markup(label)

    def apply(self, *args):
        self.makeatoms()
        if self.atoms is not None:
            self.gui.new_atoms(self.atoms)
            return True
        else:
            oops(_("No valid atoms."),
                 _("You have not (yet) specified a consistent "
                   "set of parameters."))
            return False

    def ok(self, *args):
        if self.apply():
            self.destroy()
