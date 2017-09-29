from __future__ import unicode_literals
from ase.gui.i18n import _

import ase.gui.ui as ui
from ase.gui.widgets import Element


def txt2pos(txt):
    try:
        x, y, z = (float(x) for x in txt.split(','))
    except ValueError as ex:
        ui.error(_('Bad position'), ex)
    else:
        return x, y, z


class AddAtoms:
    def __init__(self, gui):
        # XXXXXXXXXXX still array based, not Atoms-based.  Will crash
        win = ui.Window(_('Add atoms'))
        self.element = Element(callback=self.add_absolute,
                               allow_molecule=True)
        win.add(self.element)
        self.absolute_position = ui.Entry('0,0,0')
        self.relative_position = ui.Entry('1.5,0,0')
        win.add([_('Absolute position:'),
                 self.absolute_position,
                 ui.Button(_('Add'), self.add_absolute)])
        win.add([_('Relative to average position (of selection):'),
                 self.relative_position,
                 ui.Button(_('Add'), self.add_relative)])
        self.gui = gui
        self.element.grab_focus()

    def add_absolute(self, callbackarg=None):
        pos = txt2pos(self.absolute_position.value)
        self.add(pos)

    def add_relative(self):
        rpos = txt2pos(self.relative_position.value)
        pos = self.gui.atoms.positions
        if self.gui.images.selected.any():
            pos = pos[self.gui.images.selected[:len(pos)]]
        if len(pos) == 0:
            ui.error('No atoms present')
        else:
            center = pos.mean(0)
            self.add(center + rpos)

    def add(self, pos):
        newatoms = self.element.get_atoms()
        if pos is None or newatoms is None:
            return

        newatoms.center(about=pos)
        atoms = self.gui.atoms
        atoms += newatoms

        if len(atoms) > self.gui.images.maxnatoms:
            self.gui.images.initialize(list(self.gui.images),
                                       self.gui.images.filenames)
        self.gui.images.selected[:] = False
        # 'selected' array may be longer than current atoms
        self.gui.images.selected[len(atoms) - 1] = True
        self.gui.set_frame()
        self.gui.draw()
