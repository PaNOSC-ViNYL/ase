"""colors.py - select how to color the atoms in the GUI."""
from gettext import gettext as _

import numpy as np

import ase.gui.ui as ui


class ColorWindow:
    """A window for selecting how to color the atoms."""
    def __init__(self, gui):
        self.win = ui.Window(_("Colors"))
        self.gui = gui
        self.win.add(ui.Label(_('Choose how the atoms are colored:')))
        values = ['jmol', 'tag', 'force', 'velocity', 'charge', 'magmom']
        labels = [_('By atomic number, default "jmol" colors'),
                  _('By tag'),
                  _('By force'),
                  _('By velocity'),
                  _('By charge'),
                  _('By magnetic moment')]
        radio = ui.RadioButtons(labels, values, self.toggle)
        radio.value = gui.colormode
        self.win.add(radio)
        self.deactivate(radio)
        self.radio = radio  # stored for testing puposes
        self.label = ui.Label()
        self.win.add(self.label)

    def deactivate(self, radio):
        images = self.gui.images
        if not images.T.any():
            radio['tag'].active = False
        if not np.isfinite(images.F).all():
            radio['force'].active = False
        if not np.isfinite(images.V).all():
            radio['velocity'].active = False
        if not images.q.any():
            radio['charge'].active = False
        if not images.M.any():
            radio['magmom'].active = False

    def toggle(self, value):
        if value == 'jmol':
            self.gui.set_colors()
            text = ''
        else:
            self.gui.colormode = value
            scalars = np.array([self.gui.get_color_scalars(i)
                                for i in range(len(self.gui.images))])
            mn = scalars.min()
            mx = scalars.max()
            colorscale = ['#{0:02X}AA00'.format(red)
                          for red in range(0, 240, 10)]
            self.gui.colormode_data = colorscale, mn, mx

            unit = {'tag': '',
                    'force': 'eV/Ang',
                    'velocity': '??',
                    'charge': '|e|',
                    'magmom': 'muB'}[value]
            text = '[{},{}]: [{},{}] {}'.format(
                _('Green'), _('Yellow'), mn, mx, unit)

        self.label.text = text
        self.gui.draw()

    def notify_atoms_changed(self):
        "Called by gui object when the atoms have changed."
        self.win.close()
