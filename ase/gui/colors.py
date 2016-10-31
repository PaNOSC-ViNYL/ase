"""colors.py - select how to color the atoms in the GUI."""
from gettext import gettext as _

import numpy as np

import ase.gui.ui as ui


class ColorWindow:
    """A window for selecting how to color the atoms."""
    def __init__(self, gui):
        win = ui.Window(_("Colors"))
        self.gui = gui
        win.add(ui.Label(_('Choose how the atoms are colored:')))
        values = ['jmol', 'tag', 'force']
        labels = [_('By atomic number, default "jmol" colors'),
                  _('By tag'),
                  _('By force')]
        win.add(ui.RadioButtons(labels, values, self.toggle))

        """
        self.radio_velocity = ui.RadioButton(self.radio_jmol, _('By velocity'))
        self.radio_charge = ui.RadioButton(self.radio_jmol, _('By charge'))
        self.radio_magnetic_moment = ui.RadioButton(
            self.radio_jmol, _('By magnetic moment'))
        self.radio_coordination = ui.RadioButton(
            self.radio_jmol, _('By coordination'))
        """

    def toggle(self, value):
        if value == 'jmol':
            self.gui.set_colors()
        else:
            self.gui.colormode = value
            scalars = np.array([self.gui.get_color_scalars(i)
                                for i in range(len(self.gui.images))])
            mn = scalars.min()
            mx = scalars.max()
            colorscale = ['#{0:02X}AA00'.format(red)
                          for red in range(0, 240, 10)]
            self.gui.colormode_data = colorscale, mn, mx

    def notify_atoms_changed(self):
        "Called by gui object when the atoms have changed."
        self.destroy()
