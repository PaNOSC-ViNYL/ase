from __future__ import print_function
import sys

import gtk
import ase.gui.ui as ui


class Debug(ui.Window):
    def __init__(self, gui):
        self.gui = gui
        ui.Window.__init__(self)
        self.set_title(_('Debug'))
        entry = ui.Entry(200)
        self.add(entry)
        entry.connect('activate', self.enter, entry)
        entry.show()
        self.show()

    def enter(self, widget, entry):
        import numpy as np
        print(eval(entry.get_text(),
                   {'g': self.gui,
                    'np': np}),
              file=sys.stderr)
