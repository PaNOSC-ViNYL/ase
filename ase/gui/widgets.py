from __future__ import print_function
import sys
import ase.gui.ui as ui
import re


def element():
    return [
class cancel_apply_ok:
    "Widget with Cancel, Apply and OK buttons.  The arguments are callbacks."
    def __init__(self, cancel, apply, ok):
        ui.HButtonBox.__init__(self)
        cancel_but = ui.Button('Cancel')
        cancel_but.connect('clicked', cancel)
        apply_but = ui.Button('Apply')
        apply_but.connect('clicked', apply)
        ok_but = ui.Button('OK')
        ok_but.connect('clicked', ok)
        for w in (cancel_but, apply_but, ok_but):
            self.pack_start(w, 0, 0)
            w.show()
        # self.show_all()

