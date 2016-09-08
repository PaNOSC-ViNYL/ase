
"""pybutton.py - a button for displaying Python code.

This module defines two classes, together they implement a button that
a module can use to display Python.

PyButton
--------

PyButton is a gkt.Button with the label 'Python'.  When pressed, it
opens a PyWindow displaying some Python code, or an error message if
no Python code is ready.

The script is stored in the attribute .python, it is the
responsibility of the owning object to keep this attribute up to date:
when pressing the Apply button would result in a sensible
configuration being created, the python attribute must be set to a
string creating this code.  When pressing Apply would cause an error,
the python attribute must be set to None.

PyWindow
--------

Displays the Python code.  This object is created by the PyButton
object when needed.
"""

import ase.gui.ui as ui
import time
from ase.gui.widgets import oops, pack


class PyButton:
    "A button for displaying Python code."
    def __init__(self, title):
        ui.Button.__init__(self, _("Python"))
        self.title = title
        self.python = None
        self.connect_after('clicked', self.run)

    def run(self, *args):
        "The method called when the button is clicked."
        if self.python:
            now = time.ctime()
            PyWindow(self.title, now, self.python)
        else:
            oops(_("No Python code"),
                 _("You have not (yet) specified a "
                   "consistent set of parameters."))

fr1_template = _("""
Title: %(title)s
Time: %(time)s
""")

class PyWindow:
    "A window displaying Python code."
    def __init__(self, title, time, code):
        ui.Window.__init__(self)
        self.set_title(_("ag: Python code"))
        vbox = ui.VBox()
        lbl = ui.Label(fr1_template % dict(title=title, time=time))
        lbl.set_alignment(0.0, 0.5)
        fr = ui.Frame(_("Information:"))
        fr.add(lbl)
        pack(vbox, fr)
        txtbuf = ui.TextBuffer()
        txtbuf.set_text(code)
        txtview = ui.TextView(txtbuf)
        txtview.set_editable(False)
        fr = ui.Frame(_("Python code:"))
        fr.add(txtview)
        fr.set_label_align(0.0, 0.5)
        pack(vbox, fr)
        but =  ui.Button('OK')
        but.connect('clicked', lambda x: self.destroy())
        pack(vbox, [but], end=True, bottom=True)
        self.add(vbox)
        self.show_all()
        
