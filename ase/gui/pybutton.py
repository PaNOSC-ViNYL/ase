"""A button for displaying Python code.

When pressed, it opens a window displaying some Python code, or an error
message if no Python code is ready.

"""
import functools
from gettext import gettext as _

import ase.gui.ui as ui


def pybutton(title, obj):
    return ui.Button('Python', functools.partial(pywindow, title, obj))


def pywindow(title, obj):
    code = obj.python
    if code is None:
        ui.oops(
            _('No Python code'),
            _('You have not (yet) specified a consistent set of parameters.'))
    else:
        win = ui.Window(title)
        win.add(ui.Text(code))
