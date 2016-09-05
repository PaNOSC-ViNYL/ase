from ase.gui.tk import *#Button, MainWindow, Window, Label
from ase.utils import import_module


def select_backend(name):
    import ase.gui.ui as ui
    module = import_module('ase.gui.' + name)
    ui.Button = module.Button
