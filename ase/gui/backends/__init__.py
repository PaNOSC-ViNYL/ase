from ase.utils import import_module


def select_backend(name):
    import ase.gui.backend as be
    module = import_module('ase.gui.backends.' + name)
    be.Button = module.Button
