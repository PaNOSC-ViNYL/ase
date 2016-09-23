from ase.gui.i18n import enable_localization
enable_localization()
from ase.gui.gui import GUI


def test(gui):
    yield
    gui.exit()


gui = GUI()
gui.run(test=test(gui))
