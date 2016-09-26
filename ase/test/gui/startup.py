from ase.gui.i18n import enable_localization
enable_localization()
from ase.gui.gui import GUI


def test(gui):
    yield
    nt = gui.nanotube_window()
    yield
    nt.apply()
    yield
    nt.element[1].value = '?'
    nt.apply()
    yield
    #gui.exit()


gui = GUI()
gui.run(test=test(gui))
