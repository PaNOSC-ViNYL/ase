from ase.gui.i18n import enable_localization
enable_localization()
import ase.gui.ui as ui
from ase.gui.gui import GUI



class OOPS:
    has_been_called = False

    def __call__(self, title, text):
        self.title = title
        self.has_been_called = True

    def called(self, title=None):
        result = self.has_been_called and (title is None or title == self.title)
        self.has_been_called = False
        return result


ui.oops = OOPS()


def nt(gui):
    yield
    nt = gui.nanotube_window()
    #yield
    nt.apply()
    #yield
    nt.element[1].value = '?'
    nt.apply()
    #yield
    assert ui.oops.called('No valid atoms.')
    nt.element[1].value = 'C'
    nt.ok()
    assert gui.images.natoms == 20
    #gui.exit()


gui = GUI()
gui.run(test=nt(gui))
