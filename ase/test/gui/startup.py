from ase.gui.i18n import enable_localization
import ase.gui.ui as ui
from ase import Atoms
from ase.gui.gui import GUI
from ase.gui.images import Images
from ase.calculators.singlepoint import SinglePointCalculator
enable_localization()


class OOPS:
    has_been_called = False

    def __call__(self, title, text):
        self.title = title
        self.has_been_called = True

    def called(self, title=None):
        result = self.has_been_called and (title is None or
                                           title == self.title)
        self.has_been_called = False
        return result


ui.oops = OOPS()


def nt(gui):
    yield
    nt = gui.nanotube_window()
    # yield
    nt.apply()
    # yield
    nt.element[1].value = '?'
    nt.apply()
    # yield
    assert ui.oops.called('No valid atoms.')
    nt.element[1].value = 'C'
    nt.ok()
    assert gui.images.natoms == 20
    # gui.exit()


def np(gui):
    np = gui.nanoparticle_window()
    yield


def color(gui):
    yield
    c = gui.colors_window()
    c.toggle('force')
    yield
    assert [button.active for button in c.radio.buttons] == [1, 0, 1, 0, 0, 0]


if 1:
    gui = GUI()
    gui.run(test=nt(gui))
    # gui.run(test=np(gui))

if 0:
    h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1)])
    h2.calc = SinglePointCalculator(h2, forces=[(0, 0, 0), (0, 0, 1)])
    gui = GUI(Images([h2]))
    gui.run(test=color(gui))
