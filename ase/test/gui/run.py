import argparse
from gettext import gettext as _
import numpy as np

from ase.build import molecule
from ase.test import NotAvailable

try:
    from ase.gui.i18n import enable_localization
    import ase.gui.ui as ui
except ImportError:
    raise NotAvailable

from ase import Atoms
from ase.gui.gui import GUI
from ase.calculators.singlepoint import SinglePointCalculator

enable_localization()


class OOPS:
    has_been_called = False

    def __call__(self, title, text):
        self.title = title
        self.has_been_called = True

    def called(self, title=None):
        result = self.has_been_called and (title is None or
                                           _(title) == self.title)
        self.has_been_called = False
        return result


ui.oops = OOPS()

tests = []


def test(f):
    tests.append(f.__name__)
    return f


@test
def nt(gui):
    nt = gui.nanotube_window()
    nt.apply()
    nt.element[1].value = '?'
    nt.apply()
    assert ui.oops.called('No valid atoms.')
    nt.element[1].value = 'C'
    nt.ok()
    assert gui.images.natoms == 20


@test
def nanopartickle(gui):
    gui.nanoparticle_window()


@test
def color(gui):
    a = Atoms('C10', magmoms=np.linspace(1, -1, 10))
    a.positions[:] = np.linspace(0, 9, 10)[:, None]
    a.calc = SinglePointCalculator(a, forces=a.positions)
    gui.new_atoms(a)
    c = gui.colors_window()
    c.toggle('force')
    text = c.toggle('magmom')
    assert [button.active for button in c.radio.buttons] == [1, 0, 1, 0, 0, 1]
    assert text.rsplit('[', 1)[1].startswith('-1.0,1.0]')


@test
def settings(gui):
    gui.new_atoms(molecule('H2O'))
    s = gui.settings()
    s.scale.value = 1.9
    s.scale_radii()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('tests', nargs='*')
    p.add_argument('-p', '--pause', action='store_true')
    args = p.parse_args()
    for name in args.tests or tests:
        for n in tests:
            if n.startswith(name):
                name = n
                break
        else:
            1 / 0
        print(name)
        test = globals()[name]
        gui = GUI()

        def f():
            test(gui)
            if not args.pause:
                gui.exit()
        gui.run(test=f)
