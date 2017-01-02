import numpy as np

import ase.gui.ui as ui


class Repeat:
    def __init__(self, gui):
        win = ui.Window('Repeat')
        win.add('Repeat atoms:')
        self.repeat = [ui.SpinBox(r, 1, 9, 1, self.change)
                       for r in gui.images.repeat]
        win.add(self.repeat)
        win.add(ui.Button('Set unit cell', self.set_unit_cell))
        self.gui = gui

    def change(self):
        repeat = [int(r.value) for r in self.repeat]
        self.gui.images.repeat_images(repeat)
        self.gui.repeat_colors(repeat)
        self.gui.set_coordinates()

    def set_unit_cell(self):
        self.gui.images.A *= self.gui.images.repeat.reshape((3, 1))
        self.gui.images.E *= self.gui.images.repeat.prod()
        self.gui.images.repeat = np.ones(3, int)
        for r in self.repeat:
            r.value = 1
        self.gui.set_coordinates()
