import ase.gui.ui as ui


class Constraints:
    def __init__(self, gui):
        win = ui.Window('Constraints')
        win.add([ui.Button('Constrain', self.selected),
                 'selected atoms'])
        win.add([ui.Button('Constrain', self.immobile),
                 'immobile atoms'])
        win.add([ui.Button('Unconstrain', self.unconstrain),
                 'selected atoms'])
        win.add(ui.Button('Clear constraints', self.clear))
        self.gui = gui

    def selected(self):
        self.gui.images.dynamic[self.gui.images.selected] = False
        self.gui.draw()

    def unconstrain(self):
        self.gui.images.dynamic[self.gui.images.selected] = True
        self.gui.draw()

    def immobile(self):
        self.gui.images.set_dynamic()
        self.gui.draw()

    def clear(self):
        self.gui.images.dynamic[:] = True
        self.gui.draw()
