import gtk
import ase.gui.ui as ui

from ase.gui.widgets import pack


class Constraints(ui.Window):
    def __init__(self, gui):
        ui.Window.__init__(self)
        self.set_title(_('Constraints'))
        vbox = ui.VBox()
        b = pack(vbox, [ui.Button(_('Constrain')),
                        ui.Label(_(' selected atoms'))])[0]
        b.connect('clicked', self.selected)
        b = pack(vbox, [ui.Button(_('Constrain')),
                        ui.Label(_(' immobile atoms:'))])[0]
        b.connect('clicked', self.immobile)
        b = pack(vbox, [ui.Button(_('Unconstrain')),
                        ui.Label(_(' selected atoms:'))])[0]
        b.connect('clicked', self.unconstrain)
        b = pack(vbox, ui.Button(_('Clear constraints')))
        b.connect('clicked', self.clear)
        close = pack(vbox, ui.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def selected(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = False
        self.gui.draw()

    def unconstrain(self, button):
        self.gui.images.dynamic[self.gui.images.selected] = True
        self.gui.draw()
        
    def immobile(self, button):
        self.gui.images.set_dynamic()
        self.gui.draw()

    def clear(self, button):
        self.gui.images.dynamic[:] = True
        self.gui.draw()

