#!/usr/bin/env python
import gtk
from ase.gui.widgets import pack
from ase.gui.languages import translate as _

class Settings(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title('Settings')
        vbox = gtk.VBox()

        # Scale atomic radii
        pack(vbox, gtk.Label(_('Radii scaling factor:')))
        self.scale = gtk.Adjustment(value=.89, lower=0.2, upper=2.0,
                                    step_incr=0.1, page_incr=0.5)
        pack(vbox, gtk.SpinButton(self.scale, climb_rate=0, digits=2))
        self.scale.connect('value-changed', self.change)

        # Hide selected atoms
        b = pack(vbox, [gtk.Button(_('Hide')),
                        gtk.Label(_(' selected atoms'))])[0]
        b.connect('clicked', self.hide_selected)

        # View all atoms
        b = pack(vbox, gtk.Button(_('View all atoms')))
        b.connect('clicked', self.view_all)

        # A close button
        close = pack(vbox, gtk.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())

        # Add elements and show frame
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def change(self, adjustment):
        self.gui.images.set_radii(float(self.scale.value))
        self.gui.draw()
        return True

    def hide_selected(self, button):
        self.gui.images.visible = ~self.gui.images.selected
        self.gui.draw()

    def view_all(self, button):
        self.gui.images.visible[:] = True
        self.gui.draw()
