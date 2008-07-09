#!/usr/bin/env python
import gtk
from ase.gui.widgets import pack

class Settings(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.set_title('Settings')
        vbox = gtk.VBox()
        pack(vbox, gtk.Label('Radii scaling factor:'))
        self.scale = gtk.Adjustment(value=.89, lower=0.2, upper=2.0,
                                    step_incr=0.1, page_incr=0.5)
        pack(vbox, gtk.SpinButton(self.scale, climb_rate=0, digits=2))
        self.scale.connect('value-changed', self.change)
        self.add(vbox)
        vbox.show()
        self.show()
        self.set_radii = gui.images.set_radii
        self.draw = gui.draw

    def change(self, adjustment):
        self.set_radii(float(self.scale.value))
        self.draw()
        return True
