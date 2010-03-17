# encoding: utf-8

"Module for calculating energies and forces."

import gtk
from ase.gui.simulation import Simulation
from ase.gui.widgets import oops, pack

class EnergyForces(Simulation):
    def __init__(self, gui):
        Simulation.__init__(self, gui)

        self.set_default_size(-1, 300)
        vbox = gtk.VBox()
        self.packtext(vbox,
                      "Calculate potential energy and the force on all atoms")
        self.forces = gtk.CheckButton("Write forces on the atoms")
        self.forces.set_active(True)
        pack(vbox, [self.forces])
        pack(vbox, [gtk.Label("")])
        pack(vbox, [gtk.Label("Output:")])
        scrwin = gtk.ScrolledWindow()
        scrwin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.output = gtk.TextBuffer()
        txtview = gtk.TextView(self.output)
        txtview.set_editable(False)
        scrwin.add(txtview)
        scrwin.show_all()
        vbox.pack_start(scrwin, True, True, 0)
        pack(vbox, gtk.Label(""))
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()

    def run(self, *args):
        if not self.setup_atoms():
            return
        e = self.atoms.get_potential_energy()
        txt = "<b>Potential Energy:</b>\n"
        txt += "  %8.3f eV\n\n" % (e,)
        if self.forces.get_active():
            txt +="<b>Forces:</b>\n"
            forces = self.atoms.get_forces()
            for f in forces:
                txt += "  %8.3f, %8.3f, %8.3f eV/Å\n" % tuple(f)
        self.output.set_text(txt)
        print "Done"
        
