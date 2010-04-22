# encoding: utf-8

"Module for performing energy minimization."

import gtk
from ase.gui.simulation import Simulation
from ase.gui.widgets import oops, pack, AseGuiCancelException
import ase
import numpy as np

class Minimize(Simulation):
    "Window for performing energy minimization."
    minimizers = ('BFGS', 'BFGSLineSearch', 'LBFGS', 'LBFGSLineSearch', 'MDMin', 'FIRE')
    
    def __init__(self, gui):
        Simulation.__init__(self, gui)
        self.set_title("Energy minimization")
        
        vbox = gtk.VBox()
        self.packtext(vbox,
                      "Minimize the energy with respect to the positions.")
        self.packimageselection(vbox)
        pack(vbox, gtk.Label(""))
        self.algo = gtk.combo_box_new_text()
        for m in self.minimizers:
            self.algo.append_text(m)
        self.algo.set_active(0)
        pack(vbox, [gtk.Label("Algorithm: "), self.algo])
        
        self.fmax = gtk.Adjustment(0.05, 0.00, 10.0, 0.01)
        self.fmax_spin = gtk.SpinButton(self.fmax, 0, 3)
        lbl = gtk.Label()
        lbl.set_markup("Convergence criterion: F<sub>max</sub> = ")
        pack(vbox, [lbl, self.fmax_spin])

        self.steps = gtk.Adjustment(100, 1, 1000000, 1)
        self.steps_spin = gtk.SpinButton(self.steps, 0, 0)
        pack(vbox, [gtk.Label("Max. number of steps: "), self.steps_spin])
        
        pack(vbox, gtk.Label(""))
        self.status_label = gtk.Label("")
        pack(vbox, [self.status_label])
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()
        self.gui.register_vulnerable(self)

    def run(self, *args):
        "User has pressed [Run]: run the minimization."
        if not self.setup_atoms():
            return
        fmax = self.fmax.value
        steps = self.steps.value
        mininame = self.minimizers[self.algo.get_active()]
        self.begin(mode="min", algo=mininame, fmax=fmax, steps=steps)
        algo = getattr(ase, mininame)
        try:
            logger_func = self.gui.simulation['progress'].get_logger_stream
        except (KeyError, AttributeError):
            logger = None
        else:
            logger = logger_func()  # Don't catch errors in the function.

        # Display status message
        self.status_label.set_text("Running ...")
        self.status_label.modify_fg(gtk.STATE_NORMAL,
                                    gtk.gdk.color_parse('#AA0000'))
        while gtk.events_pending():
            gtk.main_iteration()

        self.prepare_store_atoms()
        minimizer = algo(self.atoms, logfile=logger)
        minimizer.attach(self.store_atoms)
        try:
            minimizer.run(fmax=fmax, steps=steps)
        except AseGuiCancelException:
            # Update display to reflect cancellation of simulation.
            self.status_label.set_text("Minimization CANCELLED after %i steps."
                                       % (self.count_steps,))
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
        except MemoryError:
            self.status_label.set_text("Out of memory, consider using LBFGS instead")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
            
        else:
            # Update display to reflect succesful end of simulation.
            self.status_label.set_text("Minimization completed in %i steps."
                                       % (self.count_steps,))
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#007700'))
            
        self.end()
        if self.count_steps:
            # Notify other windows that atoms have changed.
            # This also notifies this window!
            self.gui.notify_vulnerable()

        # Open movie window and energy graph
        if self.gui.images.nimages > 1:
            self.gui.movie()
            assert not np.isnan(self.gui.images.E[0])
            if not self.gui.plot_graphs_newatoms():
                expr = 'i, e - E[-1]'            
                self.gui.plot_graphs(expr=expr)

    def notify_atoms_changed(self):
        "When atoms have changed, check for the number of images."
        self.setupimageselection()
        
