# encoding: utf-8

"Module for homogeneous deformation and calculations of elastic constants."

import gtk
from ase.gui.simulation import Simulation
from ase.gui.minimize import MinimizeMixin
from ase.gui.energyforces import OutputFieldMixin
from ase.gui.widgets import oops, pack, AseGuiCancelException
import ase
import numpy as np

class HomogeneousDeformation(Simulation, MinimizeMixin, OutputFieldMixin):
    "Window for homogeneous deformation and elastic constants."
    
    def __init__(self, gui):
        Simulation.__init__(self, gui)
        self.set_title("Homogeneous scaling")
        vbox = gtk.VBox()
        self.packtext(vbox, "XXXX Bla bla bla.")
        self.packimageselection(vbox)
        pack(vbox, gtk.Label(""))

        # Radio buttons for choosing deformation mode.
        tbl = gtk.Table(4,3)
        for i, l in enumerate(('3D', '2D', '1D')):
            l = l + " deformation   "
            lbl = gtk.Label(l)
            tbl.attach(lbl, i, i+1, 0, 1)
        self.radio_bulk = gtk.RadioButton(None, "Bulk")
        tbl.attach(self.radio_bulk, 0, 1, 1, 2)
        self.radio_xy = gtk.RadioButton(self.radio_bulk, "xy-plane")
        tbl.attach(self.radio_xy, 1, 2, 1, 2)
        self.radio_xz = gtk.RadioButton(self.radio_bulk, "xz-plane")
        tbl.attach(self.radio_xz, 1, 2, 2, 3)
        self.radio_yz = gtk.RadioButton(self.radio_bulk, "yz-plane")
        tbl.attach(self.radio_yz, 1, 2, 3, 4)
        self.radio_x = gtk.RadioButton(self.radio_bulk, "x-axis")
        tbl.attach(self.radio_x, 2, 3, 1, 2)
        self.radio_y = gtk.RadioButton(self.radio_bulk, "y-axis")
        tbl.attach(self.radio_y, 2, 3, 2, 3)
        self.radio_z = gtk.RadioButton(self.radio_bulk, "z-axis")
        tbl.attach(self.radio_z, 2, 3, 3, 4)
        tbl.show_all()
        pack(vbox, [tbl])
        self.deformtable = [
            (self.radio_bulk, (1,1,1)),
            (self.radio_xy, (1,1,0)),
            (self.radio_xz, (1,0,1)),
            (self.radio_yz, (0,1,1)),
            (self.radio_x, (1,0,0)),
            (self.radio_y, (0,1,0)),
            (self.radio_z, (0,0,1))]
        self.deform_label = gtk.Label("")
        pack(vbox, [self.deform_label])
        self.choose_possible_deformations()

        # Parameters for the deformation
        self.max_scale = gtk.Adjustment(0.010, 0.001, 10.0, 0.001)
        max_scale_spin = gtk.SpinButton(self.max_scale, 10.0, 3)
        pack(vbox, [gtk.Label("Maximal scale factor: "), max_scale_spin])
        self.scale_offset = gtk.Adjustment(0.0, -10.0, 10.0, 0.001)
        scale_offset_spin = gtk.SpinButton(self.scale_offset, 10.0, 3)
        pack(vbox, [gtk.Label("Scale offset: "), scale_offset_spin])
        self.nsteps = gtk.Adjustment(5, 3, 100, 1)
        nsteps_spin = gtk.SpinButton(self.nsteps, 1, 0)
        pack(vbox, [gtk.Label("Number of steps: "), nsteps_spin])
        pack(vbox, gtk.Label(""))
        
        # Atomic relaxations
        frame = gtk.Frame("Atomic relaxations:")
        pack(vbox, frame)
        vbox2 = gtk.VBox()
        vbox2.show()
        frame.add(vbox2)
        self.radio_relax_on = gtk.RadioButton(None, "On   ")
        self.radio_relax_off = gtk.RadioButton(self.radio_relax_on, "Off")
        self.radio_relax_off.set_active(True)
        pack(vbox2, [self.radio_relax_on, self.radio_relax_off])
        self.make_minimize_gui(vbox2)
        for r in (self.radio_relax_on, self.radio_relax_off):
            r.connect("toggled", self.relax_toggled)
        self.relax_toggled()
        pack(vbox, gtk.Label(""))
        
        # Results
        pack(vbox, [gtk.Label("Results:")])
        self.radio_results_optimal = gtk.RadioButton(
            None, "Load optimal configuration (XXX broken!)")
        self.radio_results_all =  gtk.RadioButton(
            self.radio_results_optimal, "Load all configurations")
        self.radio_results_optimal.set_active(True)
        pack(vbox, [self.radio_results_optimal])
        pack(vbox, [self.radio_results_all])

        # Output field
        self.makeoutputfield(vbox)
        pack(vbox, gtk.Label(""))

        # Status field
        self.status_label = gtk.Label("")
        pack(vbox, [self.status_label])
        
        # Run buttons etc.
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()
        self.gui.register_vulnerable(self)

    def choose_possible_deformations(self):
        """Turn on sensible radio buttons.

        Only radio buttons corresponding to deformations in directions
        with periodic boundary conditions should be turned on.
        """
        if self.setup_atoms():
            pbc = self.atoms.get_pbc()
        else:
            pbc = [False, False, False]
        for radio, requirement in self.deformtable:
            ok = True
            for i in range(3):
                if requirement[i] and not pbc[i]:
                    ok = False
            radio.set_sensitive(ok)

    def relax_toggled(self, *args):
        "Turn minimization widgets on or off."
        state = self.radio_relax_on.get_active()
        for widget in (self.algo, self.fmax_spin, self.steps_spin):
            widget.set_sensitive(state)
        
    def notify_atoms_changed(self):
        "When atoms have changed, check for the number of images."
        self.setupimageselection()
        self.choose_possible_deformations()

    def get_deformation_axes(self):
        """Return which axes the user wants to deform along."""
        for but, deform in self.deformtable:
            if but.get_active() and but.get_sensitive():
                return np.array(deform)
        # No deformation chosen!
        oops("No deformation chosen: Please choose a deformation mode.")
        return False
            
    def run(self, *args):
        """Make the deformation."""
        self.output.set_text("")
        if not self.setup_atoms():
            return
        deform_axes = self.get_deformation_axes()
        if deform_axes is False:
            return  #Nothing to do!

        # Prepare progress bar
        if self.radio_relax_on.get_active():
            fmax = self.fmax.value
            mininame = self.minimizers[self.algo.get_active()]
            self.begin(mode="scale/min", algo=mininame, fmax=fmax,
                       steps=self.steps.value, scalesteps=self.nsteps.value)
        else:
            self.begin(mode="scale", scalesteps=self.nsteps.value)
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

        # Do the scaling
        scale = self.max_scale.value
        steps = np.linspace(-scale, scale, self.nsteps.value)
        steps += self.scale_offset.value
        undef_cell = self.atoms.get_cell()
        results = []
        txt = "Strain\t\tEnergy [eV]\n"
        # If we load all configurations, prepare it.
        if self.radio_results_all.get_active():
            self.prepare_store_atoms()

        try:
            # Now, do the deformation
            for i, d in enumerate(steps):
                deformation = np.diag(1.0 + d * deform_axes)
                self.atoms.set_cell(np.dot(undef_cell, deformation),
                                    scale_atoms=True)
                if self.gui.simulation.has_key('progress'):
                    self.gui.simulation['progress'].set_scale_progress(i)
                if self.radio_relax_on.get_active():
                    algo = getattr(ase.optimize, mininame)
                    minimizer = algo(self.atoms, logfile=logger)
                    minimizer.run(fmax=fmax, steps=self.steps.value)
                e = self.atoms.get_potential_energy()
                results.append((d, e))
                txt = txt + ("%.3f\t\t%.3f\n" % (d, e))
                self.output.set_text(txt)
                if self.radio_results_all.get_active():
                    self.store_atoms()
        except AseGuiCancelException:
            # Update display to reflect cancellation of simulation.
            self.status_label.set_text("Cancellation CANCELLED.")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
        except MemoryError:
            self.status_label.set_text("Out of memory, consider using LBFGS instead")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
            
        else:
            # Update display to reflect succesful end of simulation.
            self.status_label.set_text("Calculation completed.")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#007700'))
                     
        self.end()    
        if results:
            self.activate_output()
            self.gui.notify_vulnerable()
            
        # If we store all configurations: Open movie window and energy graph
        if self.gui.images.nimages > 1:
            self.gui.movie()
            assert not np.isnan(self.gui.images.E[0])
            if not self.gui.plot_graphs_newatoms():
                expr = 'i, e - E[-1]'            
                self.gui.plot_graphs(expr=expr)

        
        
    
