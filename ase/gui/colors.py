from __future__ import print_function
# encoding: utf-8
"""colors.py - select how to color the atoms in the GUI."""


import gtk
from gettext import gettext as _
from ase.gui.widgets import pack, cancel_apply_ok, oops
import ase
from ase.data.colors import jmol_colors
import numpy as np
import colorsys

named_colors = ('Green', 'Yellow', 'Blue', 'Red', 'Orange', 'Cyan',
                'Magenta', 'Black', 'White', 'Grey', 'Violet', 'Brown',
                'Navy')


class ColorWindow(gtk.Window):
    "A window for selecting how to color the atoms."
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui
        self.colormode = gui.colormode
        self.actual_colordata = None
        self.colordata = {}
        self.set_title(_("Colors"))
        vbox = gtk.VBox()
        self.add(vbox)
        vbox.show()
        # The main layout consists of two columns, the leftmost split in an upper and lower part.
        self.maintable = gtk.Table(2,2)
        pack(vbox, self.maintable)
        self.methodbox = gtk.VBox()
        self.methodbox.show()
        self.maintable.attach(self.methodbox, 0, 1, 0, 1)
        self.scalebox = gtk.VBox()
        self.scalebox.show()
        self.maintable.attach(self.scalebox, 0, 1, 1, 2)
        self.colorbox = gtk.Frame()
        self.colorbox.show()
        self.maintable.attach(self.colorbox, 1, 2, 0, 2, gtk.EXPAND)
        # Upper left: Choose how the atoms are colored.
        lbl = gtk.Label(_("Choose how the atoms are colored:"))
        pack(self.methodbox, [lbl])
        self.radio_jmol = gtk.RadioButton(None, _('By atomic number, default "jmol" colors'))
        self.radio_atno = gtk.RadioButton(self.radio_jmol,
                                          _('By atomic number, user specified'))
        self.radio_tag = gtk.RadioButton(self.radio_jmol, _('By tag'))
        self.radio_force = gtk.RadioButton(self.radio_jmol, _('By force'))
        self.radio_velocity = gtk.RadioButton(self.radio_jmol, _('By velocity'))
        self.radio_charge = gtk.RadioButton(self.radio_jmol, _('By charge'))
        self.radio_magnetic_moment = gtk.RadioButton(
            self.radio_jmol, _('By magnetic moment'))
        self.radio_coordination = gtk.RadioButton(
            self.radio_jmol, _('By coordination'))
        self.radio_manual = gtk.RadioButton(self.radio_jmol, _('Manually specified'))
        self.radio_same = gtk.RadioButton(self.radio_jmol, _('All the same color'))
        self.force_box = gtk.VBox()
        self.velocity_box = gtk.VBox()
        self.charge_box = gtk.VBox()
        self.magnetic_moment_box = gtk.VBox()
        for widget in (self.radio_jmol, self.radio_atno, self.radio_tag,
                       self.radio_force, self.force_box, 
                       self.radio_velocity, self.velocity_box, 
                       self.radio_charge, self.charge_box,
                       self.radio_magnetic_moment,
                       self.magnetic_moment_box,
                       self.radio_coordination,
                       self.radio_manual, self.radio_same):
            pack(self.methodbox, [widget])
            if isinstance(widget, gtk.RadioButton):
                widget.connect('toggled', self.method_radio_changed)
        # Now fill in the box for additional information in case the force is used.
        self.force_label = gtk.Label(_("This should not be displayed in forces!"))
        pack(self.force_box, [self.force_label])
        self.min = gtk.Adjustment(0.0, 0.0, 100.0, 0.05)
        self.max = gtk.Adjustment(0.0, 0.0, 100.0, 0.05)
        self.steps = gtk.Adjustment(10, 2, 500, 1)
        force_apply = gtk.Button(_('Update'))
        force_apply.connect('clicked', self.set_min_max_colors, 'force')
        pack(self.force_box, [gtk.Label(_('Min: ')),
                              gtk.SpinButton(self.min, 1.0, 2),
                              gtk.Label(_('  Max: ')),
                              gtk.SpinButton(self.max, 1.0, 2),
                              gtk.Label(_('  Steps: ')),
                              gtk.SpinButton(self.steps, 1, 0),
                              gtk.Label('  '),
                              force_apply])
        self.force_box.hide()
        # Now fill in the box for additional information in case the velocity is used.
        self.velocity_label = gtk.Label("This should not be displayed!")
        pack(self.velocity_box, [self.velocity_label])
        velocity_apply = gtk.Button(_('Update'))
        velocity_apply.connect('clicked', self.set_min_max_colors, 'velocity')
        pack(self.velocity_box, [gtk.Label(_('Min: ')),
                                 gtk.SpinButton(self.min, 1.0, 3),
                                 gtk.Label(_('  Max: ')),
                                 gtk.SpinButton(self.max, 1.0, 3),
                                 gtk.Label(_('  Steps: ')),
                                 gtk.SpinButton(self.steps, 1, 0),
                                 gtk.Label('  '),
                                 velocity_apply])
        self.velocity_box.hide()
        # Now fill in the box for additional information in case
        # the charge is used.
        self.charge_label = gtk.Label(_("This should not be displayed!"))
        pack(self.charge_box, [self.charge_label])
        charge_apply = gtk.Button(_('Update'))
        charge_apply.connect('clicked', self.set_min_max_colors, 'charge')
        pack(self.charge_box, [gtk.Label(_('Min: ')),
                              gtk.SpinButton(self.min, 10.0, 2),
                              gtk.Label(_('  Max: ')),
                              gtk.SpinButton(self.max, 10.0, 2),
                              gtk.Label(_('  Steps: ')),
                              gtk.SpinButton(self.steps, 1, 0),
                              gtk.Label('  '),
                              charge_apply])
        self.charge_box.hide()
        # Now fill in the box for additional information in case
        # the magnetic moment is used.
        self.magnetic_moment_label = gtk.Label(_(
            "This should not be displayed!"))
        pack(self.magnetic_moment_box, [self.magnetic_moment_label])
        magnetic_moment_apply = gtk.Button(_('Update'))
        magnetic_moment_apply.connect('clicked', self.set_min_max_colors,
                                      'magnetic moment')
        pack(self.magnetic_moment_box, [gtk.Label(_('Min: ')),
                                        gtk.SpinButton(self.min, 10.0, 2),
                                        gtk.Label(_('  Max: ')),
                                        gtk.SpinButton(self.max, 10.0, 2),
                                        gtk.Label(_('  Steps: ')),
                                        gtk.SpinButton(self.steps, 1, 0),
                                        gtk.Label('  '),
                                        magnetic_moment_apply])
        self.magnetic_moment_box.hide()
        # Lower left: Create a color scale
        pack(self.scalebox, gtk.Label(""))
        lbl = gtk.Label(_('Create a color scale:'))
        pack(self.scalebox, [lbl])
        color_scales = (
            _('Black - white'),
            _('Black - red - yellow - white'),
            _('Black - green - white'),
            _('Black - blue - cyan'),
            _('Blue - white - red'),
            _('Hue'),
            _('Named colors')
            )
        self.scaletype_created = None
        self.scaletype = gtk.combo_box_new_text()
        for s in color_scales:
            self.scaletype.append_text(s)
        self.createscale = gtk.Button(_("Create"))
        pack(self.scalebox, [self.scaletype, self.createscale])
        self.createscale.connect('clicked', self.create_color_scale)
        # The actually colors are specified in a box possibly with scrollbars
        self.colorwin = gtk.ScrolledWindow()
        self.colorwin.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        self.colorwin.show()
        self.colorbox.add(self.colorwin)
        self.colorwin.add_with_viewport(gtk.VBox()) # Dummy contents
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [buts], end=True, bottom=True)
        # Make the initial setup of the colors
        self.color_errors = {}
        self.init_colors_from_gui()
        self.show()
        gui.register_vulnerable(self)

    def notify_atoms_changed(self):
        "Called by gui object when the atoms have changed."
        self.destroy()
  
    def init_colors_from_gui(self):
        cm = self.gui.colormode
        # Disallow methods if corresponding data is not available
        if not self.gui.images.T.any():
            self.radio_tag.set_sensitive(False)
            if self.radio_tag.get_active() or cm == 'tag':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_tag.set_sensitive(True)
        if np.isnan(self.gui.images.F).any() or not self.gui.images.F.any():
            self.radio_force.set_sensitive(False)
            if self.radio_force.get_active() or cm == 'force':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_force.set_sensitive(True)
        if np.isnan(self.gui.images.V).any() or not self.gui.images.V.any():
            self.radio_velocity.set_sensitive(False)
            if self.radio_velocity.get_active() or cm == 'velocity':
                self.radio_jmol.set_active(True)
                return
        else:
            self.radio_velocity.set_sensitive(True)
        if not self.gui.images.q.any():
            self.radio_charge.set_sensitive(False)
        else:
            self.radio_charge.set_sensitive(True)
        if not self.gui.images.M.any():
            self.radio_magnetic_moment.set_sensitive(False)
        else:
            self.radio_magnetic_moment.set_sensitive(True)
        self.radio_manual.set_sensitive(self.gui.images.natoms <= 1000)
        # Now check what the current color mode is
        if cm == 'jmol':
            self.radio_jmol.set_active(True)
            self.set_jmol_colors()
        elif cm == 'atno':
            self.radio_atno.set_active(True)
        elif cm == 'tags':
            self.radio_tag.set_active(True)
        elif cm == 'force':
            self.radio_force.set_active(True)
        elif cm == 'velocity':
            self.radio_velocity.set_active(True)
        elif cm == 'charge':
            self.radio_charge.set_active(True)
        elif cm == 'magnetic moment':
            self.radio_magnetic_moment.set_active(True)
        elif cm == 'coordination':
            self.radio_coordination.set_active(True)
        elif cm == 'manual':
            self.radio_manual.set_active(True)
        elif cm == 'same':
            self.radio_same.set_active(True)
            
    def method_radio_changed(self, widget=None):
        "Called when a radio button is changed."
        self.scaletype_created = None
        self.scaletype.set_active(-1)
        if not widget.get_active():
            # Ignore most events when a button is turned off.
            if widget is self.radio_force:
                self.force_box.hide()
            if widget is self.radio_velocity:
                self.velocity_box.hide()
            if widget is self.radio_charge:
                self.charge_box.hide()
            return
        if widget is self.radio_jmol:
            self.set_jmol_colors()
        elif widget is self.radio_atno:
            self.set_atno_colors()
        elif widget is self.radio_tag:
            self.set_tag_colors()
        elif widget is self.radio_force:
            self.show_force_stuff()
            self.set_min_max_colors(None, 'force')
        elif widget is self.radio_velocity:
            self.show_velocity_stuff()
            self.set_min_max_colors(None, 'velocity')
        elif widget is self.radio_charge:
            self.show_charge_stuff()
            self.set_min_max_colors(None, 'charge')
        elif widget is self.radio_magnetic_moment:
            self.show_magnetic_moment_stuff()
            self.set_min_max_colors(None, 'magnetic moment')
        elif widget is self.radio_coordination:
            self.set_coordination_colors()
        elif widget is self.radio_manual:
            self.set_manual_colors()
        elif widget is self.radio_same:
            self.set_same_color()
        else:
            raise RuntimeError('Unknown widget in method_radio_changed')
            
    def make_jmol_colors(self):
        "Set the colors to the default jmol colors"
        self.colordata_z = []
        hasfound = {}
        for z in self.gui.images.Z:
            if z not in hasfound:
                hasfound[z] = True
                self.colordata_z.append([z, jmol_colors[z]])

    def set_jmol_colors(self):
        "We use the immutable jmol colors."
        self.make_jmol_colors()
        self.set_atno_colors()
        for entry in self.color_entries:
            entry.set_sensitive(False)
        self.colormode = 'jmol'
        
    def set_atno_colors(self):
        "We use user-specified per-element colors."
        if not hasattr(self, 'colordata_z'):
            # No initial colors.  Use jmol colors
            self.make_jmol_colors()
        self.actual_colordata = self.colordata_z
        self.color_labels = ["%i (%s):" % (z, ase.data.chemical_symbols[z])
                             for z, col in self.colordata_z]
        self.make_colorwin()
        self.colormode = 'atno'

    def set_tag_colors(self):
        "We use per-tag colors."
        # Find which tags are in use
        tags = self.gui.images.T
        existingtags = range(tags.min(), tags.max() + 1)
        if not hasattr(self, 'colordata_tags') or len(self.colordata_tags) != len(existingtags):
            colors = self.get_named_colors(len(existingtags))
            self.colordata_tags = [[x, y] for x, y in
                                   zip(existingtags, colors)]
        self.actual_colordata = self.colordata_tags
        self.color_labels = [str(x)+':' for x, y in self.colordata_tags]
        self.make_colorwin()
        self.colormode = 'tags'

    def set_same_color(self):
        "All atoms have the same color"
        if not hasattr(self, 'colordata_same'):
            try:
                self.colordata_same = self.actual_colordata[0:1]
            except AttributeError:
                self.colordata_same = self.get_named_colors(1)
        self.actual_colordata = self.colordata_same
        self.actual_colordata[0][0] = 0
        self.color_labels = ['all:']
        self.make_colorwin()
        self.colormode = 'same'

    def set_min_max_colors(self, widget, mode):
        borders = np.linspace(self.min.value, self.max.value, self.steps.value,
                              endpoint=False)
        if self.scaletype_created is None:
            colors = self.new_color_scale([[0, [1,1,1]],
                                           [1, [0,0,1]]], len(borders))
        elif (mode not in  self.colordata or
              len(self.colordata[mode]) != len(borders)):
            colors = self.get_color_scale(len(borders), self.scaletype_created)
        else:
            colors = [y for x, y in self.colordata[mode]]
        self.colordata[mode] = [[x, y] for x, y in zip(borders, colors)]
        self.actual_colordata = self.colordata[mode]
        self.color_labels = ["%.2f:" % x for x, y in self.colordata[mode]]
        self.make_colorwin()
        self.colormode = mode
        factor = self.steps.value / (self.max.value - self.min.value)
        self.colormode_data = (self.min.value, factor)

    def set_coordination_colors(self, *args):
        "Use coordination as basis for the colors."
        if not hasattr(self.gui, 'coordination'):
            self.gui.toggle_show_bonds(None)
        if not hasattr(self, 'colordata_coordination'):
            colors = self.get_named_colors(len(named_colors))
            self.colordata_coordination = [[x, y] for x, y in
                                           enumerate(colors)]
        self.actual_colordata = self.colordata_coordination
        self.color_labels = [(str(x) + ':')
                             for x, y in self.colordata_coordination]
        self.make_colorwin()
        self.colormode = 'coordination'

    def set_manual_colors(self):
        "Set colors of all atoms from the last selection."
        # We cannot directly make np.arrays of the colors, as they may
        # be sequences of the same length, causing creation of a 2D
        # array of characters/numbers instead of a 1D array of
        # objects.
        colors = np.array([None] * self.gui.images.natoms)
        if self.colormode in ['atno', 'jmol', 'tags']:
            maxval = max([x for x, y in self.actual_colordata])
            oldcolors = np.array([None] * (maxval+1))
            for x, y in self.actual_colordata:
                oldcolors[x] = y
            if self.colormode == 'tags':
                colors[:] = oldcolors[self.gui.images.T[self.gui.frame]]
            else:
                colors[:] = oldcolors[self.gui.images.Z]
        elif self.colormode == 'force':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            F = self.gui.images.F[self.gui.frame]
            F = np.sqrt((F * F).sum(axis=-1))
            nF = (F - self.colormode_force_data[0]) * self.colormode_force_data[1]
            nF = np.clip(nF.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nF]
        elif self.colormode == 'velocity':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            V = self.gui.images.V[self.gui.frame]
            V = np.sqrt((V * V).sum(axis=-1))
            nV = (V - self.colormode_velocity_data[0]) * self.colormode_velocity_data[1]
            nV = np.clip(nV.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nV]
        elif self.colormode == 'charge':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            q = self.gui.images.q[self.gui.frame]
            nq = ((q - self.colormode_charge_data[0]) *
                  self.colormode_charge_data[1])
            nq = np.clip(nq.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nq]
        elif self.colormode == 'magnetic moment':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
            q = self.gui.images.q[self.gui.frame]
            nq = ((q - self.colormode_magnetic_moment_data[0]) *
                  self.colormode_magnetic_moment_data[1])
            nq = np.clip(nq.astype(int), 0, len(oldcolors)-1)
            colors[:] = oldcolors[nq]
        elif self.colormode == 'coordination':
            oldcolors = np.array([None] * len(self.actual_colordata))
            oldcolors[:] = [y for x, y in self.actual_colordata]
        elif self.colormode == 'same':
            oldcolor = self.actual_colordata[0][1]
            if len(colors) == len(oldcolor):
                # Direct assignment would be e.g. one letter per atom. :-(
                colors[:] = [oldcolor] * len(colors)
            else:
                colors[:] = oldcolor
        elif self.colormode == 'manual':
            if self.actual_colordata is None:   # import colors from gui, if they don't exist already
                colors = [y for x,y in self.gui.colordata]

        self.color_labels = ["%d:" % i for i in range(len(colors))]
        self.actual_colordata = [[i, x] for i, x in enumerate(colors)]
        self.make_colorwin()
        self.colormode = 'manual'

    def get_min_max_text(self, mode, vmin, vmax, min_frame, max_frame):
        nimages = self.gui.images.nimages
        txt = 'Max {0}: {1:.2f}'.format(mode, vmax)
        if nimages > 1:
            txt += '(all frames), {0:.2f} (this frame)'.format(max_frame)
        self.max.value = vmax
        if vmin is None:
            self.min.value = 0.
        else:
            txt += ', Min {0:.2f}'.format(vmin)
            if nimages > 1:
                txt += '(all frames), {0:.2f} (this frame)'.format(min_frame)
            self.min.value = vmin
        return txt

    def show_force_stuff(self):
        "Show and update widgets needed for selecting the force scale."
        self.force_box.show()
        # XXX is this projected on some axis ? XXX
        F = np.sqrt(((self.gui.images.F *
                      self.gui.images.dynamic[:,np.newaxis])**2).sum(axis=-1))
        txt = self.get_min_max_text(
            'force', None, F.max(),
            None, self.gui.images.F[self.gui.frame].max())
        self.force_label.set_text(_(txt))

    def show_velocity_stuff(self):
        "Show and update widgets needed for selecting the velocity scale."
        self.velocity_box.show()
        V = np.sqrt((self.gui.images.V * self.gui.images.V).sum(axis=-1))
        Vframe = np.sqrt((self.gui.images.V[self.gui.frame] *
                          self.gui.images.V[self.gui.frame]).sum(axis=-1))
        txt = self.get_min_max_text(
            'velocity', None, V.max(), None, Vframe.max())
        self.velocity_label.set_text(_(txt))

    def show_charge_stuff(self):
        "Show and update widgets needed for selecting the charge scale."
        self.charge_box.show()
        txt = self.get_min_max_text(
            'charge', self.gui.images.q.min(), self.gui.images.q.max(),
            self.gui.images.q[self.gui.frame].min(),
            self.gui.images.q[self.gui.frame].max())
        self.charge_label.set_text(_(txt))

    def show_magnetic_moment_stuff(self):
        "Show and update widgets needed for selecting the magn. mom. scale."
        self.magnetic_moment_box.show()
        txt = self.get_min_max_text(
            'magnetic moment', self.gui.images.M.min(), self.gui.images.M.max(),
            self.gui.images.M[self.gui.frame].min(),
            self.gui.images.M[self.gui.frame].max())
        self.magnetic_moment_label.set_text(_(txt))

    def make_colorwin(self):
        """Make the list of editable color entries.

        Uses self.actual_colordata and self.color_labels.  Produces self.color_entries.
        """
        assert len(self.actual_colordata) == len(self.color_labels)
        self.color_entries = []
        old = self.colorwin.get_child()
        self.colorwin.remove(old)
        del old
        table = gtk.Table(len(self.actual_colordata)+1, 4)
        self.colorwin.add_with_viewport(table)
        table.show()
        self.color_display = []
        for i in range(len(self.actual_colordata)):
            lbl = gtk.Label(self.color_labels[i])
            entry = gtk.Entry(max=20)
            val = self.actual_colordata[i][1]
            error = False
            if not isinstance(val, str):
                assert len(val) == 3
                intval = tuple(np.round(65535*np.array(val)).astype(int))
                val = "%.3f, %.3f, %.3f" % tuple(val)
                clr = gtk.gdk.Color(*intval)
            else:
                try:
                    clr = gtk.gdk.color_parse(val)
                except ValueError:
                    error = True
            entry.set_text(val)
            blob = gtk.EventBox()
            space = gtk.Label
            space = gtk.Label("    ")
            space.show()
            blob.add(space)
            if error:
                space.set_text(_("ERROR"))
            else:
                blob.modify_bg(gtk.STATE_NORMAL, clr)
            table.attach(lbl, 0, 1, i, i+1, yoptions=0)
            table.attach(entry, 1, 2, i, i+1, yoptions=0)
            table.attach(blob, 2, 3, i, i+1, yoptions=0)
            lbl.show()
            entry.show()
            blob.show()
            entry.connect('changed', self.entry_changed, i)
            self.color_display.append(blob)
            self.color_entries.append(entry)
            
    def entry_changed(self, widget, index):
        """The user has changed a color."""
        txt = widget.get_text()
        txtfields = txt.split(',')
        if len(txtfields) == 3:
            self.actual_colordata[index][1] = [float(x) for x in txtfields]
            val = tuple([int(65535*float(x)) for x in txtfields])
            clr = gtk.gdk.Color(*val)
        else:
            self.actual_colordata[index][1] = txt
            try:
                clr = gtk.gdk.color_parse(txt)
            except ValueError:
                # Cannot parse the color
                displ = self.color_display[index]
                displ.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse('white'))
                displ.get_child().set_text(_("ERR"))
                self.color_errors[index] = (self.color_labels[index], txt)
                return
        self.color_display[index].get_child().set_text("    ") # Clear error message
        self.color_errors.pop(index, None)
        self.color_display[index].modify_bg(gtk.STATE_NORMAL, clr)
        
    def create_color_scale(self, *args):
        if self.radio_jmol.get_active():
            self.radio_atno.set_active(1)
        n = len(self.color_entries)
        s = self.scaletype.get_active()
        scale = self.get_color_scale(n, s)
        self.scaletype_created = s
        for i in range(n):
            if isinstance(scale[i], str):
                self.color_entries[i].set_text(scale[i])
            else:
                s = "%.3f, %.3f, %.3f" % tuple(scale[i])
                self.color_entries[i].set_text(s)
            self.color_entries[i].activate()

    def get_color_scale(self, n, s):
        if s == 0:
            # Black - White
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [1, [1,1,1]]], n)
        elif s == 1:
            # Black - Red - Yellow - White (STM colors)
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.33, [1,0,0]],
                                          [0.67, [1,1,0]],
                                          [1, [1,1,1]]], n)
        elif s == 2:
            # Black - Green - White
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.5, [0,0.9,0]],
                                          [0.75, [0.2,1.0,0.2]],
                                          [1, [1,1,1]]], n)
        elif s == 3:
            # Black - Blue - Cyan
            scale = self.new_color_scale([[0, [0,0,0]],
                                          [0.5, [0,0,1]],
                                          [1, [0,1,1]]], n)
        elif s == 4:
            # Blue - White - Red
             scale = self.new_color_scale([[0, [0,0,1]],
                                          [0.5, [1,1,1]],
                                          [2, [1,0,0]]], n)
        elif s == 5:
            # Hues
            hues = np.linspace(0.0, 1.0, n, endpoint=False)
            scale = ["%.3f, %.3f, %.3f" % colorsys.hls_to_rgb(h, 0.5, 1)
                     for h in hues]
        elif s == 6:
            # Named colors
            scale = self.get_named_colors(n)
        else:
            scale = None
        return scale

    def new_color_scale(self, fixpoints, n):
        "Create a homogeneous color scale."
        x = np.array([a[0] for a in fixpoints], float)
        y = np.array([a[1] for a in fixpoints], float)
        assert y.shape[1] == 3
        res = []
        for a in np.linspace(0.0, 1.0, n, endpoint=True):
            n = x.searchsorted(a)
            if n == 0:
                v = y[0]  # Before the start
            elif n == len(x):
                v = x[-1] # After the end
            else:
                x0 = x[n-1]
                x1 = x[n]
                y0 = y[n-1]
                y1 = y[n]
                v = y0 + (y1 - y0) / (x1 - x0) * (a - x0)
            res.append(v)
        return res

    def get_named_colors(self, n):
        if n <= len(named_colors):
            return named_colors[:n]
        else:
            return named_colors + ('Black',) * (n - len(named_colors))
        
    def apply(self, *args):
        #if self.colormode in ['atno', 'jmol', 'tags']:
        # Color atoms according to an integer value number
        if self.color_errors:
            oops(_("Incorrect color specification"),
                 "%s: %s" % self.color_errors.values()[0])
            return False
        colordata = self.actual_colordata
        if self.colormode in [
                'force', 'velocity', 'charge', 'magnetic moment']:
            # Use integers instead for border values
            colordata = [[i, x[1]] for i, x in enumerate(self.actual_colordata)]
            self.gui.colormode_data = self.colormode_data
        maxval = max([x for x, y in colordata])
        self.gui.colors = [None] * (maxval + 1)
        new = self.gui.drawing_area.window.new_gc
        alloc = self.gui.colormap.alloc_color
        for z, val in colordata:
            if isinstance(val, str):
                self.gui.colors[z] = new(alloc(val))
            else:
                clr = tuple([int(65535*x) for x in val])
                assert len(clr) == 3
                self.gui.colors[z] = new(alloc(*clr))
        self.gui.colormode = self.colormode
        self.gui.colordata = colordata
        self.gui.draw()
        return True

    def cancel(self, *args):
        self.destroy()

    def ok(self, *args):
        if self.apply():
            self.destroy()
        
