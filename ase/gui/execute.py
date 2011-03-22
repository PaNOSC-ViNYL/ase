#!/usr/bin/env python
import __future__
import gtk
import os
import numpy as np
import sys

from ase.gui.languages import translate as _
from ase.gui.widgets import pack, Help

class Execute(gtk.Window):
    """ The Execute class provides an expert-user window for modification
    and evaluation of system properties with a simple one-line command structure.
    There are two types of commands, one set only applies to the global image and
    one set applies to all atoms. If the command line contains any of the atom
    commands, then it is executed separately for all atoms and for all images.
    Otherwise it is executed only once per image. 

    Please do not mix global and atom commands."""
    
    terminal_help_txt="""\
    Global commands work on all frames or only on the current frame
    - Assignment of a global variable may not reference a local one
    - use 'Current frame' switch to switch off application to all frames
    <c>e</c>:\t\ttotal energy of one frame
    <c>E</c>:\t\ttotal energy array of all frames
    <c>frame</c>:\tframe number
    <c>F</c>:\t\tall forces in one frame
    <c>fmax</c>:\tmaximal force in one frame
    <c>A</c>:\tunit cell
    <c>S</c>:\tall selected atoms (boolean array)
    <c>R</c>:\t\tall atomic positions
    <c>del S</c>:\tdelete selection
    <c>gui</c>:\tadvanced: ag window python object
    <c>img</c>:\tadvanced: ag images object
    examples: <c>frame = 1</c>, <c>A[0][1] += 4</c>, <c>e-E[-1]</c>, <c>gui.movie()</c>

    Atom commands work on each atom (or a selection) individually
    - these can use global commands on the RHS of an equation
    - use 'selected atoms only' to restrict application of command
    <c>x,y,z</c>:\tatomic coordinates
    <c>s</c>:\t\tatom is selected
    <c>f</c>:\t\tforce
    <c>Z</c>:\tatomic number
    examples: x -= A[0][0], s = (z > 5), Z = 6
    """
    
    def __init__(self, gui):
        gtk.Window.__init__(self)
        self.gui = gui
        self.set_title('Expert user mode')
        vbox = gtk.VBox()
        vbox.set_border_width(5)
        self.sw = gtk.ScrolledWindow()
        self.sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.textview = gtk.TextView()
        self.textbuffer = self.textview.get_buffer()
        self.textview.set_editable(False)
        self.textview.set_cursor_visible(False)
        self.sw.add(self.textview)
        pack(vbox, self.sw, expand=True, padding = 5)
        self.sw.set_size_request(540, 150)
        self.textview.show()
        self.add_text('Welcome to the ASE Expert user mode')
        self.cmd = gtk.Entry(60)
        self.cmd.connect('activate', self.execute)
        self.cmd.connect('key-press-event', self.update_command_buffer)
        pack(vbox, [gtk.Label('>>>'),self.cmd])
        self.cmd_buffer = getattr(gui,'expert_mode_buffer',[''])
        self.cmd_position = len(self.cmd_buffer)-1
        self.selected = gtk.CheckButton('Only selected atoms (sa)   ')
        self.selected.connect('toggled',self.selected_changed)
        self.images_only = gtk.CheckButton('Only current frame (cf)  ')
        self.images_only.connect('toggled',self.images_changed)
        pack(vbox, [self.selected, self.images_only])
        save_button = gtk.Button(stock=gtk.STOCK_SAVE)
        save_button.connect('clicked',self.save_output)
        help_button = gtk.Button(stock=gtk.STOCK_HELP)
        help_button.connect('clicked',self.terminal_help,"")
        pack(vbox, [gtk.Label('Global: Use n, N, R, A, S, D, E, frame;'
                              +' Atoms: Use a, x, y, z, f, s, Z     '),
                    help_button, save_button])
        self.add(vbox)
        vbox.show()
        self.show()
        self.cmd.grab_focus()

    def execute(self, widget=None):
        global_commands = ['n','N','R','A','S','e','E','D','F','frame']  # explicitly 'implemented' commands for use on whole system or entire single frame
        index_commands  = ['a','x','y','z','s','f','Z','d']              # commands for use on all (possibly selected) atoms

        cmd = self.cmd.get_text().strip()
        if len(cmd) == 0:
            return
        self.add_text('>>> '+cmd)
        self.cmd_buffer[-1] = cmd
        self.cmd_buffer += ['']
        setattr(self.gui,'expert_mode_buffer', self.cmd_buffer)
        self.cmd_position = len(self.cmd_buffer)-1
        self.cmd.set_text('')

        gui = self.gui
        img = gui.images

        N = img.nimages
        n = img.natoms
        S = img.selected
        D = img.dynamic[:, np.newaxis]
        E = img.E
        if self.selected.get_active():
            indices = np.where(S)[0]
        else:
            indices = range(n)

        loop_images = range(N)
        if self.images_only.get_active():
            loop_images = [self.gui.frame]

        # split off the first valid command in cmd to determine whether
        # it is global or index based
        index_based = False
        first_command = cmd.split()[0]
        for c in ['=',',','+','-','/','*',';','.','[',']','(',')','{','}']:
            if c in first_command:
                first_command = first_command[:first_command.find(c)]
        for c in index_commands:
            if c == first_command:
                index_based = True

        # check various special commands: 
        if cmd == 'del S':
            gui.delete_selected_atoms()
        elif cmd == 'sa':
            self.selected.set_active(not self.selected.get_active())
        elif cmd == 'cf':
            self.images_only.set_active(not self.images_only.get_active())
        else:
            code = compile(cmd + '\n', 'execute.py', 'single',
                           __future__.CO_FUTURE_DIVISION)
            for i in loop_images:
                R = img.P[i][indices]
                A = img.A[i]
                F = img.F[i][indices]
                e = img.E[i]
                if len(indices) > 0:
                    fmax = max(((F * D[indices])**2).sum(1)**.5)
                else:
                    fmax = None
                frame = gui.frame
            
                if not index_based:
                    try:
                        self.add_text(repr(eval(cmd)))
                    except:
                        exec code
                    gui.set_frame(frame)
                    if gui.movie_window is not None:
                        gui.movie_window.frame_number.value = frame
                else:
                    for n,a in enumerate(indices):
                        x, y, z = R[n]
                        d = D[a]
                        f = np.vdot(F[n]*d,F[n]*d)**0.5
                        s = S[a]
                        Z = img.Z[a]
                        try:
                            self.add_text(repr(eval(cmd)))
                        except:
                            exec code
                        S[a] = s
                        img.P[i][a] = x, y, z
                        img.Z[a] = Z
        gui.images.set_radii(0.89)
        gui.set_colors()
        gui.set_coordinates()

    def add_text(self,val):
        text_end = self.textbuffer.get_end_iter()
        self.textbuffer.insert(text_end,val+'\n');
        if self.sw.get_vscrollbar() is not None:
            scroll = self.sw.get_vscrollbar().get_adjustment()
            scroll.set_value(scroll.get_upper())
        
    def selected_changed(self, *args):
        if self.selected.get_active():
            self.add_text('*** Only working on selected atoms')
        else:
            self.add_text('*** Working on all atoms')

    def images_changed(self, *args):
        if self.images_only.get_active():
            self.add_text('*** Only working on current image')
        else:
            self.add_text('*** Working on all images')

    def update_command_buffer(self, entry, event, *args):
        arrow = {gtk.keysyms.Up: -1, gtk.keysyms.Down: 1}.get(event.keyval, None)
        if arrow is not None:
            self.cmd_position += arrow
            self.cmd_position = max(self.cmd_position,0)
            self.cmd_position = min(self.cmd_position,len(self.cmd_buffer)-1)
            cmd = self.cmd_buffer[self.cmd_position]
            self.cmd.set_text(cmd)
            return True
        else:
            return False

    def save_output(self, *args):
        chooser = gtk.FileChooserDialog(
            _('Save Terminal text ...'), None, gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        save = chooser.run()
        if save == gtk.RESPONSE_OK or save == gtk.RESPONSE_SAVE:
            filename = chooser.get_filename()
            text = self.textbuffer.get_text(self.textbuffer.get_start_iter(),
                                            self.textbuffer.get_end_iter())
            fd = open(filename,'w')
            fd.write(text)
            fd.close()
            chooser.destroy()
            
    def terminal_help(self,*args):
        Help(self.terminal_help_txt)

    python = execute
