def modify_atoms(self, widget, data=None):
    """Presents a dialog box where the user is able to change the
    atomic type, the magnetic moment and tags of the selected atoms.
    An item marked with X will not be changed.
    """
    if data:
        if data == 'OK':
            import ase
            symbol = self.add_entries[1].get_text()
            tag = self.add_entries[2].get_text()
            mom = self.add_entries[3].get_text()
            a = None
            if symbol != 'X':
                try:
                    a = ase.Atoms([ase.Atom(symbol)])
                except:
                    self.add_entries[1].set_text('?' + symbol)
                    return ()

            y = self.images.selected.copy()
            # and them to the molecule
            atoms = self.images.get_atoms(self.frame)
            for i in range(len(atoms)):
                if self.images.selected[i]:
                    if a:
                        atoms[i].symbol = symbol
                    try:
                        if tag != 'X':
                            atoms[i].tag = int(tag)
                    except:
                        self.add_entries[2].set_text('?' + tag)
                        return ()
                    try:
                        if mom != 'X':
                            atoms[i].magmom = float(mom)
                    except:
                        self.add_entries[3].set_text('?' + mom)
                        return ()
            self.new_atoms(atoms, init_magmom=True)

            # Updates atomic labels
            a = self.ui.get_action_groups()[0].get_action("NoLabel")
            cv = a.get_current_value()
            a.set_current_value(0)
            a.set_current_value(cv)

            # and finally select the new molecule for easy moving
            # and rotation
            self.images.selected = y
            self.draw()

        self.add_entries[0].destroy()
    if data is None and sum(self.images.selected):
        atoms = self.images.get_atoms(self.frame)
        s_tag = ''
        s_mom = ''
        s_symbol = ''
        # Get the tags, moments and symbols of the selected atomsa
        for i in range(len(atoms)):
            if self.images.selected[i]:
                if not (s_tag):
                    s_tag = str(atoms[i].tag)
                elif s_tag != str(atoms[i].tag):
                    s_tag = 'X'
                if not (s_mom):
                    s_mom = ("%2.2f" % (atoms[i].magmom))
                elif s_mom != ("%2.2f" % (atoms[i].magmom)):
                    s_mom = 'X'
                if not (s_symbol):
                    s_symbol = str(atoms[i].symbol)
                elif s_symbol != str(atoms[i].symbol):
                    s_symbol = 'X'

        self.add_entries = []
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.add_entries.append(window)
        window.set_title(_('Modify'))

        vbox = gtk.VBox(False, 0)
        window.add(vbox)
        vbox.show()
        pack = False
        for i, j in [[_('Atom'), s_symbol], [_('Tag'), s_tag],
                     [_('Moment'), s_mom]]:
            label = gtk.Label(i)
            if not pack:
                vbox.pack_start(label, True, True, 0)
            else:
                pack = True
                vbox.add(label)
            label.show()

            entry = gtk.Entry()
            entry.set_text(j)
            self.add_entries.append(entry)
            entry.set_max_length(50)
            entry.show()
            vbox.add(entry)
        button = gtk.Button(_('_OK'))
        button.connect('clicked', self.modify_atoms, 'OK')
        button.show()
        vbox.add(button)
        button = gtk.Button(_('_Cancel'))
        button.connect('clicked', self.modify_atoms, 'Cancel')
        button.show()
        vbox.add(button)

        window.show()
