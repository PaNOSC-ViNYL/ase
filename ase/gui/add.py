def add_atoms(self, widget, data=None, paste=None):
    """Presents a dialogbox to the user, that allows him to add
    atoms/molecule to the current slab or to paste the clipboard.

    The molecule/atom is rotated using the current rotation of the
    coordinate system.

    The molecule/atom can be added at a specified position - if the
    keyword auto+Z is used, the COM of the selected atoms will be used
    as COM for the moleculed. The COM is furthermore
    translated Z ang towards the user.

    If no molecules are selected, the COM of all the atoms will be used
    for the x-y components of the active coordinate system, while the
    z-direction will be chosen from the nearest atom position
    along this direction.

    Note: If this option is used, all frames except the active one are
    deleted.
    """

    if data == 'load':
        chooser = gtk.FileChooserDialog(
            _('Open ...'), None, gtk.FILE_CHOOSER_ACTION_OPEN,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN,
             gtk.RESPONSE_OK))

        chooser.set_filename(_("<<filename>>"))
        ok = chooser.run()
        if ok == gtk.RESPONSE_OK:
            filename = chooser.get_filename()
            chooser.destroy()
        else:
            chooser.destroy()
            return

    if data == 'OK' or data == 'load':
        import ase
        if data == 'load':
            molecule = filename
        else:
            molecule = self.add_entries[1].get_text()
        tag = self.add_entries[2].get_text()
        mom = self.add_entries[3].get_text()
        pos = self.add_entries[4].get_text().lower()

        if paste is not None:
            a = paste.copy()
        else:
            a = None

        if a is None:
            try:
                a = ase.Atoms([ase.Atom(molecule)])
            except:
                try:
                    import ase.build
                    a = ase.build.molecule(molecule)
                except:
                    try:
                        a = ase.io.read(molecule, -1)
                    except:
                        self.add_entries[1].set_text('?' + molecule)
                        return ()

        directions = np.transpose(self.axes)
        if a is not None:
            for i in a:
                try:
                    i.set('tag', int(tag))
                except:
                    self.add_entries[2].set_text('?' + tag)
                    return ()
                try:
                    i.magmom = float(mom)
                except:
                    self.add_entries[3].set_text('?' + mom)
                    return ()
            if self.origin_radio.get_active() and paste:
                a.translate(-paste.reference_position)
            # apply the current rotation matrix to A
            for i in a:
                i.position = np.dot(self.axes, i.position)
            # find the extent of the molecule in the local coordinate
            # system
            if self.centre_radio.get_active():
                a_cen_pos = np.array([0.0, 0.0, 0.0])
                m_cen_pos = 0.0
                for i in a.positions:
                    a_cen_pos[0] += np.dot(directions[0], i)
                    a_cen_pos[1] += np.dot(directions[1], i)
                    a_cen_pos[2] += np.dot(directions[2], i)
                    m_cen_pos = max(np.dot(-directions[2], i), m_cen_pos)

                a_cen_pos[0] /= len(a.positions)
                a_cen_pos[1] /= len(a.positions)
                a_cen_pos[2] /= len(a.positions)
                a_cen_pos[2] -= m_cen_pos
            else:
                a_cen_pos = np.array([0.0, 0.0, 0.0])

            # now find the position
            cen_pos = np.array([0.0, 0.0, 0.0])
            if sum(self.images.selected) > 0:
                for i in range(len(self.R)):
                    if self.images.selected[i]:
                        cen_pos += self.R[i]
                cen_pos /= sum(self.images.selected)
            elif len(self.R) > 0:
                px = 0.0
                py = 0.0
                pz = -1e6

                for i in range(len(self.R)):
                    px += np.dot(directions[0], self.R[i])
                    py += np.dot(directions[1], self.R[i])
                    pz = max(np.dot(directions[2], self.R[i]), pz)
                px = (px / float(len(self.R)))
                py = (py / float(len(self.R)))
                cen_pos = (directions[0] * px +
                           directions[1] * py +
                           directions[2] * pz)

            if 'auto' in pos:
                pos = pos.replace('auto', '')
                import re
                pos = re.sub('\s', '', pos)
                if '(' in pos:
                    sign = eval('%s1' % pos[0])
                    a_cen_pos -= sign * np.array(eval(pos[1:]), float)
                else:
                    a_cen_pos -= float(pos) * directions[2]
            else:
                cen_pos = np.array(eval(pos))
            for i in a:
                i.position += cen_pos - a_cen_pos

        # and them to the molecule
            atoms = self.images.get_atoms(self.frame)
            atoms = atoms + a
            self.new_atoms(atoms, init_magmom=True)

            # and finally select the new molecule for easy moving and
            # rotation
            for i in range(len(a)):
                self.images.selected[len(atoms) - i - 1] = True

            self.draw()
        self.add_entries[0].destroy()

    if data == 'Cancel':
        self.add_entries[0].destroy()

    if data is None or data == 'Paste':
        from ase.gui.widgets import pack
        molecule = ''
        tag = '0'
        mom = '0'
        pos = 'auto+1'
        self.add_entries = []
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.add_entries.append(window)
        window.set_title(_('Add atoms'))
        if data == 'Paste':
            molecule = paste.get_chemical_formula()
            window.set_title(_('Paste'))

        vbox = gtk.VBox(False, 0)
        window.add(vbox)
        vbox.show()
        packed = False
        for i, j in [[_('Insert atom or molecule'), molecule],
                     [_('Tag'), tag], [_('Moment'), mom], [_('Position'),
                                                           pos]]:

            label = gtk.Label(i)
            if not packed:
                vbox.pack_start(label, True, True, 0)
            else:
                packed = True
                vbox.add(label)
            label.show()

            entry = gtk.Entry()
            entry.set_text(j)
            self.add_entries.append(entry)
            entry.set_max_length(50)
            entry.show()
            vbox.add(entry)

        pack(vbox, [gtk.Label('atom/molecule reference:')])
        self.centre_radio = gtk.RadioButton(None, "centre ")
        self.origin_radio = gtk.RadioButton(self.centre_radio, "origin")
        pack(vbox, [self.centre_radio, self.origin_radio])
        if data == 'Paste':
            self.origin_radio.set_active(True)
            self.add_entries[1].set_sensitive(False)
        if data is None:
            button = gtk.Button(_('_Load molecule'))
            button.connect('clicked', self.add_atoms, 'load')
            button.show()
            vbox.add(button)
        button = gtk.Button(_('_OK'))
        button.connect('clicked', self.add_atoms, 'OK', paste)
        button.show()
        vbox.add(button)
        button = gtk.Button(_('_Cancel'))
        button.connect('clicked', self.add_atoms, 'Cancel')
        button.show()
        vbox.add(button)
        window.show()
