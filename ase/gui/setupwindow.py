    def update_element(self, *args):
        "Called when a new element may have been entered."
        # Assumes the element widget is self.element and that a label
        # to keep updated is self.elementinfo.  The chemical symbol is
        # placed in self.legalelement - or None if the element is
        # invalid.
        elem = self.element.get_text()
        if not elem:
            self.invalid_element(_("  No element specified!"))
            return False
        try:
            z = int(elem)
        except ValueError:
            # Probably a symbol
            try:
                z = ase.data.atomic_numbers[elem]
            except KeyError:
                self.invalid_element()
                return False
        try:
            symb = ase.data.chemical_symbols[z]
        except KeyError:
            self.invalid_element()
            return False
        name = ase.data.atomic_names[z]
        ref = ase.data.reference_states[z]
        if ref is None:
            struct = _("No crystal structure data")
        else:
            struct = ref['symmetry']
            if struct == 'fcc' or struct == 'bcc':
                struct = "%s (a=%.3f Å)" % (struct, ref['a'])

        txt = "  %s: %s, Z=%i, %s" % (name, symb, z, struct)
        self.elementinfo.set_text(txt)
        self.legal_element = symb
        return True

    def invalid_element(self, txt=_("  ERROR: Invalid element!")):
        self.legal_element = False
        self.elementinfo.set_text(txt)

