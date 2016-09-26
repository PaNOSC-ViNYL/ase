from gettext import gettext as _

import ase.data
import ase.gui.ui as ui


class Element(list):
    def __init__(self, symbol='', callback=None):
        list.__init__(self,
                      [_('Element:'),
                       ui.Entry(symbol, 3, self.enter),
                       ui.Label('', 'red')])
        self.callback = callback
        self._symbol = None
        self._Z = None

    @property
    def symbol(self):
        self.check()
        return self._symbol

    @property
    def Z(self):
        self.check()
        return self._Z

    def check(self):
        self._symbol = self[1].value
        if not self._symbol:
            self.error(_('No element specified!'))
            return
        self._Z = ase.data.atomic_numbers.get(self._symbol)
        if self._Z is None:
            try:
                self._Z = int(self._symbol)
            except ValueError:
                self.error()
                return
            self._symbol = ase.data.chemical_symbols[self._Z]
        self[2].text = ''

    def enter(self):
        self.check()
        self.callback()

    def error(self, text=_('ERROR: Invalid element!')):
        self._symbol = None
        self._Z = None
        self[2].text = text
