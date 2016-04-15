from ase.db.row import AtomsRow
from ase.io.jsonio import read_json


class Collection:
    """Collection of atomic configurations and associated data.
    
    >>> from ase.collections import g2
    >>> g2.names
    >>> g2.filename
    >>> g2['CO2']
    >>> g2.data['CO2']
    >>> ???
        
    """
    def __init__(self, name):
        """Create a collection lazily.
        
        Will read data from json file when needed.
        
        Attributes:
        
        name:
        data
        filename
        names
        """
        
        self.name = name
        self._names = []
        self._systems = {}
        self._data = {}
        self.filename = __file__[:-13] + self.name + '.json'
        
    def __getitem__(self, name):
        self._read()
        return self._systems[name]

    def __iter__(self):
        for name in self.names:
            yield self[name]
            
    def __len__(self):
        return len(self.names)
        
    def __str__(self):
        return '<{0}-collection, {1} systems: {2}, {3}, ...>'.format(
            self.name, len(self), *self.names[:2])
        
    def __repr__(self):
        return 'Collection({0!r})'.format(self.name)
        
    @property
    def names(self):
        self._read()
        return self._names
        
    @property
    def data(self):
        self._read()
        return self._data
        
    def _read(self):
        if self._names:
            return
        bigdct = read_json(self.filename)
        self._description = bigdct['description']
        for id in bigdct['ids']:
            dct = bigdct[id]
            kvp = dct['key_value_pairs']
            name = kvp['name']
            self._names.append(name)
            self._systems[name] = AtomsRow(dct).toatoms()
            del kvp['name']
            self._data[name] = kvp

