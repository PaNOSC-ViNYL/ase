from ase.db.row import AtomsRow
from ase.io.jsonio import read_json


class Collection:
    def __init__(self, name):
        self.name = name
        self._names = []
        self._systems = {}
        self._data = {}
        
    def __getitem__(self, name):
        self._read()
        return self._systems[name]

    def __iter__(self):
        for name in self.names:
            yield self[name]
            
    def __len__(self):
        return len(self.names)
        
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
        bigdct = read_json(__file__[:-13] + self.name + '.json')
        for id in bigdct['ids']:
            dct = bigdct[id]
            kvp = dct['key_value_pairs']
            name = kvp['name']
            self._names.append(name)
            self._systems[name] = AtomsRow(dct).toatoms()
            del kvp['name']
            self._data[name] = kvp

            
if 0:
    from ase.data.s22 import s22, data
    import ase.db
    from ase import Atoms
    import os
    os.environ['USER'] = 'ase'
    c = ase.db.connect('s22.json')
    for n in s22:
        d = data[n]
        a = Atoms(d.pop('symbols'), d.pop('positions'),
                  magmoms=d.pop('magmoms'))
        c.write(a, name=n, cc_energy=d['interaction energy CC'])
