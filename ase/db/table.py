from __future__ import print_function

from ase.visualize import view
from ase.data import atomic_masses, chemical_symbols
from ase.db.core import float_to_time_string, now, dict2constraint, dict2atoms

import numpy as np


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)

    
def cut(txt, length):
    if len(txt) <= length:
        return txt
    return txt[:length - 3] + '...'    


def hill(numbers):
    d = {}
    for Z in numbers:
        s = chemical_symbols[Z]
        d[s] = d.get(s, 0) + 1
    result = [(s, d.pop(s)) for s in 'CH' if s in d]
    result += [(s, d[s]) for s in sorted(d)]
    return ''.join(
        '{0}<sub>{1}</sub>'.format(symbol, n) if n > 1 else symbol
        for symbol, n in result)
    
    
def dict2forces(d):
    forces = d.get('forces')
    if forces is None:
        return None
        
    constraints = [dict2constraint(c) for c in d.get('constraints', [])]
    if constraints:
        forces = forces.copy()
        for constraint in constraints:
            constraint.adjust_forces(d.positions, forces)
            
    return forces

    
class Table:
    def __init__(self, connection, query='', limit=50, verbosity=1,
                 columns=None, sort=None):
        self.connection = connection
        self.query = ''
        self.limit = limit
        self.verbosity = verbosity
        self.sort = sort
        
        self.rows = []
        
        self.columns = ['id', 'age', 'user', 'formula', 'calculator',
                        'energy', 'fmax', 'pbc', 'volume',
                        'keywords', 'keys',
                        'charge', 'mass', 'smax', 'magmom']
        
        if columns is not None:
            if columns[0] == '+':
                columns = columns[1:]
            elif columns[0] != '-':
                self.columns = []
            for col in columns.split(','):
                if col[0] == '-':
                    self.columns.remove(col[1:])
                else:
                    self.columns.append(col.lstrip('+'))
                    
        self.columns_original = self.columns
        
        self.search(query)
        
    def search(self, query):
        self.query = query
        self.columns = list(self.columns_original)
        
        limit = self.limit if not self.sort else 0
        self.rows = [Row(d, self.columns)
                     for d in self.connection.select(
                         query, verbosity=self.verbosity, limit=limit)]

        delete = set(range(len(self.columns)))
        for row in self.rows:
            for n in delete.copy():
                if row.values[n] is not None:
                    delete.remove(n)
        delete = sorted(delete, reverse=True)
        for row in self.rows:
            for n in delete:
                del row.values[n]
        for n in delete:
            del self.columns[n]
            
        if self.sort:
            reverse = self.sort[0] == '-'
            n = self.columns.index(self.sort.lstrip('-'))
            
            def key(row):
                return row.values[n]
                
            self.rows = sorted(self.rows, key=key, reverse=reverse)
            if self.limit:
                self.rows = self.rows[:self.limit]
                
    def toggle(self, key=None, n=None):
        if n is not None:
            self.rows[n].toggle()
            return
            
        if key in self.columns:
            self.columns.remove(key)
        else:
            self.columns.append(key)
        for row in self.rows:
            row.set_columns(self.columns)
            
    def remove(self, column):
        self.columns.remove(column)
        for row in self.rows:
            row.set_columns(self.columns)
        
    def gui(self, n):
        dct = self.rows[n].dct
        view(dict2atoms(dct))
    
    def format(self, mode='ascii'):
        right = set()
        allkeys = set()
        for row in self.rows:
            numbers = row.format(self.columns, mode)
            right.update(numbers)
            allkeys.update(row.dct.key_value_pairs)
            
        right.add('age')
        self.right = [column in right for column in self.columns]
        
        self.keys = sorted(allkeys)

    def write(self):
        self.format()
        L = [[len(s) for s in row.strings]
             for row in self.rows]
        L.append([len(c) for c in self.columns])
        N = np.max(L, axis=0)
        print('Rows:', len(self.rows), end='')
        if self.limit:
            print(' (limited to first {0})'.format(self.limit))
        else:
            print()

        if len(self.rows) == 0:
            return
            
        fmt = '{0:{align}{width}}'
        print('|'.join(fmt.format(c, align='<>'[a], width=w)
                       for c, a, w in zip(self.columns, self.right, N)))
        for row in self.rows:
            print('|'.join(fmt.format(c, align='<>'[a], width=w)
                           for c, a, w in
                           zip(row.strings, self.right, N)))
        if self.keys:
            print('Keys:', ', '.join(self.keys))

                
class Row:
    def __init__(self, dct, columns):
        self.dct = dct
        self.values = None
        self.strings = None
        self.more = False
        self.set_columns(columns)
        if 'key_value_pairs' not in dct:
            dct['key_value_pairs'] = {}
        if 'keywords' not in dct:
            dct['keywords'] = []
        
    def set_columns(self, columns):
        self.values = []
        for c in columns:
            f = getattr(self, c, None)
            if f is None:
                value = getattr(self.dct, c, None)
            else:
                try:
                    value = f(self.dct)
                except (AttributeError, TypeError):
                    value = None
            self.values.append(value)
            
    def toggle(self):
        self.more = not self.more
        
    def format(self, columns, mode):
        self.strings = []
        numbers = set()
        for value, column in zip(self.values, columns):
            if column == 'formula' and mode == 'ascii':
                value = value.replace('<sub>', '').replace('</sub>', '')
            elif isinstance(value, int):
                value = str(value)
                numbers.add(column)
            elif isinstance(value, float):
                numbers.add(column)
                value = '{0:.3f}'.format(value)
            elif value is None:
                value = ''
            self.strings.append(value)
        
        if self.more:
            self.cellstring = ['({0:.3f}, {1:.3f}, {2:.3f})'.format(*axis)
                               for axis in self.dct.cell]
            
        return numbers
        
    def age(self, d):
        return float_to_time_string(now() - d.ctime)

    def formula(self, d):
        return hill(d.numbers)

    def volume(self, d):
        return abs(np.linalg.det(d.cell))

    def pbc(self, d):
        return ''.join('-P'[p] for p in d.pbc)

    def fmax(self, d):
        forces = dict2forces(d)
        return (forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return cut(','.join(d.keywords), 30)

    def keys(self, d):
        return cut(','.join(['%s=%s' % (key, cut(str(value), 11))
                             for key, value in d.key_value_pairs.items()]), 40)

    def mass(self, d):
        if 'masses' in d:
            return d.masses.sum()
        return atomic_masses[d.numbers].sum()

    def smax(self, d):
        return (d.stress**2).max()**0.5
