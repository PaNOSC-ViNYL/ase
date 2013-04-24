import sys
import argparse

import numpy as np

from ase.db import connect
from ase.atoms import Atoms


def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('name', nargs=1)
    parser.add_argument('query', nargs=1)
    'limit,explain'
    args = parser.parse_args()
    con = connect(args.name[0])
    Formatter().format(list(con.select(*args.query)))


class Formatter:
    def format(self, dcts, columns=None):
        columns = ['id', 'age', 'user', 'symbols', 'calc', 
                   'energy', 'fmax', 'pbc', 'size', 'keywords',
                   'charge', 'mass', 'fixed', 'smax', 'magmom']
        table = [columns]
        chars = [len(column) for column in columns]
        fd = sys.stdout
        for dct in dcts:
            row = []
            for i, column in enumerate(columns):
                try:
                    s = getattr(self, column)(dct)
                except KeyError:
                    s = ''
                if len(s) > chars[i]:
                    chars[i] = len(s)
                row.append(s)
            table.append(row)
        for row in table:
            fd.write('|'.join('%*s' % (n, s) for n, s in zip(chars, row)))
            fd.write('\n')
        
    def id(self, d):
        return d.id
    
    def age(self, d):
        return '%.1f' % d.timestamp

    def user(self, d):
        return d.username
    
    def symbols(self, d):
        return Atoms(d.numbers).get_chemical_formula()

    def energy(self, d):
        return '%.3f' % d.results.energy

    def size(self, d):
        return '%.3f' % abs(np.linalg.det(d.cell))

    def pbc(self, d):
        a, b, c = d.pbc
        return '%s' % (' a'[a] + ' b'[b] + ' c'[c])

    def calc(self, d):
        return '%s' % d.calculator_name

    def fmax(self, d):
        return '%.3f' % (d.results.forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return '%s' % ','.join(d.keywords +
                               ['%s=%s' % (key, value)
                                for key, value in d.key_value_pairs.items()])

    def charge(self, d):
        return '%.3f' % d.results.charge

    def mass(self, d):
        return '%.3f' % d.masses.sum()

    def fixed(self, d):
        c = d.constraints
        if c is None:
            return ''
        if len(c) > 1:
            return '?'

    def smax(self, d):
        return '%.3f' % (d.results.stress**2).max()**0.5

    def magmom(self, d):
        return '%.3f' % d.results.magmom
