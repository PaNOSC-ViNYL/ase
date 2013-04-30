import sys
import argparse

import numpy as np

from ase.db import connect
from ase.atoms import Atoms


def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    add = parser.add_argument
    add('name', nargs=1)
    add('selection', nargs='?')
    add('-n', '--count', action='store_true')
    add('-c', '--columns', help='short/long+row-row')
    add('-l', '--limit', type=int)
    add('-o', '--offset', type=int)
    add('--explain', action='store_true')
    add('-y', '--yes', action='store_true')
    add('-i', '--insert-into')
    add('--set-keyword')
    add('--set-key-value-pair')
    add('--remove', action='store_true')
    add('-v', '--verbose', action='store_true')
    add('-q', '--quiet', action='store_true')

    args = parser.parse_args()

    verbosity = 1 - args.quiet + args.verbose

    con = connect(args.name[0])

    if verbosity == 2:
        print(args)

    if args.selection:
        selection = (args.selection,)
    else:
        selection = ()

    rows = con.select(*selection, limit=args.limit, offset=args.offset,
                       explain=args.explain, set_keywords=args.set_keyword,
                       count=args.count, verbosity=verbosity)
    if args.count:
        n = rows.next()[0]
        print('%d row%s' % (n, 's'[:int(n!=1)]))
        return

    if args.explain:
        for row in rows:
            print('%d %d %d %s' % row)
        return

    if args.insert_into:
        return    

    if args.remove:
        ids = [dct['id'] for dct in rows]
        if ids and not args.yes:
            if raw_input('Remove %d rows? (yes/no): ').lower() != 'yes':
                return
        for id in ids:
            con.remove(id)
        return

    f = Formatter(columns=args.columns)
    f.format(list(rows))


def cut(txt, length):
    if len(txt) <= length:
        return txt
    return txt[:length - 3] + '...'


class Formatter:
    def __init__(self, columns):
        pass

    def format(self, dcts, columns=None):
        columns = ['id', 'age', 'user', 'symbols', 'calc',
                   'energy', 'fmax', 'pbc', 'size', 'keywords',
                   'charge', 'mass', 'fixed', 'smax', 'magmom']
        table = [columns]
        widths = [0 for column in columns]
        fd = sys.stdout
        for dct in dcts:
            row = []
            for i, column in enumerate(columns):
                try:
                    s = getattr(self, column)(dct)
                except KeyError:
                    s = ''
                if len(s) > widths[i]:
                    widths[i] = len(s)
                row.append(s)
            table.append(row)
        widths = [w and max(w, len(column))
                  for w, column in zip(widths, columns)]
        for row in table:
            fd.write('|'.join('%*s' % (w, s)
                              for w, s in zip(widths, row) if w > 0))
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
        dims = d.pbc.sum()
        if dims == 0:
            return ''
        if dims == 1:
            return '%.3f' % np.linalg.norm(d.cell[d.pbc][0])
        if dims == 2:
            return '%.3f' % np.linalg.norm(np.cross(*d.cell[d.pbc]))
        return '%.3f' % abs(np.linalg.det(d.cell))

    def pbc(self, d):
        a, b, c = d.pbc
        return '%d%d%d' % tuple(d.pbc)

    def calc(self, d):
        return '%s' % d.calculator_name

    def fmax(self, d):
        return '%.3f' % (d.results.forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return '%s' % ','.join(d.keywords +
                               ['%s=%s' % (key, cut(value, 8))
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
