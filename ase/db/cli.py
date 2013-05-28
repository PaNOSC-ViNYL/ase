import sys
import argparse
from time import time

import numpy as np

from ase.db import connect
from ase.db.core import float_to_time_string, T0, YEAR, dict2atoms
from ase.atoms import Atoms


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)


def run(args=sys.argv[1:]):
    if isinstance(args, str):
        args = args.split(' ')
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    add = parser.add_argument
    add('name', nargs=1)
    add('selection', nargs='?')
    add('-n', '--count', action='store_true')
    add('-c', '--columns', help='short/long+row-row')
    add('--explain', action='store_true')
    add('-y', '--yes', action='store_true')
    add('-i', '--insert-into')
    add('-k', '--add-keywords')
    add('-K', '--add-key-value-pairs')
    add('--delete-keywords')
    add('--delete-key-value-pairs')
    add('--delete', action='store_true')
    add('-v', '--verbose', action='store_true')
    add('-q', '--quiet', action='store_true')

    args = parser.parse_args(args)

    verbosity = 1 - args.quiet + args.verbose

    con = connect(args.name[0])

    if verbosity == 2:
        print(args)

    rows = con.select(args.selection, explain=args.explain,
                      verbosity=verbosity)

    if args.count:
        n = 0
        for row in rows:
            n += 1
        print('%s' % plural(n, 'row'))
        return n

    if args.explain:
        for row in rows:
            print('%d %d %d %s' % row)
        return

    if args.add_keywords:
        add_keywords = args.add_keywords.split(',')
    else:
        add_keywords = []

    add_key_value_pairs = {}
    if args.add_key_value_pairs:
        for pair in args.add_key_value_pairs.split(','):
            key, value = pair.split('=')
            try:
                value = int(value)
                value = float(value)
            except ValueError:
                pass
            add_key_value_pairs[key] = value

    if args.insert_into:
        con2 = connect(args.insert_into)
        n = 0
        ids = []
        for dct in rows:
            for keyword in add_keywords:
                if keyword not in dct.keywords:
                    dct.keywords.append(keyword)
                    n += 1
            rollback = True
            if 1:#try:
                id = con2.write(dct.id, dct, timestamp=dct.timestamp)
                rollback = False
            if 0:#finally:
                if rollback:
                    con2.delete(ids)
                    return
            ids.append(id)
        print('Added %s' % plural(n, 'keyword'))
        print('Inserted %s' % plural(len(ids), 'row'))
        return ids  

    if add_keywords or add_key_value_pairs:
        ids = [dct['id'] for dct in rows]
        m, n = con.update(ids, add_keywords, add_key_value_pairs)
        print('Added %s and %s (%s updated)' %
              (plural(m, 'keyword'),
               plural(n, 'key-value pair'),
               plural(len(add_key_value_pairs) * len(ids) - n, 'pair')))
        return ids

    if args.delete:
        ids = [dct['id'] for dct in rows]
        if ids and not args.yes:
            msg = 'Delete %s? (yes/no): ' % plural(len(ids), 'row')
            if raw_input(msg).lower() != 'yes':
                return
        con.delete(ids)
        print('Deleted %s' % plural(len(ids), 'row'))
        return ids

    dcts = list(rows)
    if len(dcts) > 0:
        f = Formatter(columns=args.columns)
        return f.format(dcts)

    return []


def cut(txt, length):
    if len(txt) <= length:
        return txt
    return txt[:length - 3] + '...'


class Formatter:
    def __init__(self, columns):
        pass

    def format(self, dcts, columns=None):
        columns = ['id', 'age', 'user', 'formula', 'calc',
                   'energy', 'fmax', 'pbc', 'size', 'keywords', 'keyvals',
                   'charge', 'mass', 'fixed', 'smax', 'magmom']
        table = [columns]
        widths = [0 for column in columns]
        signs = [1 for column in columns]  # left or right adjust
        ids = []
        fd = sys.stdout
        for dct in dcts:
            row = []
            for i, column in enumerate(columns):
                try:
                    s = getattr(self, column)(dct)
                except AttributeError:
                    s = ''
                if isinstance(s, int):
                    s = '%d' % s 
                elif isinstance(s, float):
                    s = '%.3f' % s
                else:
                    signs[i] = -1
                if len(s) > widths[i]:
                    widths[i] = len(s)
                row.append(s)
            table.append(row)
            ids.append(dct.id)
        widths = [w and max(w, len(column))
                  for w, column in zip(widths, columns)]
        for row in table:
            fd.write('|'.join('%*s' % (w * sign, s)
                              for w, sign, s in zip(widths, signs, row)
                              if w > 0))
            fd.write('\n')
        return ids
        
    def id(self, d):
        return d.id
    
    def age(self, d):
        return float_to_time_string((time() - T0) / YEAR - d.timestamp)

    def user(self, d):
        return d.username
    
    def formula(self, d):
        return Atoms(d.numbers).get_chemical_formula()

    def energy(self, d):
        return d.results.energy

    def size(self, d):
        dims = d.pbc.sum()
        if dims == 0:
            return ''
        if dims == 1:
            return np.linalg.norm(d.cell[d.pbc][0])
        if dims == 2:
            return np.linalg.norm(np.cross(*d.cell[d.pbc]))
        return abs(np.linalg.det(d.cell))

    def pbc(self, d):
        a, b, c = d.pbc
        return '%d%d%d' % tuple(d.pbc)

    def calc(self, d):
        return d.calculator_name

    def fmax(self, d):
        return (d.results.forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return ','.join(d.keywords)

    def keyvals(self, d):
        return ','.join(['%s=%s' % (key, cut(value, 8))
                         for key, value in d.key_value_pairs.items()])

    def charge(self, d):
        return d.results.charge

    def mass(self, d):
        return d.masses.sum()

    def fixed(self, d):
        c = d.constraints
        if c is None:
            return ''
        if len(c) > 1:
            return '?'

    def smax(self, d):
        return (d.results.stress**2).max()**0.5

    def magmom(self, d):
        return d.results.magmom
