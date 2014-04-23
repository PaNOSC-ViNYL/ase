from __future__ import print_function
import sys
import optparse

import ase.io
from ase.db import connect
from ase.visualize import view
from ase.data import atomic_masses, chemical_symbols
from ase.calculators.calculator import get_calculator
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
    
    
description = """Selecton is a comma-separated list of
selections where each selection is of the type "ID", "keyword" or
"key=value".  Instead of "=", one can also use "<", "<=", ">=", ">"
and  "!=" (these must be protected from the shell by using quotes).
Special keys: id, user, calculator, age, natoms, energy, magmom,
and charge.  Chemical symbols can also be used to select number of
specific atomic species (H, He, Li, ...)."""

examples = ['calculator=nwchem',
            'age<1d',
            'natoms=1',
            'user=alice',
            '2.2<bandgap<4.1',
            'Cu>=10']


def main(args=sys.argv[1:]):
    if isinstance(args, str):
        args = args.split(' ')
    parser = optparse.OptionParser(
        usage='Usage: %prog db-name [selection] [options]',
        description=description,
        epilog='Selection examples: ' + ', '.join(examples) + '.')
    
    add = parser.add_option
    add('-v', '--verbose', action='store_true', default=False)
    add('-q', '--quiet', action='store_true', default=False)
    add('-n', '--count', action='store_true',
        help='Count number of selected rows.')
    add('-l', '--long', action='store_true',
        help='Long description of selected row')
    add('-i', '--insert-into', metavar='db-name',
        help='Insert selected rows into another database.')
    add('-a', '--add-from-file', metavar='[type:]filename',
        help='Add results from file.')
    add('-k', '--add-keywords', metavar='word1,word2,...',
        help='Add keywords to selected rows.  Keywords can only contain ' +
        'letters, numbers and the underscore character and the first ' +
        'character can not be a number.')
    add('-K', '--add-key-value-pairs', metavar='key1=val1,key2=val2,...',
        help='Add key-value pairs to selected rows.  Values must be numbers ' +
        'or strings and keys must follow the same rules as keywords.')
    add('--limit', type=int, default=500, metavar='N',
        help='Show only first N rows (default is 500 rows).  Use --limit=0 ' +
        'to show all.')
    add('--delete', action='store_true',
        help='Delete selected rows.')
    add('--delete-keywords', metavar='key1=word1,word2,...',
        help='Delete keywords for selected rows.')
    add('--delete-key-value-pairs', metavar='key1=val1,key2=val2,...',
        help='Delete key-value pairs for selected rows.')
    add('-y', '--yes', action='store_true',
        help='Say yes.')
    add('--explain', action='store_true',
        help='Explain query plan.')
    add('-c', '--columns', metavar='col1,col2,...',
        help='Specify columns to show.  Precede the column specification ' +
        'with a "+" in order to add columns to the default set of columns.  ' +
        'Precede by a "-" to remove columns.')
    add('-s', '--sort', metavar='column',
        help='Sort rows using column.  Default is to sort after ID.')
    add('-p', '--python-expression', metavar='expression',
        help='Examples: "id,energy", "id,mykey".')
    add('-w', '--open-web-browser', action='store_true',
        help='Open results in web-browser.')

    opts, args = parser.parse_args(args)

    if not args:
        parser.error('No database given')

    verbosity = 1 - opts.quiet + opts.verbose

    try:
        run(opts, args, verbosity)
    except Exception as x:
        if verbosity < 2:
            print('{0}: {1}'.format(x.__class__.__name__, x.message))
            sys.exit(1)
        else:
            raise

        
def run(opts, args, verbosity):
    filename = args.pop(0)
    query = ','.join(args)

    if query.isdigit():
        query = int(query)
    
    if opts.add_keywords:
        add_keywords = opts.add_keywords.split(',')
    else:
        add_keywords = []

    add_key_value_pairs = {}
    if opts.add_key_value_pairs:
        for pair in opts.add_key_value_pairs.split(','):
            key, value = pair.split('=')
            for type in [int, float]:
                try:
                    value = type(value)
                except ValueError:
                    pass
                else:
                    break
            add_key_value_pairs[key] = value

    con = connect(filename)
    
    if opts.add_from_file:
        filename = opts.add_from_file
        if ':' in filename:
            calculator_name, filename = filename.split(':')
            atoms = get_calculator(calculator_name)(filename).get_atoms()
        else:
            atoms = ase.io.read(filename)
        con.write(atoms, add_keywords, key_value_pairs=add_key_value_pairs)
        print('Added {0} from {1}'.format(atoms.get_chemical_formula(),
                                          filename))
        return
        
    if opts.count:
        n = 0
        for dct in con.select(query):
            n += 1
        print('%s' % plural(n, 'row'))
        return

    if opts.explain:
        for dct in con.select(query, explain=True,
                              verbosity=verbosity, limit=opts.limit):
            print('%d %d %d %s' % dct['explain'])
        return

    if opts.insert_into:
        con2 = connect(opts.insert_into)
        nkw = 0
        nkvp = 0
        nrows = 0
        for dct in con.select(query):
            keywords = dct.get('keywords', [])
            for keyword in add_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    nkw += 1

            kvp = dct.get('key_value_pairs', {})
            nkvp = -len(kvp)
            kvp.update(add_key_value_pairs)
            nkvp += len(kvp)
            con2.write(dct, keywords, data=dct.get('data'), **kvp)
            nrows += 1
            
        print('Added %s and %s (%s updated)' %
              (plural(nkw, 'keyword'),
               plural(nkvp, 'key-value pair'),
               plural(len(add_key_value_pairs) * nrows - nkvp, 'pair')))
        print('Inserted %s' % plural(nrows, 'row'))
        return

    if add_keywords or add_key_value_pairs:
        ids = [dct['id'] for dct in con.select(query)]
        nkw, nkv = con.update(ids, add_keywords, **add_key_value_pairs)
        print('Added %s and %s (%s updated)' %
              (plural(nkw, 'keyword'),
               plural(nkv, 'key-value pair'),
               plural(len(add_key_value_pairs) * len(ids) - nkv, 'pair')))
        return

    if opts.delete:
        ids = [dct['id'] for dct in con.select(query)]
        if ids and not opts.yes:
            msg = 'Delete %s? (yes/no): ' % plural(len(ids), 'row')
            if raw_input(msg).lower() != 'yes':
                return
        con.delete(ids)
        print('Deleted %s' % plural(len(ids), 'row'))
        return

    if opts.python_expression:
        for dct in con.select(query):
            row = eval(opts.python_expression, dct)
            if not isinstance(row, (list, tuple, np.ndarray)):
                row = [row]
            print(', '.join(str(x) for x in row))
        return

    if opts.long:
        dct = con.get(query)
        summary = Summary(dct)
        if opts.open_web_browser:
            from ase.db.web import run
            run(summary=dct)
        else:
            summary.write()
    else:
        rows = Rows(con, query, opts.limit, verbosity, opts.columns,
                    opts.sort)
        if opts.open_web_browser:
            from ase.db.web import run
            run(rows)
        else:
            rows.write()


class Rows:
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
        view(dict2atoms(self.rows[n].dct))
    
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

        
class Summary:
    def __init__(self, dct):
        self.dct = dct
        
        self.cell = [['{0:.3f}'.format(a) for a in axis] for axis in dct.cell]
        
        forces = dict2forces(dct)
        if forces is None:
            fmax = None
            self.forces = None
        else:
            fmax = (forces**2).sum(1).max()**0.5
            N = len(forces)
            self.forces = []
            for n, f in enumerate(forces):
                if n < 5 or n >= N - 5:
                    f = tuple('{0:10.3f}'.format(x) for x in f)
                    symbol = chemical_symbols[dct.numbers[n]]
                    self.forces.append((n, symbol) + f)
                elif n == 5:
                    self.forces.append((' ...', '',
                                        '       ...',
                                        '       ...',
                                        '       ...'))
                    
        self.stress = dct.get('stress')
        if self.stress is not None:
            self.stress = ', '.join('{0:.3f}'.format(s) for s in self.stress)
            
        if 'masses' in dct:
            mass = dct.masses.sum()
        else:
            mass = atomic_masses[dct.numbers].sum()
            
        table = [
            ('id', dct.id),
            ('age', float_to_time_string(now() - dct.ctime, True)),
            ('formula', hill(dct.numbers)),
            ('user', dct.user),
            ('calculator', dct.get('calculator')),
            ('energy [eV]', dct.get('energy')),
            ('fmax [eV/Ang]', fmax),
            ('charge [|e|]', dct.get('charge')),
            ('mass [au]', mass),
            ('unique id', dct.unique_id),
            ('volume [Ang^3]', abs(np.linalg.det(dct.cell)))]
        self.table = [(name, value) for name, value in table
                      if value is not None]

        if 'key_value_pairs' in dct:
            self.key_value_pairs = sorted(dct.key_value_pairs.items())
        else:
            self.key_value_pairs = None

        if 'keywords' in dct:
            self.keywords = sorted(dct.keywords)
        else:
            self.keywords = None
            
        self.dipole = dct.get('dipole')
        if self.dipole is not None:
            self.dipole = ', '.join('{0:.3f}'.format(d) for d in self.dipole)
        
        self.data = dct.get('data')
        if self.data:
            self.data = ', '.join(self.data.keys())
            
        self.constraints = dct.get('constraints')
        if self.constraints:
            self.constraints = ', '.join(d[name] for d in self.constraints)
        
    def write(self):
        dct = self.dct
        
        width = max(len(name) for name, value in self.table)
        for name, value in self.table:
            print('{0:{width}}|{1}'.format(name, value, width=width))

        print('\nUnit cell in Ang:')
        print('axis|periodic|     x     |     y     |    z')
        c = 1
        for p, axis in zip(dct.pbc, self.cell):
            print('   {0}|     {1}|{2[0]:>11}|{2[1]:>11}|{2[2]:>11}'.format(
                c, [' no', 'yes'][p], axis))
            c += 1
            
        if self.key_value_pairs:
            print('\nKey-value pairs:')
            width = max(len(key) for key, value in self.key_value_pairs)
            for key, value in self.key_value_pairs:
                print('{0:{width}}|{1}'.format(key, value, width=width))
                
        if self.keywords:
            print('\nKeywords:', ', '.join(self.keywords))
                
        if self.forces:
            print('\nForces in ev/Ang:')
            for f in self.forces:
                print('{0:4}|{1:2}|{2}|{3}|{4}'.format(*f))

        if self.stress:
            print('\nStress tensor (xx, yy, zz, zy, zx, yx) in eV/Ang^3:')
            print('   ', self.stress)

        if self.dipole:
            print('\nDipole moment in e*Ang: ({0})'.format(self.dipole))
        
        if self.constraints:
            print('\nConstraints:', self.constraints)
            
        if self.data:
            print('\nData:', self.data)
