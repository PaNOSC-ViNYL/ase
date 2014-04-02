from __future__ import print_function
import sys
import optparse

import ase.io
from ase.db import connect
from ase.atoms import Atoms
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
        help='Long description of selected row (or first row selected)')
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
        help="""Examples: "d.id", "d.mykey", where """ +
        '"d" is a dictionary representing a row.')
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
    expressions = ','.join(args)
        
    if verbosity == 2:
        print('Options:')
        for k, v in opts.__dict__.items():
            print('    {0:24}{1}'.format(k + ':', v))
        print('Arguments:', ', '.join(args))

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
        for dct in con.select(expressions):
            n += 1
        print('%s' % plural(n, 'row'))
        return

    if opts.explain:
        for dct in con.select(expressions, explain=True,
                              verbosity=verbosity, limit=opts.limit):
            print('%d %d %d %s' % dct['explain'])
        return

    if opts.insert_into:
        con2 = connect(opts.insert_into)
        nkw = 0
        nkvp = 0
        nrows = 0
        for dct in con.select(expressions, limit=opts.limit):
            keywords = dct.get('keywords', [])
            for keyword in add_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    nkw += 1

            kvp = dct.get('key_value_pairs', {})
            nkvp = -len(kvp)
            kvp.update(add_key_value_pairs)
            nkvp += len(kvp)
            con2.write(dct, keywords, **kvp)
            nrows += 1
            
        print('Added %s and %s (%s updated)' %
              (plural(nkw, 'keyword'),
               plural(nkvp, 'key-value pair'),
               plural(len(add_key_value_pairs) * nrows - nkvp, 'pair')))
        print('Inserted %s' % plural(nrows, 'row'))
        return

    if add_keywords or add_key_value_pairs:
        ids = [dct['id'] for dct in con.select(expressions, limit=opts.limit)]
        nkw, nkv = con.update(ids, add_keywords, **add_key_value_pairs)
        print('Added %s and %s (%s updated)' %
              (plural(nkw, 'keyword'),
               plural(nkv, 'key-value pair'),
               plural(len(add_key_value_pairs) * len(ids) - nkv, 'pair')))
        return

    if opts.delete:
        ids = [dct['id'] for dct in con.select(expressions, limit=opts.limit)]
        if ids and not opts.yes:
            msg = 'Delete %s? (yes/no): ' % plural(len(ids), 'row')
            if raw_input(msg).lower() != 'yes':
                return
        con.delete(ids)
        print('Deleted %s' % plural(len(ids), 'row'))
        return

    if opts.python_expression:
        for n, dct in enumerate(con.select(expressions, limit=opts.limit)):
            row = eval(opts.python_expression, {'d': dct, 'n': n})
            if not isinstance(row, (list, tuple, np.ndarray)):
                row = [row]
            print(', '.join(str(x) for x in row))
        return

    if opts.long:
        dct = con.get(expressions)
        summary = Summary(dct)
        if opts.open_web_browser:
            from ase.db.web import run
            run(summary=dct)
        else:
            summary.write()
    else:
        rows = Rows(con, expressions, opts.limit, verbosity, opts.columns)
        if opts.open_web_browser:
            from ase.db.web import run
            run(rows)
        else:
            rows.write()


class Rows:
    def __init__(self, connection, query='', limit=50, verbosity=1,
                 columns=None):
        self.connection = connection
        self.query = ''
        self.rows = []
        self.time = 0.0
        
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
                    
        self.columns_original = list(self.columns)
        
        self.search(query)
        
    def search(self, query):
        self.query = query
        if not isinstance(query, int) and query.isdigit():
            query = int(query)
        self.columns = self.columns_original
        self.rows = [Row(d, self.columns)
                     for d in self.connection.select(query)]

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
        
        self.keys = ['id'] + sorted(allkeys)

    def write(self, fd=sys.stdout):
        self.format()
        L = [[len(s) for s in row.strings]
             for row in self.rows]
        L.append([len(c) for c in self.columns])
        N = np.max(L, axis=0)
        print('Rows:', len(self.rows), '(limited to first 500) Time:',
              self.time, file=fd)
        if len(self.rows) == 0:
            return
        fmt = '{0:{align}{width}}'
        print('|'.join(fmt.format(c, align='<>'[a], width=w)
                       for c, a, w in zip(self.columns, self.right, N)),
              file=fd)
        for row in self.rows:
            print('|'.join(fmt.format(c, align='<>'[a], width=w)
                           for c, a, w in
                           zip(row.strings, self.right, N)), file=fd)

                
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
                except AttributeError:
                    value = None
            self.values.append(value)
            
    def toggle(self):
        self.more = not self.more
        
    def format(self, columns, mode):
        self.strings = []
        numbers = set()
        for value, column in zip(self.values, columns):
            if column == 'formula':
                pass
                #'{0}<sub>{1}</sub>'.format(symbol, n) if n > 1 else symbol
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
        return ''.join(
            '{0}<sub>{1}</sub>'.format(symbol, n) if n > 1 else symbol
            for symbol, n in hill(d.numbers))

    def volume(self, d):
        return abs(np.linalg.det(d.cell))

    def pbc(self, d):
        return ''.join('-P'[p] for p in d.pbc)

    def fmax(self, d):
        forces = d.forces
        constraints = [dict2constraint(c) for c in d.constraints]
        if constraints:
            forces = forces.copy()
            for constraint in constraints:
                constraint.adjust_forces(d.positions, forces)
        return (forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return cut(','.join(d.keywords), 30)

    def keys(self, d):
        return cut(','.join(['%s=%s' % (key, cut(str(value), 8))
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
    result += [(s, d[s]) for s in sorted(d.keys())]
    return result
    
    
def long(d, verbosity=1):
    print('id:', d.id)
    print('formula:', Atoms(d.numbers).get_chemical_formula())
    print('user:', d.user)
    print('age: {0}'.format(float_to_time_string(now() - d.ctime)))
    if 'calculator' in d:
        print('calculator:', d.calculator)
    if 'energy' in d:
        print('energy: {0:.3f} eV'.format(d.energy))
    if 'forces' in d:
        print('maximum atomic force: {0:.3f} eV/Ang'.format(
                  (d.forces**2).sum(1).max()**0.5))
    if 'stress' in d:
        print('stress tensor: (xx, yy, zz, zy, zx, yx) eV/Ang^3')
        for s in d.stress:
            print('{0:10.3f}'.format(s), end='')
        print()
    if 'dipole' in d:
        print('dipole moment [e*Ang]: {0:10.3f}, {1:10.3f}, {2:10.3f}'.format(
            *d.dipole))
    print('magnetic moment:', d.get('magmom', 0))
    print('periodic boundary conditions:', d.pbc)
    print('unit cell [Ang]:')
    for axis in d.cell:
        print('{0:10.3f}{1:10.3f}{2:10.3f}'.format(*axis))
    dims = d.pbc.sum()
    if dims == 1:
        print('length: {0:.3f} Ang'.format(np.linalg.norm(d.cell[d.pbc][0])))
    elif dims == 2:
        print('area: {0:.3f} Ang^2'.format(
                  np.linalg.norm(np.cross(*d.cell[d.pbc]))))
    print('volume: {0:.3f} Ang^3'.format(abs(np.linalg.det(d.cell))))
    if 'charge' in d:
        print('charge: {0:.6f}'.format(d.charge))
    if 'masses' in d:
        m = d.masses.sum()
    else:
        m = atomic_masses[d.numbers].sum()
    print('mass: {0:.3f} au'.format(m))
    if d.get('keywords'):
        print('keywords: ', ', '.join(d.keywords))
    kvp = d.get('key_value_pairs')
    if kvp:
        print('key-value pairs:')
        for key in sorted(kvp):
            print('    {0}: {1}'.format(key, kvp[key]))
    print('unique id:', d.unique_id)
