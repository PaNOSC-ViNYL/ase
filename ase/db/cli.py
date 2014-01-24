from __future__ import print_function
import sys
import optparse

import numpy as np

import ase.io
from ase.db import connect
from ase.atoms import Atoms
from ase.data import atomic_masses
from ase.calculators.calculator import get_calculator
from ase.db.core import float_to_time_string, now


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)


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


def run(args=sys.argv[1:]):
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
    add('-p', '--python-expression', metavar='expression',
        help="""Examples: "d.key", "d.keywords['keyword']", where """ +
        '"d" is a dictionary representing a row.')

    opts, args = parser.parse_args(args)

    verbosity = 1 - opts.quiet + opts.verbose

    if not args:
        parser.error('No database given')
        
    con = connect(args.pop(0))

    if verbosity == 2:
        print(opts, args)

    if args:
        if len(args) == 1 and args[0].isdigit():
            expressions = int(args[0])
        else:
            expressions = ','.join(args)
    else:
        expressions = []
    
    if opts.count:
        opts.limit = 0
        
    if not opts.add_from_file:
        rows = con.select(expressions, explain=opts.explain,
                          verbosity=verbosity, limit=opts.limit)

    if opts.count:
        n = 0
        for row in rows:
            n += 1
        print('%s' % plural(n, 'row'))
        return

    if opts.explain:
        for row in rows:
            print('%d %d %d %s' % row['explain'])
        return

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

    if opts.insert_into:
        con2 = connect(opts.insert_into)
        nkw = 0
        nkvp = 0
        nrows = 0
        for dct in rows:
            keywords = dct.get('keywords', [])
            for keyword in add_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    nkw += 1

            kvp = dct['key_value_pairs']
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
        
    if add_keywords or add_key_value_pairs:
        ids = [dct['id'] for dct in rows]
        nkw, nkv = con.update(ids, add_keywords, **add_key_value_pairs)
        print('Added %s and %s (%s updated)' %
              (plural(nkw, 'keyword'),
               plural(nkv, 'key-value pair'),
               plural(len(add_key_value_pairs) * len(ids) - nkv, 'pair')))
        return

    if opts.delete:
        ids = [dct['id'] for dct in rows]
        if ids and not opts.yes:
            msg = 'Delete %s? (yes/no): ' % plural(len(ids), 'row')
            if raw_input(msg).lower() != 'yes':
                return
        con.delete(ids)
        print('Deleted %s' % plural(len(ids), 'row'))
        return

    if opts.python_expression:
        for n, dct in enumerate(rows):
            row = eval(opts.python_expression, {'d': dct, 'n': n})
            if not isinstance(row, (list, tuple, np.ndarray)):
                row = [row]
            print(', '.join(str(x) for x in row))
        return
        
    dcts = list(rows)
    if len(dcts) > 0:
        if opts.long:
            long(dcts[0], verbosity)
            return
        
        f = Formatter()
        f.format(dcts)


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


def cut(txt, length):
    if len(txt) <= length:
        return txt
    return txt[:length - 3] + '...'


class Formatter:
    def format(self, dcts, columns=None, sort=None):
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
                else:
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
        return float_to_time_string(now() - d.ctime)

    def user(self, d):
        return d.user
    
    def formula(self, d):
        return Atoms(d.numbers).get_chemical_formula()

    def energy(self, d):
        return d.energy

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
        return d.calculator

    def fmax(self, d):
        return (d.forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return cut(','.join(d.keywords), 30)

    def keyvals(self, d):
        return cut(','.join(['%s=%s' % (key, cut(str(value), 8))
                             for key, value in d.key_value_pairs.items()]), 40)

    def charge(self, d):
        return d.charge

    def mass(self, d):
        if 'masses' in d:
            return d.masses.sum()
        return atomic_masses[d.numbers].sum()

    def fixed(self, d):
        c = d.constraints
        if c is None:
            return ''
        if len(c) > 1:
            return '?'
        c = c[0]
        if 'mask' in c:
            return sum(c['mask'])
        return len(c['indices'])

    def smax(self, d):
        return (d.stress**2).max()**0.5

    def magmom(self, d):
        return d.magmom
