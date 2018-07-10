"""SQLite3 backend.

Versions:

1) Added 3 more columns.
2) Changed "user" to "username".
3) Now adding keys to keyword table and added an "information" table containing
   a version number.
4) Got rid of keywords.
5) Add fmax, smax, mass, volume, charge
6) Use REAL for magmom and drop possibility for non-collinear spin
7) Volume can be None
8) Added name='metadata' row to "information" table
"""

from __future__ import absolute_import, print_function
import json
import numbers
import os
import sqlite3
import sys
import functools

import numpy as np

from ase.data import atomic_numbers
from ase.db.row import AtomsRow
from ase.db.core import Database, ops, now, lock, invop, parse_selection
from ase.io.jsonio import encode, numpyfy, mydecode
from ase.parallel import parallel_function
from ase.utils import basestring

if sys.version >= '3':
    buffer = memoryview

VERSION = 8

init_statements = [
    """CREATE TABLE systems (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- ID's, timestamps and user name
    unique_id TEXT UNIQUE,
    ctime REAL,
    mtime REAL,
    username TEXT,
    numbers BLOB,  -- stuff that defines an Atoms object
    positions BLOB,
    cell BLOB,
    pbc INTEGER,
    initial_magmoms BLOB,
    initial_charges BLOB,
    masses BLOB,
    tags BLOB,
    momenta BLOB,
    constraints TEXT,  -- constraints and calculator
    calculator TEXT,
    calculator_parameters TEXT,
    energy REAL,  -- calculated properties
    free_energy REAL,
    forces BLOB,
    stress BLOB,
    dipole BLOB,
    magmoms BLOB,
    magmom REAL,
    charges BLOB,
    key_value_pairs TEXT,  -- key-value pairs and data as json
    data TEXT,
    natoms INTEGER,  -- stuff for making queries faster
    fmax REAL,
    smax REAL,
    volume REAL,
    mass REAL,
    charge REAL)""",

    """CREATE TABLE species (
    Z INTEGER,
    n INTEGER,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",

    """CREATE TABLE keys (
    key TEXT,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",

    """CREATE TABLE text_key_values (
    key TEXT,
    value TEXT,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",

    """CREATE TABLE number_key_values (
    key TEXT,
    value REAL,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",

    """CREATE TABLE information (
    name TEXT,
    value TEXT)""",

    "INSERT INTO information VALUES ('version', '{}')".format(VERSION)]

index_statements = [
    'CREATE INDEX unique_id_index ON systems(unique_id)',
    'CREATE INDEX ctime_index ON systems(ctime)',
    'CREATE INDEX username_index ON systems(username)',
    'CREATE INDEX calculator_index ON systems(calculator)',
    'CREATE INDEX species_index ON species(Z)',
    'CREATE INDEX key_index ON keys(key)',
    'CREATE INDEX text_index ON text_key_values(key)',
    'CREATE INDEX number_index ON number_key_values(key)']

all_tables = ['systems', 'species', 'keys',
              'text_key_values', 'number_key_values']


class SQLite3Database(Database, object):
    type = 'db'
    initialized = False
    _allow_reading_old_format = False
    default = 'NULL'  # used for autoincrement id
    connection = None
    version = None
    columnnames = [line.split()[0].lstrip()
                   for line in init_statements[0].splitlines()[1:]]

    def _connect(self):
        return sqlite3.connect(self.filename, timeout=600)

    def __enter__(self):
        assert self.connection is None
        self.connection = self._connect()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None

    def _initialize(self, con):
        if self.initialized:
            return

        self._metadata = {}

        cur = con.execute(
            'SELECT COUNT(*) FROM sqlite_master WHERE name="systems"')

        if cur.fetchone()[0] == 0:
            for statement in init_statements:
                con.execute(statement)
            if self.create_indices:
                for statement in index_statements:
                    con.execute(statement)
            con.commit()
            self.version = VERSION
        else:
            cur = con.execute(
                'SELECT COUNT(*) FROM sqlite_master WHERE name="user_index"')
            if cur.fetchone()[0] == 1:
                # Old version with "user" instead of "username" column
                self.version = 1
            else:
                try:
                    cur = con.execute(
                        'SELECT value FROM information WHERE name="version"')
                except sqlite3.OperationalError:
                    self.version = 2
                else:
                    self.version = int(cur.fetchone()[0])

                cur = con.execute(
                    'SELECT value FROM information WHERE name="metadata"')
                results = cur.fetchall()
                if results:
                    self._metadata = json.loads(results[0][0])

        if self.version > VERSION:
            raise IOError('Can not read new ase.db format '
                          '(version {}).  Please update to latest ASE.'
                          .format(self.version))
        if self.version < 5 and not self._allow_reading_old_format:
            raise IOError('Please convert to new format. ' +
                          'Use: python -m ase.db.convert ' + self.filename)

        self.initialized = True


    def _get_values(self, atoms, key_value_pairs=None, data=None, id=None):
        if self.type == 'postgresql':
            encode = functools.partial(_encode, pg=True)
        else:
            encode = _encode

        mtime = now()

        if not isinstance(atoms, AtomsRow):
            row = AtomsRow(atoms)
            row.ctime = mtime
            row.user = os.getenv('USER')
        else:
            row = atoms

        constraints = row._constraints
        if constraints:
            if isinstance(constraints, list):
                constraints = encode(constraints)
        else:
            constraints = None

        if self.type == 'postgresql':
            blob = functools.partial(_blob, pg=True)
        else:
            blob = _blob

        values = (row.unique_id,
                  row.ctime,
                  mtime,
                  row.user,
                  blob(row.numbers),
                  blob(row.positions),
                  blob(row.cell),
                  int(np.dot(row.pbc, [1, 2, 4])),
                  blob(row.get('initial_magmoms')),
                  blob(row.get('initial_charges')),
                  blob(row.get('masses')),
                  blob(row.get('tags')),
                  blob(row.get('momenta')),
                  constraints)

        if 'calculator' in row:
            values += (row.calculator, encode(row.calculator_parameters))
        else:
            values += (None, None)

        if not key_value_pairs and not id:
            key_value_pairs = row.key_value_pairs
        if not data and not id:
            data = row._data
        if not isinstance(data, basestring):
            data = encode(data)

        values += (row.get('energy'),
                   row.get('free_energy'),
                   blob(row.get('forces')),
                   blob(row.get('stress')),
                   blob(row.get('dipole')),
                   blob(row.get('magmoms')),
                   row.get('magmom'),
                   blob(row.get('charges')),
                   encode(key_value_pairs),
                   data,
                   len(row.numbers),
                   float_if_not_none(row.get('fmax')),
                   float_if_not_none(row.get('smax')),
                   float_if_not_none(row.get('volume')),
                   float(row.mass),
                   float(row.charge))

        count = row.count_atoms()
        return values, count, key_value_pairs

    def _write(self, atoms, key_value_pairs, data, id):
        Database._write(self, atoms, key_value_pairs, data)

        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        text_key_values = []
        number_key_values = []
        if id:
            self._delete(cur, [id], ['keys', 'text_key_values',
                                     'number_key_values', 'species'])
        values, count, key_value_pairs \
            = self._get_values(atoms, key_value_pairs, data, id)

        if id is None:
            q = self.default + ', ' + ', '.join('?' * len(values))
            cur.execute("""INSERT INTO systems VALUES ({})""".format(q), values)
            id = self.get_last_id(cur)
        else:
            q = ', '.join(name + '=?' for name in self.columnnames[1:])
            cur.execute('UPDATE systems SET {} WHERE id=?'.format(q),
                        values + (id,))

        if count:
            species = [(atomic_numbers[symbol], n, id)
                       for symbol, n in count.items()]
            cur.executemany('INSERT INTO species VALUES (?, ?, ?)',
                            species)

        text_key_values = []
        number_key_values = []
        for key, value in key_value_pairs.items():
            if isinstance(value, (numbers.Real, np.bool_)):
                number_key_values.append([key, float(value), id])
            else:
                assert isinstance(value, basestring)
                text_key_values.append([key, value, id])

        cur.executemany('INSERT INTO text_key_values VALUES (?, ?, ?)',
                        text_key_values)
        cur.executemany('INSERT INTO number_key_values VALUES (?, ?, ?)',
                        number_key_values)
        cur.executemany('INSERT INTO keys VALUES (?, ?)',
                        [(key, id) for key in key_value_pairs])

        if self.connection is None:
            con.commit()
            con.close()

        return id

    def _writemany(self, atomslist):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        assert isinstance(atomslist, list), 'Please pass a list of AtomsRows'
        assert isinstance(atomslist[0], AtomsRow), \
            'Please pass a list of AtomsRows'

        values_collect = []
        species = []
        text_key_values = []
        number_key_values = []
        keys = []
        for i, atoms in enumerate(atomslist):
            values, count, key_value_pairs = self._get_values(atoms)
            values_collect += [values]
            if count:
                species += [[atomic_numbers[symbol], n, i]
                            for symbol, n in count.items()]

            for key, value in key_value_pairs.items():
                keys.append([key, i])
                if isinstance(value, (numbers.Real, np.bool_)):
                    number_key_values.append([key, float(value), i])
                else:
                    assert isinstance(value, basestring)
                    text_key_values.append([key, value, i])

        N_rows = len(values_collect)
        statement = """INSERT INTO systems VALUES ({})"""

        last_id = self.get_last_id(cur)
        q = self.default + ', ' + ', '.join('?' * len(values[0]))
        if self.type == 'postgresql':
            statement += ' RETURNING id;'

        cur.executemany(statement.format(q), values_collect)

        if self.type == 'postgresql':
            ids = [id[0] for id in cur.fetchall()]
        else:
            last_id = self.get_last_id(cur)
            ids = range(last_id + 1 - N_rows, last_id + 1)

        # Update with id from systems
        if len(ids) == 0:
            if self.connection is None:
                con.commit()
                con.close()
            return ids
        
        for spec in species:
            spec[2] = ids[spec[2]]
            spec = tuple(spec)
        for tkv in text_key_values:
            tkv[2] = ids[tkv[2]]
            tkv = tuple(tkv)
        for nkv in number_key_values:
            nkv[2] = ids[nkv[2]]
            nkv = tuple(nkv)
        for key in keys:
            key[1] = ids[key[1]]
            key = tuple(key)

        cur.executemany('INSERT INTO species VALUES (?, ?, ?)',
                        species)
        cur.executemany('INSERT INTO text_key_values VALUES (?, ?, ?)',
                        text_key_values)
        cur.executemany('INSERT INTO number_key_values VALUES (?, ?, ?)',
                        number_key_values)
        cur.executemany('INSERT INTO keys VALUES (?, ?)',
                        keys)

        if self.connection is None:
            con.commit()
            con.close()

        return ids[-1]

    def get_last_id(self, cur):
        cur.execute('SELECT seq FROM sqlite_sequence WHERE name="systems"')
        result = cur.fetchone()
        if result is not None:
            id = result[0]
            return id
        else:
            return 0

    def _get_row(self, id):
        con = self._connect()
        self._initialize(con)
        c = con.cursor()
        if id is None:
            c.execute('SELECT COUNT(*) FROM systems')
            assert c.fetchone()[0] == 1
            c.execute('SELECT * FROM systems')
        else:
            c.execute('SELECT * FROM systems WHERE id=?', (id,))
        values = c.fetchone()

        values = self._old2new(values)
        return self._convert_tuple_to_row(values)

    def _convert_tuple_to_row(self, values):
        if self.type == 'postgresql':
            deblob = functools.partial(_deblob, pg=True)
            decode = functools.partial(_decode, pg=True)
        else:
            deblob = _deblob
            decode = _decode

        cn = self.columnnames
        values = self._old2new(values)
        dct = {'id': values[cn.index('id')],
               'unique_id': values[cn.index('unique_id')],
               'ctime': values[cn.index('ctime')],
               'mtime': values[cn.index('mtime')],
               'user': values[cn.index('username')],
               'numbers': deblob(values[cn.index('numbers')], np.int32),
               'positions': deblob(values[cn.index('positions')], shape=(-1, 3)),
               'cell': deblob(values[cn.index('cell')], shape=(3, 3))}

        if values[cn.index('pbc')] is not None:
            dct['pbc'] = (values[cn.index('pbc')] &
                          np.array([1, 2, 4])).astype(bool)
        if values[cn.index('initial_magmoms')] is not None:
            dct['initial_magmoms'] = deblob(values[cn.index('initial_magmoms')])
        if values[cn.index('initial_charges')] is not None:
            dct['initial_charges'] = deblob(values[cn.index('initial_charges')])
        if values[cn.index('masses')] is not None:
            dct['masses'] = deblob(values[cn.index('masses')])
        if values[cn.index('tags')] is not None:
            dct['tags'] = deblob(values[cn.index('tags')], np.int32)
        if values[cn.index('momenta')] is not None:
            dct['momenta'] = deblob(values[cn.index('momenta')], shape=(-1, 3))
        if values[cn.index('constraints')] is not None:
            dct['constraints'] = values[cn.index('constraints')]
        if values[cn.index('calculator')] is not None:
            dct['calculator'] = values[cn.index('calculator')]
        if values[cn.index('calculator_parameters')] is not None:
            dct['calculator_parameters'] = \
            decode(values[cn.index('calculator_parameters')])
        if values[cn.index('energy')] is not None:
            dct['energy'] = values[cn.index('energy')]
        if values[cn.index('free_energy')] is not None:
            dct['free_energy'] = values[cn.index('free_energy')]
        if values[cn.index('forces')] is not None:
            dct['forces'] = deblob(values[cn.index('forces')], shape=(-1, 3))
        if values[cn.index('stress')] is not None:
            dct['stress'] = deblob(values[cn.index('stress')])
        if values[cn.index('dipole')] is not None:
            dct['dipole'] = deblob(values[cn.index('dipole')])
        if values[cn.index('magmoms')] is not None:
            dct['magmoms'] = deblob(values[cn.index('magmoms')])
        if values[cn.index('magmom')] is not None:
            dct['magmom'] = values[cn.index('magmom')]
        if values[cn.index('charges')] is not None:
            dct['charges'] = deblob(values[cn.index('charges')])
        if values[cn.index('key_value_pairs')] != '{}':
            dct['key_value_pairs'] = decode(values[cn.index('key_value_pairs')])
        if len(values) >= cn.index('data') + 1 and \
        values[cn.index('data')] != 'null':
            dct['data'] = decode(values[cn.index('data')])
        return AtomsRow(dct)

    def _old2new(self, values):
        if self.type == 'postgresql':
            assert self.version >= 8, 'Your db-server is too old!'
        assert self.version >= 4, 'Your db-file is too old!'
        if self.version < 5:
            pass  # should be ok for reading by convert.py script
        if self.version < 6:
            m = values[23]
            if m is not None and not isinstance(m, float):
                magmom = float(_deblob(m, shape=()))
                values = values[:23] + (magmom,) + values[24:]
        return values

    def create_select_statement(self, keys, cmps,
                                sort=None, order=None, sort_table=None,
                                what='systems.*'):
        tables = ['systems']
        where = []
        args = []
        for key in keys:
            if key == 'forces':
                where.append('systems.fmax IS NOT NULL')
            elif key == 'strain':
                where.append('systems.smax IS NOT NULL')
            elif key in ['energy', 'fmax', 'smax',
                         'constraints', 'calculator']:
                where.append('systems.{} IS NOT NULL'.format(key))
            else:
                if '-' not in key:
                    q = 'systems.id in (select id from keys where key=?)'
                else:
                    key = key.replace('-', '')
                    q = 'systems.id not in (select id from keys where key=?)'
                where.append(q)
                args.append(key)

        # Special handling of "H=0" and "H<2" type of selections:
        bad = {}
        for key, op, value in cmps:
            if isinstance(key, int):
                bad[key] = bad.get(key, True) and ops[op](0, value)

        for key, op, value in cmps:
            if key in ['id', 'energy', 'magmom', 'ctime', 'user',
                       'calculator', 'natoms', 'pbc', 'unique_id',
                       'fmax', 'smax', 'volume', 'mass', 'charge']:
                if key == 'user' and self.version >= 2:
                    key = 'username'
                elif key == 'pbc':
                    assert op in ['=', '!=']
                    value = int(np.dot([x == 'T' for x in value], [1, 2, 4]))
                elif key == 'magmom':
                    assert self.version >= 6, 'Update your db-file'
                where.append('systems.{}{}?'.format(key, op))
                args.append(value)
            elif isinstance(key, int):
                if self.type == 'postgresql':
                    where.append(
                        'cardinality(array_positions(' +
                        'numbers::int[], ?)){}?'.format(op))
                    args += [key, value]
                else:
                    if bad[key]:
                        where.append('systems.id not in (select id from species ' +
                                     'where Z=? and n{}?)'.format(invop[op]))
                        args += [key, value]
                    else:
                        where.append('systems.id in (select id from species ' +
                                     'where Z=? and n{}?)'.format(op))
                        args += [key, value]

            elif self.type == 'postgresql':
                jsonop = '->'
                if isinstance(value, basestring):
                    jsonop = '->>'
                where.append("systems.key_value_pairs {} '{}'{}?"
                             .format(jsonop, key, op))
                args.append(str(value))

            elif isinstance(value, basestring):
                where.append('systems.id in (select id from text_key_values ' +
                             'where key=? and value{}?)'.format(op))
                args += [key, value]
            else:
                where.append('systems.id in (select id from number_key_values ' +
                             'where key=? and value{}?)'.format(op))
                args += [key, float(value)]

        if sort:
            if sort_table != 'systems':
                tables.append('{} AS sort_table'.format(sort_table))
                where.append('systems.id=sort_table.id AND '
                             'sort_table.key=?')
                args.append(sort)
                sort_table = 'sort_table'
                sort = 'value'

        sql = 'SELECT {} FROM\n  '.format(what) + ', '.join(tables)
        if where:
            sql += '\n  WHERE\n  ' + ' AND\n  '.join(where)
        if sort:
            # XXX use "?" instead of "{}"
            sql += '\nORDER BY {0}.{1} IS NULL, {0}.{1} {2}'.format(
                sort_table, sort, order)

        return sql, args

    def _select(self, keys, cmps, explain=False, verbosity=0,
                limit=None, offset=0, sort=None, include_data=True,
                columns='all'):
        con = self._connect()
        self._initialize(con)

        n_values = self.columnnames.index('data') + 1
        values = np.array([None for i in range(n_values)])
        values[self.columnnames.index('key_value_pairs')] = '{}'
        values[self.columnnames.index('data')] = 'null'
        if columns == 'all':
            columnindex = list(range(n_values))
        else:
            columnindex = [c for c in range(n_values)
                           if self.columnnames[c] in columns]

        if not include_data:
            if self.columnnames.index('data') in columnindex:
                columnindex.remove(self.columnnames.index('data'))

        if sort:
            if sort[0] == '-':
                order = 'DESC'
                sort = sort[1:]
            else:
                order = 'ASC'
            if sort in ['id', 'energy', 'username', 'calculator',
                        'ctime', 'mtime', 'magmom', 'pbc',
                        'fmax', 'smax', 'volume', 'mass', 'charge', 'natoms']:
                sort_table = 'systems'
            else:
                for dct in self._select(keys + [sort], cmps=[], limit=1,
                                        include_data=False,
                                        columns=['key_value_pairs']):
                    if isinstance(dct['key_value_pairs'][sort], basestring):
                        sort_table = 'text_key_values'
                    else:
                        sort_table = 'number_key_values'
                    break
                else:
                    # No rows.  Just pick a table:
                    sort_table = 'number_key_values'

        else:
            order = None
            sort_table = None

        what = ', '.join('systems.' + name
                         for name in np.array(self.columnnames)
                         [np.array(columnindex)])

        sql, args = self.create_select_statement(keys, cmps, sort, order,
                                                 sort_table, what)

        if explain:
            sql = 'EXPLAIN QUERY PLAN ' + sql

        if limit:
            sql += '\nLIMIT {0}'.format(limit)

        if offset:
            sql += '\nOFFSET {0}'.format(offset)

        if verbosity == 2:
            print(sql, args)

        cur = con.cursor()
        cur.execute(sql, args)
        if explain:
            for row in cur.fetchall():
                yield {'explain': row}
        else:
            n = 0
            for shortvalues in cur.fetchall():
                values[columnindex] = shortvalues
                yield self._convert_tuple_to_row(tuple(values))
                n += 1

            if sort and sort_table != 'systems':
                # Yield rows without sort key last:
                if limit is not None:
                    if n == limit:
                        return
                    limit -= n
                for row in self._select(keys + ['-' + sort], cmps,
                                        limit=limit, offset=offset,
                                        include_data=include_data,
                                        columns=columns):
                    yield row

    @parallel_function
    def count(self, selection=None, **kwargs):
        keys, cmps = parse_selection(selection, **kwargs)
        sql, args = self.create_select_statement(keys, cmps, what='COUNT(*)')
        con = self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute(sql, args)
        return cur.fetchone()[0]

    def analyse(self):
        con = self._connect()
        self._initialize(con)
        con.execute('ANALYZE')

    @parallel_function
    @lock
    def delete(self, ids):
        if len(ids) == 0:
            return
        con = self._connect()
        self._delete(con.cursor(), ids)
        con.commit()
        con.close()

    def _delete(self, cur, ids, tables=None):
        tables = tables or all_tables[::-1]
        for table in tables:
            cur.execute('DELETE FROM {} WHERE id in ({});'.
                        format(table, ', '.join([str(id) for id in ids])))

    @property
    def metadata(self):
        if self._metadata is None:
            self._initialize(self._connect())
        return self._metadata.copy()

    @metadata.setter
    def metadata(self, dct):
        self._metadata = dct
        con = self._connect()
        self._initialize(con)
        md = json.dumps(dct)
        cur = con.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM information WHERE name='metadata'")

        if cur.fetchone()[0]:
            cur.execute(
                "UPDATE information SET value=? WHERE name='metadata'", [md])
        else:
            cur.execute('INSERT INTO information VALUES (?, ?)',
                        ('metadata', md))
        con.commit()


def float_if_not_none(x):
    """Convert numpy.float64 to float - old db-interfaces need that."""
    if x is not None:
        return float(x)


def _blob(array, pg=False):
    """Convert array to blob/buffer object."""

    if array is None:
        return None
    if len(array) == 0:
        array = np.zeros(0)
    if array.dtype == np.int64:
        array = array.astype(np.int32)
    if pg:
        return array.tolist()
    if not np.little_endian:
        array = array.byteswap()
    return buffer(np.ascontiguousarray(array))


def _deblob(buf, dtype=float, shape=None, pg=False):
    """Convert blob/buffer object to ndarray of correct dtype and shape.

    (without creating an extra view)."""
    if buf is None:
        return None
    if pg:
        return np.array(buf, dtype=dtype)
    if len(buf) == 0:
        array = np.zeros(0, dtype)
    else:
        if len(buf) % 2 == 1:
            # old psycopg2:
            array = np.fromstring(str(buf)[1:].decode('hex'), dtype)
        else:
            array = np.frombuffer(buf, dtype)
        if not np.little_endian:
            array = array.byteswap()
    if shape is not None:
        array.shape = shape
    return array


def _encode(obj, pg=False):
    if pg:
        return encode(obj).replace('NaN', '"NaN"').replace('Infinity', '"Infinity"')
    else:
        return encode(obj)


def _decode(txt, pg=False):
    if pg:
        txt = encode(txt).replace('"NaN"', 'NaN').replace('"Infinity"', 'Infinity')
    return numpyfy(mydecode(txt))


if __name__ == '__main__':
    import sys
    from ase.db import connect
    con = connect(sys.argv[1])
    con._initialize(con._connect())
    print('Version:', con.version)
