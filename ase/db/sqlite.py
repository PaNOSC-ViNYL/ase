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

import numpy as np

from ase.data import atomic_numbers
from ase.db.row import AtomsRow
from ase.db.core import Database, ops, now, lock, invop, parse_selection
from ase.io.jsonio import encode, numpyfy, mydecode, object_hook, read_json
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


def float_if_not_none(x):
    """Convert numpy.float64 to float - old db-interfaces need that."""
    if x is not None:
        return float(x)


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

    def _write(self, atoms, key_value_pairs, data, id):
        Database._write(self, atoms, key_value_pairs, data)
        
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        pg = False
        if self.type == 'postgresql':
            pg = True

        mtime = now()

        if isinstance(atoms, list):
            assert isinstance(atoms[0], AtomsRow)
            assert not key_value_pairs,  'Please pass single row to write key_value_pairs'
            assert not data,  'Please pass single row to write data'
        else:
            atoms = [atoms]

        values_collect = []
        species = []
        text_key_values = []
        number_key_values = []
        keys = []
        for i, atoms in enumerate(atoms):           
            if not isinstance(atoms, AtomsRow):
                row = AtomsRow(atoms)
                row.ctime = mtime
                row.user = os.getenv('USER')
            else:
                row = atoms

            if id:
                self._delete(cur, [id], ['keys', 'text_key_values',
                                         'number_key_values', 'species'])

            constraints = row._constraints
            if constraints:
                if isinstance(constraints, list):
                    constraints = encode(constraints)
            else:
                constraints = None

            values = (row.unique_id,
                      row.ctime,
                      mtime,
                      row.user,
                      blob(row.numbers, pg),
                      blob(row.positions, pg),
                      blob(row.cell, pg),
                      int(np.dot(row.pbc, [1, 2, 4])),
                      blob(row.get('initial_magmoms'), pg),
                      blob(row.get('initial_charges'), pg),
                      blob(row.get('masses'), pg),
                      blob(row.get('tags'), pg),
                      blob(row.get('momenta'), pg),
                      constraints)

            if 'calculator' in row:
                values += (row.calculator, encode(row.calculator_parameters))
            else:
                values += (None, None)

            if not id:
                if not key_value_pairs or i > 0:
                    key_value_pairs = row.key_value_pairs
            if not data or i > 0:
                data = row._data
            if not isinstance(data, basestring):
                data = encode(data)
            
            values += (row.get('energy'),
                       row.get('free_energy'),
                       blob(row.get('forces'), pg),
                       blob(row.get('stress'), pg),
                       blob(row.get('dipole'), pg),
                       blob(row.get('magmoms'), pg),
                       row.get('magmom'),
                       blob(row.get('charges'), pg),
                       encode(key_value_pairs),
                       data,
                       len(row.numbers),
                       float_if_not_none(row.get('fmax')),
                       float_if_not_none(row.get('smax')),
                       float_if_not_none(row.get('volume')),
                       float(row.mass),
                       float(row.charge))      
            
            if id is None:
                values_collect += [values]
            elif row is not None:
                q = ', '.join(name + '=?' for name in self.columnnames[1:])
                cur.execute('UPDATE systems SET {} WHERE id=?'.format(q),
                            values + (id,))
                ids = [id]

            count = row.count_atoms()
            if count:
                species += [[atomic_numbers[symbol], n, i] #id
                            for symbol, n in count.items()]

            for key, value in key_value_pairs.items():
                keys.append([key, i])
                
                if isinstance(value, (numbers.Real, np.bool_)):
                    number_key_values.append([key, float(value), i])
                else:
                    assert isinstance(value, basestring)
                    text_key_values.append([key, value, i])


        N_rows = len(values_collect)
        statement = 'INSERT INTO systems VALUES ({})'
        if N_rows == 1:  # One row
            q = self.default + ', ' + ', '.join('?' * len(values[0]))
            cur.execute(statement.format(q), values)
            ids = [self.get_last_id(cur)]

        elif N_rows > 1:  # Several rows
            last_id = self.get_last_id(cur)
            q = self.default + ', ' + ', '.join('?' * len(values[0]))
            if self.type == 'postgresql':
                statement += ' returning id'

            cur.executemany(statement.format(q), values_collect)

            if self.type == 'postgresql':
                ids = cur.fetchall()
                ids = [ids[i][0] for i in range(len(ids))]
            else:
                last_id = self.get_last_id(cur)
                ids =  range(last_id + 1 - N_rows, last_id + 1)

        # Update with id from systems
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
        pg = False
        if self.type == 'postgresql':
            pg = True

        values = self._old2new(values)
        dct = {'id': values[0],
               'unique_id': values[1],
               'ctime': values[2],
               'mtime': values[3],
               'user': values[4],
               'numbers': deblob(values[5], np.int32, pg=pg),
               'positions': deblob(values[6], shape=(-1, 3), pg=pg),
               'cell': deblob(values[7], shape=(3, 3), pg=pg)}

        if values[8] is not None:
            dct['pbc'] = (values[8] & np.array([1, 2, 4])).astype(bool)
        if values[9] is not None:
            dct['initial_magmoms'] = deblob(values[9], pg=pg)
        if values[10] is not None:
            dct['initial_charges'] = deblob(values[10], pg=pg)
        if values[11] is not None:
            dct['masses'] = deblob(values[11], pg=pg)
        if values[12] is not None:
            dct['tags'] = deblob(values[12], np.int32, pg=pg)
        if values[13] is not None:
            dct['momenta'] = deblob(values[13], shape=(-1, 3), pg=pg)
        if values[14] is not None:
            dct['constraints'] = values[14]
        if values[15] is not None:
            dct['calculator'] = values[15]
        if values[16] is not None:
            dct['calculator_parameters'] = decode(values[16], pg=pg)
        if values[17] is not None:
            dct['energy'] = values[17]
        if values[18] is not None:
            dct['free_energy'] = values[18]
        if values[19] is not None:
            dct['forces'] = deblob(values[19], shape=(-1, 3), pg=pg)
        if values[20] is not None:
            dct['stress'] = deblob(values[20], pg=pg)
        if values[21] is not None:
            dct['dipole'] = deblob(values[21], pg=pg)
        if values[22] is not None:
            dct['magmoms'] = deblob(values[22], pg=pg)
        if values[23] is not None:
            dct['magmom'] = values[23]
        if values[24] is not None:
            dct['charges'] = deblob(values[24], pg=pg)
        if values[25] != '{}':
            dct['key_value_pairs'] = decode(values[25], pg=pg)
        if len(values) >= 27 and values[26] != 'null':
            dct['data'] = decode(values[26], pg=pg)
            
        return AtomsRow(dct)

    def _old2new(self, values):
        assert self.version >= 4, 'Your db-file is too old!'
        if self.version < 5:
            pass  # should be ok for reading by convert.py script
        if self.version < 6:
            m = values[23]
            if m is not None and not isinstance(m, float):
                magmom = float(deblob(m, shape=()), pg=pg)
                values = values[:23] + (magmom,) + values[24:]
        return values

    def create_select_statement(self, keys, cmps,
                                sort=None, order=None, sort_table=None,
                                what='systems.*'):
        tables = ['systems']
        where = []
        args = []

        for n, key in enumerate(keys):
            if key == 'forces':
                where.append('systems.fmax IS NOT NULL')
            elif key == 'strain':
                where.append('systems.smax IS NOT NULL')
            elif key in ['energy', 'fmax', 'smax',
                         'constraints', 'calculator']:
                where.append('systems.{} IS NOT NULL'.format(key))
            else:
                if '-' not in key:
                    tables.append('keys AS keys{}'.format(n))
                    q = 'systems.id=keys{0}.id AND keys{0}.key=?'.format(n)
                else:
                    key = key.replace('-', '')
                    q = 'systems.id=keys.id AND keys.key=?'
                    q = 'NOT EXISTS (SELECT id FROM keys WHERE {})'.format(q)
                where.append(q)
                args.append(key)

        # Special handling of "H=0" and "H<2" type of selections:
        bad = {}
        for key, op, value in cmps:
            if isinstance(key, int):
                bad[key] = bad.get(key, True) and ops[op](0, value)

        found_sort_table = False
        nspecies = 0
        ntext = 0
        nnumber = 0
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
                        where.append(
                            'NOT EXISTS (SELECT id FROM species WHERE\n' +
                            '  species.id=systems.id AND species.Z=? AND ' +
                            'species.n{}?)'.format(invop[op]))
                        args += [key, value]
                    else:
                        tables.append('species AS specie{}'.format(nspecies))
                        where.append(('systems.id=specie{0}.id AND ' +
                                      'specie{0}.Z=? AND ' +
                                      'specie{0}.n{1}?').format(nspecies, op))
                        args += [key, value]
                        nspecies += 1

            elif self.type == 'postgresql':
                jsonop = '->'
                if isinstance(value, str):
                    jsonop = '->>'
                where.append("systems.key_value_pairs {} '{}'{}?".format(jsonop, key, op))
                args.append(str(value))

            elif isinstance(value, basestring):
                tables.append('text_key_values AS text{0}'.format(ntext))
                where.append(('systems.id=text{0}.id AND ' +
                              'text{0}.key=? AND ' +
                              'text{0}.value{1}?').format(ntext, op))
                args += [key, value]
                if sort_table == 'text_key_values' and sort == key:
                    sort_table = 'text{0}'.format(ntext)
                    found_sort_table = True
                ntext += 1
            else:
                tables.append('number_key_values AS number{}'.format(nnumber))
                where.append(('systems.id=number{0}.id AND ' +
                              'number{0}.key=? AND ' +
                              'number{0}.value{1}?').format(nnumber, op))
                args += [key, float(value)]
                if sort_table == 'number_key_values' and sort == key:
                    sort_table = 'number{}'.format(nnumber)
                    found_sort_table = True
                nnumber += 1

        if sort:
            if sort_table != 'systems':
                if not found_sort_table:
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

        values = np.array([None for i in range(27)])
        values[25] = '{}'
        values[26] = 'null'

        if columns == 'all':
            columnindex = list(range(27))
        else:
            columnindex = [c for c in range(0, 27) if self.columnnames[c] in columns]

        if not include_data:
            if 26 in columnindex:
                columnindex.remove(26)
        #['id', 'ctime', 'mtime', 'username', 'numbers', 'pbc', 'charge', 'mass']

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
                for dct in self._select(keys + [sort], cmps, limit=1,
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
                         for name in np.array(self.columnnames)[np.array(columnindex)])
        #if not include_data:
        #    what = 'systems.*'
        #else:
        #    #columnlist = np.array([0, 2, 3, 4, 5, 7, 8, 11, 25])
        #    what = ', '.join('systems.' + name
        #                     for name in np.array(self.columnnames)[columnlist])#self.columnnames[:26])

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
                                        columns=['id', 'key_value_pairs']):
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
        con = self._connect()
        self._delete(con.cursor(), ids)
        con.commit()
        con.close()

    def _delete(self, cur, ids, tables=None):
        tables = tables or all_tables[::-1]
        for table in tables:
            cur.executemany('DELETE FROM {0} WHERE id=?'.format(table),
                            [(id,) for id in ids])

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


def blob(array, postgresql=False):
    """Convert array to blob/buffer object."""

    if array is None:
        return None
    if len(array) == 0:
        array = np.zeros(0)
    if array.dtype == np.int64:
        array = array.astype(np.int32)
    if postgresql:
        return array.tolist()
    if not np.little_endian:
        array.byteswap(True)
    return buffer(np.ascontiguousarray(array))


def deblob(buf, dtype=float, shape=None, pg=False):
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
            array.byteswap(True)
    if shape is not None:
        array.shape = shape
    return array


def decode(txt, pg=False):
    if pg:
        txt = encode(txt)
    return numpyfy(mydecode(txt))


if __name__ == '__main__':
    import sys
    from ase.db import connect
    con = connect(sys.argv[1])
    con._initialize(con._connect())
    print('Version:', con.version)
