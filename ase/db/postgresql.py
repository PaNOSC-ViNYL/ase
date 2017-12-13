import numpy as np
import json
import psycopg2

from ase.data import atomic_numbers
from ase.db.sqlite import VERSION #init_statements, index_statements, VERSION
from ase.db.sqlite import SQLite3Database, float_if_not_none
from ase.db.core import Database, now
from ase.db.row import AtomsRow
from ase.io.jsonio import encode, decode


init_statements = [
    """CREATE TABLE systems (
    id SERIAL PRIMARY KEY,  -- ID's, timestamps and user name
    unique_id TEXT UNIQUE,
    ctime DOUBLE PRECISION,
    mtime DOUBLE PRECISION,
    username TEXT,
    numbers INTEGER[],  -- stuff that defines an Atoms object
    positions DOUBLE PRECISION[][],
    cell DOUBLE PRECISION[][],
    pbc INTEGER,
    initial_magmoms DOUBLE PRECISION[],
    initial_charges DOUBLE PRECISION[],
    masses DOUBLE PRECISION[],
    tags INTEGER[],
    momenta DOUBLE PRECISION[],
    constraints TEXT,  -- constraints and calculator
    calculator TEXT,
    calculator_parameters TEXT,
    energy DOUBLE PRECISION,  -- calculated properties
    free_energy DOUBLE PRECISION,
    forces DOUBLE PRECISION[][],
    stress DOUBLE PRECISION[],
    dipole DOUBLE PRECISION[],
    magmoms DOUBLE PRECISION[],
    magmom DOUBLE PRECISION,
    charges DOUBLE PRECISION[],
    key_value_pairs TEXT,  -- key-value pairs and data as json
    data TEXT,
    natoms INTEGER,  -- stuff for making queries faster
    fmax DOUBLE PRECISION,
    smax DOUBLE PRECISION,
    volume DOUBLE PRECISION,
    mass DOUBLE PRECISION,
    charge DOUBLE PRECISION)""",

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
    value DOUBLE PRECISION,
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


def list_if_not_none(x):
    if x is not None:
        return x.tolist()

class Connection:
    def __init__(self, con):
        self.con = con

    def cursor(self):
        return Cursor(self.con.cursor())

    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()


class Cursor:
    def __init__(self, cur):
        self.cur = cur

    def fetchone(self):
        return self.cur.fetchone()

    def fetchall(self):
        return self.cur.fetchall()

    def execute(self, statement, *args):
        self.cur.execute(statement.replace('?', '%s'), *args)

    def executemany(self, statement, *args):
        self.cur.executemany(statement.replace('?', '%s'), *args)


class PostgreSQLDatabase(SQLite3Database):
    default = 'DEFAULT'

    def _connect(self):
        return Connection(psycopg2.connect(self.filename))

    def _initialize(self, con):
        if self.initialized:
            return

        self._metadata = {}

        cur = con.cursor()

        try:
            cur.execute('SELECT name, value FROM information')
        except psycopg2.ProgrammingError:
            # Initialize database:
            sql = ';\n'.join(init_statements)
            #for a, b in [('BLOB', 'BYTEA'),
            #             ('REAL', 'DOUBLE PRECISION'),
            #             ('INTEGER PRIMARY KEY AUTOINCREMENT',
            #              'SERIAL PRIMARY KEY')]:
            #    sql = sql.replace(a, b)
            
            con.commit()
            cur = con.cursor()
            cur.execute(sql)
            if self.create_indices:
                cur.execute(';\n'.join(index_statements))
            con.commit()
            self.version = VERSION
        else:
            for name, value in cur.fetchall():
                if name == 'version':
                    self.version = int(value)
                elif name == 'metadata':
                    self._metadata = json.loads(value)

        assert 5 < self.version <= VERSION

        self.initialized = True

    def get_last_id(self, cur):
        cur.execute('SELECT last_value FROM systems_id_seq')
        id = cur.fetchone()[0]
        return int(id)


    def _write(self, atoms, key_value_pairs, data):
        Database._write(self, atoms, key_value_pairs, data)
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        id = None

        if not isinstance(atoms, AtomsRow):
            row = AtomsRow(atoms)
            row.ctime = mtime = now()
            row.user = os.getenv('USER')
        else:
            row = atoms
            cur.execute('SELECT id FROM systems WHERE unique_id=?',
                        (row.unique_id,))
            results = cur.fetchall()
            if results:
                id = results[0][0]
                self._delete(cur, [id], ['keys', 'text_key_values',
                                         'number_key_values'])
                
            mtime = now()

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
                  row.numbers.tolist(),
                  row.positions.tolist(),
                  row.cell.tolist(),
                  int(np.dot(row.pbc, [1, 2, 4])),
                  list_if_not_none(row.get('initial_magmoms')),
                  list_if_not_none(row.get('initial_charges')),
                  list_if_not_none(row.get('masses')),
                  list_if_not_none(row.get('tags')),
                  list_if_not_none(row.get('momenta')),
                  constraints)
        if 'calculator' in row:
            values += (row.calculator, encode(row.calculator_parameters))
        else:
            values += (None, None)

        if key_value_pairs is None:
            key_value_pairs = row.key_value_pairs

        if not data:
            data = row._data
        if not isinstance(data, basestring):
            data = encode(data)
        values += (row.get('energy'),
                   row.get('free_energy'),
                   list_if_not_none(row.get('forces')),
                   list_if_not_none(row.get('stress')),
                   list_if_not_none(row.get('dipole')),
                   list_if_not_none(row.get('magmoms')),
                   row.get('magmom'),
                   list_if_not_none(row.get('charges')),
                   encode(key_value_pairs),
                   data,
                   len(row.numbers),
                   float_if_not_none(row.get('fmax')),
                   float_if_not_none(row.get('smax')),
                   float_if_not_none(row.get('volume')),
                   float(row.mass),
                   float(row.charge))


        if id is None:
            q = self.default + ', ' + ', '.join('?' * len(values))
            cur.execute('INSERT INTO systems VALUES ({})'.format(q),
                        values)
        else:
            q = ', '.join(name + '=?' for name in self.columnnames[1:])
            cur.execute('UPDATE systems SET {} WHERE id=?'.format(q),
                        values + (id,))

        if id is None:
            id = self.get_last_id(cur)

            count = row.count_atoms()
            if count:
                species = [(atomic_numbers[symbol], n, id)
                           for symbol, n in count.items()]
                cur.executemany('INSERT INTO species VALUES (?, ?, ?)',
                                species)

        text_key_values = []
        number_key_values = []
        for key, value in key_value_pairs.items():
            if isinstance(value, (float, int, np.bool_)):
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


    def _convert_tuple_to_row(self, values):
        values = self._old2new(values)
        dct = {'id': values[0],
               'unique_id': values[1],
               'ctime': values[2],
               'mtime': values[3],
               'user': values[4],
               'numbers': np.array(values[5], dtype=np.int32),
               'positions': np.array(values[6]),
               'cell': np.array(values[7]),
               'pbc': (values[8] & np.array([1, 2, 4])).astype(bool)}
        if values[9] is not None:
            dct['initial_magmoms'] = np.array(values[9])
        if values[10] is not None:
            dct['initial_charges'] = np.array(values[10])
        if values[11] is not None:
            dct['masses'] = np.array(values[11])
        if values[12] is not None:
            dct['tags'] = np.array(values[12], dtype=np.int32)
        if values[13] is not None:
            dct['momenta'] = np.array(values[13])
        if values[14] is not None:
            dct['constraints'] = values[14]
        if values[15] is not None:
            dct['calculator'] = values[15]
            dct['calculator_parameters'] = decode(values[16])
        if values[17] is not None:
            dct['energy'] = values[17]
        if values[18] is not None:
            dct['free_energy'] = values[18]
        if values[19] is not None:
            dct['forces'] = np.array(values[19])
        if values[20] is not None:
            dct['stress'] = np.array(values[20])
        if values[21] is not None:
            dct['dipole'] = np.array(values[21])
        if values[22] is not None:
            dct['magmoms'] = np.array(values[22])
        if values[23] is not None:
            dct['magmom'] = values[23]
        if values[24] is not None:
            dct['charges'] = np.array(values[24])
        if values[25] != '{}':
            dct['key_value_pairs'] = decode(values[25])
        if len(values) >= 27 and values[26] != 'null':
            dct['data'] = values[26]

        return AtomsRow(dct)


