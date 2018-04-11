import numpy as np
import json
import numbers
import os
import psycopg2

from ase.data import atomic_numbers
from ase.db.sqlite import VERSION
from ase.db.sqlite import SQLite3Database, float_if_not_none
from ase.db.core import Database, now
from ase.db.row import AtomsRow
from ase.io.jsonio import encode, decode
from ase.utils import basestring

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
    calculator_parameters JSONB,
    energy DOUBLE PRECISION,  -- calculated properties
    free_energy DOUBLE PRECISION,
    forces DOUBLE PRECISION[][],
    stress DOUBLE PRECISION[],
    dipole DOUBLE PRECISION[],
    magmoms DOUBLE PRECISION[],
    magmom DOUBLE PRECISION,
    charges DOUBLE PRECISION[],
    key_value_pairs JSONB,  -- key-value pairs and data as json
    data JSONB,
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
        self.type = 'postgresql'
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
            dct['calculator_parameters'] = values[16]
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
            dct['key_value_pairs'] = values[25]
        if len(values) >= 27 and values[26] != 'null':
            dct['data'] = values[26]

        return AtomsRow(dct)
