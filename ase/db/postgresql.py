import numpy as np
import json
import psycopg2

from ase.db.sqlite import VERSION
from ase.db.sqlite import SQLite3Database


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
