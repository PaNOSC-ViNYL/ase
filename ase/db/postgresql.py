import numpy as np
import json
import psycopg2

from ase.db.sqlite import init_statements, index_statements, VERSION
from ase.db.sqlite import VERSION
from ase.db.sqlite import SQLite3Database


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
            sql = schema_update(sql)
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


def schema_update(sql):
    for a, b in [('REAL', 'DOUBLE PRECISION'),
                 ('INTEGER PRIMARY KEY AUTOINCREMENT',
                  'SERIAL PRIMARY KEY'),
                 ('BLOB', 'DOUBLE PRECISION')]:
        sql = sql.replace(a, b)

    arrays_1D = ['numbers', 'initial_magmoms', 'initial_charges', 'masses',
                 'tags', 'momenta', 'stress', 'dipole', 'magmoms', 'charges']
    arrays_2D = ['positions', 'cell', 'forces']

    txt2jsonb = ['calculator_parameters', 'key_value_pairs', 'data']

    for column in arrays_1D:
        sql = sql.replace('{} DOUBLE PRECISION,'.format(column),
                          '{} DOUBLE PRECISION[],'.format(column))
    for column in arrays_2D:
        sql = sql.replace('{} DOUBLE PRECISION,'.format(column),
                          '{} DOUBLE PRECISION[][],'.format(column))
    for column in txt2jsonb:
        sql = sql.replace('{} TEXT,'.format(column),
                          '{} JSONB,'.format(column))

    return sql
