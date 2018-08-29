import json
from psycopg2 import connect
from psycopg2.extras import execute_values

from ase.db.sqlite import (init_statements, index_statements, VERSION,
                           SQLite3Database)

jsonb_indices = ['CREATE INDEX idxkeys ON systems USING GIN (key_value_pairs);',
                 'CREATE INDEX idxcalc ON systems USING GIN (calculator_parameters);']

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
        if len(args[0]) > 0:
            N = len(args[0][0])
        else:
            return
        if 'INSERT INTO systems' in statement:
            q = 'DEFAULT' + ', ' + ', '.join('?' * N)  # DEFAULT for id
        else:
            q = ', '.join('?' * N)
        statement = statement.replace('({})'.format(q), '%s')
        q = '({})'.format(q.replace('?', '%s'))

        execute_values(self.cur, statement.replace('?', '%s'),
                       argslist=args[0], template=q, page_size=len(args[0]))


class PostgreSQLDatabase(SQLite3Database):
    type = 'postgresql'
    default = 'DEFAULT'

    def _connect(self):
        return Connection(connect(self.filename))

    def _initialize(self, con):
        if self.initialized:
            return

        self._metadata = {}

        cur = con.cursor()
        cur.execute("show search_path;")
        schema = cur.fetchone()[0].split(', ')
        if schema[0] == '"$user"':
            schema = schema[1]
        else:
            schema = schema[0]

        cur.execute("""
        SELECT EXISTS(select * from information_schema.tables where
        table_name='information' and table_schema='{}');
        """.format(schema))

        if not cur.fetchone()[0]:  # information schema doesn't exist.
            # Initialize database:
            sql = ';\n'.join(init_statements)
            sql = schema_update(sql)
            cur.execute(sql)
            if self.create_indices:
                cur.execute(';\n'.join(index_statements))
                cur.execute(';\n'.join(jsonb_indices))
            con.commit()
            self.version = VERSION
        else:
            cur.execute('select * from information;')
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
                  'SERIAL PRIMARY KEY')]:
        sql = sql.replace(a, b)

    arrays_1D = ['numbers', 'initial_magmoms', 'initial_charges', 'masses',
                 'tags', 'momenta', 'stress', 'dipole', 'magmoms', 'charges']

    arrays_2D = ['positions', 'cell', 'forces']

    txt2jsonb = ['calculator_parameters', 'key_value_pairs', 'data']

    for column in arrays_1D:
        if column in ['numbers', 'tags']:
            dtype = 'INTEGER'
        else:
            dtype = 'DOUBLE PRECISION'
        sql = sql.replace('{} BLOB,'.format(column),
                          '{} {}[],'.format(column, dtype))
    for column in arrays_2D:
        sql = sql.replace('{} BLOB,'.format(column),
                          '{} DOUBLE PRECISION[][],'.format(column))
    for column in txt2jsonb:
        sql = sql.replace('{} TEXT,'.format(column),
                          '{} JSONB,'.format(column))

    return sql
