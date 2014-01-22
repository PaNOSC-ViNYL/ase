from __future__ import absolute_import, print_function
import sqlite3

import numpy as np

from ase.db.core import Database, ops, now
from ase.db.json import encode, decode


init_statements = """\
create table systems (
    id integer primary key autoincrement,
    unique_id text unique,
    ctime real,
    mtime real,
    user text,
    numbers blob,
    positions blob,
    cell blob,
    pbc integer,
    initial_magmoms blob,
    initial_charges blob,
    masses blob,
    tags blob,
    momenta blob,
    constraints text,
    calculator text,
    calculator_parameters text,
    energy real,
    free_energy real,
    forces blob,
    stress blob,
    dipole blob,
    magmoms blob,
    magmom blob,
    charges blob,
    data text); -- contains keywords and key_value_pairs also
create table species (
    Z integer,
    n integer,
    id text,
    foreign key (id) references systems(id));
create table keywords (
    keyword text,
    id text,
    foreign key (id) references systems(id));
create table text_key_values (
    key text,
    value text,
    id text,
    foreign key (id) references systems(id));
create table number_key_values (
    key text,
    value real,
    id text,
    foreign key (id) references systems (id))
"""

index_statements = """\
create index unique_id_index on systems(unique_id);
create index ctime_index on systems(ctime);
create index user_index on systems(user);
create index calculator_index on systems(calculator);
create index species_index on species(Z);
create index keyword_index on keywords(keyword);
create index text_index on text_key_values(key);
create index number_index on number_key_values(key)
"""

tables = ['systems', 'species', 'keywords',
          'text_key_values', 'number_key_values']


class SQLite3Database(Database):
    initialized = False
    
    def _connect(self):
        return sqlite3.connect(self.filename)

    def _initialize(self, con):
        if self.initialized:
            return
        cur = con.execute(
            'select count(*) from sqlite_master where name="systems"')
        if cur.fetchone()[0] == 0:
            for statement in init_statements.split(';'):
                con.execute(statement)
            if self.create_indices:
                for statement in index_statements.split(';'):
                    con.execute(statement)
            con.commit()
            self.initialized = True
            
    def _write(self, atoms, keywords, key_value_pairs, data):
        Database._write(self, atoms, keywords, key_value_pairs, data)
        
        con = self._connect()
        self._initialize(con)
        cur = con.cursor()
                
        id = None
        
        if isinstance(atoms, dict):
            dct = atoms
            unique_id = dct['unique_id']
            cur.execute('select id from systems where unique_id=?',
                        (unique_id,))
            rows = cur.fetchall()
            if rows:
                id = rows[0][0]
                self._delete(cur, [id])
            dct['mtime'] = now()
        else:
            dct = self.collect_data(atoms)

        if 'constraints' in dct:
            constraints = encode(dct['constraints'])
        else:
            constraints = None
            
        row = (id,
               dct['unique_id'],
               dct['ctime'],
               dct['mtime'],
               dct['user'],
               blob(dct.get('numbers')),
               blob(dct.get('positions')),
               blob(dct.get('cell')),
               int(np.dot(dct.get('pbc'), [1, 2, 4])),
               blob(dct.get('initial_magmoms')),
               blob(dct.get('initial_charges')),
               blob(dct.get('masses')),
               blob(dct.get('tags')),
               blob(dct.get('momenta')),
               constraints)

        if 'calculator' in dct:
            row += (dct['calculator'],
                    encode(dct['calculator_parameters']))
        else:
            row += (None, None)

        magmom = dct.get('magmom')
        if magmom is not None:
            # magmom can be one or three numbers (non-collinear case)
            magmom = np.array(magmom)
        row += (dct.get('energy'),
                dct.get('free_energy'),
                blob(dct.get('forces')),
                blob(dct.get('stress')),
                blob(dct.get('dipole')),
                blob(dct.get('magmoms')),
                blob(magmom),
                blob(dct.get('charges')))

        row += (encode({'data': data,
                        'keywords': keywords,
                        'key_value_pairs': key_value_pairs}),)

        q = ', '.join('?' * len(row))
        cur.execute('insert into systems values (%s)' % q, row)
        
        if id is None:
            cur.execute('select seq from sqlite_sequence where name="systems"')
            id = cur.fetchone()[0]

        count = np.bincount(dct['numbers'])
        unique_numbers = count.nonzero()[0]
        species = [(int(Z), int(count[Z]), id) for Z in unique_numbers]
        cur.executemany('insert into species values (?, ?, ?)', species)

        text_key_values = []
        number_key_values = []
        for key, value in key_value_pairs.items():
            if isinstance(value, (str, unicode)):
                text_key_values.append([key, value, id])
            elif isinstance(value, (float, int)):
                number_key_values.append([key, float(value), id])
            else:
                assert 0, value
 
        if text_key_values:
            cur.executemany('insert into text_key_values values (?, ?, ?)',
                            text_key_values)
        if number_key_values:
            cur.executemany('insert into number_key_values values (?, ?, ?)',
                            number_key_values)
        if keywords:
            cur.executemany('insert into keywords values (?, ?)',
                            [(keyword, id) for keyword in keywords])

        con.commit()
        con.close()
        return id
        
    def _get_dict(self, id):
        con = self._connect()
        c = con.cursor()
        if id is None:
            c.execute('select count(*) from systems')
            assert c.fetchone()[0] == 1
            c.execute('select * from systems')
        else:
            c.execute('select * from systems where id=?', (id,))
        row = c.fetchone()
        return self.row_to_dict(row)

    def row_to_dict(self, row):
        dct = {'id': row[0],
               'unique_id': row[1],
               'ctime': row[2],
               'mtime': row[3],
               'user': row[4],
               'numbers': deblob(row[5], np.int32),
               'positions': deblob(row[6], shape=(-1, 3)),
               'cell': deblob(row[7], shape=(3, 3)),
               'pbc': (row[8] & np.array([1, 2, 4])).astype(bool)}
        if row[9] is not None:
            dct['magmoms'] = deblob(row[9])
        if row[10] is not None:
            dct['charges'] = deblob(row[10])
        if row[11] is not None:
            dct['masses'] = deblob(row[11])
        if row[12] is not None:
            dct['tags'] = deblob(row[12], np.int32)
        if row[13] is not None:
            dct['momenta'] = deblob(row[13], shape=(-1, 3))
        if row[14] is not None:
            dct['constraints'] = decode(row[14])
        if row[15] is not None:
            dct['calculator'] = row[15]
            dct['calculator_parameters'] = decode(row[16])
        if row[17] is not None:
            dct['energy'] = row[17]
        if row[18] is not None:
            dct['free_energy'] = row[18]
        if row[19] is not None:
            dct['forces'] = deblob(row[19], shape=(-1, 3))
        if row[20] is not None:
            dct['stress'] = deblob(row[20])
        if row[21] is not None:
            dct['dipole'] = deblob(row[21])
        if row[22] is not None:
            dct['magmoms'] = deblob(row[22])
        if row[23] is not None:
            dct['magmom'] = deblob(row[23])[0]
        if row[24] is not None:
            dct['charges'] = deblob(row[24])

        extra = decode(row[25])
        for key in ['keywords', 'key_value_pairs', 'data']:
            if extra[key]:
                dct[key] = extra[key]
        return dct

    def _select(self, keywords, cmps, explain=False, verbosity=0, limit=None):
        tables = ['systems']
        where = []
        args = []
        for n, keyword in enumerate(keywords):
            tables.append('keywords as keyword{0}'.format(n))
            where.append(
                'systems.id=keyword{0}.id and keyword{0}.keyword=?'.format(n))
            args.append(keyword)
            
        # Special handling of "H=0" and "H<2" type of selections:
        bad = {}
        for key, op, value in cmps:
            if isinstance(key, int):
                bad[key] = bad.get(key, True) and ops[op](0, value)
                
        cmps2 = []
        nspecies = 0
        ntext = 0
        nnumber = 0
        for key, op, value in cmps:
            if key in ['id', 'energy', 'magmom', 'ctime', 'user',
                       'calculator']:
                where.append('systems.{0}{1}?'.format(key, op))
                args.append(value)
            elif key == 'natoms':
                cmps2.append((key, ops[op], value))
            elif isinstance(key, int):
                if bad[key]:
                    cmps2.append((key, ops[op], value))
                else:
                    tables.append('species as specie{0}'.format(nspecies))
                    where.append(('systems.id=specie{0}.id and ' +
                                  'specie{0}.Z=? and ' +
                                  'specie{0}.n{1}?').format(nspecies, op))
                    args += [key, value]
                    nspecies += 1
            elif isinstance(value, str):
                tables.append('text_key_values as text{0}'.format(ntext))
                where.append(('systems.id=text{0}.id and ' +
                              'text{0}.key=? and ' +
                              'text{0}.value{1}?').format(ntext, op))
                args += [key, value]
                ntext += 1
            else:
                tables.append('number_key_values as number{0}'.format(nnumber))
                where.append(('systems.id=number{0}.id and ' +
                              'number{0}.key=? and ' +
                              'number{0}.value{1}?').format(nnumber, op))
                args += [key, value]
                nnumber += 1
                
        sql = 'select systems.* from\n  ' + ', '.join(tables)
        if where:
            sql += '\n  where\n  ' + ' and\n  '.join(where)
            
        if explain:
            sql = 'explain query plan ' + sql
            
        if limit:
            sql += '\nlimit {0}'.format(limit)
            
        if verbosity == 2:
            print(sql, args, cmps2)
        con = self._connect()
        cur = con.cursor()
        self._initialize(con)
        cur.execute(sql, args)
        if explain:
            for row in cur.fetchall():
                yield {'explain': row}
        else:
            for row in cur.fetchall():
                if cmps2:
                    numbers = deblob(row[5], np.int32)
                    for key, op, value in cmps2:
                        if key == 'natoms':
                            if not op(len(numbers), value):
                                break
                        elif not op((numbers == key).sum(), value):
                            break
                    else:
                        yield self.row_to_dict(row)
                else:
                    yield self.row_to_dict(row)
        
    def delete(self, ids):
        con = self._connect()
        self._delete(con.cursor(), ids)
        con.commit()
        con.close()

    def _delete(self, cur, ids):
        for table in tables[::-1]:
            cur.executemany('delete from {0} where id=?'.format(table),
                            ((id,) for id in ids))


def blob(array):
    """Convert array to blob/buffer object."""

    if array is None:
        return None
    if array.dtype == np.int64:
        array = array.astype(np.int32)
    return buffer(array)


def deblob(buf, dtype=float, shape=None):
    """Convert blob/buffer object to ndarray of correct dtype and shape.

    (without creating an extra view)."""

    if buf is None:
        return None
    array = np.frombuffer(buf, dtype)
    if shape is not None:
        array.shape = shape
    return array
