from __future__ import absolute_import, print_function
import os
import copy
import datetime
from json import JSONEncoder, JSONDecoder

import numpy as np

from ase.parallel import world
from ase.db.core import NoDatabase, ops, parallel, lock


class MyEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, datetime.datetime):
            return {'__datetime__': obj.isoformat()}
        if hasattr(obj, 'todict'):
            return obj.todict()
        return JSONEncoder.default(self, obj)


encode = MyEncoder().encode


def object_hook(dct):
    if '__datetime__' in dct:
        return datetime.datetime.strptime(dct['__datetime__'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
    return dct


mydecode = JSONDecoder(object_hook=object_hook).decode


def intkey(key):
    if key[0].isdigit():
        return int(key)
    return key
    
    
def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((intkey(key), numpyfy(value))
                    for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
        else:
            if a.dtype in [bool, int, float]:
                obj = a
    return obj


def decode(txt):
    return numpyfy(mydecode(txt))

    
def read_json(name):
    if isinstance(name, str):
        fd = open(name, 'r')
    else:
        fd = name
    dct = decode(fd.read())
    fd.close()
    return dct


class JSONDatabase(NoDatabase):
    def _write(self, atoms, keywords, key_value_pairs, data):
        bigdct = {}
        ids = []
        nextid = 1

        if isinstance(self.filename, str) and os.path.isfile(self.filename):
            try:
                bigdct, ids, nextid = self._read_json()
            except SyntaxError:
                pass

        if isinstance(atoms, dict):
            dct = copy.deepcopy(atoms)
            unique_id = dct['unique_id']
            for id in ids:
                if bigdct[id]['unique_id'] == unique_id:
                    break
            else:
                id = None
        else:
            dct = self.collect_data(atoms)
            id = None

        if keywords:
            dct['keywords'] = keywords
        if key_value_pairs:
            dct['key_value_pairs'] = key_value_pairs
        if data:
            dct['data'] = data
        
        if id is None:
            id = nextid
            ids.append(id)
            nextid += 1
            
        bigdct[id] = dct
        self._write_json(bigdct, ids, nextid)

    def _read_json(self):
        bigdct = read_json(self.filename)
        return bigdct, list(bigdct['ids']), bigdct['nextid']
        
    def _write_json(self, bigdct, ids, nextid):
        if world.rank > 0:
            return
            
        if isinstance(self.filename, str):
            fd = open(self.filename, 'w')
        else:
            fd = self.filename
        print('{', end='', file=fd)
        for id in ids:
            print('"{0}":\n{1},'.format(id, encode(bigdct[id])), file=fd)
        print('"ids": {0},'.format(ids), file=fd)
        print('"nextid": {0}}}'.format(nextid), file=fd)
        fd.close()

    @lock
    @parallel
    def delete(self, ids):
        bigdct, myids, nextid = self._read_json()
        for id in ids:
            del bigdct[id]
            myids.remove(id)
        self._write_json(bigdct, myids, nextid)

    def _get_dict(self, id):
        bigdct, ids, nextid = self._read_json()
        if id is None:
            assert len(ids) == 1
            id = ids[0]
        dct = bigdct[id]
        dct['id'] = id
        return dct

    def _select(self, keywords, cmps, explain, verbosity, limit):
        if explain:
            return
        bigdct, ids, nextid = self._read_json()
        cmps = [(key, ops[op], val) for key, op, val in cmps]
        n = 0
        for id in ids:
            if n == limit:
                return
            dct = bigdct[id]
            for keyword in keywords:
                if 'keywords' not in dct or keyword not in dct['keywords']:
                    break
            else:
                for key, op, val in cmps:
                    value = get_value(id, dct, key)
                    if not op(value, val):
                        break
                else:
                    dct['id'] = id
                    n += 1
                    yield dct

    @lock
    @parallel
    def update(self, ids, add_keywords=[], **add_key_value_pairs):
        bigdct, myids, nextid = self._read_json()
        m = 0
        n = 0
        for id in ids:
            dct = bigdct[id]
            if add_keywords:
                keywords = dct.setdefault('keywords', [])
                for keyword in add_keywords:
                    if keyword not in keywords:
                        keywords.append(keyword)
                        m += 1
            if add_key_value_pairs:
                key_value_pairs = dct.setdefault('key_value_pairs', {})
                n -= len(key_value_pairs)
                key_value_pairs.update(add_key_value_pairs)
                n += len(key_value_pairs)
        self._write_json(bigdct, myids, nextid)
        return m, n


def get_value(id, dct, key):
    pairs = dct.get('key_value_pairs')
    if pairs is None:
        value = None
    else:
        value = pairs.get(key)
    if value is not None:
        return value
    if key in ['energy', 'magmom', 'timestamp', 'username']:
        return dct.get(key)
    if key == 'calculator':
        return dct.get('calculator_name')
    if isinstance(key, int):
        return np.equal(dct['numbers'], key).sum()
    if key == 'natoms':
        return len(dct['numbers'])
    if key == 'id':
        return id
