#from __future__ import absolute_import  # PY24
import os
import copy

import numpy as np

if 1:
    def encode(obj):
        if isinstance(obj, (str, unicode)):
            return '"' + obj + '"'
        if isinstance(obj, (bool, np.bool_)):
            return repr(obj).lower()
        if isinstance(obj, (int, float)):
            return repr(obj)
        if isinstance(obj, dict):
            return '{' + ', '.join(['"' + key + '": ' + encode(value)
                                    for key, value in obj.items()]) + '}'
        if isinstance(obj, complex):  # no complex in json
            return 'null'
        if obj is None:
            return 'null'
        if isinstance(obj, (list, tuple, np.ndarray)):
            return '[' + ','.join([encode(value) for value in obj]) + ']'
        return encode(obj.todict())

    def loads(txt):
        return eval(txt, {'false': False, 'true': True, 'null': None})
else:
    from json import JSONEncoder, loads
    class NDArrayEncoder(JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return JSONEncoder.default(self, obj)
    encode = NDArrayEncoder().encode
        
from ase.parallel import world
from ase.db import IdCollisionError
from ase.db.core import NoDatabase, ops, parallel, lock


def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((key, numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
        else:
            if a.dtype in [bool, int, float]:
                obj = a
    return obj


def write_json(name, results):
    if world.rank == 0:
        fd = open(name, 'w')
        fd.write(encode(results))
        fd.close()


def read_json(name):
    fd = open(name, 'r')
    results = loads(fd.read())
    fd.close()
    world.barrier()
    return numpyfy(results)


class JSONDatabase(NoDatabase):
    def _write(self, id, atoms, keywords, key_value_pairs, data, replace):
        bigdct = {}
        if os.path.isfile(self.filename):
            try:
                bigdct = read_json(self.filename)
            except SyntaxError:
                pass
            else:
                if not replace and id in bigdct:
                    raise IdCollisionError

        if isinstance(atoms, dict):
            dct = copy.deepcopy(atoms)
            unique_id = dct['unique_id']
            for d in bigdct.values():
                if d['unique_id'] == unique_id:
                    id = d['id']
                    break
        else:
            dct = self.collect_data(atoms)

        if id is None:
            nrows = len(bigdct)
            while id is None:
                id = self.create_random_id(nrows)
                if id in bigdct:
                    id = None

        dct['id'] = id

        if keywords:
            dct['keywords'] = keywords
        if key_value_pairs:
            dct['key_value_pairs'] = key_value_pairs
        if data:
            dct['data'] = data

        bigdct[id] = dct
        write_json(self.filename, bigdct)

    @lock
    @parallel
    def delete(self, ids):
        bigdct = read_json(self.filename)
        for id in ids:
            del bigdct[id]
        write_json(self.filename, bigdct)

    def _get_dict(self, id):
        bigdct = read_json(self.filename)
        if id in [-1, 0]:
            assert len(bigdct) == 1
            id = bigdct.keys()[0]
        return bigdct[id]

    def _select(self, keywords, cmps, limit, offset,
                explain=False, verbosity=1):
        if explain:
            return
        bigdct = read_json(self.filename)
        cmps = [(key, ops[op], val) for key, op, val in cmps]
        offset = offset or 0
        ids = bigdct.keys()[offset:offset + limit]
        for id in ids:
            dct = bigdct[id]
            for keyword in keywords:
                if 'keywords' not in dct or keyword not in dct['keywords']:
                    break
            else:
                for key, op, val in cmps:
                    value = get_value(dct, key)
                    if value is None or not op(value, val):
                        break
                else:
                    yield dct

    @lock
    @parallel
    def update(self, ids, set_keywords):
        bigdct = read_json(self.filename)
        n = 0
        for id in ids:
            keywords = bigdct[id].get('keywords', [])
            for keyword in set_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    n += 1
        if n > 0:
            write_json(self.filename, bigdct)
        return n


def get_value(dct, key):
    pairs = dct.get('key_value_pairs')
    if pairs is None:
        value = None
    else:
        value = pairs.get(key)
    if value is not None:
        return value
    if key in ['energy', 'magmom']:
        return dct.get('results', {}).get(key)
    if key in ['id', 'timestamp', 'username', 'calculator']:
        return dct.get(key)
    if isinstance(key, int):
        return (dct['numbers'] == key).sum()
    if key == 'natoms':
        return len(dct['numbers'])
