#from __future__ import absolute_import
import os
import numpy as np

if 1:
    def encode(obj):
        if isinstance(obj, str):
            return '"' + obj + '"'
        if isinstance(obj, (int, float)):
            return repr(obj)
        if isinstance(obj, (bool, np.bool_)):
            return repr(obj).lower()
        if isinstance(obj, dict):
            return '{' + ', '.join(['"' + key + '": ' + encode(value)
                                    for key, value in obj.items()]) + '}'
        return '[' + ','.join([encode(value) for value in obj]) + ']'
    def loads(txt):
        return eval(txt, {'false': False, 'true': True})
else:
    from json import JSONEncoder, loads
    class NDArrayEncoder(JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return JSONEncoder.default(self, obj)
    encode = NDArrayEncoder().encode
        
from ase.parallel import world
from ase.db import KeyCollisionError
from ase.db.core import NoDatabase, dict2atoms


def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((key, numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list):
        try:
            a = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
        else:
            if a.dtype != object:
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
    def _write(self, id, atoms, extra, replace):
        if os.path.isfile(self.filename):
            data = read_json(self.filename)
            if not replace and id in data:
                raise KeyCollisionError
        else:
            data = {}

        dct = self.collect_data(atoms)
        dct['timestamp'] = dct['timestamp'].isoformat(' ')
        dct['extra'] = extra
        data[id] = dct
        write_json(self.filename, data)

    def get_dict(self, id):
        dct = read_json(self.filename)
        if id in [-1, 0]:
            assert len(dct) == 1
            id = dct.keys()[0]
        return dct[id]
