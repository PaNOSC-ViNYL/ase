import os
import json

import numpy as np

from ase.parallel import world
from ase.db import KeyCollisionError
from ase.db.core import NoDatabase, dict2atoms


class NDArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


encode = NDArrayEncoder().encode


def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((key, numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list):
        try:
            obj = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
    return obj


def write_json(name, results):
    if world.rank == 0:
        fd = open(name, 'w')
        fd.write(encode(results))
        fd.close()


def read_json(name):
    fd = open(name, 'r')
    results = json.loads(fd.read())
    fd.close()
    world.barrier()
    return numpyfy(results)


class JSONDatabase(NoDatabase):
    def __init__(self, filename, use_lock_file=False):
        NoDatabase.__init__(self, use_lock_file)
        self.filename = filename
        
    def _write(self, name, atoms, extra, replace):
        if name is None:
            name = self.create_random_key(atoms)

        if os.path.isfile(self.filename):
            dct = read_json(self.filename)
            if not replace and name in dct:
                raise KeyCollisionError
        else:
            dct = {}

        data = self.collect_data(atoms)
        data['extra'] = extra
        dct[name] = data
        write_json(self.filename, dct)

    def _get(self, names, attach_calculator=False):
        dct = read_json(self.filename)
        if names == [-1] or names == [0]:
            assert len(dct) == 1
            names = dct.keys()
        elif names == [slice(None, None, None)]:
            names = dct.keys()
        return [(dict2atoms(dct[name], attach_calculator), dct[name]['extra'])
                for name in names]
