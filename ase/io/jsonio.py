import datetime
import json

import numpy as np


class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, datetime.datetime):
            return {'__datetime__': obj.isoformat()}
        if hasattr(obj, 'todict'):
            return obj.todict()
        return json.JSONEncoder.default(self, obj)


encode = MyEncoder().encode


def object_hook(dct):
    if '__datetime__' in dct:
        return datetime.datetime.strptime(dct['__datetime__'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
    return dct


mydecode = json.JSONDecoder(object_hook=object_hook).decode


def intkey(key):
    try:
        return int(key)
    except ValueError:
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
