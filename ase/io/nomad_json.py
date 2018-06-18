import json
from ase.nomad import read, dict2images
from ase.utils import basestring

def read_nomad_json(fd, index):
    # wth, we should not be passing index like this!
    from ase.io.formats import string2index
    if isinstance(index, basestring):
        index = string2index(index)
    d = json.load(fd)
    images = dict2images(d)
    return list(images)[index]
