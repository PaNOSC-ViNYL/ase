import json
from ase.nomad import dict2images
from ase.utils import basestring

def read_nomad_json(fd, index, silent=False):
    # wth, we should not be passing index like this!
    from ase.io.formats import string2index
    if isinstance(index, basestring):
        index = string2index(index)
    d = json.load(fd)
    images = dict2images(d, silent=silent)
    return list(images)[index]
