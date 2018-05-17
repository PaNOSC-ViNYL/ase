import json
from ase.nomad import read


def read_nomad_json(fd, index):
    # wth, we should not be passing index like this!
    from ase.io.formats import string2index
    i = string2index(index)
    print(i)
    d = json.load(fd)
    images = dict2images(d)
    return images[i]
