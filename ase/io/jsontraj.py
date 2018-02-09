from __future__ import absolute_import
import json

from ase.db.row import AtomsRow, atoms2dict
from ase.io.jsonio import encode, decode


def read_jsontraj(fd, index):  # XXX ignores index
    for line in fd:
        obj = json.loads(line)
        row = AtomsRow(obj)
        atoms = row.toatoms()
        yield atoms


def write_jsontraj(fd, images):
    for image in images:
        obj = atoms2dict(image)
        obj.pop('unique_id')
        json = encode(obj)
        fd.write(json)
        fd.write('\n')


def main():
    import io
    from ase.io import iread, write
    from ase.build import bulk, molecule

    images1 = [molecule('H2O'), molecule('O2'),
               bulk('Si'), bulk('Fe')]

    buf = io.StringIO()
    write(buf, images1, format='jsontraj')
    buf.seek(0)
    print(buf.getvalue())

    images2 = list(iread(buf, format='jsontraj'))

    for img1, img2 in zip(images1, images2):
        assert img1 == img2

if __name__ == '__main__':
    main()
