# creates: io.csv
from __future__ import print_function
from ase.io.formats import all_formats, get_ioformat
with open('io.csv', 'w') as fd:
    print('format, description, capabilities', file=fd)
    for format in sorted(all_formats):
        io = get_ioformat(format)
        c = ''
        if io.read:
            c = 'R'
        if io.write:
            c += 'W'
        if not io.single:
            c += '+'
        print('{0}, {1}, {2}'.format(format, all_formats[format][0], c),
              file=fd)
