#!/usr/bin/env python
import sys

def convert(filename):
    lines = open(filename).readlines()
    t1 = ''.join(lines)

    first = True
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('from ASE'):
            if first:
                lines[i] = 'from ase import *\n'
                first = False
            else:
                lines[i] = ''

    t = ''.join(lines)

    for old, new in [('GetCartesianPositions', 'get_positions'),
                     ('SetCartesianPositions', 'set_positions'),
                     ('SetUnitCell', 'set_cell'),
                     ('GetUnitCell', 'get_cell'),
                     ('GetBoundaryConditions', 'get_pbc'),
                     ('GetCartesianForces', 'get_forces'),
                     ('GetCartesianVelocities', 'get_velocities'),
                     ('SetCartesianVelocities', 'set_velocities'),
                     ('GetCartesianMomenta', 'get_momenta'),
                     ('SetCartesianMomenta', 'set_momenta'),
                     ('ListOfAtoms', 'Atoms'),
                     ('periodic', 'pbc'),
                     ('.Converge(', '.run('),
                     ('Numeric', 'numpy'):
                     ('numpyal', 'Numerical')]:
        t = t.replace(old, new)

    t2 = ''
    while 1:
        i = t.find('.')
        if i == -1:
            i = t.find('def ')
            if i == -1:
                break
        t2 += t[:i + 4]
        t = t[i + 4:]
        if t[0].isupper() and t[1].islower():
            j = t.find('(')
            if j != -1 and t[2: j].isalpha():
                for k in range(j):
                    if t[k].isupper() and k > 0:
                        t2 += '_'
                    t2 += t[k].lower()
                t = t[j:]

    t2 += t

    if t2 != t1:
        print filename, len(t1) - len(t2)
        open(filename + '.bak', 'w').write(t1)
        open(filename, 'w').write(t2)

for filename in sys.argv[1:]:
    convert(filename)
