#!/usr/bin/env python
import sys
lines = open(sys.argv[1]).readlines()
open(sys.argv[1] + '.bak', 'w').write(''.join(lines))
f = open(sys.argv[1], 'w')

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
                 ('GetUnitCell', 'get_cell'),
                 ('GetBoundaryConditions', 'get_pbc'),
                 ('GetCartesianForces', 'get_forces'),
                 ('ListOfAtoms', 'Atoms'),
                 ('periodic', 'pbc')]:
    t = t.replace(old, new)
    
while 1:
    i = t.find('.')
    if i == -1:
        break
    f.write(t[:i + 1])
    t = t[i + 1:]
    if t[0].isupper() and t[1].islower():
        j = t.find('(')
        if j != -1 and t[2: j].isalpha():
            for k in range(j):
                if t[k].isupper() and k > 0:
                    f.write('_')
                f.write(t[k].lower())
            t = t[j:]

f.write(t)
