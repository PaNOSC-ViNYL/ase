import numpy as npy

from ase.atoms import Atom, Atoms

def read_dacapo_text(fileobj):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    i = lines.index(' Structure:             A1           A2            A3\n')
    cell = npy.array([[float(w) for w in line.split()[2:5]]
                      for line in lines[i + 1:i + 4]]).transpose()
    i = lines.index(' Structure:  >>         Ionic positions/velocities ' +
                    'in cartesian coordinates       <<\n')
    atoms = []
    for line in lines[i + 4:]:
        words = line.split()
        if len(words) != 9:
            break
        Z, x, y, z = words[2:6]
        atoms.append(Atom(int(Z), [float(x), float(y), float(z)]))
    return Atoms(atoms, cell=cell.tolist())

