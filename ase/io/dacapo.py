import numpy as npy

from ase.calculators import SinglePointCalculator
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



def read_dacapo(filename):
    from ase.io.pupynere import NetCDFFile

    nc = NetCDFFile(filename)
    dims = nc.dimensions
    vars = nc.variables

    cell = vars['UnitCell'][-1]
    atoms = Atoms(positions=npy.dot(vars['DynamicAtomPositions'][-1], cell),
                  symbols=[(a + b).strip() 
                           for a, b in vars['DynamicAtomSpecies'][:]],
                  cell=cell,
                  magmoms=vars['InitialAtomicMagneticMoment'][:],
                  tags=vars['AtomTags'][:],
                  pbc=True)

    calc = SinglePointCalculator(vars['TotalEnergy'].getValue(),
                                 vars['DynamicAtomForces'][-1], None, atoms)
    atoms.set_calculator(calc)
        
    return atoms
