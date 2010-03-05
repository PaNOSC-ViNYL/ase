from ase.atoms import Atoms
from ase.units import Bohr


def read_dftb(filename='dftb_in.hsd'):
    """Method to read coordinates form DFTB+ input file dftb_in.hsd
    additionally read information about fixed atoms
    and periodic boundary condition
    """
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms
    import sys, string

    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()
    atoms_pos = []
    atom_symbols = []
    type_names = []
    my_pbc = False
    myconstraints=[]
    moved_atoms_found = False
    range_found = False
    moved_atoms = []

    for line in lines:
        if (line.strip().startswith('#')):
            pass
        else:
            if ('TypeNames' in line):
                col=line.split()
                for i in range(3,len(col)-1):
                    type_names.append(col[i].strip("\""))
            elif ('Periodic' in line):
                if ('Yes' in line):
                    my_pbc = True
            else:
                pass

    start_reading_coords=False
    stop_reading_coords=False
    for line in lines:
        if (line.strip().startswith('#')):
            pass
        else:
            if ('TypesAndCoordinates' in line):
                start_reading_coords=True
            if start_reading_coords:
                if ('}' in line):
                    stop_reading_coords=True
            if (start_reading_coords and not(stop_reading_coords)
            and not 'TypesAndCoordinates' in line):
                typeindexstr, x, y, z = line.split()[:4]
                typeindex=int(typeindexstr)
                symbol=type_names[typeindex-1]
                atom_symbols.append(symbol)
                atoms_pos.append([float(x), float(y), float(z)])

            
    if type(filename) == str:
        f.close()

    atoms = Atoms(positions = atoms_pos, symbols = atom_symbols, pbc = my_pbc)


    return atoms


def write_dftb(filename,atoms):
    """Method to write coordinates in DFTB+ format
    """

    import numpy as np
    from ase.constraints import FixAtoms

    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()

    if type(filename) == str:
        f.close()

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    coord = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()


    start_writing_coords = False
    stop_writing_coords = False
    i=0
    for line in lines:
        if ('TypesAndCoordinates' in line):
            start_writing_coords=True
        if (start_writing_coords and not(stop_writing_coords)):
            if ('}' in line):
                stop_writing_coords = True
        if (start_writing_coords and not(stop_writing_coords)and 
            not ('TypesAndCoordinates' in line)):
            atom_type_index = line.split()[0]
            f.write('%6s  %20.14f  %20.14f  %20.14f\n'
                    % (atom_type_index,coord[i][0],coord[i][1],coord[i][2]))
            i=i+1
        else:
            f.write(line)
