
def read_aims(filename):
    """Import FHI-aims geometry type files.

    Reads unitcell, atom positions and constraints from
    a geometry.in file.
    """

    from ase import Atoms, FixAtoms

    atoms = Atoms()
    fd = open(filename, 'r')
    lines = fd.readlines()
    fd.close()
    positions = []
    cell = []
    symbols = []
    fix = []
    i = -1
    for n, line in enumerate(lines):
        inp = line.split()
        if inp == []:
            continue
        if inp[0] == 'atom':
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            symbols.append(inp[-1])
            i += 1
        elif inp[0] == 'lattice_vector':
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            cell.append(floatvect)
        try:
            if lines[n+1].split()[0] == 'constraint_relaxation':
                fix.append(i)
        except IndexError:
            continue
    atoms = Atoms(symbols, positions)
    if len(cell)==3:
        atoms.set_cell(cell)
    if len(fix):
        atoms.set_constraint(FixAtoms(indices=fix))
    return atoms

def write_aims(filename, atoms, cell=False):
    """Method to write FHI-aims geometry files.

    Writes the atoms positions and constraints (only FixAtoms is
    supported at the moment). If cell=True also the unitcell is
    written out.
    """

    from ase.constraints import FixAtoms
    import numpy as np

    fd = open(filename, 'w')
    if cell:
        for vector in atoms.get_cell():
            fd.write('lattice_vector ')
            for i in range(3):
                fd.write('%16.16f ' % vector[i])
            fd.write('\n')
        fd.write('\n')
    fix = []
    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix = constr.index
    fix = np.array(fix)%len(atoms)
    for i, atom in enumerate(atoms):
        fd.write('atom ')
        for pos in atom.get_position():
            fd.write('%16.16f ' % pos)
        fd.write(atom.symbol)
        fd.write('\n')
        if i in fix:
            fd.write('constraint_relaxation .true.\n')



def read_energy(filename):
    for line in open(filename, 'r'):
        if line.startswith('  | Total energy corrected'):
            E = float(line.split()[-2])
    return E
