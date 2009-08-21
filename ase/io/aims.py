
def read_aims(filename):
    """Import FHI-aims geometry type files.

    Reads unitcell, atom positions and constraints from
    a geometry.in file.
    """

    from ase import Atoms, FixAtoms, FixCartesian
    import numpy as np

    atoms = Atoms()
    fd = open(filename, 'r')
    lines = fd.readlines()
    fd.close()
    positions = []
    cell = []
    symbols = []
    fix = []
    fix_cart = []
    xyz = np.array([0, 0, 0])
    i = -1
    for n, line in enumerate(lines):
        inp = line.split()
        if inp == []:
            continue
        if inp[0] == 'atom':
            if xyz.all():
                fix.append(i)
            elif xyz.any():
                print 1
                fix_cart.append(FixCartesian(i, xyz))
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            positions.append(floatvect)
            symbols.append(inp[-1])
            i += 1
            xyz = np.array([0, 0, 0])
        elif inp[0] == 'lattice_vector':
            floatvect = float(inp[1]), float(inp[2]), float(inp[3])
            cell.append(floatvect)
        if inp[0] == 'constraint_relaxation':
            if inp[1] == '.true.':
                fix.append(i)
            elif inp[1] == 'x':
                xyz[0] = 1
            elif inp[1] == 'y':
                xyz[1] = 1
            elif inp[1] == 'z':
                xyz[2] = 1
    if xyz.all():
        fix.append(i)
    elif xyz.any():
        fix_cart.append(FixCartesian(i, xyz))
    atoms = Atoms(symbols, positions)
    if len(cell)==3:
        atoms.set_cell(cell)
    if len(fix):
        atoms.set_constraint([FixAtoms(indices=fix)]+fix_cart)
    else:
        atoms.set_constraint(fix_cart)
    return atoms

def write_aims(filename, atoms, cell=None):
    """Method to write FHI-aims geometry files.

    Writes the atoms positions and constraints (only FixAtoms is
    supported at the moment). If cell=True also the unitcell is
    written out.
    """

    from ase.constraints import FixAtoms, FixCartesian
    import numpy as np

    if cell is None:
        cell = []
        for b in atoms.get_pbc():
            cell.append(b)
    elif type(cell) is bool:
        cell = [bool, bool, bool]
    fd = open(filename, 'w')
    i = 0
    for n, vector in enumerate(atoms.get_cell()):
        if cell[n]:
            i = 1
            fd.write('lattice_vector ')
            for i in range(3):
                fd.write('%16.16f ' % vector[i])
            fd.write('\n')
    if i:
        fd.write('\n')
    fix = []
    fix_cart = {}
    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fix = constr.index
            elif isinstance(constr, FixCartesian):
                fix_cart[str(constr.a)] = -constr.mask+1
    fix = np.array(fix)%len(atoms)
    for i, atom in enumerate(atoms):
        fd.write('atom ')
        for pos in atom.get_position():
            fd.write('%16.16f ' % pos)
        fd.write(atom.symbol)
        fd.write('\n')
        if i in fix:
            fd.write('constraint_relaxation .true.\n')
        else:
            try:
                xyz = fix_cart[str(i)]
                for n in range(3):
                    if xyz[n]:
                        fd.write('constraint_relaxation %s\n' % 'xyz'[n])
            except KeyError:
                continue


def read_energy(filename):
    for line in open(filename, 'r'):
        if line.startswith('  | Total energy corrected'):
            E = float(line.split()[-2])
    return E
