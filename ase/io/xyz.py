from ase.atoms import Atoms


def read_xyz(fileobj, index):
    lines = fileobj.readlines()
    natoms = int(lines[0])
    nimages = len(lines) // (natoms + 2)
    for i in range(*index.indices(nimages)):
        symbols = []
        positions = []
        n = i * (natoms + 2) + 2
        for line in lines[n:n + natoms]:
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        yield Atoms(symbols=symbols, positions=positions)


def write_xyz(fileobj, images, comment=''):
    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n%s\n' % (natoms, comment))
        for s, (x, y, z) in zip(symbols, atoms.positions):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
