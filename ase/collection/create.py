import os

import ase.db
from ase.io import read


def create_dcdft_database():
    os.environ['USER'] = 'ase'
    c = ase.db.connect('dcdft.json')
    with open('WIEN2k.txt') as fd:
        lines = fd.readlines()
    for line in lines[2:73]:
        words = line.split()
        symbol = words.pop(0)
        vol, B, Bp = (float(x) for x in words)
        filename = 'cif/' + symbol + '.cif'
        atoms = read(filename)
        M = {'Fe': 2.3,
             'Co': 1.2,
             'Ni': 0.6,
             'Cr': 1.5,
             'O': 1.5,
             'Mn': 2.0}.get(symbol)
        if M is not None:
            magmoms = [M] * len(atoms)
            if symbol in ['Cr', 'O', 'Mn']:
                magmoms[len(atoms) // 2:] = [-M] * (len(atoms) // 2)
            atoms.set_initial_magnetic_moments(magmoms)
        # c.write(atoms, name=symbol, w2k_B=B, w2k_Bp=Bp, w2k_volume=vol)
        print(symbol, vol - atoms.get_volume() / len(atoms))
        filename = 'pcif/' + symbol + '.cif'
        p = read(filename, primitive_cell=True)
        print(len(atoms), len(p),
              atoms.get_volume() / len(atoms) - p.get_volume() / len(p))

create_dcdft_database()
