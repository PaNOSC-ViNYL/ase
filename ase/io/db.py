import numpy as np
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

def read_db(filename, index):
    import cmr
    r = cmr.read(filename)
    if r.is_group():
       hashes = r.get_member_hashes()
       hash = hashes[index]
       d = r.get_xmldata(hash)
    else:
       d = r
    if type(d) == str or len(d.keys("ase")) == 0:
        raise RuntimeError('db-file: there is no ase readable data available in the specified file!')
    positions = d.get('ase_positions')
    numbers = d.get('ase_atomic_numbers')
    cell = d.get('ase_cell')
    pbc = d.get('ase_pbc')
    tags = np.array(d.get('ase_tags'))
    magmoms = np.array(d.get('ase_magnetic_moments'))
    energy = d.get('ase_potential_energy')

    forces = d.get('ase_forces')

    atoms = Atoms(positions=positions,
                  numbers=numbers,
                  cell=cell,
                  pbc=pbc)
    if tags.any():
        atoms.set_tags(tags)

    if magmoms.any():
        atoms.set_initial_magnetic_moments(magmoms)
    else:
        magmoms = None

    atoms.calc = SinglePointCalculator(energy, forces, None, magmoms,
                                           atoms)

    return atoms

