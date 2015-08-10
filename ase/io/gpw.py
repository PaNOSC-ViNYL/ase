from ase import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.calculators.singlepoint import SinglePointKPoint
from ase.units import Bohr, Hartree


def read_gpw(filename):
    import gpaw
    r = gpaw.io.open(filename, 'r')
    positions = r.get('CartesianPositions') * Bohr
    numbers = r.get('AtomicNumbers')
    cell = r.get('UnitCell') * Bohr
    pbc = r.get('BoundaryConditions')
    tags = r.get('Tags')
    magmoms = r.get('MagneticMoments')
    energy = r.get('PotentialEnergy') * Hartree

    if r.has_array('CartesianForces'):
        forces = r.get('CartesianForces') * Hartree / Bohr
    else:
        forces = None

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

    atoms.calc = SinglePointDFTCalculator(atoms, energy=energy,
                                          forces=forces, magmoms=magmoms)
    kpts = []
    if r.has_array('IBZKPoints'):
        for w, kpt, eps_n, f_n in zip(r.get('IBZKPointWeights'),
                                      r.get('IBZKPoints'),
                                      r.get('Eigenvalues'),
                                      r.get('OccupationNumbers')):
            kpts.append(SinglePointKPoint(w, kpt[0], kpt[1],
                                          eps_n[0], f_n[0]))
    atoms.calc.kpts = kpts

    yield atoms
