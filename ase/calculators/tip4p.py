""" TIP4P potential for water.
    http://dx.doi.org/10.1063/1.445869

    Requires an atoms object of OHH,OHH, ... sequence
    Correct TIP4P charges and LJ parameters set automatically.

    Virtual interaction sites implemented in the following scheme:
    Original atoms object has no virtual sites.
    When energy/forces are requested:
        - virtual sites added to temporary xatoms object
        - energy / forces calculated
        - forces redistributed from virtual sites to actual atoms object
    This means you do not get into trouble when propagating your system with MD
    while having to skip / account for massless virtual sites.

    This also means that if using for QM/MM MD with GPAW, the EmbedTIP4P class
    must be used.
    """

import numpy as np
import ase.units as unit
from ase import Atoms
import warnings

# Electrostatic constant and parameters:
k_c = 332.1 * unit.kcal / unit.mol
sigma0 = 3.15365
epsilon0 = 0.6480 * unit.kJ / unit.mol
rOH = 0.9572
thetaHOH = 104.52 / 180 * np.pi


class TIP4P:
    implemented_properties = ['energy', 'forces']

    def __init__(self):
        self.energy = None
        self.forces = None
        self.name = 'TIP4P'
        self.rc = 9.
        nm = 4
        self.nm = nm
        LJ = np.zeros((2, nm))
        LJ[0, 0] += epsilon0
        LJ[1, 0] += sigma0
        self.LJ = LJ

    def calculate(self, atoms):
        self.atoms = atoms
        self.cell = atoms.get_cell()
        self.pbc = atoms.get_pbc()
        self.numbers = atoms.get_atomic_numbers()
        self.positions = atoms.get_positions()
        N = self.nm

        natoms = len(atoms)
        nmol = natoms // self.nm

        self.energy = 0.0
        self.forces = np.zeros((natoms, 3))

        C = self.cell.diagonal()
        cell = C
        pbc = atoms.pbc
        assert (atoms.cell == np.diag(cell)).all(), 'not orthorhombic'
        # assert ((cell >= 2 * self.rc) | ~pbc).all(), 'cutoff too large'
        if not ((cell >= 2 * self.rc) | ~pbc).all():
            warnings.warn('Cutoff too large.')

        # Get dx,dy,dz from first atom of each mol to same atom of all other
        # and find min. distance. Everything moves according to this analysis.
        for a in range(nmol-1):
            D = self.positions[(a+1)*N::N] - self.positions[a*N]
            n = np.rint(D / C) * self.pbc
            q_v = self.atoms[(a+1)*N:].get_initial_charges()

            # Min. img. position list as seen for molecule !a!
            position_list = np.zeros(((nmol-1-a)*N, 3))

            for j in range(N):
                position_list[j::N] += self.positions[(a+1)*N+j::N] - n*C

            self.energy_and_forces(a, position_list, q_v, nmol)

    def energy_and_forces(self, a, position_list, q_v, nmol):
        """ The combination rules for the LJ terms follow Waldman-Hagler:
            J. Comp. Chem. 14, 1077 (1993)
        """
        N = self.nm
        LJ = self.LJ
        for i in range(N):
            D = position_list - self.positions[a*N+i]
            d = (D**2).sum(axis=1)

            # Create arrays to hold on to epsilon and sigma
            epsilon = np.zeros(N)
            sigma = np.zeros(N)

            for j in range(N):
                if self.LJ[1, i] * self.LJ[1, j] == 0:
                    epsilon[j] = 0
                else:
                    epsilon[j] = 2 * LJ[1, i]**3 * LJ[1, j]**3 \
                              * np.sqrt(LJ[0, i] * LJ[0, j]) \
                              / (LJ[1, i]**6 + LJ[1, j]**6)

                sigma[j] = ((LJ[1, i]**6 + LJ[1, j]**6) / 2)**(1./6)

            # Create list of same length as position_list
            epsilon_list = np.tile(epsilon, (nmol-1-a))
            sigma_list = np.tile(sigma, (nmol-1-a))

            self.energy += (k_c * q_v[i] * q_v / d**0.5 + 4 * epsilon_list *
                            (sigma_list**12 / d**6 -
                             sigma_list**6 / d**3)).sum()

            F = ((k_c * q_v[i] * q_v / d**0.5 + 4 * epsilon_list *
                 (12 * sigma_list**12 / d**6 - 6 * sigma_list**6 / d**3))
                 / d)[:, np.newaxis] * D

            self.forces[a*N+i] -= F.sum(axis=0)
            self.forces[(a+1)*N:] += F

    def add_virtual_sites(self, atoms):
        # Order: OHHM,OHHM,...
        # DOI: 10.1002/(SICI)1096-987X(199906)20:8
        nm = self.nm
        b = 0.15
        a = 0.5
        pos = atoms.get_positions()
        xatomspos = np.zeros((int((4./3)*len(atoms)), 3))
        molct = 0
        for w in range(0, len(atoms), 3):
            r_i = pos[w]    # O pos
            r_j = pos[w+1]  # H1 pos
            r_k = pos[w+2]  # H2 pos
            r_ij = r_j - r_i
            r_jk = r_k - r_j
            r_d = r_i + b*(r_ij + a*r_jk)/np.linalg.norm(r_ij + a*r_jk)

            xatomspos[w+0+molct] = r_i
            xatomspos[w+1+molct] = r_j
            xatomspos[w+2+molct] = r_k
            xatomspos[w+3+molct] = r_d

            molct += 1

        nummols = len(atoms)/3
        assert(nummols - int(nummols) == 0)
        nummols = int(nummols)

        # ase max recursion depth...
        mol = Atoms('OHHH')
        xatoms = mol.copy()
        for i in range(nummols-1):
            xatoms += mol
        xatoms.set_pbc(atoms.get_pbc())
        xatoms.set_cell(atoms.get_cell())
        xatoms.set_positions(xatomspos)
        charges = np.zeros(len(xatoms))
        charges[0::nm] = 0.00  # O
        charges[1::nm] = 0.52  # H1
        charges[2::nm] = 0.52  # H2
        charges[3::nm] = -1.04  # X1
        xatoms.set_initial_charges(charges)
        return xatoms

    def redistribute_forces(self, atoms):
        nm = self.nm
        f = self.forces
        b = 0.15
        a = 0.5
        pos = atoms.get_positions()
        for w in range(0, len(atoms), nm):
            r_i = pos[w]    # O pos
            r_j = pos[w+1]  # H1 pos
            r_k = pos[w+2]  # H2 pos
            r_ij = r_j - r_i
            r_jk = r_k - r_j
            r_d = r_i + b*(r_ij + a*r_jk)/np.linalg.norm(r_ij + a*r_jk)
            r_id = r_d - r_i
            gamma = b/np.linalg.norm(r_ij + a * r_jk)

            Fd = f[w+3]  # force on M
            F1 = (np.dot(r_id, Fd) / np.dot(r_id, r_id)) * r_id
            Fi = Fd - gamma*(Fd - F1)  # Force from M on O
            Fj = (1-a)*gamma*(Fd - F1)  # Force from M on H1
            Fk = a*gamma*(Fd-F1)       # Force from M on H2

            f[w] += Fi
            f[w+1] += Fj
            f[w+2] += Fk

        # remove virtual sites from force array
        f = np.delete(f, list(range(3, f.shape[0], 4)), axis=0)
        self.forces = f

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise NotImplementedError

    def update(self, atoms):
        self.realpositions = atoms.get_positions()
        xatoms = self.add_virtual_sites(atoms)
        self.calculate(xatoms)
        self.redistribute_forces(xatoms)
