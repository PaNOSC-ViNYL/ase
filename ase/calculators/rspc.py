""" Rigid Simple Point Charge / LJ Potential.
    You must use ase.constraints.FixBondLengths to fix all internal inter-
    atomic distances in every molecule, according to the geometry specified
    in the chosen force field.

    Remember to set_initial_charges() before running a calculation.

    Based on tip3p.py and qmmm.py, and only intended for testing and/or for
    QM/MM MD Calculations, where the QM performance is the major bottleneck."""

import numpy as np
import ase.units as units
from ase.calculators.calculator import Calculator
import warnings
import itertools as it

# Electrostatic constant:
k_c = units.Hartree * units.Bohr

def wrap(D, cell, pbc):
    """Wrap distances to nearest neighbor (minimum image convention)."""
    shift = np.zeros_like(D)
    for i, periodic in enumerate(pbc):
        if periodic:
            d = D[:, i]
            L = cell[i]
            shift[:, i] = (d + L / 2) % L - L / 2 - d
    return shift


class RSPC(Calculator):
    implemented_properties = ['energy', 'forces']
    nolabel = True
    pcpot = None

    def __init__(self, atoms=None, parameters=None, molecule_size=3,
                 rc=7.0, width=1.0):
        """
        rc: float
            Cutoff radius for Coulomb part. Done molecule-wise using the
            center of mass for each molecule.
            WIP: Implemented, but total generality needs further testing.
        width: float
            Width for cutoff function for Coulomb part.
        molecule_size: int
            number of atoms per molecule.
            WIP: (Counter)ions... Currently not possible.
        parameters: dict
            Mapping from pair of atoms to tuple containing epsilon and sigma.
            Example: parameters = {('O', 'O'): (eps, sigma)}
            WIP: implement this parsing of this:
                parameters = {'O': (eps, sigma), 'H', (eps, sigma),
                              'combine':'rule'}
            WIP: Avoid keyerror if  pairs are missing (i.e. find missing pairs,
                 set to eps, sigma to 0)
        """
        self.rc = rc
        self.width = width
        Calculator.__init__(self)
        self.atoms = atoms
        self.molecule_size = molecule_size
        fparams = {}
        for (Z1, Z2), (epsilon, sigma) in parameters.items():
            fparams[(Z1, Z2)] = epsilon, sigma
            fparams[(Z2, Z1)] = epsilon, sigma
        self.fparams = fparams

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        pbc = atoms.pbc
        cell = atoms.cell.diagonal()
        assert (atoms.cell == np.diag(cell)).all(), 'not orthorhombic'
        # assert ((cell >= 2 * self.rc) | ~pbc).all(), 'cutoff too large'
        if not ((cell >= 2 * self.rc) | ~pbc).all():
            warnings.warn('Cutoff too large.')
        if not atoms.constraints:  # Only warn so single point calcs will run
            warnings.warn('No constraints set.')

        apm = self.molecule_size
        mc = self.get_molcoms()
        nmol = len(mc)
        R = self.atoms.positions.reshape((-1, apm, 3))

        # Init LJ and Coloumb quantities
        lj_forces = np.zeros_like(self.atoms.positions)
        c_forces = np.zeros_like(self.atoms.positions)
        lj_energy = 0
        c_energy = 0
        cs = atoms.get_chemical_symbols()
        cs = np.array(cs, dtype='S2').reshape((-1, nmol, apm))[0]
        q_all = atoms.get_initial_charges()
        for m in range(nmol - 1):
            # get vectors between COM of mol m and m+1: the rest
            Dmm = mc[m + 1:] - mc[m]
            shift = wrap(Dmm, cell, pbc)
            Rshift = R[m + 1:] + shift[:, None, :]
            # get quantities from all atoms of mol m and all other atoms
            types = cs[m+1:]
            q = atoms[(m+1)*apm:].get_initial_charges()
            cm_cut, cm_dcutdd = self.cutoff(np.linalg.norm(Dmm + shift, axis=1))
            all_cut = np.repeat(cm_cut, apm)
            all_dcutdd = np.repeat(cm_dcutdd, apm)
            for i, (a, t) in enumerate(zip(R[m], cs[m])):
                D = Rshift.reshape(-1, apm) - a
                T = types.reshape(-1, (nmol-1-m)*apm)[0]
                d2 = (D**2).sum(1)
                d = d2**0.5
                """ WIP: From here: Turn into function that can be exchanged
                         with other nonbonded interaction models.
                """
                epssig = [self.fparams[t, tt] for tt in T]
                eps = np.array([l[0] for l in epssig])
                sig = np.array([l[-1] for l in epssig])
                c6 = (sig**2 / d2)**3
                c12 = c6**2
                e = 4 * eps * (c12 - c6)
                lj_energy += np.dot(e, all_cut)
                # this one works without cutoff:
                # f = (24 * eps * ((2 * c12 - c6) / d2))[:, np.newaxis] * D
                f = (24 * eps * (2 * c12 - c6) / d2 * all_cut -
                     e * all_dcutdd / d)[:, np.newaxis] * D
                lj_forces[m*apm+i] -= f.sum(0)  # add forces to this atom
                lj_forces[(m+1)*apm:] += f  # add forces to all other atoms
                # Coulomb:
                # c_energy += (k_c * q_all[m*apm+i] * q / d).sum()  # no cut
                e = k_c * q_all[m*apm+i] * q / d
                c_energy += np.dot(all_cut, e)
                # For some reasons this does not work:
                # f = - ((e / d2) * all_cut + e * all_dcutdd
                #       )[:, np.newaxis] * D / d[:, np.newaxis]
                # following TIP3P instead (generality problems?):
                F = (e.reshape(nmol-1-m, apm) / d2.reshape(nmol-1-m, apm) *
                     cm_cut[:, np.newaxis])[:, :, np.newaxis] * \
                    D.reshape(nmol-1-m, apm, 3)
                FOO = -(e.reshape(nmol-1-m, apm).sum(1) * cm_dcutdd /
                        d.reshape(nmol-1-m, apm)[:, 0])[:, np.newaxis] \
                    * D.reshape(nmol-1-m, apm, 3)[:, 0]
                c_forces[(m + 1) * apm::apm] += FOO
                c_forces[m * apm] -= FOO.sum(0)
                c_forces[(m + 1) * apm:] += F.reshape((-1, 3))
                c_forces[m * apm + i] -= F.sum(axis=0).sum(axis=0)

        forces = lj_forces + c_forces
        energy = lj_energy + c_energy
        self.results['energy'] = energy
        self.results['forces'] = forces

    def cutoff(self, d):
        x1 = d > self.rc - self.width
        x2 = d < self.rc
        x12 = np.logical_and(x1, x2)
        y = (d[x12] - self.rc + self.width) / self.width
        cut = np.zeros(len(d))  # cutoff function
        cut[x2] = 1.0
        cut[x12] -= y**2 * (3.0 - 2.0 * y)
        dtdd = np.zeros(len(d))
        dtdd[x12] -= 6.0 / self.width * y * (1.0 - y)
        return cut, dtdd

    def get_molcoms(self):
        apm = self.molecule_size
        nmols = len(self.atoms)/apm
        molcoms = np.zeros((nmols, 3))
        for mol in range(nmols):
            molcoms[mol] = self.atoms[mol*apm:(mol+1)*apm].get_center_of_mass()
            # molcoms[mol] = self.atoms[mol*apm].position
        return molcoms

    def check_state(self, atoms, tol=1e-15):
        system_changes = Calculator.check_state(self, atoms, tol)
        return system_changes


class CombineLJ():
    """ Different Combining Rules for combining LJ parameters
        WIP: Add more combining rules """
    def lorenz_berthelot(self, p):
        combined = {}
        for comb in it.product(p.keys(), repeat=2):
            combined[comb] = ((p[comb[0]][0] * p[comb[1]][0])**0.5,
                              (p[comb[0]][1] + p[comb[1]][1])/2)
        return combined

