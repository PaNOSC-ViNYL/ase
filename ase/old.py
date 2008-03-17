import numpy as npy
import Numeric as num
from ase.data import chemical_symbols

def npy2num(a, typecode=num.Float):
    return num.array(a, typecode)
if num.__version__ <= '23.8':
    #def npy2num(a, typecode=num.Float):
    #    return num.array(a.tolist(), typecode)
    def npy2num(a, typecode=num.Float):
        b = num.fromstring(a.tostring(), typecode)
        b.shape = a.shape
        return b


class OldASEListOfAtomsWrapper:
    def __init__(self, atoms):
        self.atoms = atoms
        self.constraints = []

    def get_positions(self):
        return npy.array(self.atoms.GetCartesianPositions())

    def get_potential_energy(self):
        return self.atoms.GetPotentialEnergy()

    def get_forces(self):
        return npy.array(self.atoms.GetCartesianForces())

    def get_stress(self):
        return npy.array(self.atoms.GetStress())

    def get_atomic_numbers(self):
        return npy.array(self.atoms.GetAtomicNumbers())

    def get_tags(self):
        return npy.array(self.atoms.GetTags())
    
    def get_momenta(self):
        return npy.array(self.atoms.GetCartesianMomenta())
    
    def get_masses(self):
        return npy.array(self.atoms.GetMasses())
    
    def get_magnetic_moments(self):
        return npy.array(self.atoms.GetMagneticMoments())
    
    def get_charges(self):
        return None
    
    def get_cell(self):
        return npy.array(self.atoms.GetUnitCell())

    def get_pbc(self):
        return npy.array(self.atoms.GetBoundaryConditions(), bool)

    def __len__(self):
        return len(self.atoms)


class OldASECalculatorWrapper:
    def __init__(self, calc, atoms=None):
        self.calc = calc
        self.atoms = calc.GetListOfAtoms()

        if self.atoms is None:
            from ASE import Atom, ListOfAtoms
            
            numbers = atoms.get_atomic_numbers()
            positions = atoms.get_positions()
            self.atoms = ListOfAtoms(
                [Atom(Z=numbers[a], position=positions[a])
                 for a in range(len(atoms))],
                cell=npy2num(atoms.get_cell()),
                periodic=tuple(atoms.get_pbc()))
            self.atoms.SetCalculator(calc)

    def get_atoms(self):
        return OldASEListOfAtomsWrapper(self.atoms)
            
    def get_potential_energy(self, atoms):
        # XXXX what about the cell?
        self.atoms.SetCartesianPositions(npy2num(atoms.get_positions()))
        return self.calc.GetPotentialEnergy()

    def get_forces(self, atoms):
        self.atoms.SetCartesianPositions(npy2num(atoms.get_positions()))
        return npy.array(self.calc.GetCartesianForces())

    def get_stress(self, atoms):
        # XXXX
        return npy.array(self.calc.GetStress())

    def get_number_of_bands(self):
        return self.calc.GetNumberOfBands()

    def get_kpoint_weights(self):
        return npy.array(self.calc.GetIBZKPointWeights())

    def get_number_of_spins(self):
        return 1 + int(self.calc.GetSpinPolarized())

    def get_eigenvalues(self, kpt=0, spin=0):
        return npy.array(self.calc.GetEigenvalues(kpt, spin))

    def get_fermi_level(self):
        return self.calc.GetFermiLevel()

    def get_number_of_grid_points(self):
        return self.get_pseudo_wave_function(0, 0, 0).shape

    def get_pseudo_wave_function(self, n=0, k=0, s=0):
        return npy.array(self.calc.GetWaveFunctionArray(n, k, s))

    def get_bz_k_points(self):
        return npy.array(self.calc.GetBZKPoints())

    def get_ibz_k_points(self):
        return npy.array(self.calc.GetIBZKPoints())

    def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        nextkpoint, G_I, spin):
        return npy.array(self.calc.GetWannierLocalizationMatrix(
            npy2num(G_I, num.Int), nbands, npy2num(dirG, num.Int), kpoint,
            nextkpoint, spin))
    
    def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        edf, spin):
        # Use initial guess to determine U and C
        init = self.calc.InitialWannier(initialwannier, self.atoms,
                                        npy2num(kpointgrid, num.Int))
        waves = [[self.calc.GetWaveFunction(band, kpt, spin)
                  for band in xrange(self.calc.GetNumberOfBands())]
                 for kpt in xrange(len(self.calc.GetBZKPoints()))]
        init.SetupMMatrix(waves, self.calc.GetBZKPoints())
        c, U = init.GetListOfCoefficientsAndRotationMatrices(
            (self.calc.GetNumberOfBands(), fixedstates, edf))
        U = npy.array(U)
        for ck in c:
            ck[:] = npy.array(ck)
        return c, U

                         
# Some day we will turn on this message:
if 0: 
    from os import env
    if 'NO_OLD_ASE_MESSAGE' not in env:
        print """\
Please consider converting your script to use the new ase
module - it's real simple:

  http://wiki.fysik.dtu.dk/ase/Converting_from_old_ASE
"""
