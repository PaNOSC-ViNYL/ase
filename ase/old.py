import numpy as npy


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
    
    def get_magmoms(self):
        return npy.array(self.atoms.GetMagneticMoments())
    
    def get_charges(self):
        return None
    
    def get_cell(self):
        return npy.array(self.atoms.GetUnitCell())

    def get_pbc(self):
        return npy.array(self.atoms.GetBoundaryConditions(), bool)


class ListOfAtoms:
    def __init__(self, atoms):
        self.atoms = atoms
        from Numeric import array
        self.array = array
        
    def GetCartesianPositions(self):
        return self.array(self.atoms.get_positions())
    
    def GetAtomicNumbers(self):
        return self.array(self.atoms.get_atomic_numbers())
    
    def GetUnitCell(self):
        return self.array(self.atoms.get_cell())
    
    def GetBoundaryConditions(self):
        return tuple(self.atoms.get_pbc())

    def GetTags(self):
        return self.array(self.atoms.get_tags())

    def GetMagneticMoments(self):
        magmoms = self.atoms.get_magnetic_moments()
        if magmoms is None:
            magmoms = npy.zeros(len(self))
        return self.array(magmoms)

    def GetCount(self):
        return 0
    
    def __len__(self):
        return len(self.atoms)

    def SetCalculator(self, calc):
        calc._SetListOfAtoms(self)
    
class OldASECalculatorWrapper:
    def __init__(self, calc, atoms):
        self.calc = calc
        self.atoms = ListOfAtoms(atoms)
        self.atoms.SetCalculator(calc)
        
    def get_potential_energy(self, atoms):
        return self.calc.GetPotentialEnergy()

    def get_forces(self, atoms):
        return self.calc.GetCartesianForces()

    def get_number_of_bands(self):
        return self.calc.GetNumberOfBands()

    def get_kpoint_weights(self):
        return npy.array(self.calc.GetIBZKPointWeights())

    def get_number_of_spin(self):
        return 1 + int(self.calc.GetSpinPolarized())

    def get_eigenvalues(self, k=0, s=0):
        return npy.array(self.calc.GetEigenvalues(k, s))

    def get_fermi_level(self):
        return self.calc.GetFermiLevel()

    def get_wave_function_array(self, n=0, k=0, s=0):
        return npy.array(self.calc.GetWaveFunctionArray(n, k, s))



                         
# Some day we will turn on this message:
if 0: 
    from os import env
    if 'NO_OLD_ASE_MESSAGE' not in env:
        print """\
Please consider converting your script to use the new ase
module - it's real simple:

  http://wiki.fysik.dtu.dk/ase/Converting_from_old_ASE
"""
