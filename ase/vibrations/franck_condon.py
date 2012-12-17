import numpy as np
from operator import mul
from itertools import combinations, product, chain

from ase.units import C, kg, _hbar
from ase.vibrations import Vibrations


class FranckCondon:
    def __init__(self, atoms, vibrations, minfreq=None):
        self.atoms = atoms
        self.vib = vibrations
        self.minfreq = minfreq

        # V = a * v is the combined atom and xyz-index
        self.mm05_V= np.repeat(1. / np.sqrt(atoms.get_masses()), 3) 
        self.shape = (len(self.atoms), 3)

        # from vibrations  
        self.energies = np.real(self.vib.get_energies(
                method='frederiksen')) # eV
        self.frequencies = np.real(self.vib.get_frequencies(
                method='frederiksen')) # cm^-1  
        self.modes = self.vib.modes
        
    def get_Huang_Rhys_from_forces(self, forces):
        """Evaluate Huang-Rhys factors from forces.

        The double harmonic approximation is used."""

        assert(forces.shape == self.shape) 

        # Hesse matrix
        H_VV = self.vib.H
        # sqrt of inverse mass matrix
        mm05_V = self.mm05_V
        # mass weighted Hesse matrix 
        Hm_VV = mm05_V[:, None] * H_VV * mm05_V
        # mass weighted displacements
        Fm_V = forces.flat * mm05_V
        X_V = np.linalg.solve(Hm_VV, Fm_V)
        # projection onto the modes
        modes_VV = self.modes
        d_V = np.dot(modes_VV, X_V)
		
        # Huang-Rhys factors S
        s=1.661*10**(-27)*10**(-20)*1.602*10**(-19)/(1.055*10**(-34))**2 #conversion to SI-units
#        S=s**d_V**2*self.energies/2
        S=s * d_V**2*self.energies/2
        # reshape for minfreq
        indices = np.where(self.frequencies <= self.minfreq)
        S = np.delete(S,indices)
        frequencies = np.delete(self.frequencies, indices)
		
        return S, frequencies


    def Huang_Rhys(self, dR):
        """Evaluate Huang-Rhys factors from displacements.""" 
        assert(dR.shape == self.shape)

        # mass weighted displacements
        d_V = dR.flat * self.mm05_V
        # Huang-Rhys factors S
        s = 1.e-20 / _hbar**2 / kg / C # conversion to SI-units
        S = s * d_V**2 * self.energies / 2

        # reshape for minfreq
        indices = np.where(self.frequencies <= self.minfreq)
        S = np.delete(S, indices)
        frequencies = np.delete(self.frequencies, indices)
		
        return S#, frequencies		

    def displacements_from_forces(self, forces):
        """Evaluate dR from forces in the harmonic approximation"""
        assert(forces.shape == self.shape) 

        # Hesse matrix
        H_VV = self.vib.H
        # displacements of atoms
        dR_V = np.linalg.solve(H_VV, forces.flat)

        return dR_V.reshape(self.shape)

    def get_from_forces(self, order, forces):
        dR_V = self.displacements_from_forces(forces)
        return self.get(order, dR_V)

    def get(self, order, dR_V):
        """Return FC factors up to given order.
        vibrational excitations are normalized to the 0-0 excitation

        Returns franck_condon factors and freqeuncies 
        in an n-dimensional array (n=given order),
        the n'th array contains factors and frequencies for the n'th quanta.
        Only combinations of 2 vibrations are calculated, which should 
        be enough for most of the cases."""

        assert(dR_V.shape == self.shape)
        HR, f_hr = self.Huang_Rhys(dR_V)
		
        def factorial(x):
            if x == 0:
                return 1
            else:
                return x * factorial(x - 1)
            
        # 0 quanta
        n = order
        FC_00 = np.exp(-HR)    
        FC00 = reduce(mul, FC_00)
        FC_0n = [[]*i for i in range(n-1)]
        frequencies_0n = [[]*i for i in range(n-1)]
        for i in range(1, n):
            # single vibration excitation up to given order
            FC_0n[i - 1] = HR**i * np.exp(-HR) / factorial(i)			
            frequencies_0n[i - 1] = self.frequencies * i
            FC_0n[i - 1] = FC_0n[i - 1] / FC_00 # normalize
			
        # combination of two vibrations up to given order
        FC_0nn = [x for x in combinations(chain(*FC_0n), 2)]	
        frequencies_0nn = [x for x in combinations(chain(*frequencies_0n), 2)]
        for i in range(len(FC_0nn)):
            FC_0nn[i] = FC_0nn[i][0] * FC_0nn[i][1]
            frequencies_0nn[i] = frequencies_0nn[i][0] + frequencies_0nn[i][1]
        indices1 = []	
        for i, y in enumerate(self.frequencies):
            ind = [j for j,x in enumerate(frequencies_0nn) 
                   if x / y == 2. or x / y == 3.]
            indices1.append(ind)
        indices1 = [x for x in chain(*indices1)]
        FC_0nn = np.delete(FC_0nn, indices1)	
        frequencies_0nn = np.delete(frequencies_0nn, indices1)
		
        return FC_0n, frequencies_0n, FC_0nn, frequencies_0nn
        
