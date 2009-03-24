
''' functions to return data from the old CamposASE2 code. Eventually this module should disappear.
'''

from Dacapo import *
from ASE import *
import numpy as np

def get_wf(ncfile,band,kpt=0,spin=0):

    atoms = Dacapo.ReadAtoms(ncfile)
    calc = atoms.GetCalculator()

    #this apparently does not have phase information
    wf = calc.GetWaveFunctionArray(band,kpt,spin)

    
    return wf


#Wannier function code Jacapo should return these functions for ase.
#from the new dacapo.py

def get_pseudo_wave_function(ncfile, band=0, k=0, spin=0, pad=True):

    atoms =  Dacapo.ReadAtoms(ncfile)
    calc = atoms.GetCalculator()
    kpt = calc.GetBZKPoints()[k]
    state = calc.GetElectronicStates().GetState(band=band,
                                                spin=spin,
                                                kptindex=k)

    # Get wf, without bolch phase (Phase = True doesn't do anything!)
    wave = state.GetWavefunctionOnGrid(phase=False)
    
    # Add bloch phase if this is not the Gamma point
    if np.all(kpt == 0):
        return wave
    coord = state.GetCoordinates()
    phase = coord[0] * kpt[0] + coord[1] * kpt[1] + coord[2] * kpt[2]
    return np.array(wave) * np.exp(-2.j * np.pi * phase) # sign! XXX

#return np.array(self.calc.GetWaveFunctionArray(n, k, s)) # No phase!

        
def get_wannier_localization_matrix(ncfile,
                                    nbands,
                                    dirG,
                                    kpoint,
                                    nextkpoint,
                                    G_I,
                                    spin):

    atoms =  Dacapo.ReadAtoms(ncfile)
    calc = atoms.GetCalculator()

    locmat = calc.GetWannierLocalizationMatrix(G_I=G_I.tolist(),
                                               nbands=nbands,
                                               dirG=dirG.tolist(),
                                               kpoint=kpoint,
                                               nextkpoint=nextkpoint,
                                               spin=spin)

    
    return np.array(locmat)
    
def initial_wannier(ncfile,
                    initialwannier,
                    kpointgrid,
                    fixedstates,
                    edf,
                    spin):
    atoms =  Dacapo.ReadAtoms(ncfile)
    calc = atoms.GetCalculator()
    
    # Use initial guess to determine U and C
    init = calc.InitialWannier(initialwannier, self.atoms,
                                    np2num(kpointgrid, num.Int))
    
    states = calc.GetElectronicStates()
    waves = [[state.GetWaveFunction()
              for state in states.GetStatesKPoint(k, spin)]
             for k in self.calc.GetIBZKPoints()] 
    
    init.SetupMMatrix(waves, calc.GetBZKPoints())
    c, U = init.GetListOfCoefficientsAndRotationMatrices(
        (calc.GetNumberOfBands(), fixedstates, edf))
    U = np.array(U)
    for k in range(len(c)):
        c[k] = np.array(c[k])
    return c, U

