
''' functions to return data from the old CamposASE2 code. Eventually this module should disappear.
'''

from Dacapo import *
from ASE import *
import numpy as np

def get_wf(ncfile,band,kpt=0,spin=0):

    atoms = Dacapo.ReadAtoms(ncfile)
    calc = atoms.GetCalculator()

    wf = calc.GetWaveFunctionArray(band,kpt,spin)

    
    return wf
    
