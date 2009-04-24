from gpaw import * 
from __init__ import STM

t = GPAW('Altip.gpw',txt=None)
s = GPAW('Al110.gpw',txt=None) 


stm = STM(surfacecalc=s, tipcalc=t,tipapex=[2.8638,4.0500,6.0000],
          surface_extention=5.0, tip_cutoff=[3.5,4.0,8.5], vacuum=3.0,
          surface_z = 9.432)
