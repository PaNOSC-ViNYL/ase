from ase import *
from ase.calculators.jacapo import *
import os

co = Atoms([Atom('C',[0,0,0]),
            Atom('O',[1.2,0,0])],
            cell=(6,6,6))

calc = Jacapo('Jacapo-co.nc',
              pw=350,
	      dw=350,
              nbands=6,
              kpts=(1,1,1),
              spinpol=False,
              symmetry=False,
              ft=0.01)

co.set_calculator(calc)

#if (os.system('which dacapo.run') == 0):
#    calc.calculate()
#    calc.strip()
