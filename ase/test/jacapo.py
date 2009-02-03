from ase import *
from ase.calculators.Jacapo import *

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
calc.calculate()
#print 'Status', calc.get_status()
#print 'Setting status to \'finished\''
#calc.set_status('finished')
#print 'Status', calc.get_status()
#print 'Energy = %f eV' % calc.get_potential_energy()
#print 'Forces = (eV/ang)\n', co.get_forces()
calc.strip()
