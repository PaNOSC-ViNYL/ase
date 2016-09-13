import numpy as np

import sander
from ase.units import kcal, mol
from ase.calculators.calculator import Calculator

class SANDER(Calculator):
    implemented_properties=['energy', 'forces']

    def __init__(self, atoms=None, label=None, top=None, crd=None, mm_options=None, qm_options=None, **kwargs):
        Calculator.__init__(self, label, atoms)
        if qm_options is not None:
            sander.setup(top, crd.coordinates, crd.box, mm_options, qm_options)
        else:
            sander.setup(top, crd.coordinates, crd.box, mm_options)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if system_changes:
            if 'energy' in self.results:
                del self.results['energy']
            if 'forces' in self.results:
                del self.results['forces']
        if 'energy' not in self.results:
            crd=np.reshape(atoms.get_positions(),(1,len(atoms),3))
            sander.set_positions(crd)
            e, f = sander.energy_forces()
            self.results['energy']=e.tot*kcal/mol
            self.results['forces']=np.reshape(np.array(f), (len(atoms),3))*kcal/mol
        if 'forces' not in self.results:
            crd=np.reshape(atoms.get_positions(),(1,len(atoms),3))
            sander.set_positions(crd)
            e, f = sander.energy_forces()
            self.results['energy']=e.tot*kcal/mol
            self.results['forces']=np.reshape(np.array(f), (len(atoms),3))*kcal/mol
