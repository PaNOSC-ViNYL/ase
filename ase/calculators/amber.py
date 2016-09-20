import numpy as np

import sander
from ase.units import kcal, mol
from ase.calculators.calculator import Calculator

def map(atoms, top):

    p = np.zeros((2,len(atoms)), dtype="int")

    elements = atoms.get_chemical_symbols()
    unique_elements = np.unique(atoms.get_chemical_symbols())

    for i in range(len(unique_elements)):
        idx = 0
        for j in range(len(atoms)):
            if elements[j] == unique_elements[i]:
                idx += 1
                symbol = unique_elements[i]+np.str(idx)
                for k in range(len(atoms)):
                    if top.atoms[k].name == symbol:
                        p[0,k] = j
                        p[1,j] = k
                        break
    return p

class SANDER(Calculator):
    implemented_properties=['energy', 'forces']

    def __init__(self, atoms=None, label=None, top=None, crd=None, mm_options=None, qm_options=None, permutation=None, **kwargs):
        Calculator.__init__(self, label, atoms)
        self.permutation = permutation
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
            if self.permutation is None:
                crd=np.reshape(atoms.get_positions(),(1,len(atoms),3))
            else:
                crd=np.reshape(atoms.get_positions()[self.permutation[0,:]],(1,len(atoms),3))
            sander.set_positions(crd)
            e, f = sander.energy_forces()
            self.results['energy']=e.tot*kcal/mol
            if self.permutation is None:
                self.results['forces']=np.reshape(np.array(f), (len(atoms),3))*kcal/mol
            else:
                ff = np.reshape(np.array(f), (len(atoms),3))*kcal/mol
                self.results['forces']=ff[self.permutation[1,:]]
        if 'forces' not in self.results:
            if self.permutation is None:
                crd=np.reshape(atoms.get_positions(),(1,len(atoms),3))
            else:
                crd=np.reshape(atoms.get_positions()[self.permutation[0,:]],(1,len(atoms),3))
            sander.set_positions(crd)
            e, f = sander.energy_forces()
            self.results['energy']=e.tot*kcal/mol
            if self.permutation is None:
                self.results['forces']=np.reshape(np.array(f), (len(atoms),3))*kcal/mol
            else:
                ff = np.reshape(np.array(f), (len(atoms),3))*kcal/mol
                self.results['forces']=ff[self.permutation[1,:]]
