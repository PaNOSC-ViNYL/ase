from random import randint

import numpy as np

from ase import Atoms
from ase.constraints import dict2constraint
from ase.calculators.calculator import get_calculator, all_properties
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import chemical_symbols, atomic_masses
from ase.io.jsonio import decode
from ase.utils import hill


class FancyDict(dict):
    """Dictionary with keys available as attributes also."""
    def __getattr__(self, key):
        if key not in self:
            return dict.__getattribute__(self, key)
        value = self[key]
        if isinstance(value, dict):
            return FancyDict(value)
        return value

    def __dir__(self):
        return self.keys()  # for tab-completion
        

def atoms2dict(atoms):
    dct = {
        'numbers': atoms.numbers,
        'pbc': atoms.pbc,
        'cell': atoms.cell,
        'positions': atoms.positions,
        'unique_id': '%x' % randint(16**31, 16**32 - 1)}
    if atoms.has('magmoms'):
        dct['initial_magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('charges'):
        dct['initial_charges'] = atoms.get_initial_charges()
    if atoms.has('masses'):
        dct['masses'] = atoms.get_masses()
    if atoms.has('tags'):
        dct['tags'] = atoms.get_tags()
    if atoms.has('momenta'):
        dct['momenta'] = atoms.get_momenta()
    if atoms.constraints:
        dct['constraints'] = [c.todict() for c in atoms.constraints]
    if atoms.calc is not None:
        dct['calculator'] = atoms.calc.name.lower()
        dct['calculator_parameters'] = atoms.calc.todict()
        if len(atoms.calc.check_state(atoms)) == 0:
            dct.update(atoms.calc.results)
    return dct
    
    
class AtomsRow:
    def __init__(self, dct):
        if not isinstance(dct, dict):
            dct = atoms2dict(dct)
        self.dct = dct
        
    def __getattr__(self, key):
        try:
            return self.dct[key]
        except KeyError:
            raise AttributeError
        
    def __contains__(self, key):
        return key in self.__dir__()
        
    def __dir__(self):
        keys = ['mass', 'volume', 'natoms', 'charge',
                'constraints', 'symbols', 'formula']
        keys.extend(self.dct.keys())
        if 'forces' in self.dct:
            keys += ['fmax', 'constrained_forces']
        if 'stress' in self.dct:
            keys.append('smax')
        return keys
        
    def get(self, key, default=None):
        return getattr(self, key, default)

    def count_atoms(self):
        count = {}
        for symbol in self.symbols:
            count[symbol] = count.get(symbol, 0) + 1
        return count
        
    def __getitem__(self, key):
        return getattr(self, key)
        
    @property
    def constraints(self):
        return [dict2constraint(d) for d in self.dct.get('constraints', [])]
        
    @property
    def data(self):
        if 'data' in self.dct:
            return FancyDict(decode(self.dct['data']))
        raise AttributeError
        
    @property
    def natoms(self):
        return len(self.numbers)
        
    @property
    def formula(self):
        return hill(self.numbers)
    
    @property
    def symbols(self):
        return [chemical_symbols[Z] for Z in self.numbers]
        
    @property
    def fmax(self):
        forces = self.constrained_forces
        return (forces**2).sum(1).max()**0.5
        
    @property
    def constrained_forces(self):
        forces = self.forces
        constraints = self.constraints
        if constraints:
            forces = forces.copy()
            for constraint in constraints:
                constraint.adjust_forces(self.positions, forces)
            
        return forces

    @property
    def smax(self):
        return (self.stress**2).max()**0.5

    @property
    def mass(self):
        if 'masses' in self:
            return self.masses.sum()
        return atomic_masses[self.numbers].sum()

    @property
    def volume(self):
        return abs(np.linalg.det(self.cell))

    @property
    def charge(self):
        charges = self.get('inital_charges')
        if charges is None:
            return 0.0
        return charges.sum()

    def toatoms(self, attach_calculator=False,
                add_additional_information=False):
        atoms = Atoms(self.numbers,
                      self.positions,
                      cell=self.cell,
                      pbc=self.pbc,
                      magmoms=self.get('initial_magmoms'),
                      charges=self.get('initial_charges'),
                      tags=self.get('tags'),
                      masses=self.get('masses'),
                      momenta=self.get('momenta'),
                      constraint=self.constraints)
    
        if attach_calculator:
            params = decode(self.get('calculator_parameters', '{}'))
            atoms.calc = get_calculator(self.calculator)(**params)
        else:
            results = {}
            for prop in all_properties:
                if prop in self:
                    results[prop] = self[prop]
            if results:
                atoms.calc = SinglePointCalculator(atoms, **results)
                atoms.calc.name = self.calculator

        if add_additional_information:
            atoms.info = {}
            for key in ['unique_id', 'key_value_pairs', 'data']:
                if key in self:
                    atoms.info[key] = self[key]
                    
        return atoms
