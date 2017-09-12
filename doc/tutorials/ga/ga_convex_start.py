from ase.build import fcc111
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score
from ase.calculators.emt import EMT

import random

metals = ['Cu', 'Au']
baseslab = fcc111(metals[0], size=(3, 3, 3), vacuum=5)

# The population size should be at least the number of different compositions
pop_size = 2 * len(baseslab)

# Create the references (pure slabs) manually
pure_slabs = []
refs = {}
print('References:')
for m in metals:
    slab = baseslab.copy()
    slab.numbers[:] = atomic_numbers[m]
    atoms_string = ''.join(slab.get_chemical_symbols())
    slab.set_calculator(EMT())
    e = slab.get_potential_energy()
    e_per_atom = e / len(slab)
    refs[m] = e_per_atom
    print('{0} = {1:.3f} eV/atom'.format(m, e_per_atom))

    set_raw_score(slab, 0.0)
    pure_slabs.append(slab)

db = PrepareDB('hull.db', population_size=pop_size,
               references=refs, metals=metals)

for slab in pure_slabs:
    db.add_relaxed_candidate(slab,
                             atoms_string=''.join(slab.get_chemical_symbols()))


for i in range(pop_size - 2):
    slab = baseslab.copy()

    # max_tag = 3

    noAu = random.randint(1, len(slab) - 1)
    slab.numbers[:noAu] = atomic_numbers['Au']
    random.shuffle(slab.numbers)

    atoms_string = ''.join(slab.get_chemical_symbols())

    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
