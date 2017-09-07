from ase.build import fcc111
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score
from ase.calculators.emt import EMT

import random

baseslab = fcc111('Cu', size=(3, 3, 3))

# The population size should be at least the number of different compositions
pop_size = 2 * len(baseslab)

db = PrepareDB('hull.db', population_size=pop_size)

# Add the end points manually
print('References:')
for m in ['Cu', 'Au']:
    slab = baseslab.copy()
    slab.numbers[:] = atomic_numbers[m]
    atoms_string = ''.join(slab.get_chemical_symbols())
    slab.set_calculator(EMT())
    e = slab.get_potential_energy()
    print('{0} = {1}'.format(m, e / len(slab)))

    set_raw_score(slab, 0.0)
    db.add_relaxed_candidate(slab, atoms_string=atoms_string)


for i in range(pop_size - 2):
    slab = baseslab.copy()

    # max_tag = 3

    noAu = random.randint(0, len(slab))
    slab.numbers[:noAu] = atomic_numbers['Au']
    random.shuffle(slab.numbers)

    atoms_string = ''.join(slab.get_chemical_symbols())

    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
