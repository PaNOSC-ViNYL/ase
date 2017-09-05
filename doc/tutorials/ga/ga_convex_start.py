from ase.build import fcc111
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB

import random

baseslab = fcc111('Cu', size=(3, 3, 3))

# The population size should be at least the number of different compositions
pop_size = 2 * len(baseslab)

db = PrepareDB('hull.db', population_size=pop_size)

for i in range(pop_size):
    slab = baseslab.copy()

    # max_tag = 3

    noAu = random.randint(0, len(slab))
    slab.numbers[:noAu] = atomic_numbers['Au']
    random.shuffle(slab.numbers)

    atoms_string = ''.join(slab.get_chemical_symbols())

    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
