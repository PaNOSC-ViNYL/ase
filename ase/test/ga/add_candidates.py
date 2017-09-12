from ase.ga.data import PrepareDB
from ase.ga.data import DataConnection
from ase.build import fcc111
from ase.ga import set_raw_score

import os

db_file = 'gadb.db'
if os.path.isfile(db_file):
    os.remove(db_file)

db = PrepareDB(db_file)

slab1 = fcc111('Ag', size=(2, 2, 2))
db.add_unrelaxed_candidate(slab1)

slab2 = fcc111('Cu', size=(2, 2, 2))
set_raw_score(slab2, 4)
db.add_relaxed_candidate(slab2)
assert slab2.info['confid'] == 3

db = DataConnection(db_file)
assert db.get_number_of_unrelaxed_candidates() == 1

slab3 = db.get_an_unrelaxed_candidate()
slab3[0].symbol = 'Au'
db.add_unrelaxed_step(slab3, 'mutated: Parent 2')

try:
    db.add_relaxed_step(slab3)
except AssertionError:
    pass
set_raw_score(slab3, 3)
db.add_relaxed_step(slab3)

from ase.ga.offspring_creator import OffspringCreator
slab4 = OffspringCreator.initialize_individual(slab1,
                                               fcc111('Au', size=(2, 2, 2)))
set_raw_score(slab4, 67)
db.add_relaxed_step(slab4)
assert slab4.info['confid'] == 6

os.remove(db_file)
