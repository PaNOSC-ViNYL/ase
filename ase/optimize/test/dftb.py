import ase.db
from ase.optimize.test.test import (test_optimizer, all_optimizers,
                                    get_optimizer)
from ase.calculators.dftb import Dftb


db1 = ase.db.connect('systems.db')
db = ase.db.connect('results-dftb.db')


def dftb(txt):
    return Dftb(#kpts=2.0,
                label=txt)


systems = [row.toatoms() for row in db1.select()]

for opt in all_optimizers:
    optimizer = get_optimizer(opt)
    test_optimizer(systems, optimizer, dftb, 'dftb-', db)
