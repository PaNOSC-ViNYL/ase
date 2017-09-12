from ase.ga.population import RankFitnessPopulation
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.slab_operators import (CutSpliceSlabCrossover,
                                   RandomSlabPermutation,
                                   SymmetrySlabPermutation,
                                   RandomCompositionMutation)
from ase.ga import set_raw_score

from ase.calculators.emt import EMT

# Connect to the database containing all candidates
db = DataConnection('hull.db')

# Retrieve saved parameters
pop_size = db.get_param('population_size')
refs = db.get_param('references')
metals = db.get_param('metals')


def get_mixing_energy(atoms):
    atoms.set_calculator(EMT())
    e = atoms.get_potential_energy()
    syms = atoms.get_chemical_symbols()
    for m in set(syms):
        e -= syms.count(m) * refs[m]
    return e


def get_comp(atoms):
    return atoms.get_chemical_formula()


# Specify the number of generations this script will run
num_cands = 500

# Specify the procreation operators for the algorithm
# Try and play with the mutation operators that move to nearby
# places in the periodic table
oclist = ([3, 1, 1, 1], [CutSpliceSlabCrossover(),
                         RandomSlabPermutation(), SymmetrySlabPermutation(),
                         RandomCompositionMutation()])
operation_selector = OperationSelector(*oclist)

# Pass parameters to the population instance
pop = RankFitnessPopulation(data_connection=db,
                            population_size=pop_size,
                            variable_function=get_comp)

# Relax the starting population
while db.get_number_of_unrelaxed_candidates() > 0:
    a = db.get_an_unrelaxed_candidate()
    set_raw_score(a, get_mixing_energy(a))
    db.add_relaxed_step(a)
pop.update()

for _ in range(num_cands):
    # Select parents for a new candidate
    parents = pop.get_two_candidates()

    # Select an operator and use it
    op = operation_selector.get_operator()
    print(op.descriptor)
    offspring, desc = op.get_new_individual(parents)
    # An operator could return None if an offspring cannot be formed
    # by the chosen parents
    if offspring is None:
        continue

    set_raw_score(offspring, get_mixing_energy(offspring))
    print(desc, offspring.info['key_value_pairs']
          ['raw_score'], offspring.get_chemical_formula())
    db.add_relaxed_step(offspring)

    pop.update()
