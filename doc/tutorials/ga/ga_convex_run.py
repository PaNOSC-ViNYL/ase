from ase.ga.population import RankFitnessPopulation
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.slab_operators import CutSpliceSlabCrossover, RandomSlabPermutation, SymmetrySlabPermutation
from ase.ga import set_raw_score

from ase.calculators.emt import EMT

refs = {'Cu': 0.23780976087291691, 'Au': 1.2957112575371597}


def get_mixing_energy(atoms):
    atoms.set_calculator(EMT())
    e = atoms.get_potential_energy()
    print(e)
    syms = atoms.get_chemical_symbols()
    for m in set(syms):
        e -= syms.count(m) * refs[m]
    print(e)
    return e


def get_comp(atoms):
    return atoms.get_chemical_formula()


# Specify the number of generations this script will run
num_gens = 40

db = DataConnection('hull.db')

# Retrieve saved parameters
pop_size = db.get_param('population_size')

# Specify the procreation operators for the algorithm
# Try and play with the mutation operators that move to nearby
# places in the periodic table
oclist = ([1, 1, 1], [CutSpliceSlabCrossover(),
                      RandomSlabPermutation(), SymmetrySlabPermutation()])
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
