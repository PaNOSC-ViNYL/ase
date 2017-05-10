from __future__ import print_function
from ase.data import s22
from ase.calculators.tip3p import TIP3P, epsilon0, sigma0
from ase.calculators.qmmm import EIQMMM, LJInteractions, Embedding
from gpaw import GPAW

# Create system
atoms = s22.create_s22_system('Water_dimer')
atoms.center(vacuum=4.0)

# Make QM atoms selection of first water molecule:
qm_idx = range(3)

# Set up interaction & embedding object
interaction = LJInteractions({('O', 'O'): (epsilon0, sigma0)})
embedding = Embedding(rc=0.02)  # Short range analytical potential cutoff 

# Set up calculator
atoms.calc = EIQMMM(qm_idx,
                    GPAW(txt='qm.out'),
                    TIP3P(),
                    interaction,
                    embedding=embedding,
                    vacuum=None,  # if None, QM cell = MM cell
                    output='qmmm.log')

print(atoms.get_potential_energy())
