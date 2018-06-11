from ase.vibrations.resonant_raman import ResonantRaman

from gpaw.cluster import Cluster
from gpaw import GPAW, FermiDirac
from gpaw.lrtddft import LrTDDFT

h = 0.25
atoms = Cluster('relaxed.traj')
atoms.minimal_box(3.5, h=h)

# relax the molecule
calc = GPAW(h=h, occupations=FermiDirac(width=0.1),
            eigensolver='cg', symmetry={'point_group': False},
            nbands=10, convergence={'eigenstates':1.e-5, 
                                    'bands':4})
atoms.set_calculator(calc)

# use only the 4 converged states for linear response calculation
rr = ResonantRaman(atoms, LrTDDFT, exkwargs={'jend':3})
rr.run()
