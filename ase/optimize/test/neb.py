#/usr/bin/env python
#PBS -l nodes=4:ppn=8
#PBS -l walltime=13:00:00
import development
from ase import QuasiNewton, FixAtoms, EMT, NEB
from ase.lattice.surface import fcc100, add_adsorbate
from ase.optimize.test import run_test
from gpaw import GPAW
from gpaw import Mixer
from gpaw.poisson import PoissonSolver

name = 'neb'

def get_atoms():
    # 2x2-Al(001) surface with 3 layers and an
    # Au atom adsorbed in a hollow site:
    slab = fcc100('Al', size=(2, 2, 3))
    add_adsorbate(slab, 'Au', 1.7, 'hollow')
    slab.center(axis=2, vacuum=4.0)

    # Fix second and third layers:
    mask = [atom.tag > 1 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    # Use EMT potential:
    slab.set_calculator(EMT())

    # Initial state:
    qn = QuasiNewton(slab, logfile=None)
    qn.run(fmax=0.05)
    initial = slab.copy()

    # Final state:
    slab[-1].x += slab.get_cell()[0, 0] / 2
    qn = QuasiNewton(slab, logfile=None)
    qn.run(fmax=0.05)
    final = slab.copy()

    # Setup a NEB calculation
    constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

    images = [initial]
    for i in range(3):
        image = initial.copy()
        image.set_constraint(constraint)
        images.append(image)

    images.append(final)

    neb = NEB(images)
    neb.interpolate()

    def set_calculator(calc):
        for image in neb.images:
            image.set_calculator(calc) 
    neb.set_calculator = set_calculator

    return neb

def get_calculator_emt():
    return EMT()

def get_calculator_gpaw():
    calc = GPAW(h=0.2,
                mode = 'lcao',
                basis = 'szp',
                nbands=-5,
                xc='LDA',
                width=0.1,
                mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0),
                poissonsolver=PoissonSolver(nn='M', relax='GS'),
                convergence={'energy':1e-4, 'bands':-3},
                stencils=(3,3),
                txt='neb.txt')
    return calc

run_test(get_atoms, get_calculator_gpaw, name + '-gpaw', 'GPAW (lcao)')
run_test(get_atoms, get_calculator_emt, name + '-emt', 'EMT')
