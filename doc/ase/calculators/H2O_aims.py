from ase import *
import os

water = Atoms([Atom('H',(1,0,0)),Atom('O',(0,0,0)),Atom('H',(0,1,0))])

calc=Aims(xc='pbe',
          sc_accuracy_etot=1e-6,
          sc_accuracy_eev=1e-3,
          sc_accuracy_rho=1e-6,
          sc_accuracy_forces=1e-4,
          species_dir='/home/hanke/codes/fhi-aims/fhi-aims.workshop/species_defaults/light/',
          run_command='aims.workshop.serial.x')

water.set_calculator(calc)
dynamics = QuasiNewton(water,trajectory='square_water.traj')
dynamics.run(fmax=0.01)

view(water)
