from ase import *
import os

# calculating H2O from scratch using ASE and FHI-aims

# definitions and variables
species='~/codes/fhi-aims/fhi-aims.workshop/species_defaults/light/' # location of the species defaults
run_command = 'aims.workshop.serial.x' # aims run command 
outputname = 'output.water.relax'      # output file 

# make a rectangular water
water = Atoms([Atom('H',(1,0,0)),Atom('O',(0,0,0)),Atom('H',(0,1,0))])

# make an FHI-aims calculator and set it for the water molecule
calc=Aims(xc='pbe',
          sc_accuracy_etot=1e-6,
          sc_accuracy_eev=1e-3,
          sc_accuracy_rho=1e-6,
          sc_accuracy_forces=1e-4,
          species_dir='/home/hanke/codes/fhi-aims/fhi-aims.workshop/species_defaults/light/',
          run_command='aims.workshop.serial.x')

# relax using ASE procedure ... possibly not very efficient but a good start for the moment
water.set_calculator(calc)
dynamics = QuasiNewton(water,trajectory='square_water.traj')
dynamics.run(fmax=0.01)

# read data and display
view(water)
