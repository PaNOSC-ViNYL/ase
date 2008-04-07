"""Test the Langevin dynamics.

Tests Langevin dynamics using the EMT Copper potential.

For time reasons, long term temperature average is only calculated with asap.
"""

import time
from ase import *
#from Scientific.Functions.LeastSquares import leastSquaresFit

try:
    import Asap
except ImportError:
    useasap = False
else:
    useasap = True

if useasap:
    nequil = 1000
    nequilprint = 25
    nsteps = 20000
    nprint = 250
    tolerance = 0.01
else:
    nsteps = 2000
    nequil = 500
    nequilprint = 25
    nprint = 25
    tolerance = 0.05
nminor = 25
timestep = 0.5

# Set up atoms in a regular simple-cubic lattice.
a = 3.5
atoms = Atoms('Cu4', a * array([[0,0,0],
                                [0.5, 0.5, 0],
                                [0.5, 0, 0.5],
                                [0, 0.5, 0.5]]),
              cell=(a, a, a))
atoms *= (1, 5, 5)

if useasap:
    atoms.set_calculator(ASAP())
else:
    atoms.set_calculator(EMT())
    print """
WARNING: You have chosen to use the Python-implemented PairPotential.
It is exceedingly slow, and will not accumulate enough statistics for
a good test of the dynamics.  For orders of magnitude better
perfomance, make sure Asap is installed!"""
    


#def targetfunc(params, x):
#    return params[0] * exp(-params[1] * x) + params[2]


def test():
    output = file('Langevin.dat', 'w')
    
    # Make a small perturbation of the momenta
    atoms.set_momenta(1e-6 * random.random([len(atoms), 3]))
    print 'Initializing ...'
    predyn = VelocityVerlet(atoms)
    predyn.run(0.5, 2500)

    temp = 0.01
    frict = 0.001
    dyn = Langevin(atoms, temp, frict)
    print ''
    print ('Testing Langevin dynamics with T = %f eV and lambda = %f' %
           (temp, frict))
    ekin = atoms.get_kinetic_energy()/len(atoms)
    print ekin
    output.write('%.8f\n' % ekin)
    temperatures = [(0, 2.0 / 3.0 * ekin)]
    a = 0.1
    b = frict
    c = temp
    print 'Equilibrating ...'
    tstart = time.time()
    for i in xrange(1,nequil+1):
        dyn.run(timestep, nminor)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        temperatures.append((i*nminor*timestep, 2.0/3.0 * ekin))
        if i % nequilprint == 0:
            #data = array(temperatures)
            #print data.shape
            (a, b, c) = 1,2,3#leastSquaresFit(targetfunc, (0.1, 2*frict, temp),
                             #           temperatures)[0]
            print '%.6f  T_inf = %.6f (goal: %f)  k = %.6f' % \
                  (ekin, c, temp, b)
        output.write('%.8f\n' % ekin)
    tequil = time.time() - tstart
    print 'This took %s minutes.' % (tequil / 60)
    output.write('&\n')
    temperatures = []
    print 'Taking data - this takes ten times longer!'
    tstart = time.time()
    for i in xrange(1,nsteps+1):
        dyn.run(timestep, nminor)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        temperatures.append(2.0/3.0 * ekin)
        if i % nprint == 0:
            tnow = time.time() - tstart
            tleft = (nsteps-i) * tnow / i
            print '%.6f    (time left: %.1f minutes)' % (ekin, tleft/60)
        output.write('%.8f\n' % ekin)
    output.write('&\n')
    temperatures = array(temperatures)
    mean = sum(temperatures) / len(temperatures)
    print 'Mean temperature:', mean, 'eV'
    print
    print 'This test is statistical, and may in rare cases fail due to a'
    print 'statistical fluctuation.'
    print
    print 'Mean temperature:', mean, temp, tolerance
            
    output.close()

if __name__ == '__main__':
    test()
