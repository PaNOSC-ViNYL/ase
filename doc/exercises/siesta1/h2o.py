from ase import *
from gpaw import Calculator

h2o = Atoms(symbols='H2O',
            positions=[( 0.776070, 0.590459, 0.00000),
                       (-0.776070, 0.590459, 0.00000),
                       (0.000000,  -0.007702,  -0.000001)],
            pbc=(1,1,1))

h2o.center(vacuum=6.)

traj = PickleTrajectory('h2o.traj', 'w', h2o)


e_shifts = [0.01,0.1,0.2,0.3,0.4,0.5]

for e_s in e_shifts:
    calc = Siesta('h2o',basis='SZ',meshcutoff=200.0*Ry,mix=0.5,pulay=4)
    calc.set_fdf('PAO.EnergyShift',e_s * eV)
    calc.set_fdf('PAO.SplitNorm',0.15)
    calc.set_fdf('DM.UseSaveDM','Y')
    h2o.set_calculator(calc)
    dyn = QuasiNewton(h2o)
    dyn.run(fmax=0.02)
    E = h2o.get_potential_energy()
    traj.write()
    print "%.2f %.4f" % (e_s,E)


