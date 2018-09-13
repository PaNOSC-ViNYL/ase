from __future__ import print_function
from ase.test import NotAvailable
from ase.build import bulk
from ase.dft.bandgap import bandgap
from ase.calculators.openmx import OpenMX

kpts = (4, 4, 4)


def run(name):
    calc = OpenMX(label=name + '_bandgap', xc='PBE', kpts=kpts)
    si = bulk('Si', crystalstructure='diamond', a=5.43)
    si.calc = calc
    si.get_potential_energy()
    print(name, bandgap(si.calc))
    del si.calc
    # test spin-polarization
    calc = OpenMX(label=name + '_bandgap_spinpol', xc='PBE', kpts=kpts)
    si.set_initial_magnetic_moments([-0.1, 0.1])
    # this should not be necessary in the new ase interface standard ...
    if si.get_initial_magnetic_moments().any():  # spin-polarization
        if name == 'aims':
            calc.set(spin='collinear')
        if name == 'elk':
            calc.set(spinpol=True)
    si.set_calculator(calc)
    si.get_potential_energy()
    print(name, bandgap(si.calc))


# gpaw does not conform to the new ase interface standard:
# https://trac.fysik.dtu.dk/projects/gpaw/ticket/268
names = ['openmx']  # , 'gpaw']
for name in names:
    try:
        run(name)
    except NotAvailable:
        pass
