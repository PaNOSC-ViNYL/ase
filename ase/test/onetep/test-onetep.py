from ase.build import molecule
from ase.calculators import onetep
from os.path import isfile

mol = molecule('H2O')
mol.center(8)
calc = onetep.Onetep(label="water")
#Tests conducted with the JTH PAW data set.
#http://www.abinit.org/downloads/PAW2
if not (isfile("H.abinit") and isfile("O.abinit")):
    raise Exception( """You must supply PAW data sets for
        hydrogen and oxygen to run this test.
        Please see http://www.abinit.org/downloads/PAW2
        for suitable data. ONETEP takes PAW data sets in the
        abinit format. I need H.abinit and O.abinit""")
calc.set_pseudos([('H',"H.abinit"),('O',"O.abinit")])
calc.set(PAW=True,xc="PBE")
mol.set_calculator(calc)
energy =  mol.get_total_energy()
ref_energy = -471.576999443
if abs(energy-ref_energy) < 1e-5:
    print "Test passed."
else:
    print "Test failed. Expected %f, got %f" % (energy,ref_energy)
