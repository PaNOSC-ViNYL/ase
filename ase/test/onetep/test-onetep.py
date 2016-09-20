from __future__ import print_function
from ase.build import molecule
from ase.calculators import onetep
from os.path import isfile, dirname, abspath, join
from ase.calculators.calculator import equal

def main():
    mol = molecule('H2O')
    mol.center(8)
    calc = onetep.Onetep(label="water")
    # Tests conducted with the JTH PAW data set.
    # http://www.abinit.org/downloads/PAW2
    prefix = dirname(abspath(__file__))
    h_path = join(prefix, "H.abinit")
    o_path = join(prefix, "O.abinit")
    if not (isfile(h_path) and isfile(o_path)):
        raise Exception("""You must supply PAW data sets for
            hydrogen and oxygen to run this test.
            Please see http://www.abinit.org/downloads/PAW2
            for suitable data. ONETEP takes PAW data sets in the
            abinit format. I need H.abinit and O.abinit""")
    calc.set_pseudos([('H', h_path), ('O', o_path)])
    calc.set(paw=True, xc="PBE")
    mol.set_calculator(calc)
    energy = mol.get_total_energy()
    ref_energy = -471.576999443
    if equal(energy, ref_energy, 1e-8):
        print("Passed: Expected energy: %.9f, got %.9f" % (ref_energy, energy))
    else:
        raise Exception("Test failed. Expected energy %.9f, got %.9f"
                        % (ref_energy, energy))

main()
