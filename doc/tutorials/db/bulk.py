from ase.build import bulk
from ase.calculators.emt import EMT
from ase.eos import calculate_eos
from ase.db import connect

db = connect('bulk.db')
for symb in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
    atoms = bulk(symb, 'fcc')
    atoms.calc = EMT()
    eos = calculate_eos(atoms)
    v, e, B = eos.fit()
    atoms.cell *= (v / atoms.get_volume())**(1 / 3)
    atoms.get_potential_energy()
    db.write(atoms, bm=B)
