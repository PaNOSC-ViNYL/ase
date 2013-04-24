from ase.cli import run
import ase.units as units
from ase.db import connect

# warning! parameters are not converged - only an illustration!

run('-x bcc -a 3.6 Li eos -c jacapo -d li.json ' +
    '-p kpts=1.5,symmetry=True,deletenc=True')
run('li.json eos -c jacapo -d li.json ' +
    '-p kpts=1.5,symmetry=True,deletenc=True')

db = connect('li.json')
assert abs(db.get_dict('Li').data.bulk_modulus / kJ * 1e24 - 8.93) < 0.02
assert abs(db.get_dict('li.json').data.bulk_modulus / kJ * 1e24 - 14.73) < 0.02
