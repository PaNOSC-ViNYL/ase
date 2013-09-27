from ase.test import cli, require
import ase.units as units
from ase.db import connect

# warning! parameters are not converged - only an illustration!

require('jacapo')
cli("""ase-build -x bcc -a 3.6 Li | \
ase-run jacapo ... -d li.json -p kpts=1.5,symmetry=True,deletenc=True""")

db = connect('li.json')
assert abs(db.get_dict('Li').data.bulk_modulus / kJ * 1e24 - 8.93) < 0.02
