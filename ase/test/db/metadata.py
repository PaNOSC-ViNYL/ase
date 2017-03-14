from ase.db import connect
from ase import Atoms
db = connect('hmm.json')
db.write(Atoms('H'), answer=42)
db.metadata = {'columns': ['formula', 'answer']}
