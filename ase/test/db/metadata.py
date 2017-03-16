from ase.db import connect
from ase import Atoms
db = connect('md.json')
db.write(Atoms('H'), answer=42)
db.metadata = {'columns': ['formula', 'answer']}
