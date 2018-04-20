from ase.db import connect
from ase.calculators.emt import EMT


def create_database():
    db = connect('systems.db', append=False)
    from ase.optimize.test.H2 import atoms
    atoms.calc = EMT()
    db.write(atoms)


if __name__ == '__main__':
    create_database()
