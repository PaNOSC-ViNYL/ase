import os
import ase

tests = [
    'H2.py',
    #'C2_Cu100.py', # I don't understand this system
    'CO_Au111.py',
    'N2Ru-relax.py',
    'nanoparticle.py', #SLOW
]

for test in tests:
    filename = ase.__path__[0] + '/optimize/test/' + test
    execfile(filename, {})

