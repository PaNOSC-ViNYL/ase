from ase.db import connect
from ase import Atoms
db = connect('md.json')

# Data for a plot:
plot = {'a': [0, 1, 2],
        'b': [1.2, 1.1, 1.0],
        'abplot': {'title': 'Example',
                   'data': [{'x': 'a',
                             'y': 'b',
                             'label': 'label1',
                             'style': 'o-g'}],
                   'xlabel': 'blah-blah [eV]',
                   'ylabel': 'Answers'}}

db.write(Atoms('H'), answer=42, kind='atom')
db.write(Atoms('H2O'), answer=117, kind='molecule', data=plot)
db.metadata = {
    'default_columns': ['formula', 'answer'],
    'special_keys': {'kind': ('COMBO', ['atom', 'molecule'])},
    'key_descriptions': {
        'kind': ('Type', 'Type of system', 'string'),
        'answer': ('Answer', 'Answer to question', 'int')},
    'summary_sections': [
        ['Basic Properties',
         ['Item', ['energy', 'fmax', 'charge', 'mass', 'magmom', 'volume']],
         ['Things', ['answer', 'kind']],
         ['STRUCTUREPLOT'],
         ['AXIS'],
         ['PLOT', 'abplot']],
        ['Calculation Details',
         ['Calculator Setting', ['calculator']],
         ['FORCES']]]}
