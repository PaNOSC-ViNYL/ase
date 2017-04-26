from ase.db.web import creates

title = 'TEST'

default_columns = ['formula', 'answer', 'kind']

special_keys = ['kind', 'foo']

key_descriptions = {
    'kind': ('Type', 'Type of system', 'string', ''),
    'answer': ('Answer', 'Answer to question', 'int', 'eV')}


@creates('xy.png')
def xy(row):
    import matplotlib.pyplot as plt
    ax = plt.figure().add_subplot(111)
    ax.plot([0, 1, 2, 3, 0, 1])
    plt.savefig('xy.png')


stuff = ('Stuff', ['energy', 'fmax', 'charge', 'mass', 'magmom', 'volume'])
things = ('Things', ['answer', 'kind'])
calc = ('Calculator Setting', ['calculator'])

layout = [
    ('Basic Properties',
     [[stuff, things],
      ['ATOMS', 'CELL']]),
    ('Calculation Details',
     [[calc, 'FORCES'],
      ['xy.png', None],
      ['abc.png', None]])]
