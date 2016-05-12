from ase import Atoms
from ase.db import connect
import ase.db.app as app
c = connect('test.db', append=False)
plot = {'title': 'A test',
        'data': [{'label': 't1', 'x': [0, 1, 2], 'y': [1, 2, 0]},
                 {'label': 't2', 'style': 'y--',
                  'x': [0, 1, 2], 'y': [[2, 2, 1], [1, 1, 1]]}],
        'xlabel': 'x',
        'ylabel': 'y'}
c.write(Atoms('H2O'), foo='bar', data={'plots': [plot]})
app.db = c
app.app.testing = True
d = app.app.test_client().get('/')
print(d)
d = app.app.test_client().get('/id/1')
print(d)
d = app.app.test_client().get('/plot/0-1.png')
print(d)
