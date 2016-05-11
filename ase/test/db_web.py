from ase import Atoms
from ase.db import connect
import ase.db.app as app
c = connect('test.db', append=False)
plot = {'tag': 'test',
        'text': 'A test',
        'names': ['t1'],
        't1': ([0, 1, 2], [[1, 2, 0]]),
        'xlabel': 'x',
        'ylabel': 'y'}
c.write(Atoms('H2O'), foo='bar', data={'plots': [plot]})
app.db = c
app.app.testing = True
d = app.app.test_client().get('/')
print(d)
#d = app.app.test_client().get('/id/1').data
d = app.app.test_client().get('/plot/test-1.png')#.data
print(d)
