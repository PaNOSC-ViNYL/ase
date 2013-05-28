from ase.cli import run as cli
from ase.db.cli import run as db
for name in ['x.json', 'x.sqlite']:#, 'postgres://localhost']:
    cli('H H2O O O2 run -d %s' % name)
    cli('H2 optimize -d %s' % name)
    cli('O2 optimize -d %s' % name)
    cli('-x fcc Cu eos -d %s' % name)
    ids = db('%s' % name)
    ids.sort()
    print ids
    db('%s id=H --delete --yes' % name)
    db('%s H>0 -k emt' % name)
    ids = db('%s' % name)
