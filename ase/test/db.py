from ase.test import cli
from ase.db import connect

cmd = """
ase-build H | ase-run emt -d x.json &&
ase-build H2O | ase-run emt -d x.json &&
ase-build O2 | ase-run emt -d x.json &&
ase-build H2 | ase-run emt -f 0.02 -d x.json &&
ase-build O2 | ase-run emt -f 0.02 -d x.json &&
ase-build -x fcc Cu | ase-run emt -E 5 -d x.json &&
ase-db x.json id=H --delete --yes &&
ase-db x.json "H>0" -k hydro"""

for name in ['x.json', 'x.db']:#, 'postgres://localhost']:
    cli(cmd.replace('x.json', name))
    con = connect(name)
    assert len(list(con.select())) == 4
    assert len(list(con.select('hydro'))) == 2

