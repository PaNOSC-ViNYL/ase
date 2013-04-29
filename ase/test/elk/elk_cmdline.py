from ase.cli import run

# warning! parameters are not converged - only an illustration!
atoms = run('-x fcc -a 4.04 Al run -c elk -p ' +
            'tasks=0,kpts=1.5,rgkmax=5.0,tforce=True,' +
            'smearing=(fermi-dirac,0.05)')
