from ase.cli import run

atoms = run('-x fcc -a 4.04 Al run -c abinit ' +
            '-p xc=PBE,kpts=3.0,ecut=340,toldfe=1e-5,chksymbreak=0')
