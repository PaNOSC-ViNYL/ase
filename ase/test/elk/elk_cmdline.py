from ase.cli import run

run('-x fcc -a 4.04 Al run -c elk -p kpts=3.0,xc=PBE,rgkmax=5.0,tforce=True')
