from ase.lattice.cubic import FaceCenteredCubic
atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,0], [0,0,1]],
                          size=(3,3,3), symbol='Cu', pbc=(1,1,1))
atoms.rattle(stdev=0.1,seed=42)
