import numpy as npy

def distance(atoms, a1, a2):
    return npy.linalg.norm(atoms.positions[a2] - atoms.positions[a1])

    
