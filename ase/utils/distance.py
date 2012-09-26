import numpy as np

def distance(s1, s2):
    """Get the distance between two structures."""

    s1 = s1.copy()
    s2 = s2.copy()
    for s in [s1, s2]:
        s.translate(-s.get_center_of_mass())
    s2pos = 1. * s2.get_positions()
    
    def align(struct, xaxis='x', yaxis='y'):
        Is, Vs = struct.get_moments_of_inertia(True)
        IV = zip(Is, Vs)
        IV.sort(cmp=lambda x, y: cmp(y[0], x[0]))
        struct.rotate(IV[0][1], xaxis)
        
        Is, Vs = struct.get_moments_of_inertia(True)
        IV = zip(Is, Vs)
        IV.sort(cmp=lambda x, y: cmp(x[0], y[0]))
        struct.rotate(IV[1][1], yaxis)
#        print struct.get_moments_of_inertia(True)

    align(s1)

    def dd(s1, s2):
        return np.linalg.norm(s1.get_positions() - s2.get_positions())

    dists = []
    # principles 
    for x, y in zip(['x','-x','x','-x'], ['y','y','-y','-y']):
        s2.set_positions(s2pos)
        align(s2, x, y)
        dists.append(dd(s1, s2))
   
    return min(dists)
 
