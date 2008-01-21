class NeighborList:
    def __init__(self, cutoffs, skin=0.3):
        self.cutoffs = dict([(Z, rc + skin) for Z, rc in cutoffs.items()])
        self.skin = skin
        self.nupdates = 0

    def update(self, atoms):
        if (self.positions is None or
            ((atoms.positions - self.positions)**2).max() > self.skin**2):
            self.build(atoms)

    def build(self, atoms):
        self.positions = atoms.get_positions()
        pbc = atoms.get_pbc()
        cell = atoms.get_cell()

        rcmax = max(cutoffs.values())
        
        icell = npy.linalg.inv(cell)
        scaled = npy.dot(self.positions, icell)
        scaled0 = scaled.copy()

        N = []
        for i in range(3):
            if pbc[i]:
                scaled0[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(npy.dot(v, v))
                N.append(int(rcmax / h) + 1)
            else:
                N.append(0)

        offsets = npy.empty((len(atoms), 3), int)
        (scaled - scaled0).round(out=offsets)
        positions0 = npy.dot(scaled0, cell)

        self.neighbors = [npy.empty(0, int) for a in range(natoms)]
        self.displacements = [npy.empty((0, 3), int) for a in range(natoms)]
        for n1 in range(-N[0], N[0] + 1):
            for n2 in range(-N[1], N[1] + 1):
                for n3 in range(-N[2], N[2] + 1):
                    displacement = npy.dot((n1, n2, n3), cell)
                    for a1 in range(natoms):
                        d = positions0 + displacement - positions0[a1]
                        for a2 in range(natoms):
                            mask = 
    def get_neighbors(self, a):
        return self.nbors[a], displ
