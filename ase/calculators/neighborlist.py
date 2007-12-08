class NeighborList:
    def __init__(self, atoms, cutoffs, skin=0.3):
        positions = atoms.get_positions()
        self.pbc = atoms.get_pbc()
        self.cell = atoms.get_cell()

        rcmax = max(cutoffs.values())
        
        icell = npy.linalg.inv(self.cell)
        scaled = npy.dot(self.positions, icell)
        N = []
        for i in range(3):
            if self.pbc[i]:
                scaled[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(npy.dot(v, v))
                N.append(int(rcmax / h) + 1)
            else:
                N.append(0)

        R = npy.dot(scaled, self.cell)
