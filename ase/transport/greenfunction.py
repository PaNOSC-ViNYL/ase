import numpy as npy

class GreenFunction:
    """Equilibrium retarded Green function."""
    
    def __init__(self, H, S=None, selfenergies=[], eta=1e-4):
        self.H = H
        self.S = S
        self.selfenergies = selfenergies
        self.eta = eta
        self.energy = None
        self.Ginv = npy.empty(H.shape, complex)

    def __call__(self, energy, inverse=False):
        if energy != self.energy:
            self.energy = energy
            z = energy + self.eta * 1.j
            if self.S is None:
                self.Ginv[:] = 0.0
                self.Ginv.flat[::len(self.S)+1] = z
            else:
                self.Ginv[:] = z
                self.Ginv *= self.S
            self.Ginv -= self.H
            for selfenergy in self.selfenergies:
                self.Ginv -= selfenergy(energy)
        if inverse:
            return self.Ginv
        else:
            return npy.linalg.inv(self.Ginv)
    def calculate(self, energy, sigma):
        ginv = energy * self.S - self.H - sigma
        return npy.linalg.inv(ginv)
    def set_hs(self, H, S):
        self.H = H
        self.S = S
    def dos(self, energy):
        """Total density of states -1/pi Im(Tr(GS))"""
        if self.S is None:
            return -self(energy).imag.trace() / npy.pi
        else:
            GS = npy.linalg.solve(self(energy, inverse=True), self.S)
            return -GS.imag.trace() / npy.pi
        
    def pdos(self, energy):
        """Projected density of states -1/pi Im(SGS/S)"""
        if self.S is None:
            return -self(energy).imag.diagonal() / npy.pi
        else:
            S = self.S
            SGS = npy.dot(S, npy.linalg.solve(self(energy, inverse=True), S))
            return -(SGS.diagonal() / S.diagonal()).imag / npy.pi
