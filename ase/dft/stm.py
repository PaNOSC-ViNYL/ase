import numpy as np


class STM:
    def __init__(self, atoms, symmetries=None):
        calc = atoms.get_calculator()
        self.nbands = calc.get_number_of_bands()
        self.weights = calc.get_k_point_weights()
        self.nkpts = len(self.weights)
        self.nspins = calc.get_number_of_spins()
        self.eigs = np.array([[calc.get_eigenvalues(k, s)
                               for k in range(self.nkpts)]
                              for s in range(self.nspins)])
        self.eigs -= calc.get_fermi_level()
        self.calc = calc
        self.cell = atoms.get_cell()
        assert not self.cell[2, :2].any() and not self.cell[:2, 2].any()
        self.ldos = None
        self.bias = None
        self.symmetries = symmetries or []
                               
    def calculate_ldos(self, bias):
        if self.ldos is not None and bias == self.bias:
            return

        if bias < 0:
            emin = bias
            emax = 0.0
        else:
            emin = 0
            emax = bias

        ldos = 0.0
        for s in range(self.nspins):
            for k in range(self.nkpts):
                for n in range(self.nbands):
                    e = self.eigs[s, k, n]
                    if emin < e < emax:
                        psi = self.calc.get_pseudo_wave_function(n, k, s)
                        ldos += self.weights[k] * (psi * np.conj(psi)).real

        if 0 in self.symmetries:
            # (x,y) -> (-x,y)
            ldos[1:] += ldos[:0:-1].copy()
            ldos[1:] *= 0.5

        if 1 in self.symmetries:
            # (x,y) -> (x,-y)
            ldos[:, 1:] += ldos[:, :0:-1].copy()
            ldos[:, 1:] *= 0.5
            
        if 2 in self.symmetries:
            # (x,y) -> (y,x)
            ldos += ldos.transpose((1, 0, 2)).copy()
            ldos *= 0.5
            
        self.ldos = ldos
        self.bias = bias

    #def save_ldos(self, filename='ldos.pckl'):
        
    def get_averaged_current(self, bias, z):
        self.calculate_ldos(bias)
        nz = self.ldos.shape[2]

        # Find grid point:
        n = z / self.cell[2, 2] * nz
        dn = n - np.floor(n)
        n = int(n) % nz

        # Average and do linear interpolation:
        return ((1 - dn) * self.ldos[:, :, n].mean() +
                dn * self.ldos[:, :, (n + 1) % nz].mean())
    
    def scan(self, bias, current):
        self.calculate_ldos(bias)

        L = self.cell[2, 2]
        nz = self.ldos.shape[2]
        h = L / nz

        ldos = self.ldos.reshape((-1, nz))

        heights = np.empty(ldos.shape[0])
        for i, a in enumerate(ldos):
            heights[i] = find_height(a, current, h)

        heights.shape = self.ldos.shape[:2]
        return heights
    
    def linescan(self, bias, current, p1, p2, npoints=50):
        heights = self.scan(bias, current)

        p1 = np.asarray(p1)
        p2 = np.asarray(p2)
        d = p2 - p1
        s = np.dot(d, d)**0.5

        cell = self.cell[:2, :2]
        shape = np.array(self.ldos.shape[:2], float)
        M = np.linalg.inv(cell)
        line = np.empty(npoints)
        for i in range(npoints):
            p = p1 + i * d / (npoints - 1)
            q = np.dot(p, M) * shape
            qi = q.astype(int)
            n0, n1 = qi
            f = q - qi
            g = 1 - f
            z = (g[0] * g[1] * heights[n0, n1] +
                 f[0] * g[1] * heights[n0 + 1, n1] +
                 g[0] * f[1] * heights[n0, n1 + 1] +
                 f[0] * f[1] * heights[n0 + 1, n1 + 1])
            line[i] = z
        return np.linspace(0, s, npoints), line

    def cube(self, filename, atoms=None):
        pass


def find_height(ldos, current, h):
    n = len(ldos) - 1
    while n >= 0:
        if ldos[n] > current:
            break
        n -= 1
    else:
        raise RuntimeError('Tip crash!')

    c2, c1 = ldos[n:n + 2]
    return (n + 1 - (current - c1) / (c2 - c1)) * h
