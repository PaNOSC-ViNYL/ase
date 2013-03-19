import numpy as np
from ase.atoms import Atoms
from ase.parallel import rank

class BEEF_Ensemble:
    """BEEF type ensemble error estimation"""
    def __init__(self, atoms=None, contribs=None, xc=None):
        if atoms is None:
            assert contribs is not None
            assert xc is not None
        else:
            if isinstance(atoms, Atoms):
                calc = atoms.get_calculator()
            else:
                calc = atoms
            xc = calc.get_xc_functional()
        self.calc = calc
        self.contribs = contribs
        self.xc = xc
        if self.xc in ['BEEF-vdW', 'BEEF-1']:
            self.beef_type = 'beefvdw'
        elif self.xc == 'mBEEF':
            self.beef_type = 'mbeef'
        else:
            raise NotImplementedError('No ensemble for xc = %s' % self.xc)

    def get_ensemble_energies(self, size=2000, seed=0):
        """Returns an array of ensemble total energies"""
        if rank == 0:
            print '\n'
            print '%s ensemble started' % self.xc

        if self.contribs is None:
            self.contribs = self.calc.get_nonselfconsistent_energies(self.beef_type)
        if self.beef_type == 'beefvdw':
            assert len(self.contribs) == 32
            coefs = self.get_beefvdw_ensemble_coefs(size, seed)
        elif self.beef_type == 'mbeef':
            assert len(self.contribs) == 64
            coefs = self.get_mbeef_ensemble_coefs(size, seed)
        de = np.dot(coefs, self.contribs)

        if rank == 0:
            print '%s ensemble finished' % self.xc
            print '\n'
        return de

    def get_beefvdw_ensemble_coefs(self, size, seed):
        """Pertubation coefficients of the BEEF-vdW ensemble"""
        from pars_beefvdw import uiOmega as omega
        assert np.shape(omega) == (31, 31)

        Wo, Vo = np.linalg.eig(omega)
        np.random.seed(seed)
        RandV = np.random.randn(31, size)

        for j in range(size):
            v = RandV[:,j]
            coefs_i = (np.dot(np.dot(Vo, np.diag(np.sqrt(Wo))), v)[:])
            if j == 0:
                ensemble_coefs = coefs_i
            else:
                ensemble_coefs = np.vstack((ensemble_coefs, coefs_i))
        PBEc_ens = -ensemble_coefs[:, 30]
        return (np.vstack((ensemble_coefs.T, PBEc_ens))).T

    def get_mbeef_ensemble_coefs(self, size, seed):
        """Pertubation coefficients of mBEEF ensemble"""
        from pars_mbeef import uiOmega as omega
        assert np.shape(omega) == (64, 64)

        mu, sigma = 0, 1.0 # mean zero and standard deviation one
        Wo, Vo = np.linalg.eig(omega)
        Wo = Wo.real
        for i, Woi in enumerate(Wo):
            if np.abs(Woi) < 1.e-15 or Woi < 0.0:
                Wo[i] = 0.

        np.random.seed(seed)
        rand = np.array(np.random.normal(mu, sigma, (len(Wo), size)))
        coefs = (np.sqrt(2.)*np.dot(np.dot(Vo, np.diag(np.sqrt(Wo))), rand)[:]).T
        return coefs
