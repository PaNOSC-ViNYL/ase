import numpy as np
import os
from ase.atoms import Atoms
from ase.parallel import rank

class BEEF_Ensemble:
    """BEEF type ensemble error estimation"""
    def __init__(self, atoms=None, e=None, contribs=None, xc=None):
        if (atoms is not None or contribs is not None or xc is not None):
            if atoms is None:
                assert e is not None
                assert contribs is not None
                assert xc is not None
            else:
                if isinstance(atoms, Atoms):
                    calc = atoms.get_calculator()
                else:
                    calc = atoms
                xc = calc.get_xc_functional()
            self.calc = calc
            self.e = e
            self.contribs = contribs
            self.xc = xc
            self.done = False
            if self.xc in ['BEEF-vdW', 'BEEF-1']:
                self.beef_type = 'beefvdw'
            elif self.xc == 'mBEEF':
                self.beef_type = 'mbeef'
            else:
                raise NotImplementedError('No ensemble for xc = %s' % self.xc)

    def get_ensemble_energies(self, size=2000, seed=0, write=False, fname=None):
        """Returns an array of ensemble total energies"""
        if write:
            isinstance(fname, str)
        if rank == 0:
            print '\n'
            print '%s ensemble started' % self.xc

        if self.contribs is None:
            self.contribs = self.calc.get_nonselfconsistent_energies(self.beef_type)
            self.e = self.calc.get_potential_energy()
        if self.beef_type == 'beefvdw':
            assert len(self.contribs) == 32
            coefs = self.get_beefvdw_ensemble_coefs(size, seed)
        elif self.beef_type == 'mbeef':
            assert len(self.contribs) == 64
            coefs = self.get_mbeef_ensemble_coefs(size, seed)
        de = np.dot(coefs, self.contribs)
        self.done = True

        if rank == 0:
            print '%s ensemble finished' % self.xc
            print '\n'

        if write:
            self.write(fname, de, seed)
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
        """Pertubation coefficients of the mBEEF ensemble"""
        from pars_mbeef import uiOmega as omega
        assert np.shape(omega) == (64, 64)

        Wo, Vo = np.linalg.eig(omega)
        np.random.seed(seed)
        mu, sigma = 0.0, 1.0
        rand = np.array(np.random.normal(mu, sigma, (len(Wo), size)))
        return (np.sqrt(2.)*np.dot(np.dot(Vo, np.diag(np.sqrt(Wo))), rand)[:]).T

    def write(self, fname, de, seed):
        """Write ensemble data file"""
        import cPickle as pickle
        if fname[-4:] != '.ens':
            fname += '.ens'
        assert self.done
        if rank is 0:
            if os.path.isfile(fname):
                os.rename(fname, fname + '.old')
            f = open(fname, 'w')
            obj = [self.e, de, self.contribs, seed, self.xc]
            pickle.dump(obj, f)
            f.close()

    def read(self, fname, all=False):
        import cPickle as pickle
        isinstance(fname, str)
        if fname[-4:] != '.ens':
            fname += '.ens'
        assert os.path.isfile(fname)
        f = open(fname, 'r')
        e, de, contribs, seed, xc = pickle.load(f)
        f.close()
        if all:
            return e, de, contribs, seed, xc
        else:
            return e, de
