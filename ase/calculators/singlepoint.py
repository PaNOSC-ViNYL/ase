import numpy as np

from ase.calculators.calculator import Calculator, all_properties


class SinglePointCalculator(Calculator):
    """Special calculator for a single configuration.

    Used to remember the energy, force and stress for a given
    configuration.  If the positions, atomic numbers, unit cell, or
    boundary conditions are changed, then asking for
    energy/forces/stress will raise an exception."""
    
    name = 'unknown'
    
    def __init__(self, *args, **results):
        """Save energy, forces, stress, ... for the current configuration."""
        if args and isinstance(args[0], float):
            # Old interface:
            assert not results
            for key, value in zip(['energy', 'forces', 'stress', 'magmoms'],
                                  args):
                if value is not None:
                    results[key] = value
            atoms = args[-1]
        else:
            if args:
                atoms = args[0]
            else:
                atoms = results.pop('atoms')
            
        Calculator.__init__(self)
        self.results = {}
        for property, value in results.items():
            assert property in all_properties
            if value is None:
                continue
            if property in ['energy', 'magmom']:
                self.results[property] = value
            else:
                self.results[property] = np.array(value, float)
        self.atoms = atoms.copy()

    def get_property(self, name, atoms=None, allow_calculation=True):
        if name not in self.results or self.check_state(atoms):
            if allow_calculation:
                raise NotImplementedError(
                    'The property "{0}" is not available.'.format(name))
            return None

        result = self.results[name]
        if isinstance(result, np.ndarray):
            result = result.copy()
        return result

    
class SinglePointKPoint:
    def __init__(self, weight, s, k, eps_n=[], f_n=[]):
        self.weight = weight
        self.s = s  # spin index
        self.k = k  # k-point index
        self.eps_n = eps_n
        self.f_n = f_n


class SinglePointDFTCalculator(SinglePointCalculator):
    def __init__(self, *args, **results):
        self.bz_kpts = results.pop('bz_kpts', None)
        self.ibz_kpts = results.pop('ibz_kpts', None)
        if args and isinstance(args[0], float):
            # Old interface:
            assert not results
            for key, value in zip(['energy', 'forces', 'stress', 'magmoms'],
                                  args):
                if value is not None:
                    results[key] = value
            atoms = args[4]
            if len(args) > 5:
                eFermi = args[5]
                if len(args) > 6:
                    energies = args[6]
        else:
            if args:
                atoms = args[0]
            else:
                atoms = results.pop('atoms')
            eFermi = results.pop('eFermi', None)
            energies = results.pop('energies', None)
            eref = results.pop('Eref', None)

        SinglePointCalculator.__init__(self, atoms, **results)

        if eFermi is not None:
            self.eFermi = eFermi
        if energies is not None:
            self.energies = energies
        if eref is not None:
            self.eref = eref
        self.kpts = None

    def get_fermi_level(self):
        """Return the Fermi-level(s)."""
        return self.eFermi

    def get_bz_k_points(self):
        """Return the k-points."""
        return self.bz_kpts

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return len(self.kpts)
        return None

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        nos = self.get_number_of_spins()
        if nos is not None:
            return nos == 2
        return None
    
    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone."""
        return self.ibz_kpts

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.f_n
        return None

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.eps_n
        return None

    def get_homo_lumo(self):
        """Return HOMO and LUMO energies."""
        if self.kpts is None:
            raise RuntimeError('No kpts')
        eHs = []
        eLs = []
        for kpt in self.kpts:
            eH, eL = self.get_homo_lumo_by_spin(kpt.s)
            eHs.append(eH)
            eLs.append(eL)
        return np.array(eHs).max(), np.array(eLs).min()
        
    def get_homo_lumo_by_spin(self, spin=0):
        """Return HOMO and LUMO energies for a give spin."""
        if self.kpts is None:
            raise RuntimeError('No kpts')
        for kpt in self.kpts:
            if kpt.s == spin:
                eH = -1.e32
                eL = 1.e32
                for e, f in zip(kpt.eps_n, kpt.f_n):
                    if e <= self.eFermi:
                        eH = max(eH, e)
                    else:
                        eL = min(eL, e)
                return eH, eL
        raise RuntimeError('No kpt with spin {0}'.format(spin))
