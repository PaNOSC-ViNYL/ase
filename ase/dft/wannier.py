""" Maximally localized Wannier Functions

    Find the set of maximally localized Wannier functions
    using the spread functional of Marzari and Vanderbilt
    (PRB 56, 1997 page 12847). 
"""
import numpy as npy
from time import time
from math import sqrt, pi
import cPickle as pickle
from ase.parallel import paropen


def dag(a):
    """Return Hermitian conjugate of input"""
    return npy.conj(a.T)


def gram_schmidt(U, order=None):
    """Orthogonalize according to the Gram-Schmidt procedure.
    
    U is an NxM matrix containing M vectors as its columns.
    These will be orthogonalized using the order specified in the list 'order'

    Warning: The values of the original matrix is changed, return is just for
    convenience.
    """
    def project(a, b):
        """Return the projection of b onto a."""
        return a * (npy.dot(a.conj(), b) / npy.dot(a.conj(), a))

    if order is None:
        order = range(U.shape[1])

    for i, m in enumerate(order):
        col = U[:, m]
        for i2 in range(i):
            col -= project(U[:, order[i2]], col)
        col /= npy.sqrt(npy.dot(col.conj(), col))
    return U


def get_kpoint_dimensions(kpts):
    """Returns number of k-points along each axis of input Monkhorst pack.
    The set of k-points must not have been symmetry reduced.
    """
    nkpts = len(kpts)
    if nkpts == 1:
        return npy.ones(3, int)
    
    tol = 1e-5
    Nk_c = npy.zeros(3, int)
    for c in range(3):
        # Sort kpoints in ascending order along current axis
        slist = npy.argsort(kpts[:, c])
        skpts = npy.take(kpts, slist, axis=0)

        # Determine increment between kpoints along current axis
        DeltaK = max([skpts[n + 1, c] - skpts[n, c] for n in range(nkpts - 1)])

        # Determine number of kpoints as inverse of distance between kpoints
        if DeltaK > tol: Nk_c[c] = int(round(1. / DeltaK))
        else: Nk_c[c] = 1
    return Nk_c


def calculate_weights(cell_cc):
    """ Weights are used for non-cubic cells, see PRB **61**, 10040"""
    alldirs_dc = npy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                            [1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=int)
    g = npy.dot(cell_cc, cell_cc.T)
    # NOTE: Only first 3 of following 6 weights are presently used:
    w = npy.zeros(6)              
    w[0] = g[0, 0] - g[0, 1] - g[0, 2]
    w[1] = g[1, 1] - g[0, 1] - g[1, 2]
    w[2] = g[2, 2] - g[0, 2] - g[1, 2]
    w[3] = g[0, 1]
    w[4] = g[0, 2]
    w[5] = g[1, 2]
    # Make sure that first 3 Gdir vectors are included - 
    # these are used to calculate Wanniercenters.
    Gdir_dc = alldirs_dc[:3]
    weight_d = w[:3]
    for d in range(3, 6):
        if abs(w[d]) > 1e-5:
            Gdir_dc = npy.concatenate(Gdir_dc, alldirs_dc[d])
            weight_d = npy.concatenate(weight_d, w[d])
    weight_d /= max(abs(weight_d))
    return weight_d, Gdir_dc


def random_orthogonal_matrix(dim, seed=None, real=False):
    """ Generate a random orthogonal matrix"""
    if seed is not None:
        npy.random.seed(seed)

    H = npy.random.rand(dim, dim)
    npy.add(dag(H), H, H)
    npy.multiply(.5, H, H)

    if real:
        return gram_schmidt(H)
    else: 
        val, vec = npy.linalg.eig(H)
        return npy.dot(vec * npy.exp(1.j * val), dag(vec))


def steepest_descent(func, step=0.005, tolerance=1.0e-6, **kwargs):
    fvalueold = 0.
    fvalue = fvalueold + 10
    count=0
    while abs((fvalue - fvalueold) / fvalue) > tolerance:
        fvalueold = fvalue
        dF = func.get_gradients()
        func.step(dF * step, **kwargs)
        fvalue = func.get_functional_value()
        count += 1
        print 'SteepestDescent: iter=%s, value=%s' % (count, fvalue)


def md_min(func, step=0.25, tolerance=1.0e-6, verbose=True, **kwargs):
        t = - time()
        fvalueold = 0.
        fvalue = fvalueold + 10
        count = 0
        V = npy.zeros(func.get_gradients().shape, dtype=complex)
        while abs((fvalue - fvalueold) / fvalue) > tolerance:
            fvalueold = fvalue
            dF = func.get_gradients()
            V *= (dF * V.conj()).real > 0
	    V += step * dF
            func.step(V, **kwargs)
	    fvalue = func.get_functional_value()
            if fvalue < fvalueold:
		step *= 0.5
	    count += 1
            if verbose:
                print 'MDmin: iter=%s, step=%s, value=%s' % (count,step,fvalue)
        t += time()
        print '%d iterations in %0.2f seconds (%0.2f ms/iter), endstep = %s'%(
            count, t, t * 1000. / count, step)


class Wannier:
    def __init__(self,
                 numberofwannier,
                 calc,
                 numberofbands=None,
                 occupationenergy=0,
                 numberoffixedstates=None,
                 spin=0,
                 dtype=None):

        # Bloch phase sign convention
        sign = -1
        try:
            from Dacapo import Dacapo
        except ImportError:
            pass
        else:
            if isinstance(calc, Dacapo):
                sign = +1
            
        self.nwannier = numberofwannier
        self.spin = spin
        self.kpt_kc = sign * calc.get_ibz_k_points()
        assert len(calc.get_bz_k_points()) == len(self.kpt_kc)
        
        self.kptgrid = get_kpoint_dimensions(self.kpt_kc)
        self.Nk = len(self.kpt_kc)
        self.unitcell_cc = calc.get_atoms().get_cell()
        self.largeunitcell_cc = (self.unitcell_cc.T * self.kptgrid).T
        self.weight_d, self.Gdir_dc = calculate_weights(self.largeunitcell_cc)
        self.Ndir = len(self.weight_d) # Number of directions
        # Dacapo's initial wannier makes complex rotations even with 1 k-point
        if dtype is None:
            if sign == -1 and self.Nk == 1:
                dtype = float
            else:
                dtype = complex
        self.dtype = dtype
        
        if numberofbands is not None:
            self.nbands = numberofbands
        else:
            self.nbands = calc.get_number_of_bands()

        if numberoffixedstates is not None:
            self.fixedstates_k = numberoffixedstates
            self.edf_k = (npy.array(self.Nk * [numberofwannier])
                          - numberoffixedstates)
        else:
            # Setting number of fixed states and EDF from specified energy.
            # All states below this energy (relative to Fermi level) are fixed.
            occupationenergy += calc.get_fermi_level()
            fixedstates = []
            extradof = []
            for k in range(self.Nk):
                count = 0
                for eig in calc.get_eigenvalues(k, spin):
                    if eig < occupationenergy and count < self.nwannier:
                        count += 1
                fixedstates.append(count)
                extradof.append(self.nwannier - count)
            self.fixedstates_k = fixedstates
            self.edf_k = extradof
            print "Wannier: Fixed states            :", fixedstates
            print "Wannier: Extra degrees of freedom:", extradof

        # Set the list of neighboring k-points
        self.kklst_dk = npy.zeros((self.Ndir, self.Nk), dtype=int)
        for d in range(self.Ndir):
            for k in range(self.Nk):
                k1, k0_c = self.get_next_kpoint(d, k)
                self.kklst_dk[d, k] = k1

        # Set the inverse list of neighboring k-points
        self.invkklst_dk = npy.zeros((self.Ndir, self.Nk), dtype=int)
        for d in range(self.Ndir):
            for k1 in range(self.Nk):
                k = self.kklst_dk[d].tolist().index(k1)
                self.invkklst_dk[d, k1] = k


    def initialize(self, calc, first=True, initialwannier=None,
                   seed=None, bloch=False):
        Nw = self.nwannier
        Nb = self.nbands

        if first:
            self.Z_dkww = npy.empty((self.Ndir, self.Nk, Nw, Nw), self.dtype)
            self.Z_dknn = npy.empty((self.Ndir, self.Nk, Nb, Nb), self.dtype)
            self.V_knw = npy.zeros((self.Nk, Nb, Nw), self.dtype)
            for d in range(self.Ndir):
                for k in range(self.Nk):                
                    # get next kpoint K1 and reciprocal lattice vector K0
                    # given kpoint K that fulfills the criteria : K1-K+K0=Gdir
                    k1, k0_c = self.get_next_kpoint(d, k)
                    self.Z_dknn[d, k] = calc.get_wannier_localization_matrix(
                        nbands=Nb, dirG=self.Gdir_dc[d], kpoint=k,
                        nextkpoint=k1, G_I=k0_c, spin=self.spin)

        # Initialize rotation and coefficient matrices
        if bloch is True:
            # Set U and C to pick the lowest Bloch states
            self.U_kww = npy.zeros((self.Nk, Nw, Nw), self.dtype)
            self.C_kul = []
            for U, N, L in zip(self.U_kww, self.fixedstates_k, self.edf_k):
                U[:] = npy.identity(Nw, self.dtype)
                if L > 0:
                    self.C_kul.append(
                        npy.identity(Nb - N, self.dtype)[:, :L])
                else:
                    self.C_kul.append([])
        elif initialwannier is not None:
            # Use initial guess to determine U and C
            self.C_kul, self.U_kww = calc.initial_wannier(
                initialwannier, self.kptgrid, self.fixedstates_k,
                self.edf_k, self.spin)
        else:
            # Set U and C to random (orthogonal) matrices
            self.U_kww = npy.zeros((self.Nk, Nw, Nw), self.dtype)
            self.C_kul = []
            real = (self.dtype == float)
            for U, N, L in zip(self.U_kww, self.fixedstates_k, self.edf_k):
                U[:] = random_orthogonal_matrix(Nw, seed, real)
                if L > 0:
                    self.C_kul.append(random_orthogonal_matrix(
                        Nb - N, seed=seed, real=real)[:, :L])
                else:
                    self.C_kul.append(npy.array([]))        
        self.update()

    def get_k0(self, c): 
        """ find distance to next kpoint, along one of the three reciprocal
        lattice vectors """
        # For a gamma-point calculation
        if self.Nk == 1:
            return 0.0

        # make a sorted list of the kpoint values in this direction
        slist = npy.argsort(self.kpt_kc[:, c])
        skpoints_kc = npy.take(self.kpt_kc, slist, axis=0)

        # find k0
        k0 = max([skpoints_kc[n + 1, c] - skpoints_kc[n, c]
                    for n in range(self.Nk - 1)])
        return k0

    def get_next_kpoint(self, dir, k):
        """find the 'next' kpoint from kpt in the dir direction and the
        corresponding vector for the resulting pair. 

        This can also be written as K1-K-Gdir+K0=0
        Example: kpoints = (-0.375,-0.125,0.125,0.375), dir=0
        Gdir = [0.25,0,0]
        K=0.375, K1= -0.375 : -0.375-0.375-0.25 => K0=[1,0,0]

        For a gamma point calculation K1=K=0,  K0 = [1,0,0] for dir=0
        """
        tol = 1e-4

        if self.Nk == 1:
            return k, self.Gdir_dc[dir]

        # setup dist vector to next kpoint
        kdist_c = npy.zeros(3)
        for c in range(3):
            if self.Gdir_dc[dir, c] > 0:
                kdist_c[c] = self.get_k0(c)

        if max(kdist_c) < tol:
            return k, self.Gdir_dc[dir]

        # now search for K1-K-Gdir+K0=0
        # reciprocal lattice vectors we are going to search for,
        # including K0=0
        alldir_dc = npy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],
                               [1,1,0],[1,0,1],[0,1,1]], int)
        for k0_c in alldir_dc:
            for k1 in range(self.Nk):
                if npy.linalg.norm(
                    self.kpt_kc[k1] - self.kpt_kc[k] - kdist_c + k0_c) < tol:
                    return k1, k0_c

        print 'Wannier: Did not find matching kpoint for kpt=', self.kpt_kc[k]
        print 'Probably non-uniform k-point grid'
        raise NotImplementedError

    def update(self):
        # Update large rotation matrix V (from rotation U and coeff C)
        for k, N in enumerate(self.fixedstates_k):
            self.V_knw[k, :N] = self.U_kww[k,:N]
            if N < self.nwannier:
                self.V_knw[k, N:] = npy.dot(self.C_kul[k], self.U_kww[k, N:])
##             else:
##                 self.V_knw[k, N:] = 0.0

        # Calculate the Zk matrix from the large rotation matrix:
        # Zk = V^d[k] Zbloch V[k1]
        for d in range(self.Ndir):
            for k in range(self.Nk):
                k1 = self.kklst_dk[d, k]
                self.Z_dkww[d, k] = npy.dot(dag(self.V_knw[k]), npy.dot(
                    self.Z_dknn[d, k], self.V_knw[k1]))

        # Update the new Z matrix
        self.Z_dww = self.Z_dkww.sum(axis=1) / self.Nk

    def get_centers(self, scaled=False):
        """Calculate the Wannier centers

        ::
        
          pos =  L / 2pi * phase(diag(Z))
        """
        coord_wc = npy.angle(self.Z_dww[:3].diagonal(0, 1, 2)).T / (2 * pi)
        if not scaled:
            coord_wc = npy.dot(coord_wc, self.largeunitcell_cc)
        return coord_wc

    def get_radii(self):
        """Calculate the Wannier radii

        radius = sum abs(L/2pi ln abs diag(Z * Z^d))
        """
        return npy.dot(self.largeunitcell_cc.diagonal() / (2 * pi),
                       abs(npy.log(abs(self.Z_dww[:3].diagonal(0, 1, 2))**2)))
    
    def get_radii2(self):
        """Calculate the Wannier radii

        ::
          
                        --  /  L  \ 2       2
          radius**2 = - >   | --- |   ln |Z| 
                        --d \ 2pi /
        """
        return -npy.dot(self.largeunitcell_cc.diagonal()**2 / (2 * pi)**2,
                        npy.log(abs(self.Z_dww[:3].diagonal(0, 1, 2))**2))

    def get_spectral_weight_of_wannier_function(self, w):
        return abs(self.V_knw[:, :, w])**2 / self.Nk

    def get_wannier_function_dos(self, calc, n, energies, width):
        spec_kn = self.get_spectral_weight_of_wannier_function(n)
        dos = npy.zeros(len(energies))
        for k, spec_n in enumerate(spec_kn):
            eig_n = calc.get_eigenvalues(k=kpt, s=self.spin)
            for weight, eig in zip(spec_n, eig):
                dos += weight * self.get_gaussian(energies, eig, width)
        return dos

    def get_gaussian(self, energies, center, width):
        gauss = npy.zeros(len(energies))
        for i in range(len(energies)):
            exponent = min(20, ((energies[i] - center)/width)**2)
            gauss[i] = npy.exp(-exponent) / (sqrt(pi) * width)
        return gauss

    def get_max_spread(self, directions=[0,1,2]):
        """Returns the index of the most delocalized Wannier function
        together with the value of the spread functional"""
        d = npy.zeros(self.nwannier)
        for dir in directions:
            d[dir] = npy.abs(self.Z_dww[dir].diagonal())**2 *self.weight_d[dir]
        index = npy.argsort(d)[0]
        print 'Index:', index
        print 'Spread:', d[index]           

    def dump(self, file):
        pickle.dump((self.Z_dknn, self.U_kww, self.C_kul), paropen(file, 'w'))

    def load(self, file):
        Nw = self.nwannier
        Nb = self.nbands
        self.Z_dkww = npy.empty((self.Ndir, self.Nk, Nw, Nw), self.dtype)
        self.V_knw = npy.zeros((self.Nk, Nb, Nw), self.dtype)
        self.Z_dknn, self.U_kww, self.C_kul = pickle.load(paropen(file))
        self.update()

    def translate(self, w, R):
        """Translate the w'th Wannier function

        The distance vector R = [n1, n2, n3], is in units of the basis
        vectors of the small cell.
        """
        for kpt_c, U_ww in zip(self.kpt_kc, self.U_kww):
            U_ww[:, w] *= npy.exp(2.0j * pi * npy.dot(npy.array(R), kpt_c))
        self.update()

    def translate_to_cell(self, w, cell):
        """Translate the w'th Wannier function to specified cell"""
        scaled_c = npy.angle(self.Z_dww[:3, w, w]) * self.kptgrid / (2 *pi)
        trans = npy.array(cell) - npy.floor(scaled_c)
        self.translate(w, trans)

    def translate_all_to_cell(self, cell=[0, 0, 0]):
        """Translate all Wannier functions to specified cell"""
        scaled_wc = npy.angle(self.Z_dww[:3].diagonal(0, 1, 2)).T  * \
                    self.kptgrid / (2 * pi)
        trans_wc =  npy.array(cell)[None] - npy.floor(scaled_wc)
        for kpt_c, U_ww in zip(self.kpt_kc, self.U_kww):
            U_ww *= npy.exp(2.0j * pi * npy.dot(trans_wc, kpt_c))
        self.update()

    def get_hopping(self, R, calc):
        """Returns the matrix H(R)_nm=<0,n|H|R,m>.

        ::
        
                                1   _   -ik.R 
          H(R) = <0,n|H|R,m> = --- >_  e      H(k)
                                Nk  k         

        where R is the cell-distance (in units of the basis vectors of
        the small cell) and n,m are index of the Wannier functions."""

        H_ww = npy.zeros([self.nwannier, self.nwannier], self.dtype)
        for k, kpt_c in enumerate(self.kpt_kc):
            phase = npy.exp(-2.j * pi * npy.dot(npy.array(R), kpt_c))
            H_ww += self.get_hamiltonian(calc, k) * phase
        return H_ww / self.Nk

    def get_hamiltonian(self, calc, k=0):
        """Get Hamiltonian at existing k-vector of index k

        ::
        
                  dag
          H(k) = V    diag(eps )  V
                  k           k    k
        """
        eps_n = calc.get_eigenvalues(kpt=k, spin=self.spin)
        return npy.dot(dag(self.V_knw[k]) * eps_n, self.V_knw[k])

    def dist(self, R):
        Nw = self.Nwannier
        cen = self.get_center()
        r1 = npy.reshape(cen.repeat(Nwannier),[Nw, Nw, 3])
        r2 = cen.copy()
        for i in range(3):
            r2 += self.unitcell[i] * R[i]

        r2 = npy.swapaxes(npy.reshape(r2.repeat(Nw), [Nw, Nw, 3]), 0, 1)
        return npy.sqrt(npy.sum((r1-r2)**2, axis=-1))

    def get_hamiltonian_kpoint(self, kpt_c, calc):
        """Get Hamiltonian at some new arbitrary k-vector

        ::
        
                  _   ik.R 
          H(k) = >_  e     H(R)
                  R         
        """
        max = (self.kptgrid - 1) / 2
        max += max > 0
        N1, N2, N3 = max
        if not hasattr(self, hop_rww):
            self.R_rc = []
            self.hop_rww = []
            for n1 in xrange(-N1, N1 + 1):
                for n2 in xrange(-N2, N2 + 1):
                    for n3 in xrange(-N3, N3 + 1):
                        R = [n1, n2, n3]
                        if self.get_min_distance(R) < self.couplingradius: #XXX
                            self.R_rc.append(R)
                            self.hop_rww.append(self.get_hopping(R, calc))
                            
        Hk = npy.zeros([self.nwannier, self.nwannier], self.dtype)
        for R, hop_ww in zip(self.R_rc, self.hop_rww):
            filter = self.distances(R) < self.couplingradius # XXX
            phase = npy.exp(+2.j * pi * npy.dot(npy.array(R), kpt_c))
            Hk += hop_ww * filter * phase
        return npy.sort(npy.linalg.eigvals(Hk).real)

    def get_function(self, calc, index, repeat=None):
        """Index can be either a single WF or a coordinate vector
        in terms of the WFs."""

        # Default size of plotting cell is the one corresponding to k-points.
        if repeat is None:
            repeat = self.kptgrid
        N1, N2, N3 = repeat

        dim = calc.get_number_of_grid_points()
        largedim = dim * [N1, N2, N3]
        
        wanniergrid = npy.zeros(largedim, dtype=self.dtype)
        for k, kpt_c in enumerate(self.kpt_kc):
            # The coordinate vector of wannier functions
            if type(index) == int:
                vec_n = self.V_knw[k, :, index]
            else:   
                vec_n = npy.dot(self.V_knw[k], index)

            wan_G = npy.zeros(dim, self.dtype)
            for n, coeff in enumerate(vec_n):
                wan_G += coeff * calc.get_pseudo_wave_function(n, k, self.spin,
                                                               pad=True)

            # Distribute the small wavefunction over large cell:
            for n1 in xrange(N1):
                for n2 in xrange(N2):
                    for n3 in xrange(N3): # sign?
                        e = npy.exp(-2.j * pi * npy.dot([n1, n2, n3], kpt_c))
                        wanniergrid[n1 * dim[0]:(n1 + 1) * dim[0],
                                    n2 * dim[1]:(n2 + 1) * dim[1],
                                    n3 * dim[2]:(n3 + 1) * dim[2]] += e * wan_G

        # Normalization
        wanniergrid /= npy.sqrt(self.Nk)
        return wanniergrid

    def write_cube(self, calc, index, fname, repeat=None):
        from ase.io.cube import write_cube

        # Default size of plotting cell is the one corresponding to k-points.
        if repeat is None:
            repeat = self.kptgrid
        atoms = calc.get_atoms() * repeat
        write_cube(fname, atoms, data=self.get_function(calc, index, repeat))

    def localize(self, step=0.25, tolerance=1.0e-08, verbose=True,
                 updaterot=True, updatecoeff=True):
        print 'Localize with step =', step, 'and tolerance =', tolerance
##         maxi = steepest_descent(self, step, tolerance, updaterot=updaterot,
##                                 updatecoeff=updatecoeff)
        maxi = md_min(self, step, tolerance, verbose=verbose,
                      updaterot=updaterot, updatecoeff=updatecoeff)

    def get_functional_value(self): 
        """Calculate the value of the spread functional.

        ::

          Tr[|ZI|^2]=sum(I)sum(n) w_i|Z_(i)_nn|^2,

        where w_i are weights."""
        a_d = npy.sum(npy.abs(self.Z_dww.diagonal(0, 1, 2))**2, axis=1)
        return npy.dot(a_d, self.weight_d).real

    def get_gradients(self):
        """Determine gradient of the spread functional

        The gradient for a rotation A_kij is::
        
           dU = dRho/dA_{k,i,j} = sum(I) sum(k')
                   + Z_jj Z_kk',ij^* - Z_ii Z_k'k,ij^*
                   - Z_ii^* Z_kk',ji + Z_jj^* Z_k'k,ji
    
        The gradient for a change of coefficients is::
        
          dRho/da^*_{k,i,j} = sum(I) [[(Z_0)_{k} V_{k'} diag(Z^*) +
                                       (Z_0_{k''})^d V_{k''} diag(Z)] *
                                       U_k^d]_{N+i,N+j}

        where diag(Z) is a square,diagonal matrix with Z_nn in the diagonal, 
        k' = k + dk and k = k'' + dk.

        The extra degrees of freedom chould be kept orthonormal to the fixed
        space, thus we introduce lagrange multipliers, and minimize instead::

            Rho_L=Rho- sum_{k,n,m} lambda_{k,nm} <c_{kn}|c_{km}>

        for this reason the coefficient gradients should be multiplied
        by (1 - c c^d).
        """
        Nb = self.nbands
        Nw = self.nwannier
        
        dU = []
        dC = []
        for k in xrange(self.Nk):
            M = self.fixedstates_k[k]
            L = self.edf_k[k]
            U_ww = self.U_kww[k]
            C_ul = self.C_kul[k]
            Utemp_ww = npy.zeros([Nw, Nw], self.dtype)
            Ctemp_nw = npy.zeros((Nb, Nw), self.dtype)

            for d, weight in enumerate(self.weight_d):
                if abs(weight) < 1.0e-6:
                    continue

                Z_knn = self.Z_dknn[d]
                diagZ_w = self.Z_dww[d].diagonal()
                Zii_ww = npy.repeat(diagZ_w, Nw).reshape(Nw, Nw)
                k1 = self.kklst_dk[d, k]
                k2 = self.invkklst_dk[d, k]
                V_knw = self.V_knw
                Z_kww = self.Z_dkww[d]
                
                if L > 0:
                    Ctemp_nw += weight * npy.dot(
                        npy.dot(Z_knn[k], V_knw[k1]) * diagZ_w.conj() +
                        npy.dot(dag(Z_knn[k2]), V_knw[k2]) * diagZ_w,
                        dag(U_ww))

                temp = Zii_ww.T * Z_kww[k].conj() - Zii_ww * Z_kww[k2].conj()
                Utemp_ww += weight * (temp - dag(temp))
            dU.append(Utemp_ww.ravel())
            if L > 0:
                # Ctemp now has same dimension as V, the gradient is in the
                # lower-right (Nb-M) x L block
                Ctemp_ul = Ctemp_nw[M:, M:]
                G_ul = Ctemp_ul - npy.dot(npy.dot(C_ul, dag(C_ul)), Ctemp_ul)
                dC.append(G_ul.ravel())

        return npy.concatenate(dU + dC)
                        
    def step(self, dX, updaterot=True, updatecoeff=True):
        """ dX is (A, dC) where U->Uexp(-A) and C->C+dC """
        Nb = self.nbands
        Nw = self.nwannier
        Nk = self.Nk
        N_k = self.fixedstates_k
        L_k = self.edf_k
        if updaterot:
            A_kww = npy.reshape(dX[:Nk * Nw**2], (Nk, Nw, Nw))
            for U, A in zip(self.U_kww, A_kww):
                H = npy.conj(A * 1.j).copy()
                epsilon, Z = npy.linalg.eigh(H)
                # Z contains the eigenvectors as COLUMNS.
                # Since H=iA, dU = exp(-A)=exp(iH) = ZDZ^d
                dU = npy.dot(Z * npy.exp(1.j * epsilon), dag(Z))
                U[:] = npy.dot(U, dU)

        if updatecoeff:
            start = 0
            for C, N, L in zip(self.C_kul, N_k, L_k):
                Ncoeff = L * (Nb - N)
                deltaC = dX[Nk * Nw**2 + start: Nk * Nw**2 + start + Ncoeff]
                C[:] = gram_schmidt(C + deltaC.reshape(Nb - N, L))
                start += Ncoeff

        self.update()
