"""
Implementation of the Precon abstract base class and subclasses
"""
import time

import numpy as np

from ase.constraints import Filter, UnitCellFilter, FixAtoms, full_3x3_to_Voigt_6_strain
from ase.utils import sum128, dot128, norm128
from ase.geometry import undo_pbc_jumps

import ase.units as units
THz = 1e12*1./units.s

from ase.optimize.precon import (get_neighbours, estimate_nearest_neighbour_distance,
                                 logger, have_matscipy)

try:
    from scipy import sparse, rand
    from scipy.linalg import eigh, solve
    from scipy.optimize import fmin_l_bfgs_b
    from scipy.sparse.linalg import spsolve
    have_scipy = True
except ImportError:
    have_scipy = False

try:
    from pyamg import smoothed_aggregation_solver
    have_pyamg = True
except ImportError:
    have_pyamg = False

class Precon(object):
    def __init__(self, r_cut=None, r_NN=None,
                 mu=None, mu_c=None,
                 dim=3, c_stab=0.1, force_stab=False,
                 sparse=True,
                 recalc_mu=False, array_convention="C",
                 use_pyamg=True, solve_tol=1e-8,
                 apply_positions=True, apply_cell=True):

        """Initialise a preconditioner object based on passed parameters.

        Args:
            r_cut: float. This is a cut-off radius. The preconditioner matrix
                will be created by considering pairs of atoms that are within a
                distance r_cut of each other. For a regular lattice, this is
                usually taken somewhere between the first- and second-nearest
                neighbour distance. If r_cut is not provided, default is
                2 * r_NN (see below)
            r_NN: nearest neighbour distance. If not provided, this is calculated
                from input structure.
            mu: float
                energy scale for position degreees of freedom. If `None`, mu
                is precomputed using finite difference derivatives.
            mu_c: float
                energy scale for cell degreees of freedom. Also precomputed if None.
            dim: int; dimensions of the problem
            c_stab: float. The diagonal of the preconditioner matrix will have
                a stabilisation constant added, which will be the value of
                c_stab times mu.
            force_stab:
                If True, always add the stabilisation to diagnonal, regardless
                of the presence of fixed atoms.
            sparse: True if the preconditioner matrix should be returned as a
                scipy.sparse matrix (rather than a dense numpy array)
            recalc_mu: if True, the value of mu will be recalculated every time
                self.make_precon is called. This can be overridden in specific
                cases with recalc_mu argument in self.make_precon. If recalc_mu
                is set to True here, the value passed for mu will be
                irrelevant unless recalc_mu is set False the first time
                make_precon is called.
            array_convention: Either "C" or "F" for Fortran; this will change
                the preconditioner to reflect the ordering of the indices in
                the vector it will operate on. The C convention assumes the
                vector will be arranged atom-by-atom (ie [x1, y1, z1, x2, ...])
                while the F convention assumes it will be arranged component
                by component (ie [x1, x2, ..., y1, y2, ...]).
            use_pyamg: use PyAMG to solve P x = y, if available.
            solve_tol: tolerance used for PyAMG sparse linear solver, if available.
            apply_positions: if True, apply preconditioner to position DoF
            apply_cell: if True, apply preconditioner to cell DoF

        Raises:
            ValueError for problem with arguments

        """

        self.r_NN = r_NN
        self.r_cut = r_cut
        self.mu = mu
        self.mu_c = mu_c
        self.c_stab = c_stab
        self.force_stab = force_stab
        self.sparse = sparse
        self.array_convention = array_convention
        self.recalc_mu = recalc_mu
        self.P = None

        global have_pyamg
        if use_pyamg and not have_pyamg:
            use_pyamg = False
            logger.warning('use_pyamg=True but PyAMG cannot be imported!'
                            'falling back on direct inversion of preconditioner, '
                            'may be slow for large systems')
            print('use_pyamg=True but PyAMG cannot be imported!'
                   'falling back on direct inversion of preconditioner, '
                   'may be slow for large systems')

        self.use_pyamg = use_pyamg
        self.solve_tol = solve_tol
        self.apply_positions = apply_positions
        self.apply_cell = apply_cell

        if dim < 1:
            raise ValueError("Dimension must be at least 1")
        self.dim=dim

        global have_matscipy
        if not have_matscipy:
            logger.warning("Unable to import Matscipy. Neighbour list "
                           "calculations may be very slow.")
            print("WARNING: Unable to import Matscipy. Reverting to the "
                  "neighbourlist module in the ASE; this may make "
                  "preconditioning very slow.")

        global have_scipy
        if not have_scipy:
            logger.warning("Unable to import scipy. Sparse option will not be "
                           "available.")
            if self.sparse:
                print("WARNING: Unable to import scipy. Preconditioner will "
                      "be unable to use sparse matrices. For large problems "
                      "the program will run much slower or may run out of "
                      "memory altogether.")
            self.sparse = False

    def make_precon(self, atoms, recalc_mu=None):
        """Create a preconditioner matrix based on the passed set of atoms.

        Creates a general-purpose preconditioner for use with optimization
        algorithms, based on examining distances between pairs of atoms in the
        lattice. The matrix will be stored in the attribute self.P and
        returned.

        Args:
            atoms: the Atoms object used to create the preconditioner.
                Can also
            recalc_mu: if True, self.mu (and self.mu_c for variable cell)
                will be recalculated by calling self.estimate_mu(atoms)
                before the preconditioner matrix is created. If False, self.mu
                will be calculated only if it does not currently have a value
                (ie, the first time this function is called).

        Returns:
            A two-element tuple:
                P: A d*N x d*N matrix (where N is the number of atoms, d is the
                    value of self.dim), which will either be a dense numpy
                    matrix or a sparse scipy csr_matrix, depending on how the
                    precon object's 'sparse' flag is set. BE AWARE that using
                    numpy.dot() with sparse matrices will result in
                    errors/incorrect results - use the .dot method directly
                    on the matrix instead.
        """

        if self.r_NN is None:
            self.r_NN = estimate_nearest_neighbour_distance(atoms)

        if self.r_cut is None:
            # This is the first time this function has been called, and no
            # cutoff radius has been specified, so calculate it automatically.
            self.r_cut = 2.0 * self.r_NN
        elif self.r_cut < self.r_NN:
            warning = ('WARNING: r_cut (%.2f) < r_NN (%.2f), '
                       'increasing to 1.1*r_NN = %.2f' % (self.r_cut, self.r_NN,
                                                          1.1*self.r_NN))
            logger.info(warning)
            print(warning)
            self.r_cut = 1.1*self.r_NN

        if recalc_mu is None:
            # The caller has not specified whether or not to recalculate mu,
            # so the Precon's setting is used.
            recalc_mu = self.recalc_mu

        if self.mu is None:
            # Regardless of what the caller has specified, if we don't
            # currently have a value of mu, then we need one.
            recalc_mu = True

        if recalc_mu:
            self.estimate_mu(atoms)

        positions = atoms.get_positions()
        if self.P is not None:
            real_atoms = atoms
            if isinstance(atoms, Filter):
                real_atoms = atoms.atoms
            displacement = undo_pbc_jumps(real_atoms)
            max_abs_displacement = abs(displacement).max()
            logger.info('max(abs(displacements)) = %.2f A (%.2f r_NN)',
                max_abs_displacement, max_abs_displacement/self.r_NN)
            if max_abs_displacement < 0.5*self.r_NN:
                return self.P

        start_time = time.time()

        # Create the preconditioner:
        if self.sparse:
            self._make_sparse_precon(atoms, force_stab=self.force_stab)
        else:
            self._make_dense_precon(atoms, force_stab=self.force_stab)

        logger.info("--- Precon created in %s seconds ---",
                    time.time()-start_time)
        return self.P

    def _make_dense_precon(self, atoms, initial_assembly=False, force_stab=False):
        """Create a dense preconditioner matrix based on the passed atoms.

        Creates a general-purpose preconditioner for use with optimization
        algorithms, based on examining distances between pairs of atoms in the
        lattice. The matrix will be stored in the attribute self.P and
        returned. Note that this function will use self.mu, whatever it is.

        Args:
            atoms: the Atoms object used to create the preconditioner.

        Returns:
            A two-dimensional numpy array which is a d*N by d*N matrix (where
            N is the number of atoms, and d is the value of self.dim).

        """

        N = len(atoms)

        # csc_P is smaller than P, and holds the coefficients which will
        # eventually make up self.P, the complete preconditioning matrix.

        if self.apply_positions:
            csc_P = np.zeros((N, N))
            i_list, j_list, rij_list, fixed_atoms = get_neighbours(atoms, self.r_cut)

            if force_stab or len(fixed_atoms) == 0:
                for i in range(N):
                    csc_P[i,i] += self.mu * self.c_stab
        else:
            csc_P = np.eye(N)

        # P_ii is mu_c for cell DoF
        if isinstance(atoms, Filter):
            if self.apply_cell:
                csc_P[N-3,N-3] = self.mu_c
                csc_P[N-2,N-2] = self.mu_c
                csc_P[N-1,N-1] = self.mu_c
            else:
                csc_P[N-3,N-3] = 1.0
                csc_P[N-2,N-2] = 1.0
                csc_P[N-1,N-1] = 1.0

        if self.apply_positions:
            for i, j, rij in zip(i_list, j_list, rij_list):
                if i == j:
                    continue
                Cij = self.get_coeff(rij)
                csc_P[i,i] -= Cij
                csc_P[i,j] = Cij

            if not initial_assembly:
                for i in fixed_atoms:
                    csc_P[:,i] = 0.0
                    csc_P[i,:] = 0.0
                    csc_P[i,i] = 1.0

        self.csc_P = csc_P
        self.P = np.zeros((self.dim*N, self.dim*N))

        if self.dim == 1:
            self.P = csc_P

        # Now create the complete self.P matrix based on csc_P.
        elif self.array_convention == "F":
            # Build a matrix using csc_P as blocks.
            self.P = np.zeros((self.dim*N, self.dim*N))
            for i in range(self.dim):
                j = i*N
                k = (i+1)*N
                self.P[j:k, j:k] = csc_P
        else:
            # array_convention assumed to be "C". Here each element in csc_P is
            # 'expanded' into a scalar matrix.
            for i, j in itertools.product(range(N), range(N)):
                k = self.dim*i
                l = self.dim*j
                np.fill_diagonal(self.P[k:k+self.dim, l:l+self.dim],
                                 csc_P[i, j])

        return self.P

    def _make_sparse_precon(self, atoms, initial_assembly=False, force_stab=False):
        """Create a sparse preconditioner matrix based on the passed atoms.

        Creates a general-purpose preconditioner for use with optimization
        algorithms, based on examining distances between pairs of atoms in the
        lattice. The matrix will be stored in the attribute self.P and
        returned. Note that this function will use self.mu, whatever it is.

        Args:
            atoms: the Atoms object used to create the preconditioner.

        Returns:
            A scipy.sparse.csr_matrix object, representing a d*N by d*N matrix
            (where N is the number of atoms, and d is the value of self.dim).
            BE AWARE that using numpy.dot() with this object will result in
            errors/incorrect results - use the .dot method directly on the
            sparse matrix instead.

        """
        logger.info('creating sparse precon: initial_assembly=%r, force_stab=%r, apply_positions=%r, apply_cell=%r',
                    initial_assembly, force_stab, self.apply_positions, self.apply_cell)

        N = len(atoms)
        diag_i = np.arange(N, dtype=int)
        start_time = time.time()
        if self.apply_positions:
            # compute neighbour list
            i, j, rij, fixed_atoms = get_neighbours(atoms, self.r_cut)
            logger.info('--- neighbour list created in %s s ---' % (time.time() - start_time))

            # compute entries in triplet format: without the constraints
            start_time = time.time()
            coeff = self.get_coeff(rij)
            diag_coeff = np.bincount(i, -coeff, minlength=N).astype(np.float64)
            if force_stab or len(fixed_atoms) == 0:
                logger.info('adding stabilisation to preconditioner')
                diag_coeff += self.mu*self.c_stab
        else:
            diag_coeff = np.ones(N)

        # precon is mu_c*identity for cell DoF
        if isinstance(atoms, Filter):
            if self.apply_cell:
                diag_coeff[-3] = self.mu_c
                diag_coeff[-2] = self.mu_c
                diag_coeff[-1] = self.mu_c
            else:
                diag_coeff[-3] = 1.0
                diag_coeff[-2] = 1.0
                diag_coeff[-1] = 1.0
        logger.info('--- computed triplet format in %s s ---' % (time.time() - start_time))

        if self.apply_positions and not initial_assembly:
            # apply the constraints
            start_time = time.time()
            mask = np.ones(N)
            mask[fixed_atoms] = 0.0
            coeff *= mask[i] * mask[j]
            diag_coeff[fixed_atoms] = 1.0
            logger.info('--- applied fixed_atoms in %s s ---' % (time.time() - start_time))

        if self.apply_positions:
            # remove zeros
            start_time = time.time()
            inz = np.nonzero(coeff)
            i = np.hstack((i[inz], diag_i))
            j = np.hstack((j[inz], diag_i))
            coeff = np.hstack((coeff[inz], diag_coeff))
            logger.info('--- remove zeros in %s s ---' % (time.time() - start_time))
        else:
            i = diag_i
            j = diag_i
            coeff = diag_coeff

        # create the matrix
        start_time = time.time()
        csc_P = sparse.csc_matrix((coeff, (i,j)), shape=(N,N))
        logger.info('--- created CSC matrix in %s s ---' % (time.time() - start_time))

        self.csc_P = csc_P

        start_time = time.time()
        if self.dim == 1:
            self.P = csc_P
        elif self.array_convention == "F":
            csc_P = csc_P.tocsr()
            self.P = csc_P
            for i in range(self.dim - 1):
                self.P = sparse.block_diag((self.P, csc_P)).tocsr()
        else:
            # convert back to triplet and read the arrays
            csc_P = csc_P.tocoo()
            i = csc_P.row * self.dim
            j = csc_P.col * self.dim
            z = csc_P.data

            # N-dimensionalise, interlaced coordinates
            I = np.hstack([i+d for d in range(self.dim)])
            J = np.hstack([j+d for d in range(self.dim)])
            Z = np.hstack([z for d in range(self.dim)])
            self.P = sparse.csc_matrix((Z, (I,J)),
                    shape=(self.dim*N,self.dim*N))
            self.P = self.P.tocsr()
        logger.info('--- N-dim precon created in %s s ---' % (time.time() - start_time))

        # Create solver
        if self.use_pyamg and have_pyamg:
            start_time = time.time()
            self.ml = smoothed_aggregation_solver(self.P, B=None,
                    strength=('symmetric', {'theta': 0.0}),
                    smooth=('jacobi', {'filter': True, 'weighting': 'local'}),
                    improve_candidates=[('block_gauss_seidel', {'sweep': 'symmetric', 'iterations': 4}),
                    None, None, None, None, None, None, None, None, None, None, None, None, None, None],
                    aggregate="standard",
                    presmoother=('block_gauss_seidel', {'sweep': 'symmetric', 'iterations': 1}),
                    postsmoother=('block_gauss_seidel', {'sweep': 'symmetric', 'iterations': 1}),
                    max_levels=15,
                    max_coarse=300,
                    coarse_solver="pinv")
            logger.info('--- multi grid solver created in %s s ---' % (time.time() - start_time))

        return self.P

    def dot(self, x, y):
        """
        Return the preconditioned dot product <P x, y>

        Uses 128-bit floating point math for vector dot products
        """
        return dot128(self.P.dot(x), y)

    def solve(self, x):
        """
        Solve the (sparse) linear system P x = y and return y
        """
        start_time = time.time()
        if self.sparse and self.use_pyamg and have_pyamg:
            y = self.ml.solve(x, x0=rand(self.P.shape[0]),
                              tol=self.solve_tol,
                              accel="cg",
                              maxiter=300,
                              cycle="W")
        elif self.sparse:
            y = spsolve(self.P, x)
        else:
            y = np.linalg.solve(self.P, x)
        logger.info("--- Precon applied in %s seconds ---",
                    time.time()-start_time)
        return y

    def half(self, x):
        """
        Return the product of P^{1/2} x
        """
        start_time = time.time()

        if self.sparse and self.use_pyamg and have_pyamg:
            raise NotImplementedError("Cholesky decomposition is not implemented for PyAMG")
        if self.sparse:
            raise NotImplementedError("Cholesky decomposition is not implemented for sparse")
        else:
            y = np.dot(np.linalg.cholesky(self.P), x)

        logger.info("--- Precon applied in %s seconds ---",
                    time.time()-start_time)

        return y

    def get_coeff(self, r):
        raise NotImplementedError("Must be overridden by subclasses")

    def estimate_mu(self, atoms, H=None):
        """
        Estimate optimal preconditioner coefficient \mu

        \mu is estimated from a numerical solution of

            [dE(p+v) -  dE(p)] \cdot v = \mu < P1 v, v >

        with perturbation

            v(x,y,z) = H (sin(x / Lx), sin(y / Ly), sin(z / Lz))

        After the optimal \mu is found, self.mu will be set to its value.

        If `atoms` is an instance of Filter an additional \mu_c
        will be computed for the cell degrees of freedom .

        Args:
            atoms: Atoms object for initial system

            H: 3x3 array or None
                Magnitude of deformation to apply. Default is 1e-2*rNN*np.eye(3)

        Returns:
            mu   : float
            mu_c : float or None
        """

        if self.dim != 3:
            raise ValueError("Automatic calculation of mu only possible for "
                             "three-dimensional preconditioners. Try setting "
                             "mu manually instead.")

        if self.r_NN is None:
            self.r_NN = estimate_nearest_neighbour_distance(atoms)

        # deformation matrix, default is diagonal
        if H is None:
            H = 1e-2*self.r_NN*np.eye(3)

        # compute perturbation
        p = atoms.get_positions()
        Lx, Ly, Lz = [ p[:, i].max() - p[:, i].min() for i in range(3) ]
        logger.debug('estimate_mu(): Lx=%.1f Ly=%.1f Lz=%.1f',
                     Lx, Ly, Lz)

        x, y, z = p.T
        # sine_vr = [np.sin(x/Lx), np.sin(y/Ly), np.sin(z/Lz)], but we need
        # to take into account the possibility that one of Lx/Ly/Lz is zero.
        sine_vr = [x, y, z]

        for i, L in enumerate([Lx, Ly, Lz]):
            if L == 0:
                logger.warning("Cell length L[%d] == 0. Setting H[%d,%d] = 0." % (i,i,i))
                H[i,i] = 0.0
            else:
                sine_vr[i] = np.sin(sine_vr[i]/L)

        v = np.dot(H, sine_vr).T

        natoms = len(atoms)
        if isinstance(atoms, Filter):
            natoms = len(atoms.atoms)
            eps = H / self.r_NN
            v[natoms:, :] = eps

        v1 = v.reshape(-1)

        # compute LHS
        dE_p = -atoms.get_forces().reshape(-1)
        atoms_v = atoms.copy()
        atoms_v.set_calculator(atoms.get_calculator())
        if isinstance(atoms, Filter):
            atoms_v = atoms.__class__(atoms_v)
            if hasattr(atoms, 'constant_volume'):
                atoms_v.constant_volume = atoms.constant_volume
        atoms_v.set_positions(p + v)
        dE_p_plus_v = -atoms_v.get_forces().reshape(-1)

        # compute left hand side
        LHS = (dE_p_plus_v - dE_p) * v1

        # assemble P with \mu = 1
        self.mu = 1.0
        self.mu_c = 1.0
        if self.sparse:
            P1 = self._make_sparse_precon(atoms, initial_assembly=True)
        else:
            P1 = self._make_dense_precon(atoms, initial_assembly=True)

        # compute right hand side
        RHS = P1.dot(v1) * v1

        # use partial sums to compute separate mu for positions and cell DoFs
        self.mu = float(sum128(LHS[:3*natoms])/sum128(RHS[:3*natoms]))
        if self.mu < 1.0:
            logger.info('mu (%.3f) < 1.0, capping at mu=1.0', self.mu)
            self.mu = 1.0

        if isinstance(atoms, Filter):
            self.mu_c = float(sum128(LHS[3*natoms:])/sum128(RHS[3*natoms:]))
            if self.mu_c < 1.0:
                logger.info('mu_c (%.3f) < 1.0, capping at mu_c=1.0', self.mu_c)
                self.mu_c = 1.0

        logger.info('estimate_mu(): mu=%r, mu_c=%r', self.mu, self.mu_c)

        self.P = None # force a rebuild with new mu (there may be fixed atoms)
        return (self.mu, self.mu_c)


class Pfrommer(object):
    """
    Use initial guess for inverse Hessian from Pfrommer et al. as a simple preconditioner

    J. Comput. Phys. vol 131 p233-240 (1997)
    """

    def __init__(self, bulk_modulus=500*units.GPa, phonon_frequency=50*THz,
                 apply_positions=True, apply_cell=True):
        """
        Default bulk modulus is 500 GPa and default phonon frequency is 50 THz
        """

        self.bulk_modulus = bulk_modulus
        self.phonon_frequency = phonon_frequency
        self.apply_positions = apply_positions
        self.apply_cell = apply_cell
        self.H0 = None

    def make_precon(self, atoms):
        if self.H0 is not None:
            # only build H0 on first call
            return NotImplemented

        ndof = 3*len(atoms)
        variable_cell = False
        if isinstance(atoms, Filter):
            variable_cell = True
            atoms = atoms.atoms

        # position DoF
        g0 = np.dot(atoms.cell, atoms.cell.T)
        omega = self.phonon_frequency
        mass = atoms.get_masses().mean()
        block = np.eye(3)/(mass*omega**2)
        blocks = [block]*len(atoms)

        # cell DoF
        if variable_cell:
            coeff = 1.0
            if self.apply_cell:
                coeff = 1.0/(3*self.bulk_modulus)
            blocks.append(np.diag([coeff]*9))

        self.H0 = sparse.block_diag(blocks, format='csr')
        return NotImplemented

    def dot(self, x, y):
        """
        Return the preconditioned dot product <P x, y>

        Uses 128-bit floating point math for vector dot products
        """
        raise NotImplementedError

    def solve(self, x):
        """
        Solve the (sparse) linear system P x = y and return y
        """
        y = self.H0.dot(x)
        return y

    def half(self, x):
        """
        Return the product of P^{1/2} x
        """
        raise NotImplementedError


class C1(Precon):
    """Creates matrix by inserting a constant whenever r_ij is less than r_cut.
    """

    def __init__(self, r_cut=None, mu=None, mu_c=None, dim=3, c_stab=0.1,
                 force_stab=False, sparse=True,
                 recalc_mu=False, array_convention="C",
                 use_pyamg=True, solve_tol=1e-9,
                 apply_positions=True, apply_cell=True):
        Precon.__init__(self, r_cut=r_cut, mu=mu, mu_c=mu_c,
                        dim=dim, c_stab=c_stab,
                        force_stab=force_stab,
                        sparse=sparse, recalc_mu=recalc_mu,
                        array_convention=array_convention,
                        use_pyamg=use_pyamg, solve_tol=solve_tol,
                        apply_positions=apply_positions,
                        apply_cell=apply_cell)

    def get_coeff(self, r):
        return -self.mu*np.ones_like(r)


class Exp(Precon):
    """Creates matrix with values decreasing exponentially with distance.
    """

    def __init__(self, A=3.0, r_cut=None, r_NN=None, mu=None, mu_c=None, dim=3, c_stab=0.1,
                 force_stab=False, sparse=True, recalc_mu=False, array_convention="C",
                 use_pyamg=True, solve_tol=1e-9,
                 apply_positions=True, apply_cell=True):
        """Initialise an Exp preconditioner with given parameters.

        Args:
            r_cut, mu, c_stab, dim, sparse, recalc_mu, array_convention: see
                precon.__init__()
            A: coefficient in exp(-A*r/r_NN). Default is A=3.0.
        """
        Precon.__init__(self, r_cut=r_cut, r_NN=r_NN,
                        mu=mu, mu_c=mu_c, dim=dim, c_stab=c_stab,
                        force_stab=force_stab,
                        sparse=sparse, recalc_mu=recalc_mu,
                        array_convention=array_convention,
                        use_pyamg=use_pyamg,
                        solve_tol=solve_tol,
                        apply_positions=apply_positions,
                        apply_cell=apply_cell)

        self.A = A

    def get_coeff(self, r):
        return -self.mu * np.exp(-self.A*(r/self.r_NN - 1))
