#__docformat__ = "restructuredtext en"
# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

import numpy as np
from numpy import atleast_1d, eye, mgrid, argmin, zeros, shape, empty, \
     squeeze, vectorize, asarray, absolute, sqrt, Inf, asfarray, isinf
from ase.utils.linesearch import LineSearch
from ase.optimize import Optimizer

# These have been copied from Numeric's MLab.py
# I don't think they made the transition to scipy_core

# Modified from scipy_optimize
abs = absolute
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max
__version__="0.1"

class FMIN_BFGS(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=1.,
                 trajectory=None, c1=1e-4, c2=0.9, alpha=1e-2):
        """Minimize a function using the BFGS algorithm.

        Notes:

            Optimize the function, f, whose gradient is given by fprime
            using the quasi-Newton method of Broyden, Fletcher, Goldfarb,
            and Shanno (BFGS) See Wright, and Nocedal 'Numerical
            Optimization', 1999, pg. 198.

        *See Also*:

          scikits.openopt : SciKit which offers a unified syntax to call
                            this and other solvers.

        """
        self.maxstep = maxstep
        self.alpha = alpha
        self.H = None
        self.c1 = c1
        self.c2 = c2
        self.force_calls = 0
        self.function_calls = 0
        self.r0 = None
        self.g0 = None
        self.e0 = None
        self.load_restart = False
        self.task = 'START'

        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

    def read(self):
        self.r0, self.g0, self.e0, self.task, self.H = self.load()
        self.load_restart = True    

    def step(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        r = r.reshape(-1)
        g = -f.reshape(-1) / self.alpha
        self.update(r, g, self.r0, self.g0)
        e = atoms.get_potential_energy() / self.alpha
        if not self.e0:
            self.e0 = e + 5000

        p = -np.dot(self.H,g)
        ls = LineSearch()
        alpha_k, fc, gc, e, self.e0, fkp1 = \
           ls._line_search(self.func, self.fprime, r, p, g, e, self.e0,
                           maxstep=self.maxstep, c1=self.c1,
                           c2=self.c2)
        #if alpha_k is None:  # line search failed try different one.
        #    alpha_k, fc, gc, e, e0, gfkp1 = \
        #             line_search(self.func, self.fprime,r,p,g,
        #                         e,self.e0)
        if alpha_k is None:
            # This line search also failed to find a better solution.
            raise ValueError('Warning: Desired error not necessarily' +
                             'achieved due to precision loss')

        print 'al',alpha_k,fc,gc,e,self.e0,fkp1
        dr = alpha_k * p
        atoms.set_positions((r+dr).reshape(len(atoms),-1))
        self.r0 = r
        self.g0 = g
        self.dump((self.r0, self.g0, self.e0, self.task, self.H))

    def update(self, r, g, r0, g0):
        self.I = eye(len(self.atoms) * 3, dtype=int)
        if self.H is None:
            self.H = eye(3 * len(self.atoms))
            #self.H = eye(3 * len(self.atoms)) / self.alpha
            return

        dr = r - r0
        dg = g - g0 

        try: # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (np.dot(dg,dr))
        except ZeroDivisionError:
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"
        if isinf(rhok): # this is patch for np
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"
        A1 = self.I - dr[:, np.newaxis] * dg[np.newaxis, :] * rhok
        A2 = self.I - dg[:, np.newaxis] * dr[np.newaxis, :] * rhok
        self.H = np.dot(A1, np.dot(self.H, A2)) + rhok * dr[:, np.newaxis] \
                 * dr[np.newaxis, :]

    def func(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.function_calls += 1
        # Scale the problem as SciPy uses I as initial Hessian.
        return self.atoms.get_potential_energy() / self.alpha
    
    def fprime(self, x):
        """Gradient of the objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.force_calls += 1
        # Remember that forces are minus the gradient!
        # Scale the problem as SciPy uses I as initial Hessian.
        print 'h0',self.alpha
        return - self.atoms.get_forces().reshape(-1) / self.alpha

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import PickleTrajectory
            traj = PickleTrajectory(traj, 'r')
        self.H = None
        atoms = traj[0]
        r0 = atoms.get_positions().ravel()
        f0 = atoms.get_forces().ravel()
        for atoms in traj:
            r = atoms.get_positions().ravel()
            f = atoms.get_forces().ravel()
            self.update(r, f, r0, f0)
            r0 = r
            f0 = f
        self.r0 = r0
        self.f0 = f0

def wrap_function(function, args):
    ncalls = [0]
    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)
    return ncalls, function_wrapper

def _cubicmin(a,fa,fpa,b,fb,c,fc):
    # finds the minimizer for a cubic polynomial that goes through the
    #  points (a,fa), (b,fb), and (c,fc) with derivative at a of fpa.
    #
    # if no minimizer can be found return None
    #
    # f(x) = A *(x-a)^3 + B*(x-a)^2 + C*(x-a) + D

    C = fpa
    D = fa
    db = b-a
    dc = c-a
    if (db == 0) or (dc == 0) or (b==c): return None
    denom = (db*dc)**2 * (db-dc)
    d1 = empty((2,2))
    d1[0,0] = dc**2
    d1[0,1] = -db**2
    d1[1,0] = -dc**3
    d1[1,1] = db**3
    [A,B] = np.dot(d1,asarray([fb-fa-C*db,fc-fa-C*dc]).flatten())
    A /= denom
    B /= denom
    radical = B*B-3*A*C
    if radical < 0:  return None
    if (A == 0): return None
    xmin = a + (-B + sqrt(radical))/(3*A)
    return xmin

def _quadmin(a,fa,fpa,b,fb):
    # finds the minimizer for a quadratic polynomial that goes through
    #  the points (a,fa), (b,fb) with derivative at a of fpa
    # f(x) = B*(x-a)^2 + C*(x-a) + D
    D = fa
    C = fpa
    db = b-a*1.0
    if (db==0): return None
    B = (fb-D-C*db)/(db*db)
    if (B <= 0): return None
    xmin = a  - C / (2.0*B)
    return xmin

def zoom(a_lo, a_hi, phi_lo, phi_hi, derphi_lo,
         phi, derphi, phi0, derphi0, c1, c2):
    maxiter = 10
    i = 0
    delta1 = 0.2  # cubic interpolant check
    delta2 = 0.1  # quadratic interpolant check
    phi_rec = phi0
    a_rec = 0
    while 1:
        # interpolate to find a trial step length between a_lo and a_hi
        # Need to choose interpolation here.  Use cubic interpolation and then 
        #if the result is within delta * dalpha or outside of the interval 
        #bounded by a_lo or a_hi then use quadratic interpolation, if the 
        #result is still too close, then use bisection

        dalpha = a_hi-a_lo;
        if dalpha < 0: a,b = a_hi,a_lo
        else: a,b = a_lo, a_hi

        # minimizer of cubic interpolant
        #    (uses phi_lo, derphi_lo, phi_hi, and the most recent value of phi)
        #      if the result is too close to the end points (or out of the 
        #         interval) then use quadratic interpolation with phi_lo, 
        #         derphi_lo and phi_hi
        #      if the result is stil too close to the end points (or out of 
        #         the interval) then use bisection

        if (i > 0):
            cchk = delta1*dalpha
            a_j = _cubicmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, a_rec, 
                            phi_rec)
        if (i==0) or (a_j is None) or (a_j > b-cchk) or (a_j < a+cchk):
            qchk = delta2*dalpha
            a_j = _quadmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi)
            if (a_j is None) or (a_j > b-qchk) or (a_j < a+qchk):
                a_j = a_lo + 0.5*dalpha
#                print "Using bisection."
#            else: print "Using quadratic."
#        else: print "Using cubic."

        # Check new value of a_j

        phi_aj = phi(a_j)
        if (phi_aj > phi0 + c1*a_j*derphi0) or (phi_aj >= phi_lo):
            phi_rec = phi_hi
            a_rec = a_hi
            a_hi = a_j
            phi_hi = phi_aj
        else:
            derphi_aj = derphi(a_j)
            if abs(derphi_aj) <= -c2*derphi0:
                a_star = a_j
                val_star = phi_aj
                valprime_star = derphi_aj
                break
            if derphi_aj*(a_hi - a_lo) >= 0:
                phi_rec = phi_hi
                a_rec = a_hi
                a_hi = a_lo
                phi_hi = phi_lo
            else:
                phi_rec = phi_lo
                a_rec = a_lo
            a_lo = a_j
            phi_lo = phi_aj
            derphi_lo = derphi_aj
        i += 1
        if (i > maxiter):
            a_star = a_j
            val_star = phi_aj
            valprime_star = None
            break
    return a_star, val_star, valprime_star

def line_search(f, myfprime, xk, pk, gfk, old_fval, old_old_fval,
                args=(), c1=1e-4, c2=0.9, amax=50):
    """Find alpha that satisfies strong Wolfe conditions.

    Parameters:

        f : callable f(x,*args)
            Objective function.
        myfprime : callable f'(x,*args)
            Objective function gradient (can be None).
        xk : ndarray
            Starting point.
        pk : ndarray
            Search direction.
        gfk : ndarray
            Gradient value for x=xk (xk being the current parameter
            estimate).
        args : tuple
            Additional arguments passed to objective function.
        c1 : float
            Parameter for Armijo condition rule.
        c2 : float
            Parameter for curvature condition rule.

    Returns:

        alpha0 : float
            Alpha for which ``x_new = x0 + alpha * pk``.
        fc : int
            Number of function evaluations made.
        gc : int
            Number of gradient evaluations made.

    Notes:

        Uses the line search algorithm to enforce strong Wolfe
        conditions.  See Wright and Nocedal, 'Numerical Optimization',
        1999, pg. 59-60.

        For the zoom phase it uses an algorithm by [...].

    """

    global _ls_fc, _ls_gc, _ls_ingfk
    _ls_fc = 0
    _ls_gc = 0
    _ls_ingfk = None
    def phi(alpha):
        global _ls_fc
        _ls_fc += 1
        return f(xk+alpha*pk,*args)

    if isinstance(myfprime,type(())):
        def phiprime(alpha):
            global _ls_fc, _ls_ingfk
            _ls_fc += len(xk)+1
            eps = myfprime[1]
            fprime = myfprime[0]
            newargs = (f,eps) + args
            _ls_ingfk = fprime(xk+alpha*pk,*newargs)  # store for later use
            return np.dot(_ls_ingfk,pk)
    else:
        fprime = myfprime
        def phiprime(alpha):
            global _ls_gc, _ls_ingfk
            _ls_gc += 1
            _ls_ingfk = fprime(xk+alpha*pk,*args)  # store for later use
            return np.dot(_ls_ingfk,pk)

    print 'ls2'

    alpha0 = 0
    phi0 = old_fval
    derphi0 = np.dot(gfk,pk)

    alpha1 = pymin(1.,1.01*2*(phi0-old_old_fval)/derphi0)

    if alpha1 == 0:
        # This shouldn't happen. Perhaps the increment has slipped below
        # machine precision?  For now, set the return variables skip the
        # useless while loop, and raise warnflag=2 due to possible imprecision.
        alpha_star = None
        fval_star = old_fval
        old_fval = old_old_fval
        fprime_star = None

    phi_a1 = phi(alpha1)
    #derphi_a1 = phiprime(alpha1)  evaluated below

    phi_a0 = phi0
    derphi_a0 = derphi0

    i = 1
    maxiter = 10
    while 1:         # bracketing phase
        if alpha1 == 0:
            break
        if (phi_a1 > phi0 + c1*alpha1*derphi0) or \
           ((phi_a1 >= phi_a0) and (i > 1)):
            alpha_star, fval_star, fprime_star = \
                        zoom(alpha0, alpha1, phi_a0,
                             phi_a1, derphi_a0, phi, phiprime,
                             phi0, derphi0, c1, c2)
            break

        derphi_a1 = phiprime(alpha1)
        if (abs(derphi_a1) <= -c2*derphi0):
            alpha_star = alpha1
            fval_star = phi_a1
            fprime_star = derphi_a1
            break

        if (derphi_a1 >= 0):
            alpha_star, fval_star, fprime_star = \
                        zoom(alpha1, alpha0, phi_a1,
                             phi_a0, derphi_a1, phi, phiprime,
                             phi0, derphi0, c1, c2)
            break

        alpha2 = 2 * alpha1   # increase by factor of two on each iteration
        i = i + 1
        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi(alpha1)
        derphi_a0 = derphi_a1

        # stopping test if lower function not found
        if (i > maxiter):
            alpha_star = alpha1
            fval_star = phi_a1
            fprime_star = None
            break

    if fprime_star is not None:
        # fprime_star is a number (derphi) -- so use the most recently
        # calculated gradient used in computing it derphi = gfk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        fprime_star = _ls_ingfk

    return alpha_star, _ls_fc, _ls_gc, fval_star, old_fval, fprime_star

