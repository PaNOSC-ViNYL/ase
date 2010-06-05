# -*- coding: utf-8 -*-
import numpy as np

from numpy import exp, log, sqrt

eV = 1.0
_e = 1.60217733e-19          # elementary charge (C)
J = 1.0 / _e
kJ = 1000.0 * J

ang = angstrom = Ang = Angstrom = 1.0
meter = m = 1e10 * Angstrom

Pascal = Pa = J / m**3  # J/m^3 or N/m^2
GPa = 1e9*Pascal
    
from Scientific.Functions.LeastSquares import leastSquaresFit

class EquationOfState:
    """Fit equation of state for bulk systems.

    The following equations of state are available:

    poly
        A simply third order polynomial fit
        
    taylor
        A third order Taylor series expansion about the minimum volume

    murnaghan
        PRB 28, 5480 (1983)

    birch
        Intermetallic compounds: Principles and Practice, Vol I: Principles. pages 195-210

    birchmurnaghan
        PRB 70, 224107

    pouriertarantola
        PRB 70, 224107

    vinet
        PRB 70, 224107

    antonschmidt
        Intermetallics 11, 23-32 (2003)

    Use::

       eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
       v0, e0, B = eos.fit()
       eos.plot()

    """
    def __init__(self, volumes, energies,eos='poly'):
        '''
        eos is one of:
           poly (original eos)
           taylor
           murnaghan
           birch
           birchmurnaghan
           pouriertarantola
           vinet
           antonschmidt
        '''
        
        self.v = np.array(volumes)
        self.e = np.array(energies)
        self.eos_string = eos
        self.eos = eval('self.%s' % eos.lower())
        
    def poly(self,parameters,V):
        'polynomial fit'

        c0 = parameters[0]
        c1 = parameters[1]
        c2 = parameters[2]
        c3 = parameters[3]

        E = c0 + c1*V + c2*V**2 + c3*V**3
        return E
        
    def taylor(self,parameters,V):
        'Taylor Expansion up to 3rd order about V0'
        E0 = parameters[0]
        beta = parameters[1]
        alpha = parameters[2]
        V0 = parameters[3]

        E = E0 + beta/2.*(V-V0)**2/V0 + alpha/6.*(V-V0)**3/V0
        return E

    def murnaghan(self,parameters,V):
        'From PRB 28,5480 (1983'
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]

        E = E0 + B0*V/BP*(((V0/V)**BP)/(BP-1)+1) - V0*B0/(BP-1)
        return E

    def birch(self,parameters,V):
        '''
        From Intermetallic compounds: Principles and Practice, Vol. I: Principles
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
        paper downloaded from Web

        case where n=0
        '''
        E0=parameters[0]
        B0=parameters[1]
        BP=parameters[2]
        V0=parameters[3]

        E = (E0
             + 9.0/8.0*B0*V0*((V0/V)**(2.0/3.0) - 1.0)**2
             + 9.0/16.0*B0*V0*(BP-4.)*((V0/V)**(2.0/3.0) - 1.0)**3)
        return E

    def birchmurnaghan(self,parameters,V):
        'BirchMurnaghan equation from PRB 70, 224107'
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]

        eta = (V/V0)**(1./3.)
        E = E0 + 9.*B0*V0/16.*(eta**2-1)**2*(6 + BP*(eta**2-1.) - 4.*eta**2)
        return E

    def pouriertarantola(self,parameters,V):
        'Pourier-Tarantola equation from PRB 70, 224107'
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]

        eta = (V/V0)**(1./3.)
        squiggle = -3.*log(eta)

        E = E0 + B0*V0*squiggle**2/6.*(3. + squiggle*(BP - 2))
        return E

    def vinet(self,parameters,V):
        'Vinet equation from PRB 70, 224107'
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]

        eta = (V/V0)**(1./3.)

        E = (E0 + 2.*B0*V0/(BP-1.)**2
             * (2. - (5. +3.*BP*(eta-1.)-3.*eta)*exp(-3.*(BP-1.)*(eta-1.)/2.)))
        return E
    
    def antonschmidt(self,parameters,V):
        '''From Intermetallics 11, 23-32 (2003)

        Einf should be E_infinity, i.e. infinite separation, but
        according to the paper it does not provide a good estimate
        of the cohesive energy. They derive this equation from an
        empirical formula for the volume dependence of pressure,

        E(vol) = E_inf + int(P dV) from V=vol to V=infinity

        but the equation breaks down at large volumes, so E_inf
        is not that meaningful

        n should be about -2 according to the paper.

        I find this equation does not fit volumetric data as well
        as the other equtions do.
        '''
        Einf = parameters[0]
        B = parameters[1]
        n = parameters[2]
        V0 = parameters[3]

        E = B*V0/(n+1.) * (V/V0)**(n+1.)*(log(V/V0)-(1./(n+1.))) + Einf
        return E

    def parabola(self,parameters,x):
        '''
        parabola polynomial function

        this function is used to fit the data to get good guesses for
        the equation of state fits

        a 4th order polynomial fit to get good guesses for
        was not a good idea because for noisy data the fit is too wiggly
        2nd order seems to be sufficient, and guarentees a single minimum'''
        a=parameters[0]
        b=parameters[1]
        c=parameters[2]

        return a + b*x + c*x**2
    
    def fit(self):
        """Calculate volume, energy, and bulk modulus.

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          from ase.utils.eos import *
          eos = EquationOfState(vols,energies)
          v0, e0, B = eos.fit()
          print B / GPa, 'GPa'
        """
        parabola_parameters,chisq = leastSquaresFit(self.parabola,
                                                    parameters=[min(self.e),1,1],
                                                    data=zip(self.v,self.e))

        ## Here I just make sure the minimum is bracketed by the volumes
        ## this if for the solver
        minvol = min(self.v)
        maxvol = max(self.v)

        # the minimum of the parabola is at dE/dV = 0, or 2*c V +b =0
        c = parabola_parameters[2]
        b = parabola_parameters[1]
        parabola_vmin = -b/2/c

        if not (minvol < parabola_vmin and parabola_vmin < maxvol):
            print 'Warning the minimum volume of a fitted parabola is not in your volumes. You may not have a minimum in your dataset'

        # evaluate the parabola at the minimum to estimate the groundstate energy
        E0 = self.parabola(parabola_parameters,parabola_vmin)
        # estimate the bulk modulus from Vo*E''.  E'' = 2*c
        B0 = 2*c*parabola_vmin

        if self.eos == self.antonschmidt:
            BP = -2
        else:
            BP = 4

        initial_guess = [E0, B0, BP, parabola_vmin]

        # now fit the equation of state
        self.eos_parameters, chisq = leastSquaresFit(self.eos,
                                                initial_guess,
                                                data=zip(self.v,self.e))

        if self.eos_string == 'poly':
            c0,c1,c2,c3 = self.eos_parameters
            # find minimum E in E = c0 + c1*V + c2*V**2 + c3*V**3
            # dE/dV = c1+ 2*c2*V + 3*c3*V**2 = 0
            # solve by quadratic formula with the positive root

            a = 3*c3
            b = 2*c2
            c = c1
            
            self.v0 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
            self.e0 = self.eos(self.eos_parameters,self.v0)
            self.B = (2*c2 + 6*c3*self.v0)*self.v0    
        else:
            self.e0 = self.eos_parameters[0]
            self.B = self.eos_parameters[1]
            self.v0 = self.eos_parameters[3]

        return self.v0, self.e0, self.B

    def plot(self, filename=None, show=None):
        """Plot fitted energy curve.

        Uses Matplotlib to plot the energy curve.  Use *show=True* to
        show the figure and *filename='abc.png'* or
        *filename='abc.eps'* to save the figure to a file."""
        
        #import matplotlib.pyplot as plt
        import pylab as plt

        if self.v0 is None:
            self.fit()
            
        if filename is None and show is None:
            show = True

        x = 3.95
        f = plt.figure(figsize=(x * 2.5**0.5, x))
        f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
        plt.plot(self.v, self.e, 'o')
        x = np.linspace(min(self.v), max(self.v), 100)
        y = self.eos(self.eos_parameters,x)
        
        plt.plot(x, y, '-r')
        plt.xlabel(u'volume [Å$^3$]')
        plt.ylabel(u'energy [eV]')
        plt.title(u'E: %.3f eV, V: %.3f Å$^3$, B: %.3f GPa' %
                  (self.e0, self.v0, self.B / GPa))

        plt.text(self.v0,max(self.e),'EOS: %s' % self.eos_string)

        if show:
            plt.show()
        if filename is not None:
            f.savefig(filename)

        return f
