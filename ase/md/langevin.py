"""Langevin dynamics class."""


import numpy as np
from numpy.random import standard_normal
from ase.md.md import MolecularDynamics
from ase.parallel import world


class Langevin(MolecularDynamics):
    """Langevin (constant N, V, T) molecular dynamics.

    Usage: Langevin(atoms, dt, temperature, friction)

    atoms
        The list of atoms.
        
    dt
        The time step.

    temperature
        The desired temperature, in energy units.

    friction
        A friction coefficient, typically 1e-4 to 1e-2.

    fixcm
        If True, the position and momentum of the center of mass is
        kept unperturbed.  Default: True.

    The temperature and friction are normally scalars, but in principle one
    quantity per atom could be specified by giving an array.

    RATTLE constraints are implemented, see:
    E. V.-Eijnden, and G. Ciccotti, Chem. Phys. Lett. 429, 310 (2006)    

    A single step amounts to:

        x(n+1) = x(n) + dt*v(n) + A(n)
        v(n+1) = v(n) + 0.5*dt*(f(x(n+1))+f(x(n))) 
                 - dt*y*v(n) + dt**0.5*o*xi(n) + y*A(n)

        where: 
        A(n) = 0.5*dt**2(f(x(n))-y*v(n)) 
               + o*dt**3/2(0.5*xi(n)-(2*3**0.5)**-1*eta(n))

        y is the friction coeff, o(sigma) is (2*kB*T*m_i*y)**1/2

        xi and eta are random variables with mean 0 and covariance.

        However, to allow for the possibility of constraints we 
        rewrite the equations the following way:

        x(n+1) = x(n) + dt*p(n)
        v(n+1) = p(n) - 0.5*dt*y*v(n) - y*A(n) 
                 - o*dt**0.5*(2*3**0.5)**-1*eta(n) 
                 + 0.5*dt*f(n+1) + 0.5*dt**0.5*o*xi(n)

        where:
        p(n) = v(n) + A(n)*dt**-1.

    This dynamics accesses the atoms using Cartesian coordinates."""
    
    # Helps Asap doing the right thing.  Increment when changing stuff:
    _lgv_version = 3
    
    def __init__(self, atoms, timestep, temperature, friction, fixcm=True,
                 trajectory=None, logfile=None, loginterval=1,
                 communicator=world):
        MolecularDynamics.__init__(self, atoms, timestep, trajectory,
                                   logfile, loginterval)
        self.temp = temperature
        self.frict = friction
        self.fixcm = fixcm  # will the center of mass be held fixed?
        self.communicator = communicator
        self.updatevars()
        
    def set_temperature(self, temperature):
        self.temp = temperature
        self.updatevars()

    def set_friction(self, friction):
        self.frict = friction
        self.updatevars()

    def set_timestep(self, timestep):
        self.dt = timestep
        self.updatevars()

    def updatevars(self):
        dt = self.dt

        dt = self.dt
        T = self.temp
        fr = self.frict
        masses = self.masses
        sigma = np.sqrt(2*T*fr/masses)
        c1 = 0.5*dt**2
        c2 = c1 * fr
        c3 = sigma*dt*dt**0.5/2.0
        c4 = sigma*dt*dt**0.5/(2.0*np.sqrt(3))
        v1 = 0.5*dt
        v2 = c3/dt
        v3 = c4/dt

        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.fr = fr

        # Works in parallel Asap, #GLOBAL number of atoms:
        self.natoms = self.atoms.get_number_of_atoms() 

    def step(self, f):
        atoms = self.atoms
        natoms = self.natoms

        v = atoms.get_velocities()

        xi = standard_normal(size=(natoms, 3))
        eta = standard_normal(size=(natoms, 3))

        if self.communicator is not None:
            self.communicator.broadcast(xi, 0)
            self.communicator.broadcast(eta, 0)

        # Begin calculating A
        A = self.c1*f/self.masses - self.c2*v + self.c3*xi \
            + self.c4*eta

        # Make V and A/dt
        A2 = A/self.dt
        V = v + A2
        x = atoms.get_positions()

        if self.fixcm:
            old_cm = atoms.get_center_of_mass()
            m = atoms.get_masses()

        # Step: x^n -> x^(n+1) - this applies constraints if any.
        atoms.set_positions(x + self.dt*V)
    
        if self.fixcm:
            new_cm = atoms.get_center_of_mass()
            d = old_cm-new_cm
            atoms.translate(d)

        # recalc vels after RATTLE constraints are applied 
        V = (self.atoms.get_positions() - x) / self.dt
        f = atoms.get_forces(md=True)

        # Update the velocities 
        V += self.v2*xi + self.v1*f/self.masses - self.fr*A \
             - self.fr*self.v1*v - self.v3*eta

        # Second part of RATTLE taken care of here
        atoms.set_momenta(V*self.masses)

        return f
