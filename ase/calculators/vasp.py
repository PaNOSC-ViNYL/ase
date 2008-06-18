"""This module defines an ASE interface to VASP.

http://cms.mpi.univie.ac.at/vasp/
"""

# Copyright (C) 2008 CSC - Scientific Computing Ltd.

import os
from os.path import join, isfile, islink

import sys
import numpy as npy

from ase.data import chemical_symbols
from ase.data import atomic_numbers
from ase.units import Bohr, Hartree



class _DownTheDrain:
    """Definition of a stream that throws away all output."""

    def write(self, string):
        pass

    def flush(self):
        pass

devnull = _DownTheDrain()


class Vasp:
    """Class for doing Vasp calculations.

    The default parameters are very close to those that the VASP
    Fortran code would use.  These are the exceptions: XXX ???

    Use the set_inp method to set extra INPUT parameters::

      calc.set_inp('nstep', 30)

    """
    def __init__(self, **kwargs):
        """Construct Vasp-calculator object.

        """

        # Parameters that can be set in INCAR. The values which are None
        # are not written and default parameters of VASP are used for them.
        self.incar_parameters = {
            'NG':         None,    # FFT mesh for wavefunctions
            'NGF':        None,    # FFT mesh for charges
            'NBANDS': 	  None,    # number of bands included in the calculation
            'NBLK': None, #	blocking for some BLAS calls (Sec. 6.5)
            'SYSTEM': None, #	name of System
            'NWRITE': None, #	verbosity write-flag (how much is written)
            'ISTART': None, #	startjob: 0-new 1-cont 2-samecut
            'ICHARG': None, #	charge: 1-file 2-atom 10-const
            'ISPIN': None, #	spin polarized calculation (2-yes 1-no)
            'MAGMOM': None, # 	initial mag moment / atom
            'INIWAV': None, # 	initial electr wf. : 0-lowe 1-rand
            'ENCUT': None, # 	energy cutoff in eV
            'PREC': None, # 	precession: medium, high or low
            'NELM': None, #
            'NELMIN': None, #
            'NELMDL': None, # 	nr. of electronic steps
            'EDIFF': None, # 	stopping-criterion for electronic upd.
            'EDIFFG': None, # 	stopping-criterion for ionic upd.
            'NSW': None, # 	number of steps for ionic upd.
            # 'NBLOCK' and KBLOCK 	inner block; outer block
            'IBRION': None, # 	ionic relaxation: 0-MD 1-quasi-New 2-CG
            'ISIF': None, # 	calculate stress and what to relax
            'IWAVPR': None, # 	prediction of wf.: 0-non 1-charg 2-wave 3-comb
            'ISYM': None, # 	symmetry: 0-nonsym 1-usesym
            'SYMPREC': None, # 	precession in symmetry routines
            'LCORR': None, # 	Harris-correction to forces
            'POTIM': None, # 	time-step for ion-motion (fs)
            'TEBEG': None, #
            'TEEND': None, # 	temperature during run
            'SMASS': None, # 	Nose mass-parameter (am)
            # 'NPACO' and APACO 	distance and nr. of slots for P.C.
            'POMASS': None, # 	mass of ions in am
            'ZVAL': None, # 	ionic valence
            'RWIGS': None, # 	Wigner-Seitz radii
            'NELECT': None, # 	total number of electrons
            'NUPDOWN': None, # 	fix spin moment to specified value
            'EMIN': None, #
            'EMAX': None, # 	energy-range for DOSCAR file
            'ISMEAR': None, # 	part. occupancies: -5 Blochl -4-tet -1-fermi 0-gaus >0 MP
            'SIGMA': None, # 	broadening in eV -4-tet -1-fermi 0-gaus
            'ALGO': None, #	algorithm: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)
            'IALGO': None, # 	algorithm: use only 8 (CG) or 48 (RMM-DIIS)
            'LREAL': None, # 	non-local projectors in real space
            'ROPT': None, # 	number of grid points for non-local proj in real space
            'GGA': None, # 	xc-type: PW PB LM or 91
            'VOSKOWN': None, # 	use Vosko, Wilk, Nusair interpolation
            'DIPOL': None, # 	center of cell for dipol
            'AMIX': None, #
            'BMIX': None, # 	tags for mixing
            # 'WEIMIN, EBREAK, DEPER 	special control tags
            'TIME': None, # 	special control tag
            'LWAVE': None, #
            'LCHARG': None, #
            'LVTOT': None, # 	create WAVECAR/CHGCAR/LOCPOT
            'LELF': None, # 	create ELFCAR
            'LORBIT': None, # 	create PROOUT
            'NPAR': None, # 	parallelization over bands
            'LSCALAPACK': None, # 	switch off scaLAPACK
            'LSCALU': None, # 	switch of LU decomposition
            'LASYNC': None, # 	overlap communcation with calculations
            }

        self.input_parameters = {
            'xc': 'LDA', # exchange correlation potential
            'kpts': (1,1,1), # k-points
            'setups': None, # Special setups (e.g pv, sv, ...)
            'txt':'-', # Where the send information
            'workdir': None # define a working dir for calculator
            }

        self.converged = False
        self.set(**kwargs)
        
    def write_poscar(self):
        """Writes the atomic positions to the VASP's POSCAR format"""

        poscar = open('POSCAR', 'w')
        poscar.write('POSCAR created by Atomic Simulation Environment\n')
        poscar.write('%.5f\n' % 1.0)
        for v in self.cell:
            poscar.write('%.14f %.14f %.14f\n' %  tuple(v))

        for nspecies in self.symbols.values():
            poscar.write('%i ' % nspecies)
        poscar.write('\n')

        poscar.write('Cartesian\n')
        for symbol in self.symbols.keys():
            mask = [asymbol == symbol for asymbol in self.all_symbols]
            ind = npy.asarray(mask, bool)
            for pos in self.positions[ind]:
                poscar.write('%.14f %.14f %.14f\n' %  tuple(pos))

        poscar.close()

    def write_incar(self):
        p = self.incar_parameters
        incar = open('INCAR', 'w')
        incar.write('# INCAR created by Atomic Simulation Environment\n')
        for key, val in p.items():
            if val is not None:
                # special cases:
                if key == 'setups':
                    self.setups_changed = True
                elif key == 'NG':
                    incar.write('NGX = %i\n' % val[0])
                    incar.write('NGY = %i\n' % val[1])
                    incar.write('NGZ = %i\n' % val[2])
                elif key == 'NGF':
                    incar.write('NGFX = %i\n' % val[0])
                    incar.write('NGFY = %i\n' % val[1])
                    incar.write('NGFZ = %i\n' % val[2])
                elif key == 'MAGMOM':
                    if len(val) != self.natoms:
                        raise RuntimeError('Incorrect number of magnetic moments')
                    incar.write('MAGMOM = ')
                    [incar.write('%.4f ' % mom) for mom in val]
                    incar.write('\n')
                else:
                    incar.write('%s = %s\n' % (key, str(val)))
        incar.close()

    def write_potcar(self):
        potcar = open('POTCAR', 'w')
        for filename in self.ppp_list:
            source = os.popen('gunzip -c ' + filename, 'r').read()
            potcar.write(source)
        potcar.close()

    def write_kpoints(self):
        p = self.input_parameters
        kpoints = open('KPOINTS', 'w')
        kpoints.write('KPOINTS created by Atomic Simulation Environemnt\n')
        kpoints.write('0\n')
        kpoints.write('Monkhorst-Pack\n')
        [kpoints.write('%i ' % kpt) for kpt in p['kpts']]
        kpoints.write('\n')
        kpoints.close()

    def set(self, **kwargs):
        p1 = self.input_parameters
        p2 = self.incar_parameters

        for key in kwargs:
            if p1.has_key(key):
                p1[key] = kwargs[key]
            elif p2.has_key(key):
                p2[key] = kwargs[key]
            else:
                raise TypeError('Unknown keyword argument:' + key)

        self.converged = False

    def update(self, atoms):
        if (not self.converged or self.setups_changed or 
            self.natoms != len(atoms) or
            (self.all_symbols != atoms.get_chemical_symbols())):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        """Initialize a VASP calculation

        Constructs the POTCAR file. User should specify the PATH
        to the pseudopotentials in VASP_PP_PATH environment variable"""

        p = self.input_parameters

        self.all_symbols = atoms.get_chemical_symbols()
        self.natoms = len(atoms)

        # Determine the number of atomic species
        self.symbols = {}
        for atom in atoms:
            symbol = atom.get_symbol()
            if not self.symbols.has_key(symbol):
                self.symbols[symbol] = 1
            else:
                self.symbols[symbol] += 1

        if 'VASP_PP_PATH' in os.environ:
            pppaths = os.environ['VASP_PP_PATH'].split(':')
        else:
            pppaths = []

        self.ppp_list = []

        # if self.xc != 'LDA':
        #     xcname = 'GGA'
        # else:
        #     xcname = 'LDA'

        for symbol in self.symbols:
            try:
                name = symbol + '_' + p['setups'][symbol]
            except (TypeError, KeyError):
                name = symbol
            name += '/POTCAR.Z'
            found = False
            for path in pppaths:
                filename = join(path, name)
                if isfile(filename) or islink(filename):
                    found = True
                    self.ppp_list.append(filename)
                    break
            if not found:
                raise RuntimeError('No pseudopotential for %s!' % symbol)

        self.converged = False
        self.setups_changed = False

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree
        else:
            # Energy extrapolated to zero Kelvin:
            return  self.etotal_zero

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        # self.update(atoms)
        # return self.stress.copy()
        raise NotImplementedError

    def calculate(self, atoms):
        p = self.input_parameters

        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        curdir = os.getcwd()
        if p['workdir'] is not None:
            if not os.access(p['workdir'], os.F_OK):
                os.mkdir(p['workdir'])
            os.chdir(p['workdir'])
            
        self.write_incar()
        self.write_poscar()
        self.write_potcar()
        self.write_kpoints()

        vasp = os.environ['VASP_SCRIPT']
        # locals = {}
        # execfile(vasp, {}, locals)
        # exitcode = locals['exitcode']
        stderr = sys.stderr

        if p['txt'] is None:
            sys.stderr = devnull
        elif p['txt'] == '-':
            pass
        elif isinstance(p['txt'], str):
            sys.stderr = open(p['txt'], 'w')                    
        exitcode = os.system(vasp)
        sys.stderr = stderr
        if exitcode != 0:
            raise RuntimeError('Vasp exited with exit code: %d.  ' % exitcode)

        self.read()

        os.chdir(curdir)

        self.converged = True

    def read(self):
        """Read the total energy, forces etc. from OUTCAR"""

        text = open('OUTCAR', 'r').read().lower()
        lines = iter(text.split('\n'))


        self.efree = None
        for line in lines:
            # Energy:
            if line.startswith('  free energy    toten'):
                self.efree = float(line.split()[-2])
            if line.startswith('  energy without entropy'):
                self.etotal_zero = float(line.split()[-1])

            # Forces
            if line.rfind('total-force') > -1:
                print "Found a force"
                lines.next()
                forces = []
                for i in range(self.natoms):
                    forces.append(npy.array([float(f) for f in lines.next().split()[3:6]]))
                self.forces = npy.array(forces)

        if self.efree is None:
            raise RuntimeError('No total energy in OUTCAR')

