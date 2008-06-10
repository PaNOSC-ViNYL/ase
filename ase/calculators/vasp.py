"""This module defines an ASE interface to VASP.

http://cms.mpi.univie.ac.at/vasp/
"""

# Copyright 2008 Jussi Enkovaara, CSC - Scientific Computing Ltd.

import os
from os.path import join, isfile, islink

import numpy as npy

from ase.data import chemical_symbols
from ase.data import atomic_numbers
from ase.units import Bohr, Hartree


class Vasp:
    """Class for doing Vasp calculations.

    The default parameters are very close to those that the VASP
    Fortran code would use.  These are the exceptions::

    Use the set_inp method to set extra INPUT parameters::

      calc.set_inp('nstep', 30)

    """
    def __init__(self, **kwargs):
        """Construct Vasp-calculator object.

        """

        # Parameters that can be set in INCAR. The values which are None
        # are not written and default parameters of VASP are used for them.
        self.input_parameters = {
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

        self.converged = False
        self.set(**kwargs)
        
    def write_poscar(self, atoms):
        """Writes the atomic positions to the VASP's POSCAR format"""

        poscar = open('POSCAR', 'w')
        poscar.write('POSCAR created by Atomic Simulation Environment\n')
        poscar.write('%.5f\n' % 1.0)
        cell = atoms.get_cell()
        for v in cell:
            poscar.write('%.14f %.14f %.14f\n' %  tuple(v))
        poscar.write('Cartesian\n')
        # Determine the number of atomic species
        symbols = {}
        for atom in atoms:
            symbol = atom.get_symbol()
            if not symbols.has_key(symbol):
                symbols[symbol] = 1
            else:
                symbols[symbol] += 1

        for symbol in symbols.keys():
            poscar.write('%i ' % symbols[symbol])
        poscar.write('\n')
        for symbol in symbols.keys():
            mask = [atom.symbol == symbol for atom in atoms]
            ind = npy.asarray(mask, bool)
            for atom in atoms[ind]:
                poscar.write('%.14f %.14f %.14f\n' %  tuple(atom.position))

        poscar.close()

    def write_incar(self):
        incar = open('INCAR', 'w')
        incar.close()
    

    def set(self, **kwargs):
        p = self.input_parameters
        
        for key in kwargs:
            if not p.has_key(key):
                raise TypeError('Unknown keyword argument:' + key)

        p.update(kwargs)
