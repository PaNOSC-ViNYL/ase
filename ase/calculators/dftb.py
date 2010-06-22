"""This module defines an ASE interface to DftbPlus

http://http://www.dftb-plus.info//
http://www.dftb.org/

-Markus Kaukonen markus.kaukonen@iki.fi
"""


import numpy as np
import os, string

#from ase.data import chemical_symbols
from ase.units import Hartree, Bohr


class Dftb:
    """Class for doing DFTB+ calculations.
    """
    def __init__(self, label='dftb', write_dftb=False,
                 charge=0.0, include_dispersion=False,
                 do_spin_polarized=False, 
                 unpaired_electrons=0.0,
                 fermi_temperature=0.0, scc=False):
        """Construct DFTB-calculator object.


        For example:
        calc = Dftb(label='dftb',write_dftb=True,include_dispersion=True )

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.txt, ...).
            Default is 'dftb'.
        write_dftb: boolean
            True: a minimal input file (name of which is always 'dftb_in.hsd')
            is written based on values given here.
            False: input file for dftb+ is not written. User must have
            generated file 'dftb_in.hsd' in the working directory.
            Use write_dftb=False to use your own 'dftb_in.hsd'-file.
        charge: float
            Total charge of the system.
        include_dispersion: boolean
            True: Default dispersion parameters are written in the 
            file 'dftb_in.hsd' (requires that also write_dftb_input_file==True)
            False: dispersion parameters are not written here.
        do_spin_polarized: boolean
            True: Spin polarized calculation
            (requires that also write_dftb_input_file==True)
            False: Spin unpolarized calculation
        unpaired_electrons: float
            Number of spin unpaired electrons in the system.
            Relevant only if do_spin_polarized==True
        fermi_temperature: float
            Fermi temperature for electrons.
        scc: boolean
            True: Do charge self consistent dftb+
            False: No SCC, charges on atoms are not iterated
            
        Input file for DFTB+ file is 'dftb_in.hsd'. Keywords in it are
        written here or read from an existing file. The atom positions
        in file 'dftb_in.hsd' are updated during ASE geometry
        optimization.
        """


        if not(write_dftb):
            if os.path.isfile('dftb_in.hsd'):
                f = open('dftb_in.hsd')
            else:
                print 'Input file for DFTB+ dftb_in.hsd is missing'
                raise RuntimeError, \
                    'Provide it or set write_dftb=True '
            #lines = f.readlines()
            f.close()

        self.label = label
        self.write_dftb = write_dftb
        self.charge = charge
        self.include_dispersion = include_dispersion
        self.do_spin_polarized = do_spin_polarized
        self.unpaired_electrons = unpaired_electrons
        #if (do_spin_polarized):
        #    print 'Sorry, generation of file "dftb_in.hsd" with spin '
        #    print 'polarization is not inplemented for DFTB+'
        #    print 'Set write_dftb=False and'
        #    raise RuntimeError, \
        #        'Generate file "dftb_in.hsd by hand"'
        
        self.etotal = 0.0
        self.cell = None
        self.fermi_temperature = fermi_temperature
        self.scc = scc 
        self.converged = False

        #dftb has no stress
        self.stress = np.empty((3, 3))
        

    def update(self, atoms):
        """Energy and forces are calculated when atoms have moved
        by calling self.calculate
        """
        if (not self.converged or
            len(self.typenumber) != len(atoms)):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        #from ase.io.dftb import read_dftb

        atomtypes = atoms.get_chemical_symbols()
        self.allspecies = []
        self.typenumber = []
        self.max_angular_momentum = []
        for species in (atomtypes):
            if species not in self.allspecies:
                self.allspecies.append(species)
        for species in (atomtypes):
            myindex = 1 + self.allspecies.index(species)
            self.typenumber.append(myindex)
        for i in self.allspecies:
            if i == 'H':
                self.max_angular_momentum.append('s')
            elif i in ['C','N','O']:
                self.max_angular_momentum.append('p')
            elif i in ['Si','S','Fe','Ni']:
                self.max_angular_momentum.append('d')
            else:
                print 'anglular momentum is not imlemented in ASE-DFTB+'
                print 'for species '+i
                raise RuntimeError('Use option write_dftb=False')
        self.converged = False
        
        #write DFTB input file if desired 
        self.positions = atoms.get_positions().copy()
        if self.write_dftb:
            self.write_dftb_input_file(atoms)
        
        
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.etotal

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):
        """Total DFTB energy is calculated (to file 'energy'
        also forces are calculated (to file 'gradient')
        """
        
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.write_dftb:
            self.write_dftb_input_file(atoms)
        else:
            #write current coordinates to file 'dftb_in.hsd' for DFTB+
            self.change_atom_positions_dftb(atoms)


        #DFTB energy&forces calculation
        if (os.environ.has_key('DFTB_COMMAND') and
            (os.environ.has_key('DFTB_PREFIX'))):
            dftb = os.environ['DFTB_COMMAND']
            exitcode = os.system(dftb+ '> dftb.output' )
        elif not(os.environ.has_key('DFTB_COMMAND')):
            raise RuntimeError('Please set DFTB_COMMAND environment variable')
        elif not(os.environ.has_key('DFTB_PREFIX')):
            print 'Path for DFTB+ slater koster files is missing'    
            raise RuntimeError('Please set DFTB_PREFIX environment variable')
        else:
            pass
        if exitcode != 0:
            raise RuntimeError('Dftb exited with exit code: %d.  ' % exitcode)

        self.read_energy()

        #DFTB atomic forces calculation, to be read in file detailed.out
        #os.system(self.dftb_program_forces +'> output.forces.dummy')

        self.read_forces(atoms)

        self.converged = True

        
    def read_energy(self):
        """Read Energy from DFTB energy file."""
        text = open('dftb.output', 'r').read().lower()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            if 'total energy' in line:
                energy_tmp = float(line.split()[2])
                #print 'energy_tmp', energy_tmp
        self.etotal = energy_tmp * Hartree


    def read_forces(self, atoms):
        """Read Forces from DFTB+ detailed.out output file"""

        myfile = open('detailed.out','r')
        line = myfile.readline()
        line = myfile.readline()
        tmpforces = np.array([[0, 0, 0]])
        while line:
            if 'Total Forces' in line:
                for i, dummy in enumerate(atoms):
                    line = myfile.readline()
                    line2 = string.replace(line, 'D', 'E')
                    tmp = np.array\
                        ([[float(f) for f in line2.split()[0:3]]])
                    tmpforces = np.concatenate((tmpforces, tmp))  
            line = myfile.readline()
            

        self.forces = (np.delete(tmpforces, np.s_[0:1], axis=0))*Hartree/Bohr

        #print 'forces', self.forces

    def read(self):
        """Dummy stress for dftb"""
        self.stress = np.empty((3, 3))

    def write_dftb_input_file(self, atoms):
        """Write input parameters to DFTB+ input file 'dftb_in.hsd'."""
        import sys
        fh = open('dftb_in.hsd', 'w')
        # geometry
        fh.write('Geometry = {\n')
        fh.write('TypeNames = {')
        for i in self.allspecies:
            fh.write(' "'+i+'"')
        fh.write(' }\n')
        fh.write('TypesAndCoordinates [Angstrom] = {\n')
        self.positions = atoms.get_positions().copy()
        for i, pos in zip(self.typenumber, self.positions):
            fh.write('%6d ' % (i))
            fh.write('%20.14f %20.14f %20.14f' %  tuple(pos))
            fh.write('\n')
        fh.write(' }\n')

        #is it periodic and when is write also lattice vectors
        periodic = atoms.get_pbc().any()
        if periodic:
            fh.write('Periodic = Yes\n')
        else:
            fh.write('Periodic = No\n')
        if periodic:
            cell = atoms.get_cell().copy()
            fh.write('LatticeVectors [Angstrom] = {\n')
            for v in cell:
                fh.write('%20.14f %20.14f %20.14f \n' %  tuple(v))
            fh.write('  }\n')

        #end of geometry session
        fh.write('}\n')

        #zero step CG relaxation to get forces and energies
        # these are dummies because ASE takes care of these things
        fh.write('\n') 
        fh.write('Driver = ConjugateGradient {\n')
        fh.write('MovedAtoms = Range { 1 -1 }\n')
        fh.write('  MaxForceComponent = 1.0e-4\n')
        fh.write('  MaxSteps = 0\n')
        fh.write('  OutputPrefix = '+self.label+ '\n')
        fh.write('}\n')

        #Hamiltonian
        fh.write('\n') 
        fh.write('Hamiltonian = DFTB { # DFTB Hamiltonian\n')
        if (self.scc):
            fh.write('  SCC = Yes')
            fh.write(' # Use self consistent charges\n')               
            fh.write('  SCCTolerance = 1.0e-5')
            fh.write(' # Tolerance for charge consistence\n')            
            fh.write('  MaxSCCIterations = 1000')
            fh.write(' # Nr. of maximal SCC iterations\n')          
            fh.write('  Mixer = Broyden {') 
            fh.write(' # Broyden mixer for charge mixing\n')          
            fh.write('    MixingParameter = 0.2')  
            fh.write(' # Mixing parameter\n')
            fh.write('  }\n')
        else:
            fh.write('  SCC = No # NO self consistent charges\n')
        fh.write('  SlaterKosterFiles = Type2FileNames {')
        fh.write(' # File names from two atom type names\n')
        sk_prefix = os.environ['DFTB_PREFIX']
        fh.write('    Prefix = "'+sk_prefix+'"')
        fh.write(' # Path as prefix\n')
        fh.write('    Separator = "-"')
        fh.write(' # Dash between type names\n')
        fh.write('    Suffix = ".skf"')
        fh.write(' # Suffix after second type name\n')
        fh.write('  }\n')
        fh.write('  MaxAngularMomentum = {')
        fh.write(' # Maximal l-value of the various species\n')
        for i, j in zip(self.allspecies, self.max_angular_momentum):
            fh.write('   '+i+' = "'+j+'"\n')
        fh.write('  }\n')
        fh.write('  Charge = ')
        fh.write('%10.6f' % (self.charge))
        fh.write(' # System neutral\n')
        if self.do_spin_polarized:
            fh.write('  SpinPolarisation = Colinear {\n') 
            fh.write('  UnpairedElectrons = '+str(self.unpaired_electrons)+'\n')
            fh.write('  } \n')
            fh.write('  SpinConstants = {\n') 
            for i in self.allspecies:
                if i == 'H':
                    fh.write('   H={\n') 
                    fh.write('    # Wss\n') 
                    fh.write('    -0.072\n') 
                    fh.write('    }\n')
                elif i == 'C': 
                    fh.write('   C={\n') 
                    fh.write('    # Wss Wsp Wps Wpp\n') 
                    fh.write('    -0.031 -0.025 -0.025 -0.023\n') 
                    fh.write('    }\n') 
                elif i == 'N': 
                    fh.write('   N={\n') 
                    fh.write('    # Wss Wsp Wps Wpp\n')
                    fh.write('    -0.033 -0.027 -0.027 -0.026\n') 
                    fh.write('     }\n') 
                elif i == 'O':
                    fh.write('   O={\n') 
                    fh.write('    # Wss Wsp Wps Wpp\n') 
                    fh.write('    -0.035 -0.030 -0.030 -0.028\n') 
                    fh.write('     }\n') 
                elif (i == 'Si' or i == 'SI'):
                    fh.write('   Si={\n') 
                    fh.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n')
                    fh.write('    -0.020 -0.015 0.000 -0.015 -0.014 0.000 0.002 0.002 -0.032\n')
                    fh.write('    }\n')
                elif (i == 'S'):
                    fh.write('   S={\n') 
                    fh.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    fh.write('    -0.021 -0.017 0.000 -0.017 -0.016 0.000 0.000 0.000 -0.080\n')
                    fh.write('    }\n')
                elif (i == 'Fe' or i == 'FE'):
                    fh.write('   Fe={\n') 
                    fh.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    fh.write('    -0.016 -0.012 -0.003 -0.012 -0.029 -0.001 -0.003 -0.001 -0.015\n')
                    fh.write('    }\n')
                elif (i == 'Ni' or i == 'NI'):
                    fh.write('   Ni={\n') 
                    fh.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    fh.write('    -0.016 -0.012 -0.003 -0.012 -0.022 -0.001 -0.003 -0.001 -0.018\n')
                    fh.write('    }\n')
                else:
                    print 'Missing spin polarisation parameters for species'+i
                    raise RuntimeError, \
                        'Run spin unpolarised calculation'
                fh.write('   }\n') 
        else:
            fh.write('  SpinPolarisation = {}')
            fh.write(' # No spin polarisation\n')
        fh.write('  Filling = Fermi {\n')
        fh.write('    Temperature [Kelvin] = ')
        fh.write('%10.6f\n' % (self.fermi_temperature))
        fh.write('  }\n')
        if periodic:
            fh.write('# gamma only\n')
            fh.write('KPointsAndWeights = { \n')
            fh.write('  0.000000    0.000000    0.000000       1.000000 \n')
            fh.write('}\n')

        #Dispersion parameters
        if (self.include_dispersion):
            fh.write('Dispersion = SlaterKirkwood {\n')
            fh.write(' PolarRadiusCharge = HybridDependentPol {\n')
            fh.write('\n')
            fh.write('  C={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.8\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('      1.382 1.382 1.382 1.064 1.064 1.064 3.8 3.8 3.8 3.8 3.8 3.8 2.5\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write('\n')
            fh.write('  N={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.8\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('      1.030 1.030 1.090 1.090 1.090 1.090 3.8 3.8 3.8 3.8 3.8 3.8 2.82\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write('\n')
            fh.write('  O={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.8\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('    # All polarisabilities and radii set the same\n')
            fh.write('      0.560 0.560 0.000 0.000 0.000 0.000 3.8 3.8 3.8 3.8 3.8 3.8 3.15\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write('\n')


            fh.write('  H={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.4\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('    # Different polarisabilities depending on the hybridisation\n')
            fh.write('      0.386 0.396 0.000 0.000 0.000 0.000 3.5 3.5 3.5 3.5 3.5 3.5 0.8\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write('\n')

            fh.write('  P={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.9\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('    # Different polarisabilities depending on the hybridisation\n')
            fh.write('      1.600 1.600 1.600 1.600 1.600 1.600 4.7 4.7 4.7 4.7 4.7 4.7 4.50\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write('\n')

            fh.write('  S={\n')
            fh.write('    CovalentRadius [Angstrom] = 0.9\n')
            fh.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            fh.write('    # Different polarisabilities depending on the hybridisation\n')
            fh.write('      3.000 3.000 3.000 3.000 3.000 3.000 4.7 4.7 4.7 4.7 4.7 4.7 4.80\n')
            fh.write('    }\n')
            fh.write('  }\n')
            fh.write(' }\n')        
            fh.write('}\n') 
        fh.write('}\n') 
        fh.write('Options = {}\n')
        fh.write('ParserOptions = {\n')
        fh.write('  ParserVersion = 3\n')
        fh.write('}\n')
 
        fh.close()

    def change_atom_positions_dftb(self, atoms):
        """Write coordinates in DFTB+ input file dftb_in.hsd
        """

        filename = 'dftb_in.hsd'
        if isinstance(filename, str):
            myfile = open(filename)

        lines = myfile.readlines()

        if type(filename) == str:
            myfile.close()

        if isinstance(filename, str):
            myfile = open(filename, 'w')
        else: # Assume it's a 'file-like object'
            myfile = filename

        coord = atoms.get_positions()

        start_writing_coords = False
        stop_writing_coords = False
        i = 0
        for line in lines:
            if ('TypesAndCoordinates' in line):
                start_writing_coords = True
            if (start_writing_coords and not(stop_writing_coords)):
                if ('}' in line):
                    stop_writing_coords = True
            if (start_writing_coords and not(stop_writing_coords)and 
                not ('TypesAndCoordinates' in line)):
                atom_type_index = line.split()[0]
                myfile.write('%6s  %20.14f  %20.14f  %20.14f\n'
                        % (atom_type_index,coord[i][0],coord[i][1],coord[i][2]))
                i = i + 1
            else:
                myfile.write(line)
