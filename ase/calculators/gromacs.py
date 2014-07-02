"""This module defines an ASE interface to GROMACS.

http://www.gromacs.org/
"""

import os
from glob import glob
from os.path import join, isfile, islink

import numpy as np

from ase.data import atomic_numbers
from ase.data import chemical_symbols
from ase.io.gromacs import read_gromacs
from ase.calculators.calculator import FileIOCalculator, Parameters, \
    ReadError

def do_clean():
    """ remove gromacs backup files """
    files = glob('#*')
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass

class Gromacs(FileIOCalculator):
    """Class for doing GROMACS calculations.
    Before running a gromacs calculation you must prepare the input files
    separately (pdb2gmx and grompp for instance.)

    Input parameters for gromacs runs (the .mdp file)
    are given in self.params and can be set when initializing the calculator 
    or by method set_own.
    for Example 
    CALC_MM_RELAX = Gromacs()
    CALC_MM_RELAX.set_own_params('integrator', 'steep', 'use steepest descent')

    Run command line arguments for gromacs related programs:
    pdb2gmx, grompp, mdrun, g_energy, g_traj
    These can be given as:
    CALC_MM_RELAX = Gromacs()
    CALC_MM_RELAX.set_own_params_runs(
        'force_field','oplsaa','oplsaa force field')         


    """

    implemented_properties = ['energy', 'forces']
    command = 'mdrun < PREFIX.files > PREFIX.log'

    default_parameters = dict(
        index_filename = 'index.ndx',
        define = '-DFLEXIBLE',
        integrator = 'cg',
        nsteps = '10000',
        nstfout = '10',
        nstlog = '10',
        nstenergy = '10',
        nstlist = '10',
        ns_type = 'grid',
        pbc = 'xyz',
        rlist = '1.15',
        coulombtype = 'PME-Switch',
        rcoulomb = '0.8',
        vdwtype = 'shift',
        rvdw = '0.8',
        rvdw_switch = '0.75',
        DispCorr = 'Ener')

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gromacs', atoms=None,
                 do_qmmm = False, freeze_qm = False, clean=True, **kwargs):
        """Construct GROMACS-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'gromacs'.

        do_qmmm : boolean
            Is gromacs used as mm calculator for a qm/mm calculation

        freeze_qm : boolean
            In qm/mm are the qm atoms kept fixed at their initial positions

        clean :     boolean
            Remove gromacs backup files

        Examples
        ========



        """
        from glob import glob
        #self.label = label
        self.do_qmmm = do_qmmm
        self.freeze_qm = freeze_qm
        self.clean = clean
        self.params_doc = {}
        #add comments for gromacs input file
        self.params_doc['define'] = \
            'flexible/ rigid water'
        self.params_doc['integrator'] = \
            'md: molecular dynamics(Leapfrog), \n' + \
            '; md-vv: molecular dynamics(Velocity Verlet), \n' + \
            '; steep: steepest descent minimization, \n' + \
            '; cg: conjugate cradient minimization \n'

        self.positions = None
        self.atoms = None
        # storage for energy and forces
        self.energy = None
        self.forces = None
        
        self.set(**kwargs)

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.params_runs = {}
        self.params_runs['init_structure'] = self.label+'.pdb'
        self.params_runs['water'] = 'tip3p'
        self.params_runs['force_field'] = 'oplsaa'
        self.params_runs['extra_mdrun_parameters'] = ' -nt 1 '
        self.params_runs['extra_pdb2gmx_parameters'] = ' '
        self.params_runs['extra_grompp_parameters'] = ' '
        self.params_runs['extra_editconf_parameters'] = ' '
        self.params_runs['extra_genbox_parameters'] = ' '

        #these below are required by qm/mm
        self.topology_filename = self.label+'.top'
        self.force_field = self.params_runs.get('force_field')
        self.name = 'Gromacs'

        # clean up gromacs backups
        if self.clean: 
            do_clean()

        if self.do_qmmm:
            self.parameters['integrator'] = 'md'
            self.parameters['nsteps'] = '0'
            self.write_gromacs_qmmm_files()            



    def generate_g96file(self):
        """ from current coordinates (self.structure_file)
            write a structure file in .g96 format
        """
        import os.path
        #generate structure file in g96 format 
        write_gromos(self.label+'.g96', self.atoms)

    def run_editconf(self):
        """ run gromacs program editconf, typically to set a simulation box 
        writing to the input structure"""
        command = 'editconf' + ' '
        os.system(command + \
                      ' -f ' + self.label+'.g96' + \
                      ' -o ' + self.label+'.g96' + \
                      ' ' + self.params_runs.get('extra_editconf_parameters') +\
                      ' > /dev/null 2>&1')        

    def run_genbox(self):
        """ run gromacs program genbox, typically to solvate the system 
        writing to the input structure
        as extra parameter you need to define the file containing the solvent

        for instance
        CALC_MM_RELAX = Gromacs()
        CALC_MM_RELAX.set_own_params_runs(
            'extra_genbox_parameters','-cs spc216.gro')
        """
        command = 'genbox' + ' '
        os.system(command + \
                      ' -cp ' + self.label+'.g96' + \
                      ' -o ' + self.label+'.g96' + \
                      ' -p ' + self.label+'.top' + \
                      ' ' + self.params_runs.get('extra_genbox_parameters') +\
                      ' > /dev/null 2>&1')        



    def run(self):
        """ runs a gromacs-mdrun with the 
        current atom-configuration """
        from ase.io.gromos import read_gromos

        # clean up gromacs backups
        if self.clean: 
            do_clean()

        command = 'mdrun'
        if self.do_qmmm:
            os.system(command \
                          + ' -s ' + self.label + '.tpr' \
                          + ' -o ' + self.label + '.trr ' \
                          + ' -e ' + self.label + '.edr ' \
                          + ' -g ' + self.label + '.log -ffout ' \
                          + ' -rerun ' + self.label + '.g96 ' \
                      ' ' + self.params_runs.get('extra_mdrun_parameters') +\
                          + ' > mm.log 2>&1')
        else:
            os.system(command \
                          + ' -s ' + self.label + '.tpr ' \
                          + ' -o ' + self.label + '.trr ' \
                          + ' -e ' + self.label + '.edr ' \
                          + ' -g ' + self.label + '.log -ffout ' \
                          + ' -c ' + self.label + '.g96 ' \
                          + self.params_runs.get('extra_mdrun_parameters') \
                          + '  > MM.log 2>&1')
            atoms = read_gromos(self.label+'.g96')
            self.atoms = atoms.copy()


    def generate_topology_and_g96file(self):
        """ from coordinates (self.label.+'pdb')
            and gromacs run input file (self.label + '.mdp)
            generate topology (self.label+'top')
            and structure file in .g96 format (self.label + '.g96')
        """
        import os.path
        from ase.io.gromos import read_gromos
        #generate structure and topology files 
        # In case of predefinded topology file this is not done
        command = 'pdb2gmx' + ' '
        os.system(command + \
                      ' -f ' + self.params_runs.get('init_structure') + \
                      ' -o ' + self.label+'.g96' + \
                      ' -p ' + self.label+'.top' + \
                      ' -ff ' + self.params_runs.get('force_field') + \
                      ' -water ' + self.params_runs.get('water') + \
                      ' ' + self.params_runs.get('extra_pdb2gmx_parameters') +\
                      ' > /dev/null 2>&1')
#                      ' > debug.log 2>&1')

#        print command + \
#                      ' -f ' + self.params_runs.get('init_structure') + \
#                      ' -o ' + self.label+'.g96' + \
#                      ' -p ' + self.label+'.top' + \
#                      ' -ff ' + self.params_runs.get('force_field') + \
#                      ' -water ' + self.params_runs.get('water') + \
#                      ' ' + self.params_runs.get('extra_pdb2gmx_parameters') +\
#                      ' > /dev/null 2>&1'
        atoms = read_gromos(self.label+'.g96')
        self.atoms = atoms.copy()

    def generate_gromacs_run_file(self):
        """ Generates input file for a gromacs mdrun
        based on structure file and topology file
        resulting file is self.label + '.tpr
        """

        import os.path
        #generate gromacs run input file (gromacs.tpr)
        try:
            os.remove(self.label + '.tpr')
        except:
            pass
        command = 'grompp ' 
        os.system(command + \
                      ' -f ' + self.label + '.mdp' + \
                      ' -c ' + self.label + '.g96' + \
                      ' -p ' + self.label + '.top' + \
                      ' -o ' + self.label + '.tpr -maxwarn 100' + \
                      ' ' + self.params_runs.get('extra_grompp_parameters') +\
                      ' > /dev/null 2>&1')

        print command + \
                      ' -f ' + self.label + '.mdp' + \
                      ' -c ' + self.label + '.g96' + \
                      ' -p ' + self.label + '.top' + \
                      ' -o ' + self.label + '.tpr -maxwarn 100' + \
                      ' ' + self.params_runs.get('extra_grompp_parameters') +\
                      ' > /dev/null 2>&1'



    def write_gromacs_qmmm_files(self):
        """write input files for gromacs force and energy calculations """
        filename = 'inputGenergy.txt'
        output = open(filename,'w')
        output.write('Potential  \n')
        output.write('   \n')
        output.write('   \n')
        output.close()

        filename = 'inputGtraj.txt'
        output = open(filename, 'w')
        output.write('System  \n')
        output.write('   \n')
        output.write('   \n')
        output.close()


    def set(self, **kwargs):
        """ Setting values for the parameters of the gromacs calculator """
        for key in kwargs:
                self.parameter[key] = kwargs[key]

    def set_own_params(self, key, value, docstring=""):
        """Set own gromacs parameter with doc strings."""
        self.parameters[key] = value
        self.params_doc[key] = docstring

    def set_own_params_runs(self, key, value):
        """Set own gromacs parameter for program parameters 
        Add spaces to avoid errors """
        self.params_runs[key] = ' ' + value + ' '

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        return system_changes

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms=None, properties=None, system_changes=None):
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        #print self.parameters
        f = open(self.label+'.mdp', 'w')
        for key, val in self.parameters.items():
            #print "k v", key, val
            #print self.params_doc
            if val is not None:
                if (self.params_doc.get(key) == None):
                    docstring = ''
                else:
                    docstring = self.params_doc[key]
                f.write('%-35s = %s ; %s\n' \
                            % (key, val, ';'+docstring))
        f.close()

    def read_atoms(self):
        """ read atoms from file """
        from ase.io.gromos import read_gromos
        self.atoms = read_gromos(self.label+'.g96')

    def read(self, label):
        """Read results from Gromacs's text-output file."""
        FileIOCalculator.read(self, label)
        filename = self.label + '.txt'
        if not os.path.isfile(filename):
            raise ReadError

        self.atoms = read_abinit(self.label + '.in')
        self.parameters = Parameters.read(self.label + '.ase')

        self.initialize(self.atoms)
        self.read_results()

    def read_results(self):
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        
        if ('error' in text or
            'was not enough scf cycles to converge' in text):
            raise ReadError

        for line in iter(text.split('\n')):
            if line.rfind('natom  ') > -1:
                natoms = int(line.split()[-1])

        lines = iter(text.split('\n'))
        # Stress:
        # Printed in the output in the following format [Hartree/Bohr^3]:
        # sigma(1 1)=  4.02063464E-04  sigma(3 2)=  0.00000000E+00
        # sigma(2 2)=  4.02063464E-04  sigma(3 1)=  0.00000000E+00
        # sigma(3 3)=  4.02063464E-04  sigma(2 1)=  0.00000000E+00
        for line in lines:
            if line.rfind(
                'cartesian components of stress tensor (hartree/bohr^3)') > -1:
                stress = np.empty(6)
                for i in range(3):
                    entries = lines.next().split()
                    stress[i] = float(entries[2])
                    stress[i + 3] = float(entries[5])
                self.results['stress'] = stress * Hartree / Bohr**3
                break
        else:
            raise RuntimeError

        # Energy [Hartree]:
        # Warning: Etotal could mean both electronic energy and free energy!
        for line in iter(text.split('\n')):
            if line.rfind('>>>>> internal e=') > -1:
                etotal = float(line.split('=')[-1])*Hartree
                for line1 in iter(text.split('\n')):
                    if line1.rfind('>>>>>>>>> etotal=') > -1:
                        efree = float(line1.split('=')[-1])*Hartree
                        break
                else:
                    raise RuntimeError
                break
        else:
            for line2 in iter(text.split('\n')):
                if line2.rfind('>>>>>>>>> etotal=') > -1:
                    etotal = float(line2.split('=')[-1])*Hartree
                    efree = etotal
                    break
            else:
                raise RuntimeError

        # Energy extrapolated to zero Kelvin:
        self.results['energy'] = (etotal + efree) / 2
        self.results['free_energy'] = efree

        # Forces:
        for line in lines:
            if line.rfind('cartesian forces (ev/angstrom) at end:') > -1:
                forces = []
                for i in range(natoms):
                    forces.append(np.array(
                            [float(f) for f in lines.next().split()[1:]]))
                self.results['forces'] = np.array(forces)
                break
        else:
            raise RuntimeError
        #
        self.niter = self.read_number_of_iterations()


    def initialize(self, atoms):
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)

