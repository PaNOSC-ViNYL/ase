"""This module defines an ASE interface to CRYSTAL14


http://www.crystal.unito.it/

daniele.selli@unimib.it
g.fazio3@campus.unimib.it

The file 'geom.f34' contains the input and output geometry
and it will be updated during the crystal calculations.

The keywords are given, for instance, as follows::

    BLOCK3_SPIN ='YES',
    BLOCK3_TOLDEE = 8,
    BLOCK3_ANDERSON = 'YES',
    BLOCK3_FMIXING = 95,
    ... # put other examples

"""

import os

import numpy as np

from ase.calculators.calculator import FileIOCalculator, kpts2mp


class crystal(FileIOCalculator):
    """ A crystal calculator with ase-FileIOCalculator nomenclature
    """
    if 'CRY_COMMAND' in os.environ:
        command = os.environ['CRY_COMMAND'] + ' > OUTPUT'
    else:
        command = 'crystal > OUTPUT'

    implemented_properties = ['energy', 'forces']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='cry', atoms=None, kpts=None,
                 **kwargs):
        """Construct a crystal calculator.

        """
        from ase.dft.kpoints import monkhorst_pack

        if 'CRY_BASIS' in os.environ:
            basis_dir = os.environ['CRY_BASIS']
        else:
            slako_dir = './'

        # call crystal only to run a single point calculation
        # [PUT HERE DEFAULT PARAMETERS]
  #      self.default_parameters = dict(
  #          Hamiltonian_='DFTB',
  #          Driver_='ConjugateGradient',
  #          Driver_MaxForceComponent='1E-4',
  #          Driver_MaxSteps=0,
  #          Hamiltonian_SlaterKosterFiles_='Type2FileNames',
  #          Hamiltonian_SlaterKosterFiles_Prefix=slako_dir,
  #          Hamiltonian_SlaterKosterFiles_Separator='"-"',
  #          Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
  #          Hamiltonian_MaxAngularMomentum_='')

        self.lines = None
        self.atoms = None
        self.atoms_input = None
        self.outfilename = 'cry.out'

        FileIOCalculator.__init__(self, label, atoms,
                                  **kwargs)
        self.kpts = kpts

    def write_crystal_in(self, filename):
        """ Write the input file for the crystal calculation.
            Geometry is taken always from the file 'fort.34'
        """

        outfile = open(filename, 'w')
        outfile.write('Single point + Gradient crystal calculation \n')
        outfile.write('EXTERNAL \n')
        outfile.write('OPTGEOM \n')
        outfile.write('FINALRUN \n')
        outfile.write('0 \n')
        outfile.write('MAXCYCLE \n')
        outfile.write('1 \n')
        outfile.write('END \n')
        outfile.write('END \n')


        # --------MAIN KEYWORDS-------
        previous_key = 'dummy_'
        myspace = ' '
        for key, value in sorted(self.parameters.items()):
            current_depth = key.rstrip('_').count('_')
            previous_depth = previous_key.rstrip('_').count('_')
            for my_backsclash in reversed(
                    range(previous_depth - current_depth)):
                outfile.write(3 * (1 + my_backsclash) * myspace + '} \n')
            outfile.write(3 * current_depth * myspace)
            if key.endswith('_'):
                outfile.write(key.rstrip('_').rsplit('_')[-1] +
                              ' = ' + str(value) + '{ \n')
            elif key.count('_empty') == 1:
                outfile.write(str(value) + ' \n')
            else:
                outfile.write(key.rsplit('_')[-1] + ' = ' + str(value) + ' \n')
            if self.pcpot is not None and ('DFTB' in str(value)):
                outfile.write('   ElectricField = { \n')
                outfile.write('      PointCharges = { \n')
                outfile.write(
                    '         CoordsAndCharges [Angstrom] = DirectRead { \n')
                outfile.write('            Records = ' +
                              str(len(self.pcpot.mmcharges)) + ' \n')
                outfile.write(
                    '            File = "dftb_external_charges.dat" \n')
                outfile.write('         } \n')
                outfile.write('      } \n')
                outfile.write('   } \n')
            previous_key = key
        current_depth = key.rstrip('_').count('_')
        for my_backsclash in reversed(range(current_depth)):
            outfile.write(3 * my_backsclash * myspace + '} \n')
        # output to 'results.tag' file (which has proper formatting)
        outfile.write('Options { \n')
        outfile.write('   WriteResultsTag = Yes  \n')
        outfile.write('} \n')
        outfile.write('ParserOptions { \n')
        outfile.write('   IgnoreUnprocessedNodes = Yes  \n')
        outfile.write('} \n')

        outfile.close()

    def write_input(self, atoms, properties=None, system_changes=None):
        from ase.io import write
        FileIOCalculator.write_input(
            self, atoms, properties, system_changes)
        self.write_crystal_in(os.path.join(self.directory, 'INPUT'))
        write(os.path.join(self.directory, 'fort.34'), atoms)
        # self.atoms is none until results are read out,
        # then it is set to the ones at writing input
        self.atoms_input = atoms
        self.atoms = None

    def read_results(self):
        """ all results are read from results.tag file
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        from ase.units import Hartree, Bohr

        myfile = open(os.path.join(self.directory, 'results.tag'), 'r')
        self.lines = myfile.readlines()
        myfile.close()

        self.atoms = self.atoms_input
        charges = self.read_charges()
        self.results['charges'] = charges
        energy = 0.0
        forces = None
        energy = self.read_energy()
        forces = self.read_forces()
        self.results['energy'] = energy
        self.results['forces'] = forces
        self.mmpositions = None
        # stress stuff begins
        sstring = 'stress'
        have_stress = False
        stress = list()
        for iline, line in enumerate(self.lines):
            if sstring in line:
                have_stress = True
                start = iline + 1
                end = start + 3
                for i in range(start, end):
                    cell = [float(x) for x in self.lines[i].split()]
                    stress.append(cell)
        if have_stress:
            stress = -np.array(stress) * Hartree / Bohr**3
        elif not have_stress:
            stress = np.zeros((3, 3))
        self.results['stress'] = stress
        # stress stuff ends

        # calculation was carried out with atoms written in write_input
        os.remove(os.path.join(self.directory, 'results.tag'))

    def read_energy(self):
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree
        # Energy line index
        for iline, line in enumerate(self.lines):
            estring = 'total_energy'
            if line.find(estring) >= 0:
                index_energy = iline + 1
                break
        try:
            return float(self.lines[index_energy].split()[0]) * Hartree
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""
        from ase.units import Hartree, Bohr

        # Force line indexes
        for iline, line in enumerate(self.lines):
            fstring = 'forces   '
            if line.find(fstring) >= 0:
                index_force_begin = iline + 1
                line1 = line.replace(':', ',')
                index_force_end = iline + 1 + \
                    int(line1.split(',')[-1])
                break
        try:
            gradients = []
            for j in range(index_force_begin, index_force_end):
                word = self.lines[j].split()
                gradients.append([float(word[k]) for k in range(0, 3)])
            return np.array(gradients) * Hartree / Bohr
        except:
            raise RuntimeError('Problem in reading forces')

