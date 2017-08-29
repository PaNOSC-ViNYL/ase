"""This module defines an ASE interface to CRYSTAL14

http://www.crystal.unito.it/

daniele.selli@unimib.it
g.fazio3@campus.unimib.it

The file 'geom.f34' contains the input and output geometry
and it will be updated during the crystal calculations.

The keywords are given, for instance, as follows::

    DFT = True,
    DFT_CORRELAT = 'PBE',
    DFT_EXCHANGE = 'PBE',
    DFT_MET = 'B3LYP',
    DFT_POL = 'SPIN',
    TOLDEE = 8,
    ANDERSON = 'Yes',
    FMIXING = 95,
    ... # put other examples

"""

import os

import numpy as np

from ase.calculators.calculator import FileIOCalculator


class CRYSTAL(FileIOCalculator):
    """ A crystal calculator with ase-FileIOCalculator nomenclature
    """
    if 'CRY_COMMAND' in os.environ:
        command = os.environ['CRY_COMMAND'] + ' < INPUT > OUTPUT 2>&1'
    else:
        command = 'crystal < INPUT > OUTPUT'

    implemented_properties = ['energy', 'forces', 'stress', 'charges']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='cry', atoms=None, kpts=None,
                 **kwargs):
        """Construct a crystal calculator.

        """
  #      # call crystal only to run a single point calculation
  #      # [PUT HERE DEFAULT PARAMETERS]
        self.default_parameters = dict(
            xc='hf'
            spin=False,
            guess=False,
            otherkey=[])

        self.lines = None
        self.atoms = None
        self.atoms_input = None
        self.outfilename = 'cry.out'

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms,
                                  **kwargs)
        self.kpts = kpts

    def write_crystal_in(self, filename):
        """ Write the input file for the crystal calculation.
            Geometry is taken always from the file 'fort.34'
        """

        # write BLOCK 1 of crystal input (only SP with gradients)
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

        # write BLOCK 2 of crystal input from file (basis sets)
        try:
            basisfile = open(os.path.join(self.directory, 'basis'))
        except:
            raise RuntimeError('"basis" file not found. \
Create a "basis" file with CRYSTAL basis set.')
        basis = basisfile.readlines()
        for line in basis:
            outfile.write(line)

        # write BLOCK 3 according to parameters set as input
        newline = '\n'
        # ----- write hamiltonian
        p = self.parameters
        if p.xc == 'hf':
            if p.spin:
                outfile.write('UHF \n')
            else:
                outfile.write('RHF \n')
        elif p.xc == 'mp2':
            outfile.write('MP2 \n')
        else:
            outfile.write('DFT \n')
        # Standalone keyword and LDA are given by a single string.
            if isinstance(p.xc, str):
                xc = {'LDA': 'EXCHANGE\nLDA\nCORRELAT\nVWN',
                      'PBE': 'PBEXC'}.get(p.xc, p.xc)
                outfile.write(xc+'\n')
        # Custom xc functional are given by a tuple of string
            else:
                x, c = p.xc
                outfile.write('EXCHANGE \n')
                outfile.write(x + ' \n')
                outfile.write('CORRELAT \n')
                outfile.write(c + ' \n')
            if p.spin:
                outfile.write('SPIN \n')
            outfile.write('END \n')
        # When guess=True, wf is read.
        if p.guess:
            # wf will be always there after 2nd step.
            if os.path.isfile('fort.20'):
                outfile.write('GUESSP \n')

        # ----- write other CRYSTAL keywords
        # ----- in the list otherkey = ['ANDERSON', ...] .

        for keyword in p.otherkey:
            if isinstance(keyword, str):
                outfile.write(keyword + '\n')
            else:
                for key in keyword:
                    outfile.write(key + '\n')

        if any(
        outfile.write('END \n')

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
        """ all results are read from OUTPUT file
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        from ase.units import Hartree, Bohr

        myfile = open(os.path.join(self.directory, 'OUTPUT'), 'r')
        self.lines = myfile.readlines()
        myfile.close()

        self.atoms = self.atoms_input
        energy = 0.0
        forces = None
        # Energy line index
        for iline, line in enumerate(self.lines):
            estring = 'OPT END'
            if line.find(estring) >= 0:
                index_energy = iline
                break
        try:
            energy = float(self.lines[index_energy].split()[7]) * Hartree
        except:
            raise RuntimeError('Problem in reading energy')
        self.results['energy'] = energy
        # Force line indexes
        fstring = 'CARTESIAN FORCES'
        fstring_end = 'RESULTANT FORCE'
        for iline, line in enumerate(self.lines):
            if line.find(fstring) >= 0:
                index_force_begin = iline + 2
            if line.find(fstring_end) >= 0:
                index_force_end = iline - 1
                break
        try:
            gradients = []
            for j in range(index_force_begin, index_force_end):
                word = self.lines[j].split()
                gradients.append([float(word[k+2]) for k in range(0, 3)])
            forces = np.array(gradients) * Hartree / Bohr
        except:
            raise RuntimeError('Problem in reading forces')

        self.results['forces'] = forces

        # stress stuff begins
        sstring = 'STRESS TENSOR'
        have_stress = False
        stress = list()
        for iline, line in enumerate(self.lines):
            if sstring in line:
                have_stress = True
                start = iline + 4
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

        """Get partial charges on atoms
            in case we cannot find charges they are set to None
        """
        qm_charges = []
        for n, line in enumerate(self.lines):
            if ('TOTAL ATOMIC CHARGE' in line):
                chargestart = n + 1
        lines1 = self.lines[chargestart:(chargestart + len(self.atoms)/6 + 1)]
        i = 0
        atomnum = self.atoms.get_atomic_numbers()
        for line in lines1:
            words = line.split()
            for word in words:
                qm_charges.append(-float(word)+atomnum[i])
                i = i + 1
        charges = np.array(qm_charges)
        self.results['charges'] = charges

        # calculation was carried out with
        # atoms written in write_input
        os.remove(os.path.join(self.directory, 'OUTPUT' ))
