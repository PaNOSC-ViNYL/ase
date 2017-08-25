"""This module defines an ASE interface to CRYSTAL14


http://www.crystal.unito.it/

daniele.selli@unimib.it
g.fazio3@campus.unimib.it

The file 'geom.f34' contains the input and output geometry
and it will be updated during the crystal calculations.

The keywords are given, for instance, as follows::

    DFT_SPIN ='YES',
    DFT_CORRELAT = 'B3',
    TOLDEE = 8,
    ANDERSON = 'YES',
    FMIXING = 95,
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

    def __init__(self, label='cry', atoms=None, kpts=None,
                 **kwargs):
        """Construct a crystal calculator.

        """
        from ase.dft.kpoints import monkhorst_pack

        if 'CRY_BASIS' in os.environ:
            basis_dir = os.environ['CRY_BASIS']
        else:
            basis_dir = './'

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
            if key.endswith('_'):
                outfile.write(key.rstrip('_').rsplit('_')[-1] +
                              ' = ' + str(value) + '{ \n')
            elif key.count('_empty') == 1:
                outfile.write(str(value) + ' \n')
            else:
                outfile.write(key.rsplit('_')[-1] + ' = ' + str(value) + ' \n')
            previous_key = key
        current_depth = key.rstrip('_').count('_')
        for my_backsclash in reversed(range(current_depth)):
            outfile.write(3 * my_backsclash * myspace + '} \n')

