"""This module defines an ASE interface to MOPAC."""
import os
import numpy as np

from ase.units import kcal, mol
from ase.calculators.general import Calculator

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError



class MOPAC(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole', 'magmom']
    command = 'mopac PREFIX.mop'

    default_parameters = dict(
        method='PM6',
        task='gradient',
        raw='')  # additional outside of dft block control string
    str_keys = ['functional', 'job_type', 'command']
    int_keys = ['restart', 'spin']
    bool_keys = ['OPT']
    float_keys = ['RELSCF']
    self.set(restart=0,
             spin=0,
             OPT=False,
             functional='PM6',
             job_type='NOANCI 1SCF GRADIENTS AUX(0,PRECISION=9)',
             RELSCF=0.0001)

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='mopac', atoms=None, **kwargs):
        """Construct MOPAC-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.initial_magmoms = atoms.get_initial_magnetic_moments().tolist()
        p.write(self.label + '.ase')
        del p['initial_magmoms']
        f = open(self.label + '.nw', 'w')
        if p.charge is not None:
            f.write('charge %s\n' % p.charge)
        
        # start the input
        mopac_input = ''

        #write functional and job_type
        for key in 'functional', 'job_type':
            if self.str_params[key] != None:
                mopac_input += self.str_params[key] + ' '
                
        if self.float_params['RELSCF'] != None:
            mopac_input += 'RELSCF=' + str(self.float_params['RELSCF']) + ' '
            
        #write charge
        charge = sum(atoms.get_initial_charges())
        if charge != 0:
            mopac_input += 'CHARGE=%i ' % (charge)
        
        #write spin
        spin = self.int_params['spin']
        if spin == 1.:
            mopac_input += 'DOUBLET '
        elif spin == 2.:
            mopac_input += 'TRIPLET '

        #input down
        mopac_input += '\n'
        mopac_input += 'Title: ASE job\n\n'

        f = 1
        # write coordinates
        for iat in range(len(atoms)):
            atom = atoms[iat]
            xyz = atom.position
            mopac_input += ' %2s' % atom.symbol
            # write x, y, z
            for idir in range(3):
                mopac_input += '    %16.5f %i' % (xyz[idir], f)
            mopac_input += '\n'

        if atoms.pbc.any():
            for v in atoms.get_cell():
                mopac_input += 'Tv %8.3f %8.3f %8.3f\n' % (v[0], v[1], v[2])
        
        # write input
        myfile = open(fname, 'w')
        myfile.write(mopac_input)
        myfile.close()
        f.close()

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        f = open(self.label + '.nw')
        for line in f:
            if line.startswith('geometry'):
                break
        symbols = []
        positions = []
        for line in f:
            if line.startswith('end'):
                break
            words = line.split()
            symbols.append(words[0])
            positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.atoms = Atoms(symbols, positions,
                           magmoms=self.parameters.pop('initial_magmoms'))
        self.read_results()

    def read_energy(self, fname):
        """
        Reads the ENERGY from the output file (HEAT of FORMATION in kcal / mol)
        Raises RuntimeError if no energy was found
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()

        energy = None
        for line in lines:
            if line.find('HEAT OF FORMATION') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('H.o.F. per unit cell') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('UNABLE TO ACHIEVE SELF-CONSISTENCE') != -1:
                energy = None
        if energy is None:
            raise RuntimeError('MOPAC: could not find total energy')
        
        energy *= (kcal / mol)
        return energy

    def read_forces(self, fname):
        """
        Reads the FORCES from the output file
        search string: (HEAT of FORMATION in kcal / mol / AA)
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()

        nats = len(self.atoms)
        forces = np.zeros((nats, 3), float)
        
        for i, line in enumerate(lines):
            if line.find('GRADIENT\n') != -1:
                for j in range(nats * 3):
                    gline = lines[i + j + 1]
                    forces[j / 3, j % 3] = float(gline[49:62])
                break
        
        forces *= - (kcal / mol)
        return forces
