"""This module defines an ASE interface to ABINIT.

http://www.abinit.org/
"""

import os
from os.path import join, isfile, islink

import numpy as npy

from ase.data import chemical_symbols
from ase.data import atomic_numbers
from ase.units import Bohr, Hartree


class Abinit:
    """Class for doing ABINIT calculations.

    The default parameters are very close to those that the ABINIT
    Fortran code would use.  These are the exceptions::

      calc = Abinit(label='abinit', xc='LDA', pulay=5, mix=0.1)

    Use the set_inp method to set extra INPUT parameters::

      calc.set_inp('nstep', 30)

    """
    def __init__(self, label='abinit', xc='LDA', kpts=None, nbands=1,
                 width=0.04*Hartree, ecut=None, charge=0,
                 pulay=5, mix=0.1
                 ):
        """Construct ABINIT-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'abinit'.
        xc: str
            Exchange-correlation functional.  Must be one of LDA, PBE,
            revPBE, RPBE.
        kpts: list of three int
            Monkhost-Pack sampling.
        nbands: int
            Number of bands.
            Default is 1.
        width: float
            Fermi-distribution width in eV.
            Default is 0.04 Hartree.
        ecut: float
            Planewave cutoff energy in eV.
            No default.
        charge: float
            Total charge of the system.
            Default is 0.
        pulay: int
            Number of old densities to use for Pulay mixing.
        mix: float
            Mixing parameter between zero and one for density mixing.

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=Abinit())
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()

        """

        self.label = label#################### != out
        self.xc = xc
        self.kpts = kpts
        self.nbands = nbands
        self.width = width
        self.ecut = ecut
        self.charge = charge
        self.pulay = pulay
        self.mix = mix

        self.converged = False
        self.inp = {}

    def update(self, atoms):
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            if Z not in self.species:
                self.species.append(Z)

        if 'ABINIT_PP_PATH' in os.environ:
            pppaths = os.environ['ABINIT_PP_PATH'].split(':')
        else:
            pppaths = []

        self.ppp_list = []
        if self.xc != 'LDA':
            xcname = 'GGA'
        else:
            xcname = 'LDA'
        for Z in self.species:
            symbol = chemical_symbols[abs(Z)]
            number = atomic_numbers[symbol]
            name = ('%02d' % number) + '-' + symbol + '.' + xcname + '.fhi'
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

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal + self.efree) / 2

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        self.write_files()

        self.write_inp(atoms)

        abinit = os.environ['ABINIT_SCRIPT']
        locals = {'label': self.label}
        execfile(abinit, {}, locals)
        exitcode = locals['exitcode']
        if exitcode != 0:
            raise RuntimeError(('Abinit exited with exit code: %d.  ' +
                                'Check %s.log for more information.') %
                               (exitcode, self.label))

        self.read()

        self.converged = True

    def write_files(self):
        """Write input parameters to files-file."""
        fh = open(self.label + '.files', 'w')

        import getpass
        #find a suitable default scratchdir (should be writeable!)
        username=getpass.getuser()

        if os.access("/scratch/"+username,os.W_OK):
                scratch = "/scratch/"+username
        elif os.access("/scratch/",os.W_OK):
                scratch = "/scratch/"
        else:
                if os.access(os.curdir,os.W_OK):
                        scratch = os.curdir #if no /scratch use curdir
                else:
                        raise IOError,"No suitable scratch directory and no write access to current dir"

        fh.write('%s\n' % (self.label+'.in')) # input
        fh.write('%s\n' % (self.label+'.txt')) # output
        fh.write('%s\n' % (self.label+'i')) # input
        fh.write('%s\n' % (self.label+'o')) # output
        # scratch files
        fh.write('%s\n' % (os.path.join(scratch, self.label+'.abinit')))
        # Provide the psp files
        for ppp in self.ppp_list:
            fh.write('%s\n' % (ppp)) # psp file path

        fh.close()

    def set_inp(self, key, value):
        """Set INPUT parameter."""
        self.inp[key] = value

    def write_inp(self, atoms):
        """Write input parameters to in-file."""
        fh = open(self.label + '.in', 'w')

        inp = {
            #'SystemLabel': self.label,
            #'LatticeConstant': 1.0,
            'natom': len(atoms),
            'charge': self.charge,
            'nband': self.nbands,
            #'DM.UseSaveDM': self.converged,
            #'SolutionMethod': 'diagon',
            'npulayit': self.pulay, # default 7
            'diemix': self.mix
            }

        if self.ecut is not None:
            inp['ecut'] = str(self.ecut)+' eV' # default Ha

        if self.width is not None:
            inp['tsmear'] = str(self.width)+' eV' # default Ha
            fh.write('occopt 3 # Fermi-Dirac smearing\n')

        inp['ixc'] = { # default 1
            'LDA':     7,
            'PBE':    11,
            'revPBE': 14,
            'RPBE':   15,
            'WC':     23
            }[self.xc]

        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms.any():
            inp['nsppol'] = 2
            fh.write('spinat\n')
            for n, M in enumerate(magmoms):
                if M != 0:
                    fh.write('%.14f %.14f %.14f\n' % (0, 0, M))
        else:
            inp['nsppol'] = 1
            #fh.write('spinat\n')
            #for n, M in enumerate(magmoms):
            #    if M != 0:
            #        fh.write('%.14f\n' % (M))

        inp.update(self.inp)

        for key, value in inp.items():
            if value is None:
                continue

            if isinstance(value, list):
                fh.write('%block %s\n' % key)
                for line in value:
                    fh.write(' '.join(['%s' % x for x in line]) + '\n')
                fh.write('%endblock %s\n' % key)

            unit = keys_with_units.get(inpify(key))
            if unit is None:
                fh.write('%s %s\n' % (key, value))
            else:
                if 'fs**2' in unit:
                    value /= fs**2
                elif 'fs' in unit:
                    value /= fs
                fh.write('%s %f %s\n' % (key, value, unit))

        fh.write('#Definition of the unit cell\n')
        fh.write('acell\n')
        fh.write('%.14f %.14f %.14f Angstrom\n' %  (1.0, 1.0, 1.0))
        fh.write('rprim\n')
        for v in self.cell:
            fh.write('%.14f %.14f %.14f\n' %  tuple(v))

        fh.write('chkprim 0 # Allow non-primitive cells\n')

        fh.write('#Definition of the atom types\n')
        fh.write('ntypat %d\n' % (len(self.species)))
        fh.write('znucl')
        for n, Z in enumerate(self.species):
            fh.write(' %d' % (Z))
        fh.write('\n')
        fh.write('#Enumerate different atomic species\n')
        fh.write('typat')
        fh.write('\n')
        self.types = []
        for Z in self.numbers:
            for n, Zs in enumerate(self.species):
                if Z == Zs:
                    self.types.append(n+1)
        for n, type in enumerate(self.types):
            fh.write(' %d' % (type))
            if n > 1 and ((n % 20) == 1):
                fh.write('\n')
        fh.write('\n')

        fh.write('#Definition of the atoms\n')
        fh.write('xangst\n')
        a = 0
        for pos, Z in zip(self.positions, self.numbers):
            a += 1
            fh.write('%.14f %.14f %.14f\n' %  tuple(pos))

        if self.kpts is not None:
            fh.write('kptopt %d\n' % (1))
            fh.write('ngkpt ')
            fh.write('%d %d %d\n' %  tuple(self.kpts))

        fh.write('#Definition of the SCF procedure\n')
        fh.write('toldfe 1.0d-06\n')
        fh.write('chkexit 1 # abinit.exit file in the running directory terminates after the current SCF\n')

        fh.close()

    def read(self):
        """Read results from ABINIT's text-output file."""
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'ERROR' not in text
        lines = iter(text.split('\n'))
        # some consistency ckecks
        for line in iter(text.split('\n')):
            if line.rfind('natom  ') > -1:
                natom = int(line.split()[-1])
                assert natom == len(self.numbers)
        for line in iter(text.split('\n')):
            if line.rfind('znucl  ') > -1:
                znucl = [float(Z) for Z in line.split()[1:]]
                for n, Z in enumerate(self.species):
                    assert Z == znucl[n]
        for line in iter(text.split('\n')):
            if line.rfind(' typat  ') > -1:
                typat = [float(t) for t in line.split()[1:]]
                for n, t in enumerate(self.types):
                    assert t == typat[n]

        # Stress:
        # Printed in the output in the following format [Hartree/Bohr^3]:
        # sigma(1 1)=  4.02063464E-04  sigma(3 2)=  0.00000000E+00
        # sigma(2 2)=  4.02063464E-04  sigma(3 1)=  0.00000000E+00
        # sigma(3 3)=  4.02063464E-04  sigma(2 1)=  0.00000000E+00
        for line in lines:
            if line.rfind('cartesian components of stress tensor (hartree/bohr^3)') > -1:
                self.stress = npy.empty((3, 3))
                for i in range(3):
                    entries = lines.next().split()
                    self.stress[i,i] = float(entries[2])
                    self.stress[min(3, 4-i)-1, max(1, 2-i)-1] = float(entries[5])
                    self.stress[max(1, 2-i)-1, min(3, 4-i)-1] = float(entries[5])
                self.stress = self.stress*Hartree/Bohr**3
                break
        else:
            raise RuntimeError

        # Energy [Hartree]:
        # Warning: Etotal could mean both electronic energy and free energy!
        for line in iter(text.split('\n')):
            if line.rfind('>>>>> internal e=') > -1:
                self.etotal = float(line.split('=')[-1])*Hartree
                for line1 in iter(text.split('\n')):
                    if line1.rfind('>>>>>>>>> etotal=') > -1:
                        self.efree = float(line1.split('=')[-1])*Hartree
                        break
                else:
                    raise RuntimeError
                break
        else:
            for line2 in iter(text.split('\n')):
                if line2.rfind('>>>>>>>>> etotal=') > -1:
                    self.etotal = float(line2.split('=')[-1])*Hartree
                    self.efree = self.etotal
                    break
            else:
                raise RuntimeError

        # Forces:
        for line in lines:
            if line.rfind('cartesian forces (ev/angstrom) at end:') > -1:
                forces = []
                for i in range(len(self.numbers)):
                    forces.append(npy.array([float(f) for f in lines.next().split()[1:]]))
                self.forces = npy.array(forces)
                break
        else:
            raise RuntimeError

        # Now, because (stupidly) abinit when it finds a name it uses nameA
        # and when nameA exists it uses nameB, etc.
        # we need to rename our file
        # (and loose it in case of e.g. QuasiNewton relaxation)!
        filename_save = filename + '.save'
        if islink(filename) or isfile(filename):
            os.rename(filename, filename_save)

def inpify(key):
    return key.lower().replace('_', '').replace('.', '').replace('-', '')


keys_with_units = {
    }
#keys_with_units = {
#    'paoenergyshift': 'eV',
#    'zmunitslength': 'Bohr',
#    'zmunitsangle': 'rad',
#    'zmforcetollength': 'eV/Ang',
#    'zmforcetolangle': 'eV/rad',
#    'zmmaxdispllength': 'Ang',
#    'zmmaxdisplangle': 'rad',
#    'ecut': 'eV',
#    'dmenergytolerance': 'eV',
#    'electronictemperature': 'eV',
#    'oneta': 'eV',
#    'onetaalpha': 'eV',
#    'onetabeta': 'eV',
#    'onrclwf': 'Ang',
#    'onchemicalpotentialrc': 'Ang',
#    'onchemicalpotentialtemperature': 'eV',
#    'mdmaxcgdispl': 'Ang',
#    'mdmaxforcetol': 'eV/Ang',
#    'mdmaxstresstol': 'eV/Ang**3',
#    'mdlengthtimestep': 'fs',
#    'mdinitialtemperature': 'eV',
#    'mdtargettemperature': 'eV',
#    'mdtargetpressure': 'eV/Ang**3',
#    'mdnosemass': 'eV*fs**2',
#    'mdparrinellorahmanmass': 'eV*fs**2',
#    'mdtaurelax': 'fs',
#    'mdbulkmodulus': 'eV/Ang**3',
#    'mdfcdispl': 'Ang',
#    'warningminimumatomicdistance': 'Ang',
#    'rcspatial': 'Ang',
#    'kgridcutoff': 'Ang',
#    'latticeconstant': 'Ang'}
