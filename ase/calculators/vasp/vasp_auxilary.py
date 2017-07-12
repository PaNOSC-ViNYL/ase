"""Auxilary functions for the ASE VASP calculator.
Defines the following classes:
VaspChargeDensity, VaspDos and xdat2traj

www.vasp.at
"""

import os
import numpy as np

import ase.io
from ase.calculators.singlepoint import SinglePointCalculator

from .vasp import Vasp


class VaspChargeDensity(object):
    """Class for representing VASP charge density"""

    def __init__(self, filename='CHG'):
        # Instance variables
        self.atoms = []   # List of Atoms objects
        self.chg = []     # Charge density
        self.chgdiff = []  # Charge density difference, if spin polarized
        self.aug = ''     # Augmentation charges, not parsed just a big string
        self.augdiff = ''  # Augmentation charge differece, is spin polarized

        # Note that the augmentation charge is not a list, since they
        # are needed only for CHGCAR files which store only a single
        # image.
        if filename is not None:
            self.read(filename)

    def is_spin_polarized(self):
        if len(self.chgdiff) > 0:
            return True
        return False

    def _read_chg(self, fobj, chg, volume):
        """Read charge from file object

        Utility method for reading the actual charge density (or
        charge density difference) from a file object. On input, the
        file object must be at the beginning of the charge block, on
        output the file position will be left at the end of the
        block. The chg array must be of the correct dimensions.

        """
        # VASP writes charge density as
        # WRITE(IU,FORM) (((C(NX,NY,NZ),NX=1,NGXC),NY=1,NGYZ),NZ=1,NGZC)
        # Fortran nested implied do loops; innermost index fastest
        # First, just read it in
        for zz in range(chg.shape[2]):
            for yy in range(chg.shape[1]):
                chg[:, yy, zz] = np.fromfile(fobj, count=chg.shape[0],
                                             sep=' ')
        chg /= volume

    def read(self, filename='CHG'):
        """Read CHG or CHGCAR file.

        If CHG contains charge density from multiple steps all the
        steps are read and stored in the object. By default VASP
        writes out the charge density every 10 steps.

        chgdiff is the difference between the spin up charge density
        and the spin down charge density and is thus only read for a
        spin-polarized calculation.

        aug is the PAW augmentation charges found in CHGCAR. These are
        not parsed, they are just stored as a string so that they can
        be written again to a CHGCAR format file.

        """
        import ase.io.vasp as aiv
        f = open(filename)
        self.atoms = []
        self.chg = []
        self.chgdiff = []
        self.aug = ''
        self.augdiff = ''
        while True:
            try:
                atoms = aiv.read_vasp(f)
            except (IOError, ValueError, IndexError):
                # Probably an empty line, or we tried to read the
                # augmentation occupancies in CHGCAR
                break
            f.readline()
            ngr = f.readline().split()
            ng = (int(ngr[0]), int(ngr[1]), int(ngr[2]))
            chg = np.empty(ng)
            self._read_chg(f, chg, atoms.get_volume())
            self.chg.append(chg)
            self.atoms.append(atoms)
            # Check if the file has a spin-polarized charge density part, and
            # if so, read it in.
            fl = f.tell()
            # First check if the file has an augmentation charge part (CHGCAR
            # file.)
            line1 = f.readline()
            if line1 == '':
                break
            elif line1.find('augmentation') != -1:
                augs = [line1]
                while True:
                    line2 = f.readline()
                    if line2.split() == ngr:
                        self.aug = ''.join(augs)
                        augs = []
                        chgdiff = np.empty(ng)
                        self._read_chg(f, chgdiff, atoms.get_volume())
                        self.chgdiff.append(chgdiff)
                    elif line2 == '':
                        break
                    else:
                        augs.append(line2)
                if len(self.aug) == 0:
                    self.aug = ''.join(augs)
                    augs = []
                else:
                    self.augdiff = ''.join(augs)
                    augs = []
            elif line1.split() == ngr:
                chgdiff = np.empty(ng)
                self._read_chg(f, chgdiff, atoms.get_volume())
                self.chgdiff.append(chgdiff)
            else:
                f.seek(fl)
        f.close()

    def _write_chg(self, fobj, chg, volume, format='chg'):
        """Write charge density

        Utility function similar to _read_chg but for writing.

        """
        # Make a 1D copy of chg, must take transpose to get ordering right
        chgtmp = chg.T.ravel()
        # Multiply by volume
        chgtmp = chgtmp * volume
        # Must be a tuple to pass to string conversion
        chgtmp = tuple(chgtmp)
        # CHG format - 10 columns
        if format.lower() == 'chg':
            # Write all but the last row
            for ii in range((len(chgtmp) - 1) // 10):
                fobj.write(' %#11.5G %#11.5G %#11.5G %#11.5G %#11.5G\
 %#11.5G %#11.5G %#11.5G %#11.5G %#11.5G\n' % chgtmp[ii * 10:(ii + 1) * 10]
                           )
            # If the last row contains 10 values then write them without a
            # newline
            if len(chgtmp) % 10 == 0:
                fobj.write(' %#11.5G %#11.5G %#11.5G %#11.5G %#11.5G'
                           ' %#11.5G %#11.5G %#11.5G %#11.5G %#11.5G' %
                           chgtmp[len(chgtmp) - 10:len(chgtmp)])
            # Otherwise write fewer columns without a newline
            else:
                for ii in range(len(chgtmp) % 10):
                    fobj.write((' %#11.5G')
                               % chgtmp[len(chgtmp) - len(chgtmp) % 10 + ii])
        # Other formats - 5 columns
        else:
            # Write all but the last row
            for ii in range((len(chgtmp) - 1) // 5):
                fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E\n'
                           % chgtmp[ii * 5:(ii + 1) * 5])
            # If the last row contains 5 values then write them without a
            # newline
            if len(chgtmp) % 5 == 0:
                fobj.write(' %17.10E %17.10E %17.10E %17.10E %17.10E'
                           % chgtmp[len(chgtmp) - 5:len(chgtmp)])
            # Otherwise write fewer columns without a newline
            else:
                for ii in range(len(chgtmp) % 5):
                    fobj.write((' %17.10E')
                               % chgtmp[len(chgtmp) - len(chgtmp) % 5 + ii])
        # Write a newline whatever format it is
        fobj.write('\n')
        # Clean up
        del chgtmp

    def write(self, filename='CHG', format=None):
        """Write VASP charge density in CHG format.

        filename: str
            Name of file to write to.
        format: str
            String specifying whether to write in CHGCAR or CHG
            format.

        """
        import ase.io.vasp as aiv
        if format is None:
            if filename.lower().find('chgcar') != -1:
                format = 'chgcar'
            elif filename.lower().find('chg') != -1:
                format = 'chg'
            elif len(self.chg) == 1:
                format = 'chgcar'
            else:
                format = 'chg'
        f = open(filename, 'w')
        for ii, chg in enumerate(self.chg):
            if format == 'chgcar' and ii != len(self.chg) - 1:
                continue  # Write only the last image for CHGCAR
            aiv.write_vasp(f, self.atoms[ii], direct=True, long_format=False)
            f.write('\n')
            for dim in chg.shape:
                f.write(' %4i' % dim)
            f.write('\n')
            vol = self.atoms[ii].get_volume()
            self._write_chg(f, chg, vol, format)
            if format == 'chgcar':
                f.write(self.aug)
            if self.is_spin_polarized():
                if format == 'chg':
                    f.write('\n')
                for dim in chg.shape:
                    f.write(' %4i' % dim)
                self._write_chg(f, self.chgdiff[ii], vol, format)
                if format == 'chgcar':
                    f.write('\n')
                    f.write(self.augdiff)
            if format == 'chg' and len(self.chg) > 1:
                f.write('\n')
        f.close()


class VaspDos(object):
    """Class for representing density-of-states produced by VASP

    The energies are in property self.energy

    Site-projected DOS is accesible via the self.site_dos method.

    Total and integrated DOS is accessible as numpy.ndarray's in the
    properties self.dos and self.integrated_dos. If the calculation is
    spin polarized, the arrays will be of shape (2, NDOS), else (1,
    NDOS).

    The self.efermi property contains the currently set Fermi
    level. Changing this value shifts the energies.

    """

    def __init__(self, doscar='DOSCAR', efermi=0.0):
        """Initialize"""
        self._efermi = 0.0
        self.read_doscar(doscar)
        self.efermi = efermi

        # we have determine the resort to correctly map ase atom index to the
        # POSCAR.
        self.sort = []
        self.resort = []
        if os.path.isfile('ase-sort.dat'):
            file = open('ase-sort.dat', 'r')
            lines = file.readlines()
            file.close()
            for line in lines:
                data = line.split()
                self.sort.append(int(data[0]))
                self.resort.append(int(data[1]))

    def _set_efermi(self, efermi):
        """Set the Fermi level."""
        ef = efermi - self._efermi
        self._efermi = efermi
        self._total_dos[0, :] = self._total_dos[0, :] - ef
        try:
            self._site_dos[:, 0, :] = self._site_dos[:, 0, :] - ef
        except IndexError:
            pass

    def _get_efermi(self):
        return self._efermi

    efermi = property(_get_efermi, _set_efermi, None, "Fermi energy.")

    def _get_energy(self):
        """Return the array with the energies."""
        return self._total_dos[0, :]
    energy = property(_get_energy, None, None, "Array of energies")

    def site_dos(self, atom, orbital):
        """Return an NDOSx1 array with dos for the chosen atom and orbital.

        atom: int
            Atom index
        orbital: int or str
            Which orbital to plot

        If the orbital is given as an integer:
        If spin-unpolarized calculation, no phase factors:
        s = 0, p = 1, d = 2
        Spin-polarized, no phase factors:
        s-up = 0, s-down = 1, p-up = 2, p-down = 3, d-up = 4, d-down = 5
        If phase factors have been calculated, orbitals are
        s, py, pz, px, dxy, dyz, dz2, dxz, dx2
        double in the above fashion if spin polarized.

        """
        # Correct atom index for resorting if we need to. This happens when the
        # ase-sort.dat file exists, and self.resort is not empty.
        if self.resort:
            atom = self.resort[atom]

        # Integer indexing for orbitals starts from 1 in the _site_dos array
        # since the 0th column contains the energies
        if isinstance(orbital, int):
            return self._site_dos[atom, orbital + 1, :]
        n = self._site_dos.shape[1]
        if n == 4:
            norb = {'s': 1, 'p': 2, 'd': 3}
        elif n == 7:
            norb = {'s+': 1, 's-up': 1, 's-': 2, 's-down': 2,
                    'p+': 3, 'p-up': 3, 'p-': 4, 'p-down': 4,
                    'd+': 5, 'd-up': 5, 'd-': 6, 'd-down': 6}
        elif n == 10:
            norb = {'s': 1, 'py': 2, 'pz': 3, 'px': 4,
                    'dxy': 5, 'dyz': 6, 'dz2': 7, 'dxz': 8,
                    'dx2': 9}
        elif n == 19:
            norb = {'s+': 1, 's-up': 1, 's-': 2, 's-down': 2,
                    'py+': 3, 'py-up': 3, 'py-': 4, 'py-down': 4,
                    'pz+': 5, 'pz-up': 5, 'pz-': 6, 'pz-down': 6,
                    'px+': 7, 'px-up': 7, 'px-': 8, 'px-down': 8,
                    'dxy+': 9, 'dxy-up': 9, 'dxy-': 10, 'dxy-down': 10,
                    'dyz+': 11, 'dyz-up': 11, 'dyz-': 12, 'dyz-down': 12,
                    'dz2+': 13, 'dz2-up': 13, 'dz2-': 14, 'dz2-down': 14,
                    'dxz+': 15, 'dxz-up': 15, 'dxz-': 16, 'dxz-down': 16,
                    'dx2+': 17, 'dx2-up': 17, 'dx2-': 18, 'dx2-down': 18}
        return self._site_dos[atom, norb[orbital.lower()], :]

    def _get_dos(self):
        if self._total_dos.shape[0] == 3:
            return self._total_dos[1, :]
        elif self._total_dos.shape[0] == 5:
            return self._total_dos[1:3, :]
    dos = property(_get_dos, None, None, 'Average DOS in cell')

    def _get_integrated_dos(self):
        if self._total_dos.shape[0] == 3:
            return self._total_dos[2, :]
        elif self._total_dos.shape[0] == 5:
            return self._total_dos[3:5, :]
    integrated_dos = property(_get_integrated_dos, None, None,
                              'Integrated average DOS in cell')

    def read_doscar(self, fname="DOSCAR"):
        """Read a VASP DOSCAR file"""
        f = open(fname)
        natoms = int(f.readline().split()[0])
        [f.readline() for nn in range(4)]  # Skip next 4 lines.
        # First we have a block with total and total integrated DOS
        ndos = int(f.readline().split()[2])
        dos = []
        for nd in range(ndos):
            dos.append(np.array([float(x) for x in f.readline().split()]))
        self._total_dos = np.array(dos).T
        # Next we have one block per atom, if INCAR contains the stuff
        # necessary for generating site-projected DOS
        dos = []
        for na in range(natoms):
            line = f.readline()
            if line == '':
                # No site-projected DOS
                break
            ndos = int(line.split()[2])
            line = f.readline().split()
            cdos = np.empty((ndos, len(line)))
            cdos[0] = np.array(line)
            for nd in range(1, ndos):
                line = f.readline().split()
                cdos[nd] = np.array([float(x) for x in line])
            dos.append(cdos.T)
        self._site_dos = np.array(dos)


class xdat2traj:
    def __init__(self, trajectory=None, atoms=None, poscar=None,
                 xdatcar=None, sort=None, calc=None):
        """
        trajectory is the name of the file to write the trajectory to
        poscar is the name of the poscar file to read. Default: POSCAR
        """
        if not poscar:
            self.poscar = 'POSCAR'
        else:
            self.poscar = poscar

        if not atoms:
            # This reads the atoms sorted the way VASP wants
            self.atoms = ase.io.read(self.poscar, format='vasp')
            resort_reqd = True
        else:
            # Assume if we pass atoms that it is sorted the way we want
            self.atoms = atoms
            resort_reqd = False

        if not calc:
            self.calc = Vasp()
        else:
            self.calc = calc
        if not sort:
            if not hasattr(self.calc, 'sort'):
                self.calc.sort = list(range(len(self.atoms)))
        else:
            self.calc.sort = sort
        self.calc.resort = list(range(len(self.calc.sort)))
        for n in range(len(self.calc.resort)):
            self.calc.resort[self.calc.sort[n]] = n

        if not xdatcar:
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = xdatcar

        if not trajectory:
            self.trajectory = 'out.traj'
        else:
            self.trajectory = trajectory

        self.out = ase.io.trajectory.Trajectory(self.trajectory,
                                                mode='w')

        if resort_reqd:
            self.atoms = self.atoms[self.calc.resort]
        self.energies = self.calc.read_energy(all=True)[1]
        # Forces are read with the atoms sorted using resort
        self.forces = self.calc.read_forces(self.atoms, all=True)

    def convert(self):
        lines = open(self.xdatcar).readlines()
        if len(lines[7].split()) == 0:
            del(lines[0:8])
        elif len(lines[5].split()) == 0:
            del(lines[0:6])
        elif len(lines[4].split()) == 0:
            del(lines[0:5])
        elif lines[7].split()[0] == 'Direct':
            del(lines[0:8])
        step = 0
        iatom = 0
        scaled_pos = []
        for line in lines:
            if iatom == len(self.atoms):
                if step == 0:
                    self.out.write_header(self.atoms[self.calc.resort])
                scaled_pos = np.array(scaled_pos)
                # Now resort the positions to match self.atoms
                self.atoms.set_scaled_positions(scaled_pos[self.calc.resort])

                calc = SinglePointCalculator(self.atoms,
                                             energy=self.energies[step],
                                             forces=self.forces[step])
                self.atoms.set_calculator(calc)
                self.out.write(self.atoms)
                scaled_pos = []
                iatom = 0
                step += 1
            else:
                if not line.split()[0] == 'Direct':
                    iatom += 1
                    scaled_pos.append([float(line.split()[n])
                                       for n in range(3)])

        # Write also the last image
        # I'm sure there is also more clever fix...
        if step == 0:
            self.out.write_header(self.atoms[self.calc.resort])
        scaled_pos = np.array(scaled_pos)[self.calc.resort]
        self.atoms.set_scaled_positions(scaled_pos)
        calc = SinglePointCalculator(self.atoms,
                                     energy=self.energies[step],
                                     forces=self.forces[step])
        self.atoms.set_calculator(calc)
        self.out.write(self.atoms)

        self.out.close()
