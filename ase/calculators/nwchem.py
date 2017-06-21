"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import re
import numpy as np
from ase.io import read
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class NWChem(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole', 'magmom']
    command = 'nwchem PREFIX.nw > PREFIX.out'

    default_parameters = dict(
        xc='LDA',
        smearing=None,
        charge=None,
        task='gradient',
        # Warning: nwchem centers atoms by default
        # see ase-developers/2012-March/001356.html
        geometry='nocenter noautosym',
        convergence={'energy': None,
                     'density': None,
                     'gradient': None,
                     'lshift': None,
                     # set lshift to 0.0 for nolevelshifting
                     'damp': None,
                     },
        basis='3-21G',
        basispar=None,
        ecp=None,
        so=None,
        spinorbit=False,
        odft=False,
        raw='')  # additional outside of dft block control string

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='nwchem', atoms=None, **kwargs):
        """Construct NWchem-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore unit cell and boundary conditions:
        if 'cell' in system_changes:
            system_changes.remove('cell')
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.initial_magmoms = atoms.get_initial_magnetic_moments().tolist()
        p.write(self.label + '.ase')
        del p['initial_magmoms']
        f = open(self.label + '.nw', 'w')
        if p.charge is not None:
            f.write('charge %s\n' % p.charge)
        write_nwchem(f, atoms, p.geometry)

        f.write('start\n')

        if p.basispar is not None:
            basispar = 'basis ' + p.basispar
        else:
            basispar = 'basis'

        def format_basis_set(string, tag=basispar):
            formatted = tag + '\n'
            lines = string.split('\n')
            if len(lines) > 1:
                formatted += string
            else:
                formatted += '  * library ' + string
            return formatted + '\nend\n'

        basis = format_basis_set(p.basis)
        if p.ecp is not None:
            basis += format_basis_set(p.ecp, 'ecp')
        if p.so is not None:
            basis += format_basis_set(p.so, 'so')
        f.write(basis)

        if p.xc == 'RHF':
            task = 'scf'
        elif p.xc == 'MP2':
            task = 'mp2'
        else:
            if p.spinorbit:
                task = 'sodft'
            else:
                task = 'dft'
            xc = {'LDA': 'slater pw91lda',
                  'PBE': 'xpbe96 cpbe96',
                  'revPBE': 'revpbe cpbe96',
                  'RPBE': 'rpbe cpbe96'}.get(p.xc, p.xc)
            f.write('\n' + task + '\n')
            f.write('  xc ' + xc + '\n')
            for key in p.convergence:
                if p.convergence[key] is not None:
                    if key == 'lshift':
                        if p.convergence[key] <= 0.0:
                            f.write('  convergence nolevelshifting\n')
                        else:
                            f.write('  convergence %s %s\n' %
                                    (key, p.convergence[key] / Hartree))
                    else:
                        f.write('  convergence %s %s\n' %
                                (key, p.convergence[key]))
            if p.smearing is not None:
                assert p.smearing[0].lower() == 'gaussian', p.smearing
                f.write('  smear %s\n' % (p.smearing[1] / Hartree))
            if 'mult' not in p:
                # Obtain multiplicity from magnetic momenta:
                tot_magmom = atoms.get_initial_magnetic_moments().sum()
                if tot_magmom < 0:
                    mult = tot_magmom - 1  # fill minority bands
                else:
                    mult = tot_magmom + 1
            else:
                mult = p.mult
            if mult != int(mult):
                raise RuntimeError('Noninteger multiplicity not possible. ' +
                                   'Check initial magnetic moments.')
            f.write('  mult %d\n' % mult)
            if p.odft:
                f.write('  odft\n')  # open shell aka spin polarized dft
            for key in sorted(p.keys()):
                if key in ['charge', 'geometry', 'basis', 'basispar', 'ecp',
                           'so', 'xc', 'spinorbit', 'convergence', 'smearing',
                           'raw', 'mult', 'task', 'odft']:
                    continue
                f.write(u"  {0} {1}\n".format(key, p[key]))
            f.write('end\n')

        if p.raw:
            f.write(p.raw + '\n')
        f.write('\ntask ' + task + ' ' + p.task + '\n')
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

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()
        if self.parameters.task.find('optimize') > -1:
            self.read_coordinates()
            self.read_forces()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.nvector = self.read_number_of_bands()
        self.results['magmom'] = self.read_magnetic_moment()
        dipole = self.read_dipole_moment()
        if dipole is not None:
            self.results['dipole'] = dipole

    def get_ibz_k_points(self):
        return np.array([0., 0., 0.])

    def get_number_of_bands(self):
        return self.nvector

    def read_number_of_bands(self):
        nvector = 0
        for line in open(self.label + '.out'):
            if line.find('Vector ') != -1:  # count all printed vectors
                nvector += 1
        if not nvector:
            nvector = None
        return nvector

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        for line in open(self.label + '.out'):  # find last one
            if line.find('of electrons') != -1:
                nelect = float(line.split(':')[1].strip())
        return nelect

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = 0
        for line in open(self.label + '.out'):
            if line.find('d= ') != -1:  # count all iterations
                niter += 1
        if not niter:
            niter = None
        return niter

    def read_magnetic_moment(self):
        magmom = None
        for line in open(self.label + '.out'):
            if line.find('Spin multiplicity') != -1:  # last one
                magmom = float(line.split(':')[-1].strip())
                if magmom < 0:
                    magmom += 1
                else:
                    magmom -= 1
        return magmom

    def read_dipole_moment(self):
        dipolemoment = []
        for line in open(self.label + '.out'):
            for component in ['1   1 0 0',
                              '1   0 1 0',
                              '1   0 0 1']:
                if line.find(component) != -1:
                    value = float(line.split(component)[1].split()[0])
                    value = value * Bohr
                    dipolemoment.append(value)
        if len(dipolemoment) == 0:
            if len(self.atoms) == 1:
                dipolemoment = [0.0, 0.0, 0.0]
            else:
                return None
        return np.array(dipolemoment)

    def read_energy(self):
        """Read Energy from nwchem output file."""
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))

        # Energy:
        estring = 'Total '
        if self.parameters.xc == 'RHF':
            estring += 'SCF'
        elif self.parameters.xc == 'MP2':
            estring += 'MP2'
        else:
            estring += 'DFT'
        estring += ' energy'
        for line in lines:
            if line.find(estring) >= 0:
                energy = float(line.split()[-1])
        self.results['energy'] = energy * Hartree

        # Eigenstates
        spin = -1
        kpts = []
        for line in lines:
            if line.find('Molecular Orbital Analysis') >= 0:
                last_eps = -99999.0
                spin += 1
                kpts.append(KPoint(spin))
            if spin >= 0:
                if line.find('Vector') >= 0:
                    line = line.lower().replace('d', 'e')
                    line = line.replace('=', ' ')
                    word = line.split()
                    this_occ = float(word[3])
                    this_eps = float(word[5])
                    kpts[spin].f_n.append(this_occ)
                    kpts[spin].eps_n.append(this_eps)
                    if this_occ < 0.1 and this_eps < last_eps:
                        warn('HOMO above LUMO - if this is not an exicted ' +
                             'state - this might be introduced by levelshift.',
                             RuntimeWarning)
                    last_eps = this_eps
        self.kpts = kpts

    def read_forces(self):
        """Read Forces from nwchem output file."""
        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        for i, line in enumerate(lines):
            if line.find('ENERGY GRADIENTS') >= 0:
                gradients = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(5, 8)])

        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def read_coordinates(self):
        """Read updated coordinates from nwchem output file."""
        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        for i, line in enumerate(lines):
            if line.find('ENERGY GRADIENTS') >= 0:
                positions = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    positions.append([float(word[k]) for k in range(2, 5)])

        self.atoms.set_positions(np.array(positions) * Bohr)

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.array(self.kpts[spin].eps_n) * Hartree

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return self.kpts[spin].f_n

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return len(self.kpts)

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return len(self.kpts) == 2


class KITnwchem(NWChem):
    """ an enhanced nwchem calculator class for molecular orbital data """

    energy = None
    orbitals = None

    def __init__(self, **parameters):

        self.set_parameters(parameters)
        self.verify_parameters()

        NWChem.__init__(self, label=self.title, geometry=self.geometry,
                        basis=self.basis_set, charge=self.charge,
                        mult=self.multiplicity, xc=self.func,
                        convergence=self.convergence, task=self.task)
        self.reset()

    def set_parameters(self, parameters):

        self.geometry = 'nocenter noautosym noautoz units angstroms'

        # e.g. "mpirun -np 4 nwchem PREFIX.nw > PREFIX.out"
        if 'command' in parameters:
            self.command = parameters['command']

        if 'title' in parameters:
            self.title = parameters['title'].replace(" ", "_")
        else:
            self.title = 'notitle'
        if 'basis set name' in parameters:
            self.basis_set = parameters['basis set name']
        if 'charge' in parameters:
            self.charge = parameters['charge']
        else:
            self.charge = 0
        if 'multiplicity' in parameters:
            self.multiplicity = parameters['multiplicity']
        if 'dft' in parameters:
            self.dft = parameters['dft']
        if 'density functional' in parameters:
            self.func = parameters['density functional']
        else:
            self.func = 'lda'
            print "Density functional set to default: ", self.func
        if 'ri' in parameters:
            self.ri = parameters['ri']
        if 'internals' in parameters:
            self.internals = parameters['internals']
        else:
            self.internals = False
        if 'scfiterlimit' in parameters:
            self.scfiterlimit = parameters['scfiterlimit']
        else:
            self.scfiterlimit = 60
        if 'force convergence' in parameters:
            self.force_conv = parameters['force convergence']
        else:
            self.force_conv = None
        if 'energy convergence' in parameters:
            self.energy_conv = parameters['energy convergence']
        else:
            self.energy_conv = None
        self.convergence = {
            'energy': self.energy_conv,
            'density': None,
            'gradient': self.force_conv,
            'lshift': None,
            'damp': None
        }
        if 'task' in parameters:
            self.task = parameters['task']
        else:
            # note that 'gradient' is the default of the base class
            self.task = 'energy'
        if 'geometry iterations' in parameters:
            self.maxiter = parameters['geometry iterations']
        else:
            self.maxiter = None

    def verify_parameters(self):
        """ only a subset of all supported functionals, no ECPs, no check for
        availability for given elements """
        available_functionals = [
            'acm', 'b3lyp', 'beckehandh', 'pbe0', 'becke97', 'becke97-1', 'becke97-2',
            'becke97-3', 'becke97-d', 'becke98', 'hcth', 'hcth120', 'hcth147', 'hcth407',
            'becke97gga1', 'hcth407p', 'mpw91', 'mpw1k', 'xft97', 'cft97', 'ft97', 'op',
            'bop', 'pbeop', 'xpkzb99', 'cpkzb99', 'xtpss03', 'ctpss03', 'xctpssh', 'b1b95',
            'bb1k', 'mpw1b95', 'mpwb1k', 'pw6b95', 'pwb6k', 'm05', 'm05-2x', 'vs98', 'm06',
            'm06-hf', 'm06-L', 'm06-2x'
        ]
        available_basis_sets = [
            '3-21++G', '3-21++G*', '3-21G', '3-21G*', '3-21GSP', '4-22GSP', '4-31G', '6-31++G', '6-31++G*', '6-31++G**', '6-31+G*', '6-311++G(2d,2p)', '6-311++G(3df,3pd)', '6-311++G**', '6-311+G*', '6-311G', '6-311G(2df,2pd)', '6-311G*', '6-311G**', '6-31G', '6-31G-Blaudeau', '6-31G(3df,3pd)', '6-31G*', '6-31G*-Blaudeau', '6-31G**', 'Ahlrichs pVDZ', 'Ahlrichs TZV', 'Ahlrichs VDZ', 'Ahlrichs VTZ', 'ANO-RCC', 'apr-cc-pV(Q+d)Z', 'aug-cc-pCV5Z', 'aug-cc-pCVDZ', 'aug-cc-pCVQZ', 'aug-cc-pCV(T+d)Z', 'aug-cc-pCVTZ', 'aug-cc-pV(5+d)Z', 'aug-cc-pV5Z', 'aug-cc-pV(6+d)Z', 'aug-cc-pV6Z', 'aug-cc-pV(D+d)Z', 'aug-cc-pVDZ', 'aug-cc-pV(Q+d)Z', 'aug-cc-pVQZ', 'aug-cc-pV(T+d)Z', 'aug-cc-pVTZ', 'aug-cc-pVTZ-J', 'aug-cc-pwCV5Z', 'aug-cc-pwCV5Z-NR', 'aug-cc-pwCVDZ', 'aug-cc-pwCVQZ', 'aug-cc-pwCVQZ-NR', 'aug-cc-pwCVTZ', 'aug-cc-pwCVTZ-NR', 'aug-mcc-pV5Z', 'aug-mcc-pV6Z', 'aug-mcc-pV7Z', 'aug-mcc-pV8Z', 'aug-mcc-pVQZ', 'aug-mcc-pVTZ', 'aug-pc-0', 'aug-pc-1', 'aug-pc-2', 'aug-pc-3', 'aug-pc-4', 'aug-pcJ-0', 'aug-pcJ-0_2006', 'aug-pcJ-1', 'aug-pcJ-1_2006', 'aug-pcJ-2', 'aug-pcJ-2_2006', 'aug-pcJ-3', 'aug-pcJ-3_2006', 'aug-pcJ-4', 'aug-pcJ-4_2006', 'aug-pcS-0', 'aug-pcS-1', 'aug-pcS-2', 'aug-pcS-3', 'aug-pcS-4', 'aug-pV7Z', 'B2 basis set for Zn', 'Bauschlicher ANO', 'Binning/Curtiss SV', 'Binning/Curtiss SVP', 'Binning/Curtiss VTZ', 'Binning/Curtiss VTZP', 'cc-pCV5Z', 'cc-pCV6Z', 'cc-pCV6Z(old)', 'cc-pCVDZ', 'cc-pCVDZ(old)', 'cc-pCVQZ', 'cc-pCVQZ(old)', 'cc-pCVTZ', 'cc-pCVTZ(old)', 'cc-pV(5+d)Z', 'cc-pV5Z', 'cc-pV(6+d)Z', 'cc-pV6Z', 'cc-pV8Z', 'cc-pV9Z', 'cc-pV(D+d)Z', 'cc-pVDZ', 'cc-pVDZ(seg-opt)', 'cc-pV(Q+d)Z', 'cc-pVQZ', 'cc-pVQZ(seg-opt)', 'cc-pV(T+d)Z', 'cc-pVTZ', 'cc-pVTZ(seg-opt)', 'cc-pwCV5Z', 'cc-pwCV5Z Core Set', 'cc-pwCV5Z-NR', 'cc-pwCVDZ', 'cc-pwCVQZ', 'cc-pwCVQZ-NR', 'cc-pwCVTZ', 'cc-pwCVTZ-NR', 'ccemd-2', 'ccemd-3', 'ccJ-pV5Z', 'ccJ-pVDZ', 'ccJ-pVQZ', 'ccJ-pVTZ', 'coemd-2', 'coemd-3', 'coemd-4', 'coemd-ref', 'CVTZ', 'd-aug-cc-pV5Z', 'd-aug-cc-pV6Z', 'd-aug-cc-pVDZ', 'd-aug-cc-pVQZ', 'd-aug-cc-pVTZ', 'Def2-QZVP', 'Def2-QZVPD', 'Def2-SVP', 'def2-SV(P)', 'Def2-SVPD', 'Def2-TZVP', 'Def2-TZVPD', 'dhf-QZVP', 'dhf-SV(P)', 'dhf-TZVP', 'Dunning-Hay Double Rydberg', 'Dunning-Hay Rydberg', 'DZ + Double Rydberg (Dunning-Hay)', 'DZ + Rydberg (Dunning)', 'DZ (Dunning)', 'DZP + Rydberg (Dunning)', 'DZP (Dunning)', 'DZQ', 'DZVP2 (DFT Orbital)', 'DZVP (DFT Orbital)', 'Feller Misc. CVDZ', 'Feller Misc. CVQZ', 'Feller Misc. CVTZ', 'GAMESS PVTZ', 'GAMESS VTZ', 'IGLO-II', 'IGLO-III', 'jul-cc-pV(D+d)Z', 'jul-cc-pV(Q+d)Z', 'jul-cc-pV(T+d)Z', 'jun-cc-pV(D+d)Z', 'jun-cc-pV(Q+d)Z', 'jun-cc-pV(T+d)Z', 'LANL08', 'LANL08+', 'LANL08d', 'LANL08(f)', 'Lanl2-\[10s8p7d3f2g\]', 'Lanl2-\[5s4p4d2f\]', 'Lanl2-\[6s4p4d2f\]', 'Lanl2DZ+1d1f', 'Lanl2DZ+2s2p2d2f', 'LANL2TZ', 'LANL2TZ+', 'LANL2TZ(f)', 'm6-31G', 'm6-31G*', 'maug-cc-pV(D+d)Z', 'maug-cc-pVDZ', 'maug-cc-pV(Q+d)Z', 'maug-cc-pVQZ', 'maug-cc-pV(T+d)Z', 'maug-cc-pVTZ', 'may-cc-pV(Q+d)Z', 'may-cc-pV(T+d)Z', 'McLean/Chandler VTZ', 'MG3S', 'MIDI!', 'MIDI (Huzinaga)', 'MINI (Huzinaga)', 'MINI (Scaled)', 'modified LANL2DZ', 'NASA Ames ANO', 'NASA Ames ANO2', 'NASA Ames cc-pCV5Z', 'NASA Ames cc-pCVQZ', 'NASA Ames cc-pCVTZ', 'NASA Ames cc-pV5Z', 'NASA Ames cc-pVQZ', 'NASA Ames cc-pVTZ', 'Partridge Uncontracted 1', 'Partridge Uncontracted 2', 'Partridge Uncontracted 3', 'Partridge Uncontracted 4', 'pc-0', 'pc-1', 'pc-2', 'pc-3', 'pc-4', 'pcemd-2', 'pcemd-3', 'pcemd-4', 'pcJ-0', 'pcJ-0_2006', 'pcJ-1', 'pcJ-1_2006', 'pcJ-2', 'pcJ-2_2006', 'pcJ-3', 'pcJ-3_2006', 'pcJ-4', 'pcJ-4_2006', 'pcS-0', 'pcS-1', 'pcS-2', 'pcS-3', 'pcS-4', 'pSBKJC', 'Pt - mDZP', 'pV6Z', 'pV7Z', 'Roos_ANO_DZ', 'Roos_ANO_TZ', 'Roos Augmented Double Zeta ANO', 'Roos Augmented Triple Zeta ANO', 's3-21G', 's3-21G*', 's6-31G', 's6-31G*', 'Sadlej pVTZ', 'SDB-aug-cc-pVQZ', 'SDB-aug-cc-pVTZ', 'SDB-cc-pVQZ', 'SDB-cc-pVTZ', 'STO-2G', 'STO-3G', 'STO-3G*', 'STO-6G', 'SV + Double Rydberg (Dunning-Hay)', 'SV + Rydberg (Dunning-Hay)', 'SV (Dunning-Hay)', 'SVP + Rydberg (Dunning-Hay)', 'SVP (Dunning-Hay)', 'TZ (Dunning)', 'TZVP (DFT Orbital)', 'UGBS', 'un-ccemd-ref', 'un-pcemd-ref', 'Wachters+f', 'WTBS', 'Z3Pol'
        ]

        assert self.basis_set.lower() in [x.lower() for x in available_basis_sets], (
            'basis set ', self.basis_set, ' not available / not supported'
        )
        assert self.func.lower() in [x.lower() for x in available_functionals], (
            'density functional ', self.func, ' not available / not supported'
        )

    def reset(self):
        NWChem.reset(self)
        self.atoms = None
        self.results = {}
        self.initialized = False
        self.orbitals = []

    def write_input(self, atoms, properties=None, system_changes=None):
        NWChem.write_input(self, atoms, properties=properties,
                           system_changes=system_changes)
        filename = self.title + '.nw'
        buf = []
        drvset = False
        for line in open(filename, 'r'):
            buf.append(line)
            if re.search('^dft\s*$', line):
                buf.append('  print \"final vectors analysis\" \"final vectors\"\n')
                buf.append('  direct\n')
                buf.append('  noio\n')
                buf.append('  iterations ' + str(self.scfiterlimit) + '\n')
            if re.search('^end\s*$', line) and self.task == 'optimize' and not drvset:
                buf.append('\ndriver\n')
                if self.force_conv > 0.000450:
                    buf.append('  loose\n')
                if self.maxiter:
                    buf.append('  maxiter ' + str(self.maxiter) + '\n')
                buf.append('end\n\n')
                drvset = True

        outfile = open(filename, 'w')
        outfile.write(''.join(buf))
        outfile.close()

    def get_potential_energy(self, atoms):
        """ return the total energy """
        if self.task == 'optimize':
            energy = NWChem.get_potential_energy(self, atoms)
            self.read_geometry()
            atoms.set_positions(self.atoms.get_positions())
            return energy
        else:
            return NWChem.get_potential_energy(self, atoms)


    def get_mos(self, atoms):
        """ return eigenvectors and eigenvalues as matrices (numpy arrays) """
        if not self.orbitals:
            self.energy = atoms.get_potential_energy()
        orbital_coefficients = []
        eigen_energies = []
        for orbital in self.orbitals:
            orbital_coefficients.append(orbital['coefficients'])
            eigen_energies.append(orbital['energy'])
        return [
            np.array(orbital_coefficients),
            np.diag(np.array(eigen_energies)*Hartree)
        ]

    def get_fock_overlap(self):
        """ calculate overlap matrix and fock matrix: will use FC = SCe """
        [C, e] = self.get_mos(self.atoms)
        # C^-1
        C_inv = np.linalg.inv(C)
        # overlap matrix of the dimer
        overlap = np.matmul(C_inv, C_inv.T)
        # Fock/Kohn-Sham matrix
        fock = np.matmul(C_inv, np.matmul(e, C_inv.T))
        return [fock, overlap]

    def read_results(self):
        NWChem.read_results(self)
        self.read_mos()

    def read_energy(self):
        """ reads the last energy from geometry optimization output """
        NWChem.read_energy(self)
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))

        estring = 'Total '
        if self.parameters.xc == 'RHF':
            estring += 'SCF'
        elif self.parameters.xc == 'MP2':
            estring += 'MP2'
        else:
            estring += 'DFT'
        estring += ' energy'
        for line in lines:
            if line.find(estring) >= 0:
                energy = float(line.split()[-1])
                continue
        self.results['energy'] = energy*Hartree

    def read_geometry(self):
        """ reads the last geometry from geometry optimization output """
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))
        natoms = len(self.atoms)
        read_flag = False
        atoms_read = 0
        xyz_string = str(natoms)+'\n\n'
        for line in lines:
            if 'Output coordinates in angstroms' in line:
                read_flag = True
                xyz_string = str(natoms)+'\n\n'
                continue
            if read_flag:
                match = re.search('(\d+)\s+(\w+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)', line)
                if match:
                    string = match.group(2) + ' ' + match.group(4) + ' ' + match.group(5) + ' ' + match.group(6) + '\n'
                    xyz_string += string
                    atoms_read += 1
                    if int(match.group(1)) > natoms:
                        print 'error: number of xyz lines more than expected number of atoms'
            if atoms_read == natoms:
                read_flag = False
                atoms_read = 0
        new_atoms = read(StringIO(xyz_string), format='xyz')
        self.atoms.set_positions(new_atoms.get_positions())

    def read_mos(self):
        self.orbitals = []
        orbital_coefficients = []
        orbital_coefficients_line = []
        filename = self.title + '.out'
        orb_number = None
        mo_read = False
        mo_init = False
        ev_number = 0
        for line in open(filename, 'r'):
            if 'Final MO vectors' in line:
                mo_read = True
                mo_init = False
                orb_number = None
                orbital_coefficients = []
                orbital_coefficients_line = []
                continue
            if mo_read:
                regex = r'^\s*(\d+)(\s+[+-]?\d+\.\d+)(\s+[+-]?\d+\.\d+)?(\s+[+-]?\d+\.\d+)?(\s+[+-]?\d+\.\d+)?(\s+[+-]?\d+\.\d+)?(\s+[+-]?\d+\.\d+)?'
                match = re.search(regex, line)
                if match:
                    orb_number = int(match.group(1))
                    for i in range(2, 8):
                        if match.group(i):
                            orbital_coefficients_line += [float(match.group(i))]
                    if mo_init:
                        orbital_coefficients[orb_number-1] += orbital_coefficients_line
                    else:
                        orbital_coefficients.append(orbital_coefficients_line)
                    orbital_coefficients_line = []
                elif orb_number:
                    if orb_number == len(orbital_coefficients):
                        mo_init = True
                    if orb_number == len(orbital_coefficients[0]):
                        mo_read = False
                        continue
            if 'Vector' in line:
                regex = 'Vector\s+(\d+)\s+Occ=\s*([\w\.+-]+)\s+E=\s*([+-\.\w]+)'
                match = re.search(regex, line)
                if match:
                    ev_number += 1
                    # for the case of geometry optimization output
                    if int(match.group(1)) < ev_number:
                        self.orbitals = []
                        ev_number = int(match.group(1))
                    self.orbitals.append({
                        'number': int(match.group(1)),
                        'occupation': float(re.sub('[dD]', 'E', match.group(2))),
                        'energy': float(re.sub('[dD]', 'E', match.group(3)))
                        })
        # this is to transpose the nested list
        mos = [[x[i] for x in orbital_coefficients] for i in range(len(orbital_coefficients[0]))]

        for orbital in self.orbitals:
            orbital['coefficients'] = mos.pop(0)
