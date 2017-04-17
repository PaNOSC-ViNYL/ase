from __future__ import print_function
"""
This module defines an ASE interface to Turbomole: http://www.turbomole.com/
"""

import os, sys, re
from math import log10
import numpy as np
import fortranformat as ff
from subprocess import Popen, PIPE, STDOUT
from ase import Atoms
from ase.units import Hartree, Bohr
from ase.io import read, write
from ase.calculators.general import Calculator

class Turbomole(Calculator):
    name = 'Turbomole'

    implemented_properties = ['energy', 'forces', 'dipole', 'free_energy']

    available_basis_sets = [ 
        'sto-3g hondo', 'def-SV(P)', 'def2-SV(P)', 'def-TZV(P)', 'def2-TZV(P)', 
        'def-SVP', 'def2-SVP', 'def-TZVP', 'def2-TZVP'
    ]
    available_functionals = [ 
        'slater-dirac-exchange', 's-vwn', 'vwn', 's-vwn_Gaussian', 'pwlda', 
        'becke-exchange', 'b-lyp', 'b-vwn', 'lyp', 'b-p', 'pbe', 'tpss', 
        'bh-lyp', 'b3-lyp', 'b3-lyp_Gaussian', 'pbe0', 'tpssh', 'lhf', 'oep', 
        'b97-d', 'b2-plyp' 
    ]
    tm_files = [
        'control', 'coord', 'basis', 'auxbasis', 'energy', 'gradient', 'mos',
        'alpha', 'beta',  'statistics', 'GEO_OPT_CONVERGED', 'GEO_OPT_FAILED',
        'not.converged', 'nextstep', 'hessapprox', 'job.last', 'job.start',
        'optinfo', 'statistics', 'converged', 'vibspectrum', 'vib_normal_modes',
        'hessian', 'dipgrad'
    ]
    tm_tmp_files = [
        'errvec', 'fock', 'oldfock', 'dens', 'ddens', 'diff_densmat', 
        'diff_dft_density', 'diff_dft_oper', 'diff_fockmat', 'diis_errvec', 
        'diis_oldfock'
    ]
    default_parameters = {
        # general and geometry
        'title': '',
        'point group': 'c1',
        'use redundant internals': False,
        # basis set
        'use basis set library': True,          # case False not implemented
        'basis set name': 'def-SV(P)',          # use the define default basis set
        'basis set definition': None,           # not implemented
        # initial guess and occupation numbers
        'initial guess': 'eht',                 # other than 'eht' not implemented
        'total charge': 0,                      # default is neutral system
        'multiplicity': 1,                      # no default
        'uhf': None,                            # enforce UHF calculation
        'rohf': None,                           # enforce ROHF calculation
        # method
        'use dft': True,                        # dft on/off
        'density functional': 'b-p',            # func
        'grid size': 'm3',                      # gridsize
        'use resolution of identity': False,    # ri on/off
        'ri memory': 1000,                      # m <int>
        # scf parameters
        'use fermi smearing': False,            # fermi on/off
        'fermi initial temperature': 300,
        'fermi final temperature': 300,
        'fermi annealing factor': 0.95,
        'fermi homo-lumo gap criterion': 0.1,
        'fermi stopping criterion': 0.001,
        'maximum number of scf iterations': 60, # scfiterlimit
        'scf energy convergence': None,         # scfconv
        'density convergence': None,            # denconv
        'orbital shift type': None,             # not implemented
        'orbital shift energy': None,           # not implemented
        'initial damping': None,                # not implemented
        'minimal damping': None,                # not implemented
        'damping adjustment step': None,        # not implemented
        # task
        'ground state': True,
        'excited state': False,                 # not implemented
        'number of excited states': None,       # not implemented
        'optimized excited state': None,        # not implemented
        'force convergence': None,              # jobex -gcart <int>
        'energy convergence': None,             # jobex -energy <int>
        'number of geometry cycles': None,      # jobex -c <int>
        'task': 'energy',           # 'energy calculation', 'energy'
                                    # 'gradient calculation', 'gradient'
                                    # 'geometry optimization', 'optimize'
                                    # 'harmonic mode analysis', 'frequencies'
        # analysis
        # advanced
    }
    parameters = {}
    results = {}
    converged = False
    updated = False
    atoms = None
    forces = None
    e_total = None

    def __init__(self, restart=False, define_str=None, label='turbomole',
                 calculate_energy='dscf', calculate_forces='grad',
                 post_HF=False, **kwargs):
        """ calculation restart not yet implemented """

        self.restart = restart
        self.define_str = define_str
        self.label = label
        self.calculate_energy = calculate_energy
        self.calculate_forces = calculate_forces
        self.post_HF = post_HF
        self.set_parameters(kwargs)
        self.verify_parameters()
        self.reset()

    def __getitem__(self, item):
        if hasattr(self, item):
            obj = getattr(self, item)
        else:
            obj = None
        return obj

    def set_parameters(self, params):
        self.parameters = self.default_parameters
        self.parameters.update(params)
        if self.parameters['use resolution of identity']:
            self.calculate_energy = 'ridft'
            self.calculate_forces = 'rdgrad'

    def verify_parameters(self):
        """ detect wrong or not implemented parameters """

        func_list = [x.lower() for x in self.available_functionals]
        assert self.parameters['density functional'].lower() in func_list, ( 
            'density functional not available / not supported'
        )

        bas_list = [x.lower() for x in self.available_basis_sets]
        assert self.parameters['basis set name'].lower() in bas_list, (
            'basis set ', self.parameters['basis set name'], 
            ' not available / not supported'
        )
        assert self.parameters['multiplicity'], 'multiplicity not defined'

        if self.parameters['rohf']:
            raise NotImplementedError('ROHF not implemented')
        if self.parameters['initial guess'] != 'eht':
            raise NotImplementedError('Initial guess not implemented')
        if not self.parameters['use basis set library']:
            raise NotImplementedError('Explicit basis set definition')
        if self.parameters['point group'] != 'c1':
            raise NotImplementedError('Point group not impemeneted')
        if self.parameters['excited state']:
            raise NotImplementedError('Excited state not implemented')

    def reset(self):
        """ removes all turbomole input, output and scratch files, 
        and deletes results dict and the atoms object """
        self.atoms = None
        self.results = {}
        self.results['calculation parameters'] = {}
        for f in self.tm_files + self.tm_tmp_files:
            if os.path.exists(f):
                os.remove(f)
        self.initialized = False
        self.converged = False

    def set_atoms(self, atoms):
        """ Create the self.atoms object and writes the coord file. If 
        self.atoms exists a check for changes and an update of the atoms
        are performed. Note: Only positions changes are tracked in this version. 
        """
        changes = self.check_state(atoms, tol=1e-13)
        if self.atoms == atoms or 'positions' not in changes:
            # print('two atoms obj are (almost) equal')
            if (self.updated and os.path.isfile('coord')):
                self.updated = False
                a = read('coord').get_positions()
                if np.allclose(a, atoms.get_positions(), rtol=0, atol=1e-13):
                    return
            else:
                return

        changes = self.check_state(atoms, tol=1e-2)
        if 'positions' in changes:
            # print(two atoms obj are different')
            self.reset()
        else:
            # print('two atoms obj are slightly different')
            if self.parameters['use redundant internals']:
                self.reset()

        write('coord', atoms)
        self.atoms = atoms.copy()
        self.update_energy = True
        self.update_forces = True
        self.update_geometry = True
        self.update_hessian = True

    def get_define_str(self):
        define_str_tpl = (
            '\n__title__\na coord\n__inter__\n'
            'bb all __basis_set__\n*\neht\ny\n__charge_str____occ_str__'
            '__single_atom_str____norb_str____dft_str____ri_str__'
            '__scfiterlimit____fermi_str__q\n'
        )

        params = self.parameters

        if params['use redundant internals']:
            internals_str = 'ired\n*'
        else:
            internals_str = '*\nno'
        charge_str = str(params['total charge']) + '\n'

        if params['multiplicity'] == 1:
            if params['uhf']:
                occ_str = 'n\ns\n*\n'
            else:
                occ_str = 'y\n'
        elif params['multiplicity'] == 2:
            occ_str = 'y\n'
        elif params['multiplicity'] == 3:
            occ_str = 'n\nt\n*\n'
        else:
            unpaired = params['multiplicity'] - 1
            if params['use fermi smearing']:
                occ_str = 'n\nuf ' + str(unpaired) + '\n*\n'
            else:
                occ_str = 'n\nu ' + str(unpaired) + '\n*\n'

        if len(self.atoms) != 1:
            single_atom_str = ''
        else:
            single_atom_str = '\n'

        if params['multiplicity'] == 1:
            norb_str = ''
        else:
            norb_str='n\n'

        if params['use dft']:
            dft_str = 'dft\non\n*\n'
        else:
            dft_str = ''

        if params['density functional']:
            dft_str += 'dft\nfunc ' + params['density functional'] + '\n*\n'

        if params['grid size']:
            dft_str += 'dft\ngrid ' + params['grid size'] + '\n*\n'

        if params['use resolution of identity']:
            ri_str = 'ri\non\nm ' + str(params['ri memory']) + '\n*\n'
        else:
            ri_str =''

        if params['maximum number of scf iterations']:
            scfmaxiter = params['maximum number of scf iterations']
            scfiter_str = 'scf\niter\n'+str(scfmaxiter)+'\n\n'
        else:
            scfiter_str = ''
        if params['scf energy convergence']:
            conv = -log10(params['scf energy convergence']/Hartree)//1
            scfiter_str += 'scf\nconv\n'+str(int(conv))+'\n\n'

        if params['use fermi smearing']:
            fermi_str = 'scf\nfermi\n'
            if params['fermi initial temperature']:
                fermi_str+=('1\n'+str(params['fermi initial temperature'])+'\n')
            if params['fermi final temperature']:
                fermi_str+=('2\n'+str(params['fermi final temperature'])+'\n')
            if params['fermi annealing factor']:
                fermi_str+=('3\n'+str(params['fermi annealing factor'])+'\n') 
            if params['fermi homo-lumo gap criterion']:
                fermi_str+=('4\n'+str(params['fermi homo-lumo gap criterion'])+'\n') 
            if params['fermi stopping criterion']:
                fermi_str+=('5\n'+str(params['fermi stopping criterion'])+'\n') 
            fermi_str += '\n\n'
        else:
            fermi_str = ''

        define_str = define_str_tpl
        define_str = re.sub('__title__', params['title'], define_str)
        define_str = re.sub('__basis_set__', params['basis set name'], define_str)
        define_str = re.sub('__charge_str__', charge_str, define_str)
        define_str = re.sub('__occ_str__', occ_str, define_str)
        define_str = re.sub('__norb_str__', norb_str, define_str)
        define_str = re.sub('__dft_str__', dft_str, define_str)
        define_str = re.sub('__ri_str__', ri_str, define_str)
        define_str = re.sub('__single_atom_str__', single_atom_str, define_str)
        define_str = re.sub('__inter__', internals_str, define_str)
        define_str = re.sub('__scfiterlimit__', scfiter_str, define_str)
        define_str = re.sub('__fermi_str__', fermi_str, define_str)

        return define_str

    def initialize(self):
        """ prepare turbomole control file by running module 'define' """
        if self.initialized:
            return
        self.verify_parameters()
        if not self.atoms:
            raise RuntimeError('atoms missing during initialization')
        if not os.path.isfile('coord'):
            raise IOError('file coord not found')

        if self.define_str:
            assert isinstance(self.define_str, str)
            define_str = self.define_str
        else:
            define_str = self.get_define_str()
        # print(define_str)

        # run define
        try:
            command = 'define > ASE.TM.define.out'
            p = Popen(command, shell=True, stdin=PIPE, stderr=PIPE)
            error = p.communicate(input=define_str)[1]
            if 'abnormally' in error:
                raise OSError(error)
            print('TM command:  ' + command + ' successfully executed')
        except OSError('define execution failed: ') as err:
            raise err

        # add or delete data groups
        self.execute('kdg scfdump')
        if self.parameters['density convergence']:
            if len(self.read_data_group('denconv')) != 0:
                self.execute('kdg denconv')
            conv = -log10(params['density convergence'])
            self.add_data_group('denconv', str(int(conv))) 

        self.initialized = True
        self.converged = False

    def calculation_required(self, atoms, properties):
        if self.atoms != atoms:
            return True
        for prop in properties:
            if prop == 'energy' and self.e_total is None:
                return True
            elif prop == 'forces' and self.forces is None:
                return True
        return False

    def calculate(self, atoms):
        """ execute the requested job """

        if self.parameters['task'] in ['energy', 'energy calculation']:
            self.get_potential_energy(atoms)
        if self.parameters['task'] in ['gradient', 'gradient calculation']:
            self.get_forces(atoms)
        if self.parameters['task'] in ['optimize', 'geometry optimization']:
            self.relax_geometry(atoms)
        if self.parameters['task'] in ['frequencies', 'harmonic mode analysis']:
            self.harmonic_analysis(atoms)
        self.read_results()

    def execute(self, command):
        from subprocess import Popen, PIPE
        try:
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            # check the error output
            if 'abnormally' in error:
                raise OSError(error)
            print('TM command: ', command, 'successfully executed')
        except OSError as e:
            print('Execution failed:', e, file=sys.stderr)
            sys.exit(1)

    def relax_geometry(self, atoms):
        """ execute geometry optimization with script jobex """

        self.set_atoms(atoms)
        if self.converged and not self.update_geometry:
            return
        self.initialize()
        jobex_flags = ''
        if self.parameters['use resolution of identity']:
            jobex_flags += ' -ri'
        if self.parameters['force convergence']:
            conv = -log10(self.parameters['force convergence']/Hartree*Bohr)//1
            jobex_flags += ' -gcart ' + str(int(conv))
        if self.parameters['energy convergence']:
            conv = -log10(self.parameters['energy convergence']/Hartree)//1
            jobex_flags += ' -energy ' + str(int(conv))
        self.converged = False
        self.execute('jobex ' + jobex_flags + ' > ASE.TM.jobex.out')
        new_struct = read('coord')
        atoms.set_positions(new_struct.get_positions())
        self.atoms = atoms.copy()
        self.read_energy()
        if os.path.exists('GEO_OPT_CONVERGED'):
            self.update_energy = False
            self.update_forces = False
            self.update_geometry = False
            self.update_hessian = True
            self.converged = True
        else:
            raise RuntimeError('Geometry not converged!')

    def harmonic_analysis(self, atoms):
        """ execute harmonic mode analysis with module aoforce """

        self.set_atoms(atoms)
        self.initialize()
        if self.update_energy:
            self.get_potential_energy(atoms)
        if self.update_hessian:
            self.execute('aoforce > ASE.TM.aoforce.out')
            self.update_hessian = False

    def read_results(self):
        """ read all results and load them in the results entity """
        import util as abcd
        self.results['atoms'] = abcd.atoms2dict(self.atoms,plain_arrays=True)
        self.read_energy()
        self.results['total energy'] = self.e_total
        self.read_mos()
        self.read_basis_set()
        self.read_occupation_numbers()
        self.read_dipole_moment ()
        self.read_ssquare()
        self.read_calc_parameters()
        if self.parameters['task'] in ['gradient', 'optimize', 
                'gradient calculation', 'geometry optimization']:
            self.read_gradient()
            self.read_forces()
            self.results['energy gradient'] = (-self.forces).tolist()
        if self.parameters['task'] in ['frequencies', 'harmonic mode analysis']:
            self.read_hessian()
            self.read_vibrational_reduced_masses()
            self.read_normal_modes()
            self.read_vibrational_spectrum()

    """ methods to work with turbomole data groups """

    def read_data_group(self, data_group):
        """ read a turbomole data group from control file """
        args = ['sdg', data_group]
        FNULL = open(os.devnull, 'w')
        p = Popen(args, stdout=PIPE, stdin=None, stderr=FNULL)
        return p.communicate()[0].strip('\n')


    def add_data_group(self, data_group, string=None, raw=False):
        """ write a turbomole data group to control file """
        if raw:
            data = data_group
        else:
            data = '$' + data_group
            if string:
                data += ' ' + string
            data += '\n'
        f = open('control', 'rw+')
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        lines.insert(2,data)
        f.write(''.join(lines))
        f.close()

    """ data reader methods """

    def read_calc_parameters(self):
        """ put things like symmetry, number of electrons etc. here """
        if 'calculation parameters' not in self.results.keys():
            self.results['calculation parameters'] = {}
        parameters = self.results['calculation parameters']
        parameters['point group'] = str(self.read_data_group('symmetry').split()[1])
        if '$uhf' in self.read_data_group('uhf'):
            parameters['uhf'] = True
        else:
            parameters['uhf'] = False
        # Gaussian function type
        gt = self.read_data_group('pople')
        if gt == '':
            parameters['gaussian type'] = 'spherical harmonic'
        else:
            gt = gt.split()[1]
            if gt == 'AO':
                parameters['gaussian type'] = 'spherical harmonic'
            elif gt == 'CAO':
                parameters['gaussian type'] = 'cartesian'
            else:
                parameters['gaussian type'] = None

        nvibro = self.read_data_group('nvibro')
        if nvibro:
            parameters['nuclear degrees of freedom'] = int(nvibro.split()[1])

    def read_energy(self):
        """Read Energy from Turbomole energy file."""
        text = open('energy', 'r').read().lower()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            if line.startswith('$end'):
                break
            elif line.startswith('$'):
                pass
            else:
                energy_tmp = float(line.split()[1])
                if self.post_HF:
                    energy_tmp += float(line.split()[4])
        # update energy units
        self.e_total = energy_tmp * Hartree

    def read_forces(self):
        """Read Forces from Turbomole gradient file."""
        file = open('gradient', 'r')
        lines = file.readlines()
        file.close()

        forces = np.array([[0, 0, 0]])
        
        nline = len(lines)
        iline = -1
        
        for i in range(nline):
            if 'cycle' in lines[i]:
                iline = i
        
        if iline < 0:
            raise RuntimeError('Please check TURBOMOLE gradients')

        # next line
        iline += len(self.atoms) + 1
        # $end line
        nline -= 1
        # read gradients
        for i in range(iline, nline):
            line = lines[i].replace('D', 'E')
            tmp = np.array([[float(f) for f in line.split()[0:3]]])
            forces = np.concatenate((forces, tmp))
        # Note the '-' sign for turbomole, to get forces
        self.forces = (-np.delete(forces, np.s_[0:1], axis=0)) * Hartree / Bohr

    def read_occupation_numbers(self):
        """ read occupation numbers with module 'eiger' """
        if 'molecular orbitals' not in self.results.keys():
            return
        mos = self.results['molecular orbitals']
        args = ['eiger', '--all', '--pview']
        p = Popen(args, stdout=PIPE, stdin=None, stderr=None)
        stdout = p.communicate()
        lines = stdout[0].split('\n')
        for line in lines:
            regex = (
                '^\s+(\d+)\.*\s+(\w*)\s+(\d+)\s+(\S+)'
                '\s+(\d*\.*\d*)\s+([-+]?\d+\.\d*)'
            )
            match = re.search(regex, line)
            properties = {}
            if match:
                orb_index = int(match.group(3))
                if match.group(2) == 'a':
                    spin = 'alpha'
                elif match.group(2) == 'b':
                    spin = 'beta'
                else:
                    spin = None
                ar_index = next(
                    index for (index, molecular_orbital) in enumerate(mos) 
                    if (
                        molecular_orbital['index'] == orb_index 
                        and molecular_orbital['spin'] == spin
                    )
                )
                mos[ar_index]['index by energy'] = int(match.group(1))
                mos[ar_index]['irreducible representation'] = str(match.group(4))
                if match.group(5) != '':
                    mos[ar_index]['occupancy'] = float(match.group(5))
                else:
                    mos[ar_index]['occupancy'] = float(0)

    def read_mos(self):
        """ read the molecular orbital coefficients and orbital energies
        from files mos, alpha and beta """

        self.results['molecular orbitals'] = []
        mos = self.results['molecular orbitals']
        keywords = [ 'scfmo', 'uhfmo_alpha', 'uhfmo_beta' ]
        spin = [ None, 'alpha', 'beta' ]

        for index, keyword in enumerate(keywords):
            f_string = None
            f_width = 0
            mo = {}
            orbitals_coefficients_line=[]
            mo_string = self.read_data_group(keyword)
            if mo_string == '': continue
            mo_string += '\n$end'
            lines = mo_string.split('\n')
            for line in lines:
                if re.match('^\s*#',line):
                    continue
                if 'eigenvalue' in line:
                    if len(orbitals_coefficients_line) != 0:
                        mo['eigenvector'] = orbitals_coefficients_line
                        mos.append(mo)
                        mo = {}
                        orbitals_coefficients_line = []
                    match = re.search(
                        '^\s*(\d+)\s+(\S+)\s+eigenvalue=([\+\-\d\.\w]+)\s',
                        line
                    )
                    mo['index'] = int(match.group(1))
                    mo['irreducible representation'] = str(match.group(2))
                    mo['eigenvalue'] = float(re.sub('[dD]','E',match.group(3)))*Hartree
                    mo['spin'] = spin[index]
                    mo['degeneracy'] = 1
                    continue
                if keyword in line:
                    # e.g. format(4d20.14)
                    match = re.search('format\(\d+([a-zA-Z](\d+)\.\d+)\)',line)
                    if match:
                        f_string = match.group(1)
                        f_width = int(match.group(2))
                    continue
                if '$end' in line:
                    if len(orbitals_coefficients_line) != 0:
                        mo['eigenvector'] = orbitals_coefficients_line
                        mos.append(mo)
                    break
                fort_str = str(len(line.rstrip())/f_width) + f_string
                r = ff.FortranRecordReader(fort_str)
                orbitals_coefficients_line += r.read(line)

    def read_basis_set(self):
        """ read basis set, ecp not supported yet """
        self.results['basis set'] = []
        self.results['basis set formatted'] = {}
        bsf = self.read_data_group('basis')
        self.results['basis set formatted']['turbomole'] = bsf
        lines = bsf.split('\n')
        basis_set = {}
        functions = []
        function = {}
        primitives = []
        read_tag = False
        read_data = False
        for line in lines:
            if len(line) == 0:
                continue
            if '$basis' in line:
                continue
            if '$end' in line:
                break
            if re.match('^\s*#',line):
                continue
            if re.match('^\s*\*',line):
                if read_tag:
                    read_tag = False
                    read_data = True
                else:
                    if read_data:
                        # end primitives
                        function['primitive functions'] = primitives
                        function['number of primitives'] = len(primitives)
                        primitives = []
                        functions.append(function)
                        function = {}
                        # end contracted
                        basis_set['functions'] = functions
                        functions = []
                        self.results['basis set'].append(basis_set)
                        basis_set = {}
                        read_data = False
                    read_tag = True
                continue
            if read_tag:
                match = re.search('^\s*(\w+)\s+(.+)',line)
                if match:
                    basis_set['element'] = match.group(1)
                    basis_set['nickname'] = match.group(2)
                else:
                    raise RuntimeError("error reading basis set")
            else:
                match = re.search('^\s+(\d+)\s+(\w+)',line)
                if match:
                    if len(primitives) is not 0:
                        # end primitives
                        function['primitive functions'] = primitives
                        function['number of primitives'] = len(primitives)
                        primitives = []
                        functions.append(function)
                        function = {}
                        # begin contracted
                    function['shell type'] = str(match.group(2))
                    continue
                regex = (
                    '^\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
                    '\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
                )
                match = re.search(regex, line)
                if match:
                    exponent = float(match.group(1))
                    coefficient = float(match.group(3))
                    primitives.append(
                        {'exponent': exponent, 'coefficient': coefficient}
                    )

    def read_gradient (self):
        """ read all information in file 'gradient' """
        from ase import Atom
        import util as abcd
        grad_string = self.read_data_group('grad')
#       there are units errors in the ASE method: energy in Hartree and forces 
#       are multiplied by Bohr;
#       update: the bug has been fixed, try to reuse ase
#       structures = read('gradient',index=':') 
        lines = grad_string.split('\n')
        history = []
        image = {}
        gradient = []
        atoms = Atoms()
        for line in lines:
            # cycle lines
            regex = (
                '^\s*cycle =\s*(\d+)'
                '\s+SCF energy =\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
                '\s+\|dE\/dxyz\| =\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
            )
            match = re.search(regex, line)
            if match:
                if len(atoms):
                    image['optimization cycle'] = cycle
                    image['total energy'] = energy
                    image['gradient norm'] = norm
                    image['atoms'] = abcd.atoms2dict(atoms,plain_arrays=True)
                    image['energy gradient'] = gradient
                    history.append(image)
                    image = {}
                    atoms = Atoms()
                    gradient = []
                cycle = int(match.group(1))
                energy = float(match.group(2))*Hartree
                norm = float(match.group(4))*Hartree/Bohr
                continue
            # coordinate lines
            regex = (
                '^\s*([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                '\s+([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                '\s+([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                '\s+(\w+)'
            )
            match = re.search(regex, line)
            if match:
                x = float(match.group(1))*Bohr
                y = float(match.group(3))*Bohr
                z = float(match.group(5))*Bohr
                symbol = str(match.group(7))
                atoms += Atom(symbol.capitalize(), (x,y,z))
                continue
            # gradient lines
            regex = (
                '^\s*([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                '\s+([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                '\s+([-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)'
                )
            match = re.search(regex, line)
            if match:
                gradx = float(match.group(1).replace('D','E'))*Hartree/Bohr
                grady = float(match.group(3).replace('D','E'))*Hartree/Bohr
                gradz = float(match.group(5).replace('D','E'))*Hartree/Bohr
                gradient.append([ gradx, grady, gradz ])

        image['optimization cycle'] = cycle
        image['total energy'] = energy
        image['gradient norm'] = norm
        image['atoms'] = abcd.atoms2dict(atoms,plain_arrays=True)
        image['energy gradient'] = gradient
        history.append(image)
        self.results['geometry optimization history'] = history

    def read_hessian(self, noproj=False, nomw=False):
        """ Read in the hessian matrix. The keys $nomw and $noproj are planned 
        for support. """
        self.results['hessian matrix'] = {}
        self.results['hessian matrix']['array'] = []
        self.results['hessian matrix']['units'] = '?'
        self.results['hessian matrix']['projected'] = True
        self.results['hessian matrix']['mass weighted'] = True
        nvibro = int(self.read_data_group('nvibro').split()[1])
        self.results['hessian matrix']['dimension'] = nvibro
        row = []
        key = 'hessian'
        if noproj:
            key = 'npr' + key
            self.results['hessian matrix']['projected'] = False
        lines = self.read_data_group(key).split('\n')
        for line in lines:
            if key in line:
                continue
            fields = line.split()
            row.extend(fields[2:len(fields)])
            if len(row) == nvibro:
                # check whether it is mass-weighted
                float_row = [ float(element) for element in row ]
                self.results['hessian matrix']['array'].append(float_row)
                row = []

    def read_normal_modes(self, noproj=False, nomw=False):
        """ Read in vibrational normal modes. The keays $nomw and $noproj are 
        planned for support. """
        self.results['normal modes'] = {}
        self.results['normal modes']['array'] = []
        self.results['normal modes']['projected'] = True
        self.results['normal modes']['mass weighted'] = True
        self.results['normal modes']['units'] = '?'
        nvibro = int(self.read_data_group('nvibro').split()[1])
        self.results['normal modes']['dimension'] = nvibro
        row = []
        key = 'vibrational normal modes'
        if noproj:
            key = 'npr' + key
            self.results['normal modes']['projected'] = False
        lines = self.read_data_group(key).split('\n')
        for line in lines:
            if key in line:
                continue
            if '$end' in line:
                break
            fields = line.split()
            row.extend(fields[2:len(fields)])
            if len(row) == nvibro:
                # check whether it is mass-weighted
                float_row = [ float(element) for element in row ]
                self.results['normal modes']['array'].append(float_row)
                row = []

    def read_vibrational_reduced_masses(self):
        """ Read vibrational reduced masses """
        nvibro = int(self.read_data_group('nvibro').split()[1])
        self.results['vibrational reduced masses'] = []
        lines = self.read_data_group('vibrational reduced masses').split('\n')
        for line in lines:
            if '$vibrational' in line:
                continue
            if '$end' in line:
                break
            fields = [ float(element) for element in line.split() ]
            self.results['vibrational reduced masses'].extend(fields)

    def read_vibrational_spectrum(self, noproj=False):
        """ Read the vibrational spectrum """
        self.results['vibrational spectrum'] = []
        key = 'vibrational spectrum'
        if noproj:
            key = 'npr' + key
        lines = self.read_data_group(key).split('\n')
        for line in lines:
            dictionary = {}
            regex = (
                '^\s+(\d+)\s+(\S*)\s+([-+]?\d+\.\d*)'
                '\s+(\d+\.\d*)\s+(\S+)\s+(\S+)'
            )
            match = re.search(regex, line)
            if match:
                dictionary['mode number'] = int(match.group(1))
                dictionary['irreducible representation'] = str(match.group(2))
                dictionary['frequency'] = {
                    'units': 'cm^-1', 
                    'value': float(match.group(3))
                }
                dictionary['infrared intensity'] = {
                    'units': 'km/mol', 
                    'value': float(match.group(4))
                }

                if match.group(5) == 'YES':
                    dictionary['infrared active'] = True
                elif match.group(5) == 'NO':
                    dictionary['infrared active'] = False
                else:
                    dictionary['infrared active'] = None

                if match.group(6) == 'YES':
                    dictionary['Raman active'] = True
                elif match.group(6) == 'NO':
                    dictionary['Raman active'] = False
                else:
                    dictionary['Raman active'] = None

                self.results['vibrational spectrum'].append(dictionary)

    def read_ssquare(self):
        """ Read the expectation value of S^2 operator """
        s2_string = self.read_data_group('ssquare from dscf')
        if s2_string == '':
            return
        string = s2_string.split('\n')[1]
        ssquare = float(re.search('^\s*(\d+\.*\d*)', string).group(1))
        self.results['ssquare from scf calculation'] = ssquare

    def read_dipole_moment(self):
        """ Read the dipole moment """
        dip_string = self.read_data_group('dipole')
        if dip_string == '':
            return
        lines = dip_string.split('\n')
        for line in lines:
            regex = (
                '^\s+x\s+([-+]?\d+\.\d*)\s+y\s+([-+]?\d+\.\d*)'
                '\s+z\s+([-+]?\d+\.\d*)\s+a\.u\.'
            )
            match = re.search(regex, line)
            if match:
                dip_vec = [ float(match.group(c)) for c in range(1,4)]
            match = re.search('^\s+\| dipole \| =\s+(\d+\.*\d*)\s+debye', line)
            if match:
                dip_abs_val = float(match.group(1))
        self.results['electric dipole moment'] = {}
        self.results['electric dipole moment']['vector'] = {
            'array': dip_vec, 
            'units': 'a.u.'
        }
        self.results['electric dipole moment']['absolute value'] = {
            'value': dip_abs_val, 
            'units':'Debye'
        }
        self.dipole = np.array(dip_vec)*Bohr

    """ getter methods """

    def get_results(self):
        return self.results

    def get_potential_energy(self, atoms, force_consistent=True):
        # update atoms
        self.updated = self.e_total is None
        self.set_atoms(atoms)
        self.initialize()
        # if update of energy is necessary
        if self.update_energy:
            # calculate energy
            self.execute(self.calculate_energy + ' > ASE.TM.energy.out')
            # check for convergence of dscf cycle
            if os.path.isfile('dscf_problem'):
                print('Turbomole scf energy calculation did not converge')
                raise RuntimeError(
                    'Please run Turbomole define and come thereafter back')
            # read energy
            self.read_energy()

        self.update_energy = False
        return self.e_total

    def get_forces(self, atoms):
        # update atoms
        self.updated = self.forces is None
        self.set_atoms(atoms)
        # complete energy calculations
        if self.update_energy:
            self.get_potential_energy(atoms)
        # if update of forces is necessary
        if self.update_forces:
            # calculate forces
            self.execute(self.calculate_forces + ' > ASE.TM.forces.out')
            # read forces
            self.read_forces()

        self.update_forces = False
        return self.forces.copy()

    def get_dipole_moment(self, atoms):
        # this must check the state and then perform a calc if necessary
        return self.dipole


    """
    The following three functions are necessary for the atoms2dict function. 
    They cause a lot of code bloat because most of the code is identical for all 
    calculators, so why not in the base class?
    """

    def todict(self):
        dct = {}
        lst = [
            attr for attr in dir(self) if not attr.startswith('__') 
            and not callable(getattr(self,attr))
            ]
        for item in lst:
            obj = getattr(self,item)
            if type(obj) in [bool, int, float, str, None]:
                dct[item] = obj
            if isinstance(obj, (list, tuple, dict)):
                dct[item] = obj
            elif hasattr(obj, 'todict'):
                dct[item] = obj.todict()
            elif obj.__class__ == np.ndarray:
                dct[item] = obj.tolist()
            else:
                pass
        return dct


    def check_state(self, atoms, tol=1e-15):
        """Check for system changes since last calculation."""
        from ase.calculators.calculator import all_changes, equal
        if self.atoms is None:
            system_changes = all_changes[:]
        else:
            system_changes = []
            if not equal(self.atoms.positions, atoms.positions, tol):
                system_changes.append('positions')
            if not equal(self.atoms.numbers, atoms.numbers):
                system_changes.append('numbers')
            if not equal(self.atoms.cell, atoms.cell, tol):
                system_changes.append('cell')
            if not equal(self.atoms.pbc, atoms.pbc):
                system_changes.append('pbc')
            if not equal(self.atoms.get_initial_magnetic_moments(),
                         atoms.get_initial_magnetic_moments(), tol):
                system_changes.append('initial_magmoms')
            if not equal(self.atoms.get_initial_charges(),
                         atoms.get_initial_charges(), tol):
                system_changes.append('initial_charges')

        return system_changes


    def get_property(self, name, atoms=None, allow_calculation=True):
        """Returns the value of a property"""

        if name not in self.implemented_properties:
            # an ugly work around; the caller should test the raised error
            if name in ['magmom', 'magmoms', 'charges', 'stress']:
                return None
            raise NotImplementedError(name)

        if atoms is None:
            atoms = self.atoms.copy()

        persist_property = {
            'energy': 'e_total',
            'forces': 'forces',
            'dipole': 'dipole',
            'free_energy': 'e_total'
        }
        property_getter = {
            'energy': self.get_potential_energy,
            'forces': self.get_forces,
            'dipole': self.get_dipole_moment,
            'free_energy': self.get_potential_energy
            }
        getter_args = {
            'energy': [atoms],
            'forces': [atoms],
            'dipole': [atoms],
            'free_energy': [atoms, True]
            }

        if allow_calculation:
            result = property_getter[name](*getter_args[name])
        else:
            if hasattr(self, persist_property[name]):
                result = getattr(self, persist_property[name])
            else:
                result = None

        if isinstance(result, np.ndarray):
            result = result.copy()
        return result

