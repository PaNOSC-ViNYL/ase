import os

import numpy as np

from ase.units import Bohr, Hartree

elk_parameters = {
    'swidth': Hartree,
    }

class ELK:
    def __init__(self, dir='.', xc=None, kpts=None, tasks=[0], **kwargs):
        for key, value in kwargs.items():
            if key in elk_parameters:
                kwargs[key] /= elk_parameters[key]
        try:
            self.exe = os.environ['ELK']
        except KeyError:
            self.exe = 'elk'

        if xc is not None:
            if 'xctype' in kwargs:
                raise ValueError("You can't use both 'xc' and 'xctype'!")
            else:
                kwargs['xctype'] = {'LDA': 3, # PW92
                                    'PBE': 20,
                                    'REVPBE': 21,
                                    'PBESOL': 22,
                                    'WC06': 26,
                                    'AM05': 30}[xc.upper()]

        if kpts is not None:
            if 'autokpt' in kwargs:
                if kwargs['autokpt']:
                    raise ValueError("You can't use both 'kpts' and 'autokpt'!")
            if 'ngridk' in kwargs:
                raise ValueError("You can't use both 'kpts' and 'ngridk'!")
            if 'vkloff' in kwargs:
                raise ValueError("You can't use both 'kpts' and 'vkloff'!")
            else:
                kwargs['ngridk'] = kpts
                vkloff = []
                for nk in kpts:
                    if nk % 2 == 0:  # shift kpoint away from gamma point
                        vkloff.append(0.5 / nk)
                    else:
                        vkloff.append(0)
                kwargs['vkloff'] = vkloff

        kwargs['tasks'] = tasks

        if 'rmt' in kwargs:
            self.rmt = kwargs['rmt']
            assert len(self.rmt.keys()) == len(list(set(self.rmt.keys()))), 'redundant rmt definitions'
        else:
            self.rmt = None

        self.parameters = kwargs
        if 'rmt' in self.parameters:
            self.parameters.pop('rmt') # this is not an elk keyword!

        self.dir = dir
        self.energy = None

        self.converged = False

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
        self.write(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

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

        self.initialize(atoms)

        assert os.system('cd %s&& %s ' % (self.dir, self.exe)) == 0
        self.read()

        self.converged = True

    def write(self, atoms):
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        fd = open('%s/elk.in' % self.dir, 'w')
        for key, value in self.parameters.items():
            fd.write('%s\n' % key)
            if isinstance(value, bool):
                fd.write('.%s.\n\n' % ('false', 'true')[value])
            elif isinstance(value, (int, float)):
                fd.write('%s\n\n' % value)
            else:
                fd.write('%s\n\n' % ' '.join([str(x) for x in value]))

        fd.write('avec\n')
        for vec in atoms.cell:
            fd.write('%.14f %.14f %.14f\n' % tuple(vec / Bohr))
        fd.write('\n')

        species = {}
        symbols = []
        for a, symbol in enumerate(atoms.get_chemical_symbols()):
            if symbol in species:
                species[symbol].append(a)
            else:
                species[symbol] = [a]
                symbols.append(symbol)
        fd.write('atoms\n%d\n' % len(species))
        scaled = atoms.get_scaled_positions()
        for symbol in symbols:
            fd.write("'%s.in'\n" % symbol)
            fd.write('%d\n' % len(species[symbol]))
            for a in species[symbol]:
                fd.write('%.14f %.14f %.14f 0.0 0.0 0.0\n' % tuple(scaled[a]))

        customspecies = self.rmt
        if customspecies:
            # custom species definitions
            fd.write("\n")
            sfile = os.path.join(os.environ['ELK_SPECIES_PATH'], 'elk.in')
            assert os.path.exists(sfile)
            slines = open(sfile, 'r').readlines()
            # all species must be defined using species keyword
            for s in species.keys():
                if s not in customspecies.keys():
                    # use default rmt for undefined species
                    customspecies.update({s: 0.0})
            # write custom species into elk.in
            skeys = customspecies.keys()
            skeys.sort()
            for s in skeys:
                found = False
                for n, line in enumerate(slines):
                    if line.find("'" + s + "'") > -1:
                        begline = n - 1
                for n, line in enumerate(slines[begline:]):
                    if not line.strip(): # first empty line
                        endline = n
                        found = True
                        break
                assert found
                fd.write("species\n")
                # set rmt on third line
                rmt = customspecies[s]
                assert isinstance(rmt, (float,int))
                if rmt <= 0.0: # relative
                    # split needed because H is defined with comments
                    newrmt = float(slines[begline + 3].split()[0].strip()) + rmt
                else:
                    newrmt = rmt
                slines[begline + 3] = '%6s\n' % str(newrmt)
                for l in slines[begline: begline + endline]:
                    fd.write('%s' % l)
                fd.write("\n")
        else:
            # use default species
            # if sppath is present in elk.in it overwrites species blocks!
            fd.write("sppath\n'%s'\n\n" % os.environ['ELK_SPECIES_PATH'])


    def read(self):
        fd = open('%s/TOTENERGY.OUT' % self.dir, 'r')
        self.energy = float(fd.readlines()[-1]) * Hartree
        # Forces:
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            text = open(INFO_file).read().lower()
            assert 'convergence targets achieved' in text
            assert not 'reached self-consistent loops maximum' in text
            lines = iter(text.split('\n'))
            forces = []
            atomnum = 0
            for line in lines:
                if line.rfind('total force') > -1:
                    forces.append(np.array([float(f) for f in line.split(':')[1].split()]))
                    atomnum =+ 1
            self.forces = np.array(forces)
        else:
            raise RuntimeError
        # Stress
        self.stress = np.empty((3, 3))
