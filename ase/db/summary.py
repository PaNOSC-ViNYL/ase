from __future__ import print_function
import os.path as op

from ase.data import atomic_masses, chemical_symbols
from ase.db.core import float_to_time_string, now
from ase.utils import hill


class Summary:
    def __init__(self, row, meta={}, subscript=None, prefix=None, tmpdir=None):
        self.row = row

        self.cell = [['{0:.3f}'.format(a) for a in axis] for axis in row.cell]

        forces = row.get('constrained_forces')
        if forces is None:
            fmax = None
            self.forces = None
        else:
            fmax = (forces**2).sum(1).max()**0.5
            N = len(forces)
            self.forces = []
            for n, f in enumerate(forces):
                if n < 5 or n >= N - 5:
                    f = tuple('{0:10.3f}'.format(x) for x in f)
                    symbol = chemical_symbols[row.numbers[n]]
                    self.forces.append((n, symbol) + f)
                elif n == 5:
                    self.forces.append((' ...', '',
                                        '       ...',
                                        '       ...',
                                        '       ...'))

        self.stress = row.get('stress')
        if self.stress is not None:
            self.stress = ', '.join('{0:.3f}'.format(s) for s in self.stress)

        if 'masses' in row:
            mass = row.masses.sum()
        else:
            mass = atomic_masses[row.numbers].sum()

        self.formula = hill(row.numbers)
        if subscript:
            self.formula = subscript.sub(r'<sub>\1</sub>', self.formula)

        age = float_to_time_string(now() - row.ctime, True)

        table = dict((key, value)
                     for key, value in [
                         ('id', row.id),
                         ('age', age),
                         ('formula', self.formula),
                         ('user', row.user),
                         ('calculator', row.get('calculator')),
                         ('energy', row.get('energy')),
                         ('fmax', fmax),
                         ('charge', row.get('charge')),
                         ('mass', mass),
                         ('magmom', row.get('magmom')),
                         ('unique id', row.unique_id),
                         ('volume', row.get('volume'))]
                     if value is not None)

        table.update(row.key_value_pairs)

        # If meta data for summary_sections does not exists a default
        # template is generated otherwise it goes through the meta
        # data and checks if all keys are indeed present

        kd = meta.get('key_descriptions', {})

        self.layout = []
        for headline, blocks in meta['layout']:
            newblocks = []
            print(blocks)
            for block in blocks:
                if block is None:
                    pass
                elif isinstance(block, tuple):
                    title, keys = block
                    rows = []
                    for key in keys:
                        value = table.pop(key, None)
                        if value is not None:
                            desc, unit = kd.get(key, [0, key, 0, ''])[1::2]
                            rows.append((desc, value, unit))
                    block = (title, rows)
                elif block.endswith('.png'):
                    name = op.join(tmpdir, prefix + '-' + block)
                    if op.isfile(name):
                        if op.getsize(name) == 0:
                            block = None
                    else:
                        for func in meta['functions']:
                            func(prefix, tmpdir, row)

                newblocks.append(block)
            self.layout.append((headline, newblocks))

        if table:
            rows = []
            for key, value in sorted(table.items()):
                desc, unit = kd.get(key, [0, key, 0, ''])[1::2]
                rows.append((desc, value, unit))
            self.layout.append(('Other stuff', [('Things', rows)]))

        self.dipole = row.get('dipole')
        if self.dipole is not None:
            self.dipole = ', '.join('{0:.3f}'.format(d) for d in self.dipole)

        self.data = row.get('data')
        if self.data:
            self.data = ', '.join(self.data.keys())

        self.constraints = row.get('constraints')
        if self.constraints:
            self.constraints = ', '.join(d['name'] for d in self.constraints)

    def write(self):
        row = self.row

        width = max(len(name) for name, unit, value in self.table)
        print('{0:{width}}|unit  |value'.format('name', width=width))
        for name, unit, value in self.table:
            print('{0:{width}}|{1:6}|{2}'.format(name, unit, value,
                                                 width=width))

        print('\nUnit cell in Ang:')
        print('axis|periodic|          x|          y|          z')
        c = 1
        for p, axis in zip(row.pbc, self.cell):
            print('   {0}|     {1}|{2[0]:>11}|{2[1]:>11}|{2[2]:>11}'.format(
                c, [' no', 'yes'][p], axis))
            c += 1

        if self.key_value_pairs:
            print('\nKey-value pairs:')
            width = max(len(key) for key, value in self.key_value_pairs)
            for key, value in self.key_value_pairs:
                print('{0:{width}}|{1}'.format(key, value, width=width))

        if self.forces:
            print('\nForces in ev/Ang:')
            for f in self.forces:
                print('{0:4}|{1:2}|{2}|{3}|{4}'.format(*f))

        if self.stress:
            print('\nStress tensor (xx, yy, zz, zy, zx, yx) in eV/Ang^3:')
            print('   ', self.stress)

        if self.dipole:
            print('\nDipole moment in e*Ang: ({0})'.format(self.dipole))

        if self.constraints:
            print('\nConstraints:', self.constraints)

        if self.data:
            print('\nData:', self.data)
