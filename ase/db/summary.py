from __future__ import print_function
import os
import shutil
import os.path as op

from ase.data import atomic_masses, chemical_symbols
from ase.db.core import float_to_time_string, now
from ase.geometry import cell_to_cellpar
from ase.utils import formula_metal, Lock


class Summary:
    def __init__(self, row, meta={}, subscript=None, prefix='', tmpdir='.'):
        self.row = row

        self.cell = [['{:.3f}'.format(a) for a in axis] for axis in row.cell]
        par = ['{:.3f}'.format(x) for x in cell_to_cellpar(row.cell)]
        self.lengths = par[:3]
        self.angles = par[3:]

        self.stress = row.get('stress')
        if self.stress is not None:
            self.stress = ', '.join('{0:.3f}'.format(s) for s in self.stress)

        self.formula = formula_metal(row.numbers)
        if subscript:
            self.formula = subscript.sub(r'<sub>\1</sub>', self.formula)

        kd = meta.get('key_descriptions', {})
        create_layout = meta.get('layout', default_layout)
        self.layout = create_layout(row, kd)

        self.dipole = row.get('dipole')
        if self.dipole is not None:
            self.dipole = ', '.join('{0:.3f}'.format(d) for d in self.dipole)

        self.data = row.get('data')
        if self.data:
            self.data = ', '.join(self.data.keys())

        self.constraints = row.get('constraints')
        if self.constraints:
            self.constraints = ', '.join(c.__class__.__name__
                                         for c in self.constraints)

    def write(self):
        row = self.row

        print(self.formula + ':')
        for headline, columns in self.layout:
            blocks = columns[0]
            if len(columns) == 2:
                blocks += columns[1]
            print((' ' + headline + ' ').center(78, '='))
            for block in blocks:
                if block is None:
                    pass
                elif isinstance(block, tuple):
                    title, keys = block
                    print(title + ':')
                    if not keys:
                        print()
                        continue
                    width = max(len(name) for name, value, unit in keys)
                    print('{:{width}}|value'.format('name', width=width))
                    for name, value, unit in keys:
                        print('{:{width}}|{} {}'.format(name, value, unit,
                                                        width=width))
                    print()
                elif block.endswith('.png'):
                    if op.isfile(block) and op.getsize(block) > 0:
                        print(block)
                    print()
                elif block.endswith('.csv'):
                    if op.isfile(block) and op.getsize(block) > 0:
                        with open(block) as f:
                            print(f.read())
                    print()
                elif block == 'CELL':
                    print('Unit cell in Ang:')
                    print('axis|periodic|          x|          y|          z')
                    c = 1
                    fmt = '   {0}|     {1}|{2[0]:>11}|{2[1]:>11}|{2[2]:>11}'
                    for p, axis in zip(row.pbc, self.cell):
                        print(fmt.format(c, [' no', 'yes'][p], axis))
                        c += 1
                    print('Lengths: {:>10}{:>10}{:>10}'
                          .format(*self.lengths))
                    print('Angles:  {:>10}{:>10}{:>10}\n'
                          .format(*self.angles))
                elif block == 'FORCES' and self.forces is not None:
                    print('\nForces in ev/Ang:')
                    for f in self.forces:
                        print('{:4}|{:2}|{}|{}|{}'.format(*f))
                    print()

        if self.stress:
            print('Stress tensor (xx, yy, zz, zy, zx, yx) in eV/Ang^3:')
            print('   ', self.stress, '\n')

        if self.dipole:
            print('Dipole moment in e*Ang: ({})\n'.format(self.dipole))

        if self.constraints:
            print('Constraints:', self.constraints, '\n')

        if self.data:
            print('Data:', self.data, '\n')


def create_table(row, keys, title, key_descriptions, digits=3):
    table = []
    for key in keys:
        value = row.get(key)
        if value:
            if isinstance(valeu, float):
                value = '{:.{}f}'.format(value, digits)
            elif not isinstance(value, str):
                value = str(value)
            if key == 'age':
                value = float_to_time_string(now() - row.ctime, True)
                desc = 'Age'
                unit = ''
            else:
                desc, unit = key_descriptions.get(key, ['', key, ''])[1:]
            table.append((desc, value, unit))
    return ('table', {'rows': table, 'title': title})


def default_layout(row, key_descriptions):
    # types: (Row, Dict[str, Tuple[str, str, str]])
    # types -> List[Tuple[str, Dict[str, Any]]]
    keys = ['id',
            'energy', 'fmax', 'smax',
            'mass',
            'age']
    table = create_table(row, keys, 'Key Value Pairs', key_descriptions)
    layout = [('Basic properties', [[('atoms', {}),
                                     ('cell', {})],
                                    [table]])]
    keys =
    misc = create_table(row, keys, 'Items', key_descriptions)
    layout.append(('Miscellaneous', [[misc]]))
    return