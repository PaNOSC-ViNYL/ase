from ase.atoms import Atoms
from ase.data import cpk_colors, covalent_radii, chemical_symbols

def write_pov(fileobj, atoms):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    if isinstance(atoms, (list, tuple)):
        atoms = atoms[-1]
        
    fileobj.write("""\
#include "colors.inc"
#include "finish.inc"
camera {
  location <0,0,15>
  look_at <0,0,0>
}
light_source {
  <10,-10,20>
  color White*2
}
background { color White }
#default { finish { phong 0.7 } }
"""
                  )
    numbers = atoms.get_atomic_numbers()
    done = {}
    for Z in numbers:
        if Z not in done:
            symbol = chemical_symbols[Z]
            fileobj.write('#declare C_%s = ' % symbol +
                          'texture{pigment{rgb <%.2f, %.2f, %.2f>}}\n' %
                          tuple(cpk_colors[Z]))
            fileobj.write('#declare R_%s = %.2f;\n' % (symbol, covalent_radii[Z]))
            done[Z] = True

    for Z, p in zip(numbers, atoms.get_positions()):
        symbol = chemical_symbols[Z]
        fileobj.write('sphere{<%.2f, %.2f, %.2f>, ' % tuple(p) +
                      'R_%s texture{C_%s}}\n' % (symbol, symbol))
