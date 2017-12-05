"Module for displaying information about the system."

from __future__ import unicode_literals

import numpy as np
from ase.gui.i18n import _
from ase.calculators.calculator import PropertyNotImplementedError

singleimage = _('Single image loaded.')
multiimage = _('Image %d loaded (0 - %d).')
ucconst = _('Unit cell is fixed.')
ucvaries = _('Unit cell varies.')

format = _("""\
%s

%s
Number of atoms: %d.

Unit cell:
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f

%s
%s
""")

calc_format = """

%s
"""

def info(gui):
    images = gui.images
    nimg = len(images)
    atoms = gui.atoms
    natoms = len(atoms)

    if len(atoms) < 1:
        txt = _('This frame has no atoms.')
    else:
        img = gui.frame

        uc = atoms.cell
        if nimg > 1:
            equal = True
            for i in range(nimg):
                equal = equal and (uc == images[i].cell).all()
            if equal:
                uctxt = ucconst
            else:
                uctxt = ucvaries
        else:
            uctxt = ''
        if nimg == 1:
            imgtxt = singleimage
        else:
            imgtxt = multiimage % (img, nimg - 1)

        periodic = [[_('no'), _('yes')][periodic]
                    for periodic in atoms.pbc]

        # TRANSLATORS: This has the form Periodic: no, no, yes
        formula = atoms.get_chemical_formula()
        pbcstring = _('Periodic: %s, %s, %s') % tuple(periodic)
        txt = format % ((imgtxt, formula, natoms) + tuple(uc.flat) +
                        (pbcstring,) + (uctxt,))
        if atoms.number_of_lattice_vectors == 3:
            txt += _('Volume: ') + '{:8.3f}'.format(atoms.get_volume())

        # Print electronic structure information if we have a calculator
        if atoms.calc:
            try:
                energy = _('Energy: ')
                energy += '{:.3f} eV'.format(atoms.get_potential_energy())

                maxf = np.linalg.norm(atoms.get_forces(apply_constraint=True),
                                      axis=1).max()
                forces = _('Max force: ')
                forces += '{:.3f}'.format(maxf)

                try:
                    magmom = atoms.get_magnetic_moments().sum()
                except PropertyNotImplementedError:
                    # We don't always have magnetic moments in our calculator
                    magmom = 0.

                mag = _('Magmom: ')
                mag += '{:.2f}'.format(magmom)
                txt += calc_format % '\n'.join([energy, forces, mag])
            except PropertyNotImplementedError:
                # We should always have energy and forces,
                # but just in case.
                pass

    return txt
