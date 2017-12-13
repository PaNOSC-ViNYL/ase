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
            calc = atoms.calc
            if hasattr(calc, 'get_property'):
                # We now assume we're a type of Calculator class
                try:
                    # Get property should take care of non-existing results,
                    # and return None
                    energy = calc.get_property('energy', atoms=calc.atoms,
                                               allow_calculation=False)
                    forces = calc.get_property('forces', atoms=calc.atoms,
                                               allow_calculation=False)
                    magmoms = calc.get_property('magmoms', atoms=calc.atoms,
                                                allow_calculation=False)
                except PropertyNotImplementedError:
                    # This shouldn't happen with get_property(),
                    # but just in case.
                    # Perhaps we could also catch other errors, and simply pass
                    # here.
                    pass
                else:
                    # Parse the default None values
                    if energy is None:
                        energy = 0.

                    if forces is None:
                        maxf = 0.
                    else:
                        maxf = np.linalg.norm(forces, axis=1).max()

                    if magmoms is None:
                        magmom = 0.
                    else:
                        magmom = magmoms.sum()

                    # Format into string
                    energy_str = _('Energy: ')
                    energy_str += '{:.3f} eV'.format(energy)

                    forces_str = _('Max force: ')
                    forces_str += '{:.3f}'.format(maxf)

                    mag_str = _('Magmom: ')
                    mag_str += '{:.3f}'.format(magmom)
                    txt += calc_format % '\n'.join([
                        energy_str, forces_str, mag_str])

    return txt
