"Module for displaying information about the system."

from __future__ import unicode_literals

import numpy as np
from ase.gui.i18n import _
import warnings

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
        try:
            if atoms.calc:
                calc = atoms.calc

                energy = None
                forces = None
                magmoms = None

                if not calc.calculation_required(calc.atoms, ['energy']):
                    energy = atoms.get_potential_energy()

                if not calc.calculation_required(calc.atoms, ['forces']):
                    forces = atoms.get_forces()

                if not calc.calculation_required(calc.atoms, ['magmoms']):
                    magmoms = atoms.get_magnetic_moments()

                calc_strs = []
                if energy is not None:
                    energy_str = _('Energy: ')
                    energy_str += '{:.3f} eV'.format(energy)
                    calc_strs.append(energy_str)

                if forces is not None:
                    maxf = np.linalg.norm(forces, axis=1).max()
                    forces_str = _('Max force: ')
                    forces_str += '{:.3f}'.format(maxf)
                    calc_strs.append(forces_str)

                if magmoms is not None:
                    magmom = magmoms.sum()
                    mag_str = _('Magmom: ')
                    mag_str += '{:.3f}'.format(magmom)
                    calc_strs.append(mag_str)

                # Format into string
                txt += calc_format % '\n'.join(calc_strs)
        except Exception as err:
            # We don't want to kill the GUI because something
            # went wrong with a calculator object.
            msg = ('An error occured while retrieving results '
                   'from the calculator')
            warnings.warn(msg)
            print(err)

    return txt
