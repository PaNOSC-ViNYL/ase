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

        def getresult(name, get_quantity):
            # ase/io/trajectory.py line 170 does this by using
            # the get_property(prop, atoms, allow_calculation=False)
            # so that is an alternative option.
            try:
                if calc.calculation_required(calc.atoms, [name]):
                    quantity = None
                else:
                    quantity = get_quantity()
            except Exception as err:
                quantity = None
                errmsg = ('An error occured while retrieving {} '
                          'from the calculator: {}'.format(name, err))
                warnings.warn(errmsg)
            return quantity

        # Print electronic structure information if we have a calculator
        if atoms.calc:
            calc = atoms.calc
            calc_strs = []
            # SinglePointCalculators are named after the code which
            # produced the result, so this will typically list the
            # name of a code even if they are just cached results.
            from ase.calculators.singlepoint import SinglePointCalculator
            if isinstance(calc, SinglePointCalculator):
                calc_strs.append(_('Calculator: {} (cached)')
                                 .format(calc.name))
            else:
                calc_strs.append(_('Calculator: {} (attached)'))

            energy = getresult('energy', atoms.get_potential_energy)
            forces = getresult('forces', atoms.get_forces)
            magmoms = getresult('magmoms', atoms.get_magnetic_moments)

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
            if calc_strs:
                txt += calc_format % '\n'.join(calc_strs)

    return txt
