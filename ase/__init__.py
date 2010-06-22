# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Atomic Simulation Environment."""


from ase.atom import Atom
from ase.atoms import Atoms

_deprecate_things_from_ase_module = True

# Some day in the future, we will uncomment this line:
#__all__ = ['Atoms', 'Atom']  

from ase.units import *
from ase.io import read, write
from ase.io.trajectory import PickleTrajectory
from ase.dft import STM, monkhorst_pack, DOS
from ase.optimize.mdmin import MDMin
from ase.optimize.lbfgs import HessLBFGS
from ase.optimize.fire import FIRE
from ase.optimize.lbfgs import LBFGS, LBFGSLineSearch
from ase.optimize.bfgs import BFGS
from ase.optimize import QuasiNewton
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.constraints import *
from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT
from ase.calculators.siesta import Siesta
from ase.calculators.dacapo import Dacapo
from ase.calculators.vasp import Vasp
from ase.calculators.aims import Aims, AimsCube
from ase.calculators.turbomole import Turbomole
from ase.calculators.dftb import Dftb
from ase.neb import NEB, SingleCalculatorNEB
from ase.visualize import view
from ase.data import chemical_symbols, atomic_numbers, atomic_names, \
     atomic_masses, covalent_radii, reference_states
from ase.data.molecules import molecule

from math import sqrt, pi
import numpy as np



if _deprecate_things_from_ase_module:
    import types
    import ase.utils.deprecate as dep

    _locals = locals()

    for name, obj in _locals.items():
        if name.startswith('_') or name in ['Atoms', 'Atom']:
            continue

        if isinstance(obj, float):
            if name == 'pi':
                pi = dep.DeprecatedFloat(pi, 'pi', 'math')
            else:
                _locals[name] = dep.DeprecatedFloat(obj, name, 'ase.units')
        elif isinstance(obj, types.ModuleType):
            pass  # how about np?  XXX
        elif hasattr(obj, '__module__'):
            module = obj.__module__
            if module.startswith('ase.optimize'):
                module = 'ase.optimize'
            elif module.startswith('ase.md'):
                module = 'ase.md'
            elif name == 'PickleTrajectory':
                module = 'ase.io'
            _locals[name] = dep.Deprecate(obj, name, module)
        else:
            pass  # how about atomic_numbers, covalent_radii, ... ? XXX

    np = dep.DeprecatedNumpyImport()

