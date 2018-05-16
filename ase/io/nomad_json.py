import json

#try:
#    from urllib2 import urlopen
#except ImportError:
#    from urllib.request import urlopen

import numpy as np

import ase.units as units
from ase import Atoms
from ase.io import read
from ase.data import chemical_symbols


def section_system2atoms(section):
    numbers = section['atom_species']
    numbers = np.array(numbers, int)
    numbers[numbers < 0] = 0
    numbers[numbers > len(chemical_symbols)] = 0
    positions = section['atom_positions']['flatData']
    positions = np.array(positions).reshape(-1, 3) * units.m
    pbc = section.get('configuration_periodic_dimensions')
    cell = section.get('lattice_vectors')
    atoms = Atoms(numbers, positions=positions)
    atoms.info['nomad_uri'] = section['uri']
    if pbc is not None:
        assert len(pbc) == 1
        pbc = pbc[0]  # it's a list??
        pbc = pbc['flatData']
        assert len(pbc) == 3
        atoms.pbc = pbc

    # celldisp?
    if cell is not None:
        cell = cell['flatData']
        cell = np.array(cell).reshape(3, 3) * units.m
        atoms.cell = cell

    return atoms


def section_singleconfig2calc(section):
    from ase.calculators.singlepoint import SinglePointCalculator
    kwargs = {}
    # Forces, total energy, ........
    # We should be able to extract e.g. a band structure as well.
    if 'energy_free' in section:
        kwargs['free_energy'] = section['energy_free'] * units.J
    calc = SinglePointCalculator(**kwargs)
    return calc


def dict2images(d):
    assert 'section_run' in d, 'Missing section_run'
    runs = d['section_run']
    images = []
    for run in runs:
        systems = run['section_system']
        for system in systems:
            atoms = section_system2atoms(system)
            images.append(atoms)
        #configs = run['section_single_configuration_calculation']
        #for config in configs:
        #    calc = section_singleconfig2calc(config)
        #    sysref = config['single_configuration_to_calculation_system_ref']
        #    systems[sysref
        #    #print(type(config))
        #skskjdf
    return images


def read_nomad_json(fd, index):
    # wth, we should not be passing index like this!
    from ase.io.formats import string2index
    i = string2index(index)
    print(i)
    d = json.load(fd)
    images = dict2images(d)
    return images[i]
