"""This module contains functions to read from QBox output files"""

from ase import Atom, Atoms
from ase.calculators.singlepoint import SinglePointCalculator

import re
import xml.etree.ElementTree as ET


def read_qbox(file, index=-1):
    """Read data from QBox output file

    Inputs:
        file - str or fileobj, path to file or file object to read from
        index - int or slice, which frames to return
    Returns:
        list of Atoms or atoms, requested frame(s)
    """

    # Read in the output file
    tree = ET.parse(file)

    # Check whether this is a QB@all output
    is_qball = 'qb@LL' in tree.find("release").text

    # Load in atomic species
    species = dict()
    if is_qball:
        # XML for Python 3 reads comments, which means the "species" data might be split between
        #   multiple nodes on the root
        species_data = '\n'.join([x.tail for x in tree.getroot() if 'species' in x.tail])

        # Read out the species information with regular expressions
        symbols = re.findall('symbol_ = ([A-Z][a-z]?)', species_data)
        masses = re.findall('mass_ = ([0-9.]+)', species_data)
        names = re.findall('name_ = ([a-z]+)', species_data)
        numbers = re.findall('atomic_number_ = ([0-9]+)', species_data)

        # Compile them into a dictionary
        for name, symbol, mass, number in zip(names, symbols, masses, numbers):
            spec_data = dict(
                symbol=symbol,
                mass=float(mass),
                number=float(number)
            )
            species[name] = spec_data
    else:
        for spec in tree.findall('species'):
            name = spec.get('name')
            spec_data = dict(
                symbol=spec.find('symbol').text,
                mass=float(spec.find('mass').text),
                number=int(spec.find('atomic_number').text))
            species[name] = spec_data

    # Find all of the frames
    frames = tree.find("run").findall("iteration") if is_qball \
        else tree.findall("iteration")

    # If index is an int, return one frame
    if isinstance(index, int):
        return _parse_frame(frames[index], species)
    else:
        return [_parse_frame(frame, species) for frame in frames[index]]


def _parse_frame(tree, species):
    """Parse a certain frame from QBOX output

    Inputs:
        tree - ElementTree, <iteration> block from output file
        species - dict, data about species. Key is name of atom type,
            value is data about that type
    Return:
        Atoms object describing this iteration"""

    # Load in data about the system
    energy = float(tree.find("etotal").text)

    # Load in data about the cell
    unitcell = tree.find('atomset').find('unit_cell')
    cell = []
    for d in ['a', 'b', 'c']:
        cell.append([float(x) for x in unitcell.get(d).split()])

    stress_tree = tree.find('stress_tensor')
    if stress_tree is None:
        stresses = None
    else:
        stresses = [float(stress_tree.find('sigma_%s' % x).text)
                    for x in ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']]

    # Create the Atoms object
    atoms = Atoms(pbc=True, cell=cell)

    # Load in the atom information
    forces = []
    for atom in tree.find('atomset').findall('atom'):
        # Load data about the atom type
        spec = atom.get('species')
        symbol = species[spec]['symbol']
        mass = species[spec]['mass']

        # Get data about position / velocity / force
        pos = [float(x) for x in atom.find('position').text.split()]
        force = [float(x) for x in atom.find('force').text.split()]
        momentum = [float(x) * mass
                    for x in atom.find('velocity').text.split()]

        # Create the objects
        atom = Atom(symbol=symbol, mass=mass, position=pos, momentum=momentum)
        atoms += atom
        forces.append(force)

    # Create the calculator object that holds energy/forces
    calc = SinglePointCalculator(atoms,
                                 energy=energy, forces=forces, stress=stresses)
    atoms.set_calculator(calc)

    return atoms
