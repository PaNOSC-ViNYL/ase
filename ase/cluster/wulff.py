import numpy as np

import ase.cluster
from ase.cluster.cubic import FaceCenteredCubic
from ase.data import atomic_numbers, reference_states

delta = 1e-10
_debug = None

def wulff_construction(symbol, energies, size, rounding="closest",
                       symmetry=None, latticeconstant=None, debug=0):
    """Create a cluster using the Wulff construction.

    A cluster is created with approximately the number of atoms
    specified, following the Wulff construction, i.e. minimizing the
    surface energy of the cluster.

    Parameters:
    -----------

    symbol: The chemical symbol (or atomic number) of the desired element.

    energies: A list of surface energies for the surfaces.  See below.

    size: The desired number of atoms.

    rounding (optional): Specifies what should be done if no Wulff
    construction corresponds to exactly the requested number of atoms.
    Should be a string, either "above", "below" or "closest" (the
    default), meaning that the nearest cluste above or below - or the
    closest one - is created instead.

    symmetry (optional): The crystal structure.  If not specified, the
    reference crystal structure from ase.data is used. CURRENTLY, ONLY
    FCC IS SUPPORTED.

    latticeconstant (optional): The lattice constant.  If not given,
    extracted from ase.data.

    Note on surface energies.
    -------------------------

    For FCC crystals, 26 surface energies must be specified, 6 for the
    {100} surfaces, 12 for the {110} surfaces and 8 for the {111}
    surfaces.  The order is (1,0,0), (-1,0,0), (0,1,0), (0,-1,0),
    (0,0,1), (0,0,-1), (1,1,0), (-1,-1,0), (1,0,1), (-1,0,-1),
    (0,1,1), (0,-1,-1), (1,-1,0), (-1,1,0), (1,0,-1), (-1,0,1),
    (0,1,-1), (0,-1,1), (1,1,1), (-1,-1,-1), (-1,1,1), (1,-1,-1),
    (1,-1,1), (-1,1,-1), (1,1,-1), (-1,-1,1).

    """

    global _debug
    _debug = debug

    if rounding not in ["above", "below", "closest"]:
        raise ValueError("Invalid rounding: "+rounding)
    
    # Interpret symbol
    if isinstance(symbol, str):
        atomic_number = atomic_numbers[symbol]
    else:
        atomic_number = symbol

    # Interpret symmetry
    if symmetry is None:
        symmetry = reference_states[self.atomic_number]['symmetry'].lower()
    else:
        symmetry = symmetry.lower()
    if symmetry not in ["fcc"]:
        raise NotImplementedError("Crystal structure "+symmetry+
                                  " is not supported.")

    # Get info about layers and opposing layers.
    surface_names = ase.cluster.data.lattice[symmetry]['surface_names']
    surface_mapping = ase.cluster.data.lattice[symmetry]['surface_mapping']
    atoms_in_unitcell = 4   # fcc
    nsurf = len(surface_names)
    
    if len(energies) != len(surface_names):
        raise ValueError("The energies array should contain %d values."
                         % (len(surface_names),))
    # First guess a size that is not too large.
    wanted_size = (size / atoms_in_unitcell)**(1.0/3.0)
    energies_sum = np.array([energies[i] + energies[surface_mapping[i]]
                             for i in range(nsurf)])
    if energies_sum.min() <= 0.0:
        raise ValueError("Surface energies are impossible (negative sum).")
    max_e = energies_sum.max()
    factor = wanted_size / max_e
    #layers = np.array([int(round(factor * e)) for e in energies])
    atoms, layers = make_atoms(symbol, energies, factor, latticeconstant,
                               symmetry)
    if len(atoms) == 0:
        # Probably the cluster is very flat
        if debug:
            print "First try made an empty cluster, trying again."
        factor = 1 / energies_sum.min()
        atoms, layers = make_atoms(symbol, energies, factor,
                                   latticeconstant, symmetry)
        if len(atoms) == 0:
            raise RuntimeError("Failed to create a finite cluster.")

    # Second guess: scale to get closer.
    old_factor = factor
    old_layers = layers
    old_atoms = atoms
    factor *= (size / len(atoms))**(1.0/3.0)
    atoms, layers = make_atoms(symbol, energies, factor, latticeconstant,
                               symmetry)
    if len(atoms) == 0:
        print "Second guess gave an empty cluster, discarding it."
        atoms = old_atoms
        factor = old_factor
        layers = old_layers
    else:
        del old_atoms

    # Find if the cluster is too small or too large (both means perfect!)
    below = above = None
    if len(atoms) <= size:
        below = atoms
    if len(atoms) >= size:
        above = atoms

    # Now iterate towards the right cluster
    iter = 0
    while (below is None or above is None):
        if len(atoms) < size:
            # Find a larger cluster
            if debug:
                print "Making a larger cluster."
            factor = ((layers + 0.5 + delta) / energies).min()
            atoms, new_layers = make_atoms(symbol, energies, factor,
                                           latticeconstant, symmetry)
            assert (new_layers - layers).max() == 1
            assert (new_layers - layers).min() >= 0
            layers = new_layers
        else:
            # Find a smaller cluster
            if debug:
                print "Making a smaller cluster."
            factor = ((layers - 0.5 - delta) / energies).min()
            atoms, new_layers = make_atoms(symbol, energies, factor,
                                           latticeconstant, symmetry)
            assert (new_layers - layers).max() <= 0
            assert (new_layers - layers).min() == -1
            layers = new_layers
        if len(atoms) <= size:
            below = atoms
        if len(atoms) >= size:
            above = atoms
        iter += 1
        if iter == 100:
            raise RuntimeError("Runaway iteration.")
    if rounding == "below":
        return below
    elif rounding == "above":
        return above
    else:
        assert rounding == "closest"
        if (len(above) - size) < (size - len(below)):
            return above
        else:
            return below


def make_atoms(symbol, energies, factor, latticeconstant, symmetry):
    layers1 = factor * np.array(energies)
    layers = np.round(layers1).astype(int)
    #layers = np.array([int(round(factor * e)) for e in energies])
    fcc_surface_names = ase.cluster.data.lattice['fcc']['surface_names']
    atoms = FaceCenteredCubic(symbol, fcc_surface_names, layers=layers,
                              latticeconstant=latticeconstant)
    if _debug:
        print "Created a cluster with %i atoms: %s" % (len(atoms),
                                                       str(layers))
        #print layers1
    return (atoms, layers)
