"""Tools easing the work with OpenBabel"""
import openbabel as ob
from ase import Atom, Atoms

def add_bonds(mol):
    """Automatically add bonds to molecule.

    This is done by writing to an xyz-file-like string and reading it again. In
    this case OpenBabel will automatically try to add bonds.
    """
    obconv = ob.OBConversion()
    obconv.SetOutFormat("xyz")
    xyz = obconv.WriteString(mol)

    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats("xyz", "mol")

    mol = ob.OBMol()
    obconv.ReadString(mol, xyz)

    return mol

def atoms_to_obmol(atoms, bonds=None):
    """Convert an Atoms object to an OBMol object.

    Parameters
    ==========
    atoms: Atoms
    bonds: list of lists of 3xint
        Define bonds between atoms such as:
            [[begin atom index, end atom index, bond order],
             ...
            ]
        If None the OpenBabel will try to construct the bonds
        automatically.
    """
    mol = ob.OBMol()
    for atom in atoms:
        a = mol.NewAtom()
        a.SetAtomicNum(int(atom.get_atomic_number()))
        a.SetVector(atom.position[0], atom.position[1], atom.position[2])

    if bonds is None:
        mol = add_bonds(mol)
    else:
        for bond in bonds:
            mol.AddBond(bond[0] + 1, bond[1] + 1, bond[2])

    return mol

def obmol_to_atoms(mol, return_bonds=False):
    """Convert an OBMol object to an Atoms object.

    Parameters
    ==========
    mol: OBMol
    return_bonds: bool
        If True, a list of list of 3xint describing the bonds will be returned.
    """
    atoms = Atoms()
    for i in range(mol.NumAtoms()):
        obatom = mol.GetAtom(i + 1)
        atoms.append(Atom(obatom.GetAtomicNum(),
                          [obatom.GetX(),
                           obatom.GetY(),
                           obatom.GetZ()]
                         )
                    )

    if return_bonds:
        bonds = []
        for i in range(mol.NumBonds()):
            obbond = mol.GetBond(i)
            bond = [obbond.GetBeginAtomIdx() - 1,
                    obbond.GetEndAtomIdx() - 1,
                    obbond.GetBondOrder()]
            bonds.append(bond)
        return atoms, bonds
    else:
        return atoms
