
class Calculator:
    def __init__(self):
        return

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def get_name(self):
        """Return the name of the calculator (string).  """
        return self.name

    def get_version(self):
        """Return the version of the calculator (string).  """
        raise NotImplementedError

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        if force_consistent:
            return self.energy_free
        else:
            return self.energy_zero

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        if self.stress is not None:
            return self.stress
        else:
            raise NotImplementedError

    def initialize(self, atoms):
        """Prepare the input files required to
        start the program (calculator).  """
        raise NotImplementedError
