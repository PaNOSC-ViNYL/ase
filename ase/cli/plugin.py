from ase.cli.run import RunCommand


class PluginCommand(RunCommand):
    def __init__(self, calculate_function):
        self.calculate_function = calculate_function

    def calculate(self, atoms, name):
        data = self.calculate_function(atoms, name)
        if data is None:
            data = {}
        return data
