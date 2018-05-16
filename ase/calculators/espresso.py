"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""


from ase import io
from ase.calculators.calculator import FileIOCalculator


class Espresso(FileIOCalculator):
    """
    """
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms']
    command = 'pw.x -in PREFIX.pwi > PREFIX.pwo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='espresso', atoms=None, **kwargs):
        """
        All options for pw.x are copied verbatim to the input file, and put
        into the correct section. Use ``input_data`` for parameters that are
        already in a dict, all other ``kwargs`` are passed as parameters.

        Accepts all the options for pw.x as given in the QE docs, plus some
        additional options:

        input_data: dict
            A flat or nested dictionary with input parameters for pw.x
        pseudopotentials: dict
            A filename for each atomic species, e.g.
            ``{'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}``.
            A dummy name will be used if none are given.
        kspacing: float
            Generate a grid of k-points with this as the minimum distance,
            in A^-1 between them in reciprocal space. If set to None, kpts
            will be used instead.
        kpts:
            Number of kpoints in each dimension for automatic kpoint generation.
        koffset: (int, int, int)
            Offset of kpoints in each direction. Must be 0 (no offset) or
            1 (half grid offset). Setting to True is equivalent to (1, 1, 1).


        .. note::
           Set ``tprnfor=True`` and ``tstress=True`` to calculate forces and
           stresses.


        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        io.write(self.label + '.pwi', atoms, **self.parameters)

    def read_results(self):
        output = io.read(self.label + '.pwo')
        self.results = output.calc.results


class IPIServer:
    port = 27182
    def __init__(self, process_args):
        import socket
        from ase.calculators.ipi import IPIProtocol
        from subprocess import Popen
        self._closed = False
        self.serversocket = socket.socket(socket.AF_INET)
        self.serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.serversocket.bind(('localhost', self.port))
        self.serversocket.listen(1)

        self.proc = Popen(process_args, shell=True)
        self.clientsocket, self.address = self.serversocket.accept()
        self.ipi = IPIProtocol(self.clientsocket)

    def close(self):
        import socket
        # Proper way to close sockets?
        self.ipi.end()
        self.clientsocket.shutdown(socket.SHUT_RDWR)
        self.serversocket.shutdown(socket.SHUT_RDWR)
        self._closed = True
        self.proc.wait()

    def calculate(self, atoms):
        return self.ipi.calculate(atoms.positions, atoms.cell)

    def __del__(self):
        print(self.proc)
        print(self.proc.wait())
    #def __del__(self):
    #    if not self._closed:
    #        self.close()


from ase.calculators.calculator import all_changes
class IPIEspresso(Espresso):
    ipi_changes = {'positions', 'cell'}

    def __init__(self, *args, **kwargs):
        kwargs.update()
        Espresso.__init__(self, *args, **kwargs)
        self.ipi = None

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        print(system_changes)
        if any(change not in self.ipi_changes for change in system_changes):
            print('IPI=NONE')
            self.ipi = None

        if self.ipi is None:
            #self.write_input(atoms, properties=properties,
            #                 system_changes=system_changes)
            #cmd = self.command.replace('PREFIX', self.prefix)
            #print('NEW IPI SERVER')
            self.ipi = IPIServer('pw.x < pw.in --ipi localhost:27182')

        results = self.ipi.calculate(atoms)
        self.results.update(results)
