import socket
from subprocess import Popen

import numpy as np

from ase.calculators.calculator import (Calculator, all_changes,
                                        PropertyNotImplementedError)
import ase.units as units


class IPIProtocol:
    """Communication using IPI protocol."""

    statements = {'POSDATA', 'GETFORCE', 'STATUS', 'INIT', ''}
    # The statement '' means end program.
    responses = {'READY', 'HAVEDATA', 'FORCEREADY', 'NEEDINIT'}

    def __init__(self, socket, txt=None):
        self.socket = socket

        if txt is None:
            log = lambda *args: None
        else:
            def log(*args):
                print('IPI:', *args, file=txt)
        self.log = log

    def sendmsg(self, msg):
        self.log('  sendmsg', repr(msg))
        assert msg in self.statements, msg
        msg = msg.encode('ascii').ljust(12)
        self.socket.sendall(msg)

    def recvmsg(self):
        msg = self.socket.recv(12)
        msg = msg.rstrip().decode('ascii')
        assert msg in self.responses, msg
        self.log('  recvmsg', repr(msg))
        return msg

    def send(self, a, dtype):
        buf = np.asarray(a, dtype).tobytes()
        #self.log('  send {}'.format(np.array(a).ravel().tolist()))
        self.log('  send {} bytes of {}'.format(len(buf), dtype))
        self.socket.sendall(buf)

    def recv(self, shape, dtype):
        a = np.empty(shape, dtype)
        nbytes = np.dtype(dtype).itemsize * np.prod(shape)
        buf = self.socket.recv(nbytes)
        self.log('  recv {} bytes of {}'.format(len(buf), dtype))
        a.flat[:] = np.frombuffer(buf, dtype=dtype)
        #self.log('  recv {}'.format(a.ravel().tolist()))
        assert np.isfinite(a).all()
        return a

    def sendposdata(self, cell, icell, positions):
        assert cell.size == 9
        assert icell.size == 9
        assert positions.size % 3 == 0

        self.log(' sendposdata')
        self.sendmsg('POSDATA')
        self.send(cell / units.Bohr, np.float64)
        self.send(icell / units.Bohr, np.float64)
        self.send(len(positions), np.int32)
        self.send(positions / units.Bohr, np.float64)

    def sendrecv_force(self):
        self.log(' sendrecv_force')
        self.sendmsg('GETFORCE')
        msg = self.recvmsg()
        assert msg == 'FORCEREADY', msg
        e = self.recv(1, np.float64)[0]
        natoms = self.recv(1, np.int32)
        assert natoms >= 0
        forces = self.recv((int(natoms), 3), np.float64)
        virial = self.recv((3, 3), np.float64)
        nmorebytes = self.recv(1, np.int32)
        assert nmorebytes >= 0
        morebytes = self.recv(int(nmorebytes), np.byte)
        return (e * units.Ha, units.Ha / units.Bohr * forces,
                units.Ha / units.Bohr**3 * virial, morebytes)

    def status(self):
        self.log(' status')
        self.sendmsg('STATUS')
        msg = self.recvmsg()
        return msg

    def end(self):
        self.log( 'end')
        self.sendmsg('')

    def sendinit(self):
        self.log(' sendinit')
        self.sendmsg('INIT')
        self.send(0, np.int32)  # 'bead index' always zero
        # number of bits (don't they mean bytes?) in initialization string:
        self.send(-1, np.int32)
        self.send(np.empty(0), np.byte)  # initialization string

    def calculate(self, positions, cell):
        self.log('calculate')
        msg = self.status()
        if msg == 'NEEDINIT':
            self.sendinit()
            msg = self.status()
        assert msg == 'READY', msg
        icell = np.linalg.pinv(cell).transpose()
        self.sendposdata(cell, icell, positions)
        msg = self.status()
        assert msg == 'HAVEDATA', msg
        e, forces, stress, morebytes = self.sendrecv_force()
        r = dict(energy=e,
                 forces=forces,
                 stress=stress)
        if morebytes:
            r['morebytes'] = morebytes
        return r


class IPIServer:
    def __init__(self, process_args, host='localhost', port=31415):
        self.host = host
        self.port = port
        self._closed = False
        self.serversocket = socket.socket(socket.AF_INET)
        self.serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.serversocket.bind((host, self.port))
        self.serversocket.listen(1)
        print('listening for ipi connection on {}:{}'.format(host, port))

        print('run', process_args)
        self.proc = Popen(process_args, shell=True)
        print('accept')
        self.clientsocket, self.address = self.serversocket.accept()
        print('accept ok')
        self.ipi = IPIProtocol(self.clientsocket)

    def close(self):
        if self._closed:
            return

        self._closed = True

        # Proper way to close sockets?
        if hasattr(self, 'ipi') and self.ipi is not None:
            self.ipi.end()
            self.ipi = None
        if hasattr(self, 'proc'):
            exitcode = self.proc.wait()
            print('subproc exitcode {}'.format(exitcode))
        if hasattr(self, 'clientsocket'):
            self.clientsocket.shutdown(socket.SHUT_RDWR)
        if hasattr(self, 'serversocket'):
            self.serversocket.shutdown(socket.SHUT_RDWR)

    def calculate(self, atoms):
        return self.ipi.calculate(atoms.positions, atoms.cell)


class IPICalculator(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
    ipi_supported_changes = {'positions', 'cell'}

    def __init__(self, calc, host='localhost', port=31415):
        Calculator.__init__(self)
        self.calc = calc
        self.host = host
        self.port = port
        self.ipi = None

    def todict(self):
        return {'type': 'calculator',
                'name': 'ipi',
                'calculator': self.calc.todict()}

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        print('CALCULATE', system_changes)
        bad = [change for change in system_changes
               if change not in self.ipi_supported_changes]

        if self.ipi is not None and any(bad):
            raise PropertyNotImplementedError(
                'Cannot change {} through IPI protocol.  '
                'Please create new IPI calculator.'
                .format(bad if len(bad) > 1 else bad[0]))

        if self.ipi is None:
            cmd = self.calc.command
            cmd = cmd.format(host=self.host, port=self.port,
                             prefix=self.calc.prefix)
            self.calc.write_input(atoms, properties=properties,
                                  system_changes=system_changes)
            self.ipi = IPIServer(cmd, port=self.port, host=self.host)

        self.atoms = atoms.copy()
        results = self.ipi.calculate(atoms)
        self.results.update(results)

    def close(self):
        if self.ipi is not None:
            self.ipi.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()



def main():
    from ase.build import molecule
    atoms = molecule('H2O')
    atoms.center(vacuum=2.0)
    atoms.pbc = 1

    port = 27182
    serversocket = socket.socket(socket.AF_INET)
    serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    serversocket.bind(('localhost', port))
    print('listen 1')
    serversocket.listen(1)

    import subprocess

    print('popen')
    subprocess.Popen(['cp2k', 'in.cp2k'])

    print('accept')
    clientsocket, address = serversocket.accept()
    ipi = IPIProtocol(clientsocket)

    msg = ipi.status()
    assert msg == 'READY', msg
    ipi.sendposdata(atoms.cell, atoms.get_reciprocal_cell(),
                    atoms.positions)
    msg = ipi.status()
    assert msg == 'HAVEDATA', msg
    e, forces, stress, morebytes = ipi.sendrecv_force()
    print('energy')
    print(e)
    print('forces')
    print(forces)
    print('stress')
    print(stress)

    msg = ipi.status()
    assert msg == 'READY', msg
    ipi.sendmsg('')

if __name__ == '__main__':
    main()
