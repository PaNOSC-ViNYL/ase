import os
import socket
from subprocess import check_call

import numpy as np


class IPIProtocol:
    requests = {'POSDATA', 'GETFORCE', 'STATUS'}  # 'INIT' not implemented
    responses = {'READY', 'FORCEREADY'}  # 'NEEDINIT' not implemented

    def __init__(self, socket):
        self.socket = socket

    def sendheader(self, msg):
        """Send header"""
        assert msg in self.requests, msg
        msg = msg.encode('ascii').ljust(12)
        self.socket.sendall(msg)

    def recvheader(self):
        msg = self.socket.recv(12)
        msg = msg.rstrip().decode('ascii')
        assert msg in self.responses, msg
        return msg

    def send(self, a, dtype):
        b = np.asarray(a, dtype).tobytes()
        self.socket.sendall(b)

    def recv(self, shape, dtype):
        a = np.empty(shape, dtype)
        nbytes = np.dtype(dtype).itemsize * np.prod(shape)
        buf = self.socket.recv(nbytes)
        a.flat[:] = np.frombuffer(buf, dtype=dtype)
        return a

    def sendposdata(self, cell, icell, positions):
        assert cell.size == 9
        assert icell.size == 9
        assert positions.size % 3 == 0

        self.sendheader('POSDATA')
        self.send(cell, np.float64)
        self.send(icell, np.float64)
        self.send(len(positions), np.int32)
        self.send(positions, np.float64)

    def sendrecv_force(self):
        self.sendheader('GETFORCE')
        msg = self.recvheader()
        assert msg == 'FORCEREADY', msg
        e = self.recv(1, np.float64)[0]
        natoms = self.recv(1, np.int32)
        forces = self.recv((int(natoms), 3), np.float64)
        virial = self.recv((3, 3), np.float64)
        nmorebytes = self.recv(1, np.int32)
        morebytes = self.recv(int(nmorebytes), np.byte)
        return e, forces, virial, morebytes

    def status(self):
        self.sendheader('STATUS')
        msg = self.recvheader()
        return msg


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
    proc = subprocess.Popen(['cp2k', 'in.cp2k'])

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

if __name__ == '__main__':
    main()
