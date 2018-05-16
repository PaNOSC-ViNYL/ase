import socket

import numpy as np


class IPIProtocol:
    statements = {'POSDATA', 'GETFORCE', 'STATUS', 'INIT', ''}
    responses = {'READY', 'HAVEDATA', 'FORCEREADY', 'NEEDINIT'}

    def __init__(self, socket):
        self.socket = socket

    def log(self, *txt):
        print('IPI:', *txt)

    def sendmsg(self, msg):
        self.log('  sendmsg', repr(msg))
        assert msg in self.statements, msg
        # The statement '' means end program.
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
        self.log('  send {} bytes of {}'.format(len(buf), dtype))
        self.socket.sendall(buf)

    def recv(self, shape, dtype):
        a = np.empty(shape, dtype)
        nbytes = np.dtype(dtype).itemsize * np.prod(shape)
        buf = self.socket.recv(nbytes)
        self.recv('  recv {} bytes of {}'.format(len(buf), dtype))
        a.flat[:] = np.frombuffer(buf, dtype=dtype)
        assert np.isfinite(a).all()
        return a

    def sendposdata(self, cell, icell, positions):
        assert cell.size == 9
        assert icell.size == 9
        assert positions.size % 3 == 0

        self.log(' sendposdata')
        self.sendmsg('POSDATA')
        self.send(cell, np.float64)
        self.send(icell, np.float64)
        self.send(len(positions), np.int32)
        self.send(positions, np.float64)

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
        return e, forces, virial, morebytes

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
        msg = self.status()
        assert msg == 'READY', msg
        return r


def ipi(process_args):
    port = 27182
    serversocket = socket.socket(socket.AF_INET)
    serversocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    serversocket.bind(('localhost', port))
    serversocket.listen(1)

    proc = Popen(process_args)
    clientsocket, address = serversocket.accept()
    ipi = IPIProtocol(clientsocket)


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
