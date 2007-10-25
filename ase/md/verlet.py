import numpy as npy


class VelocityVerlet:
    def __init__(self, atoms):
        self.atoms = atoms
        self.callbacks = []

    def attach(self, callback):
        self.callbacks.append(callback)

    def run(self, dt=0.01, steps=10000000000000):
        self.step = 0
        for x in self.iter():
            for callback in self.callbacks:
                callback()
            self.step += 1
            if self.step == steps:
                break
            
    def iter(self):
        atoms = self.atoms
        b = self.dt * (atoms.masses**-1)[:, npy.newaxis]
        atoms.set_momenta(atoms.momenta() + 0.5 * self.dt * atoms.forces())
        while True:
            atoms.set_positions(atoms.positions() + self.dt * atoms.momenta())
            yield None
            
            atoms.set_momenta(atoms.momenta() + self.dt * atoms.forces())
