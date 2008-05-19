import numpy as npy

from ase.parallel import world, rank

class NEB:
    def __init__(self, images, k=0.1, climb=False, parallel=False):
        self.images = images
        self.k = k
        self.climb = climb
        self.parallel = parallel
        self.natoms = len(images[0])
        self.nimages = len(images)

    def interpolate(self):
        pos1 = self.images[0].get_positions()
        pos2 = self.images[-1].get_positions()
        d = (pos2 - pos1) / (self.nimages - 1.0)
        for i in range(1, self.nimages - 1):
            self.images[i].set_positions(pos1 + i * d)

    def get_positions(self):
        positions = npy.empty(((self.nimages - 2) * self.natoms, 3))
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2
        
    def get_forces(self):
        images = self.images

        forces = npy.empty(((self.nimages - 2), self.natoms, 3))

        energies = npy.empty(self.nimages - 2)

        if not self.parallel:
            # Do all images - one at a time:
            for i in range(1, self.nimages - 1):
                energies[i - 1] = images[i].get_potential_energy()
                forces[i - 1] = images[i].get_forces()
        else:
            # Parallelize over images:
            i = rank // (self.nimages - 2) + 1
            energies[i - 1] = images[i].get_potential_energy()
            forces[i - 1] = images[i].get_forces()
            for i in range(1, self.nimages - 1):
                root = (i - 1) * size // (self.nimages - 2)
                world.broadcast(energies[i:i + 1], root)
                world.broadcast(forces[i], root)

        imax = 1 + npy.argsort(energies)[-1]

        tangent1 = images[1].get_positions() - images[0].get_positions()
        for i in range(1, self.nimages - 1):
            tangent2 = images[i + 1].get_positions() - images[i].get_positions()
            if i < imax:
                tangent = tangent2
            elif i > imax:
                tangent = tangent1
            else:
                tangent = tangent1 + tangent2
                
            tt = npy.vdot(tangent, tangent)
            f = forces[i - 1]
            ft = npy.vdot(f, tangent)
            if i == imax and self.climb:
                f -= 2 * ft / tt * tangent
            else:
                f -= ft / tt * tangent
                f -= (npy.vdot(tangent1 - tangent2, tangent) *
                      self.k / tt * tangent)
                
            tangent1 = tangent2

        return forces.reshape((-1, 3))

    def get_potential_energy(self):
        return npy.nan
        """
        if not self.parallel:
            return max([image.get_potential_energy() 
                        for image in self.images[1:-1]])
        else:
            i = rank // (self.nimages - 2) + 1
            forces[i - 1] = images[i].get_forces()
            for i in range(1, self.nimages - 1):
                world.broadcast(forces[i], (i - 1) * size // (self.nimages - 2))
        """
    def __len__(self):
        return (self.nimages - 2) * self.natoms

    def writer(self, trajectory):
        return NEBTrajectoryWriter(self, trajectory).write

    def write(self, trajectory):
        NEBTrajectoryWriter(self, trajectory).write()


class NEBTrajectoryWriter:
    def __init__(self, neb, traj):
        self.neb = neb
        self.traj = traj

    def write(self):
        assert not self.neb.parallel
        for image in self.neb.images:
            self.traj.write(image)

    
