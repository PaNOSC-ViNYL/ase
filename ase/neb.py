import numpy as npy


class NEB:
    def __init__(self, images, k=0.1, climb=False):
        self.images = images
        self.k = k
        self.climb = climb
        self.natoms = len(images[0])
        self.nimages = len(images)

    def interpolate(self):
        pos1 = self.images[0].get_positions()
        pos2 = self.images[-1].get_positions()
        d = (pos2 - pos1) / (self.nimages - 1.0)
        for i in range(1, self.nimages - 1):
            self.images[i].set_positions(pos1 + i * d)

    def get_positions(self):
        positions = npy.empty((self.nimages * self.natoms, 3))
        n1 = 0
        for image in self.images:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        n1 = 0
        for image in self.images:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2
        
    def get_forces(self):
        positions = self.get_positions()
        positions.shape = (self.nimages, self.natoms, 3)
        forces = npy.empty((self.nimages * self.natoms, 3))

        forces[:self.natoms] = self.images[0].get_forces()
        imax = max([(self.images[i].get_potential_energy(), i)
                    for i in range(self.nimages)])[1]
        n1 = self.natoms
        tangent1 = positions[1] - positions[0]
        for i in range(1, self.nimages - 1):
            n2 = n1 + self.natoms
            tangent2 = positions[i + 1] - positions[i]
            if i < imax:
                tangent = tangent2
            elif i > imax:
                tangent = tangent1
            else:
                tangent = tangent1 + tangent2
                
            tt = npy.vdot(tangent, tangent)
            f = self.images[i].get_forces()
            ft = npy.vdot(f, tangent)
            if i == imax and self.climb:
                f -= 2 * ft / tt * tangent
            else:
                f -= ft / tt * tangent
                f += (npy.vdot(tangent1 - tangent2, tangent) * self.k / tt *
                      tangent)
                
            forces[n1:n2] = f
            n1 = n2
            tangent1 = tangent2

        forces[-self.natoms:] = self.images[-1].get_forces()

        return forces

    def get_potential_energy(self):
        return max([image.get_potential_energy() for image in self.images])

    def __len__(self):
        return self.nimages * self.natoms

    def writer(self, trajectory):
        return NEBTrajectoryWriter(self, trajectory).write


class NEBTrajectoryWriter:
    def __init__(self, neb, traj):
        self.neb = neb
        self.traj = traj

    def write(self):
        for image in self.neb.images:
            self.traj.write(image)

    
