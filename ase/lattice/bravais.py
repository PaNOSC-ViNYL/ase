"""Bravais.py - class for generating Bravais lattices etc.

This is a base class for numerous classes setting up pieces of crystal.
"""

import math
import numpy as np
from ase.atoms import Atoms

class Bravais:
    """Bravais lattice factory.

    This is a base class for the objects producing various lattices
    (SC, FCC, ...).
    """

    # The following methods are NOT defined here, but must be defined
    # in classes inhering from Bravais:
    #   get_lattice_constant
    #   make_crystal_basis
    # The following class attributes are NOT defined here, but must be defined
    # in classes inhering from Bravais:
    #   int_basis
    #   inverse_basis

    other = {0:(1,2), 1:(2,0), 2:(0,1)}

    # For Bravais lattices with a basis, set the basis here.  Leave as
    # None if no basis is present.
    bravais_basis = None

    # If more than one type of element appear in the crystal, give the
    # order here.  For example, if two elements appear in a 3:1 ratio,
    # bravais_basis could contain four vectors, and element_basis
    # could be (0,0,1,0) - the third atom in the basis is different
    # from the other three.  Leave as None if all atoms are of the
    # same type.
    element_basis = None
    
    def __call__(self, directions=(None,None,None), miller=(None,None,None),
                 size=(1,1,1), element=None, latticeconstant=None,
                 periodic=True, align=True, debug=0):
        "Create a lattice."
        self.size = size
        self.periodic = periodic
        self.debug = debug
        self.process_element(element)
        self.find_directions(directions, miller)
        if self.debug:
            self.print_directions_and_miller()
        self.convert_to_natural_basis()
        if self.debug:
            self.print_directions_and_miller(" (natural basis)")
        if latticeconstant is None:
            if self.element_basis is None:
                self.latticeconstant = self.get_lattice_constant()
            else:
                raise ValueError,\
                      "A lattice constant must be specified for a compound"
        else:
            self.latticeconstant = latticeconstant
        if self.debug:
            print "Expected number of atoms in unit cell:", self.calc_num_atoms()
        if self.debug >= 2:
            print "Bravais lattice basis:", self.bravais_basis
            if self.bravais_basis is not None:
                print " ... in natural basis:", self.natural_bravais_basis
        self.make_crystal_basis()
        self.make_unit_cell()
        if align:
            self.align()
        return self.make_list_of_atoms()

    def align(self):
        "Align the first axis along x-axis and the second in the x-y plane."
        degree = 180/pi
        if self.debug >= 2:
            print "Basis before alignment:"
            print self.basis
        if self.basis[0][0]**2 + self.basis[0][2]**2 < 0.01 * self.basis[0][1]**2:
            # First basis vector along y axis - rotate 90 deg along z
            t = np.array([[0, -1, 0],
                          [1, 0, 0],
                          [0, 0, 1]], Float)
            self.basis = matrixmultiply(self.basis, t)
            transf = t
            if self.debug >= 2:
                print "Rotating -90 degrees around z axis for numerical stability."
                print self.basis
        else:
            transf = identity(3, Float)
        assert abs(determinant(transf) - 1) < 1e-6
        # Rotate first basis vector into xy plane
        theta = math.atan2(self.basis[0,2], self.basis[0,0])
        t = np.array([[cos(theta), 0, -sin(theta)],
                      [         0, 1, 0          ],
                      [sin(theta), 0, cos(theta) ]])
        self.basis = matrixmultiply(self.basis, t)
        transf = matrixmultiply(transf, t)
        if self.debug >= 2:
            print "Rotating %f degrees around y axis." % (-theta*degree,)
            print self.basis
        assert abs(determinant(transf) - 1) < 1e-6
        # Rotate first basis vector to point along x axis
        theta = math.atan2(self.basis[0,1], self.basis[0,0])
        t = np.array([[cos(theta), -sin(theta), 0],
                      [sin(theta),  cos(theta), 0],
                      [         0,           0, 1]])
        self.basis = matrixmultiply(self.basis, t)
        transf = matrixmultiply(transf, t)
        if self.debug >= 2:
            print "Rotating %f degrees around z axis." % (-theta*degree,)
            print self.basis
        assert abs(determinant(transf) - 1) < 1e-6
        # Rotate second basis vector into xy plane
        theta = math.atan2(self.basis[1,2], self.basis[1,1])
        t = np.array([[1, 0, 0],
                      [0, cos(theta), -sin(theta)],
                      [0, sin(theta),  cos(theta)]])
        self.basis = matrixmultiply(self.basis, t)
        transf = matrixmultiply(transf, t)
        if self.debug >= 2:
            print "Rotating %f degrees around x axis." % (-theta*degree,)
            print self.basis
        assert abs(determinant(transf) - 1) < 1e-6
        # Now we better rotate the atoms as well
        self.atoms = matrixmultiply(self.atoms, transf)
        # ... and rotate miller_basis
        self.miller_basis = matrixmultiply(self.miller_basis, transf)
        
    def make_list_of_atoms(self):
        "Repeat the unit cell."
        nrep = self.size[0] * self.size[1] * self.size[2]
        if nrep <= 0:
            raise ValueError, "Cannot create a non-positive number of unit cells"
        # Now the unit cells must be merged.
        a2 = []
        e2 = []
        for i in xrange(self.size[0]):
            offset = self.basis[0] * i
            a2.append(self.atoms + offset[NewAxis,:])
            e2.append(self.elements)
        atoms = concatenate(a2)
        elements = concatenate(e2)
        a2 = []
        e2 = []
        for j in xrange(self.size[1]):
            offset = self.basis[1] * j
            a2.append(atoms + offset[NewAxis,:])
            e2.append(elements)
        atoms = concatenate(a2)
        elements = concatenate(e2)
        a2 = []
        e2 = []
        for k in xrange(self.size[2]):
            offset = self.basis[2] * k
            a2.append(atoms + offset[NewAxis,:])
            e2.append(elements)
        atoms = concatenate(a2)
        elements = concatenate(e2)
        del a2, e2
        assert len(atoms) == nrep * len(self.atoms)
        basis = np.array([[self.size[0],0,0],
                       [0,self.size[1],0],
                       [0,0,self.size[2]]])
        basis = matrixmultiply(basis, self.basis)
        # None should be replaced, and memory should be freed.
        lattice = Lattice(positions=atoms, cell=basis, numbers=elements,
                       pbc=self.periodic)
        lattice.millerbasis = self.miller_basis
        return lattice
        
    def process_element(self, element):
        "Extract atomic number from element"
        # The types that can be converted to Elements: integers and strings
        if self.element_basis is None:
            if isinstance(element, type("string")):
                self.atomicnumber = data.atomic_numbers[element]
            else:
                # Let's hope it is an atomic number or something compatible.
                self.atomicnumber = element
        else:
            self.atomicnumber = []
            for e in element:
                if isinstance(e, type("string")):
                    self.atomicnumber.append(data.atomic_numbers[e])
                else:
                    # Let's hope it is an atomic number or something compatible.
                    self.atomicnumber.append(e)
            assert len(self.atomicnumber) == len(self.bravais_basis)
        
    def convert_to_natural_basis(self):
        "Convert directions and miller indices to the natural basis."
        self.directions = matrixmultiply(self.directions, self.inverse_basis)
        if self.bravais_basis is not None:
            self.natural_bravais_basis = matrixmultiply(self.bravais_basis,
                                                        self.inverse_basis)
        for i in (0,1,2):
            self.directions[i] = reduceindex(self.directions[i])
        for i in (0,1,2):
            (j,k) = self.other[i]
            self.miller[i] = reduceindex(self.handedness *
                                         cross(self.directions[j],
                                                self.directions[k]))
        
    def calc_num_atoms(self):
        v = int(round(abs(determinant(self.directions))))
        if self.bravais_basis is None:
            return v
        else:
            return v * len(self.bravais_basis)
    
    def make_unit_cell(self):
        "Make the unit cell."
        # Make three loops, and find the positions in the integral
        # lattice.  Each time a position is found, the atom is placed
        # in the real unit cell by put_atom().
        self.natoms = self.calc_num_atoms()
        self.nput = 0
        self.atoms = np.zeros((self.natoms,3), Float)
        self.elements = np.zeros(self.natoms, Int)
        self.farpoint = farpoint = sum(self.directions)
        #printprogress = self.debug and (len(self.atoms) > 250)
        percent = 0
        # Find the radius of the sphere containing the whole system
        sqrad = 0
        for i in (0,1):
            for j in (0,1):
                for k in (0,1):
                    vect = (i * self.directions[0] +
                            j * self.directions[1] +
                            k * self.directions[2])
                    if dot(vect,vect) > sqrad:
                        sqrad = dot(vect,vect)
        del i,j,k
        # Loop along first crystal axis (i)
        for (istart, istep) in ((0,1), (-1,-1)):
            i = istart
            icont = True
            while icont:
                nj = 0
                for (jstart, jstep) in ((0,1), (-1,-1)):
                    j = jstart
                    jcont = True
                    while jcont:
                        nk = 0
                        for (kstart, kstep) in ((0,1), (-1,-1)):
                            k = kstart
                            #print "Starting line i=%d, j=%d, k=%d, step=(%d,%d,%d)" % (i,j,k,istep,jstep,kstep) 
                            kcont = True
                            while kcont:
                                # Now (i,j,k) loops over Z^3, except that
                                # the loops can be cut off when we get outside
                                # the unit cell.
                                point = np.array((i,j,k))
                                if self.inside(point):
                                    self.put_atom(point)
                                    nk += 1
                                    nj += 1
                                # Is k too high?
                                if dot(point,point) > sqrad:
                                    assert not self.inside(point)
                                    kcont = False
                                k += kstep
                        # Is j too high?
                        if i*i+j*j > sqrad:
                            jcont = False
                        j += jstep
                # Is i too high?
                if i*i > sqrad:
                    icont = False
                i += istep
                #if printprogress:
                #    perce = int(100*self.nput / len(self.atoms))
                #    if perce > percent + 10:
                #        print ("%d%%" % perce),
                #        percent = perce
        assert(self.nput == self.natoms)

    def inside(self, point):
        "Is a point inside the unit cell?"
        return (dot(self.miller[0], point) >= 0 and
                dot(self.miller[0], point - self.farpoint) < 0 and
                dot(self.miller[1], point) >= 0 and
                dot(self.miller[1], point - self.farpoint) < 0 and
                dot(self.miller[2], point) >= 0 and
                dot(self.miller[2], point - self.farpoint) < 0)

    def put_atom(self, point):
        "Place an atom given its integer coordinates."
        if self.bravais_basis is None:
            # No basis - just place a single atom
            pos = matrixmultiply(point, self.crystal_basis)
            if self.debug >= 2:
                print ("Placing an atom at (%d,%d,%d) ~ (%.3f, %.3f, %.3f)."
                       % (tuple(point) + tuple(pos)))
            self.atoms[self.nput] = pos
            self.elements[self.nput] = self.atomicnumber
            self.nput += 1
        else:
            for i, offset in enumerate(self.natural_bravais_basis):
                pos = matrixmultiply(point + offset, self.crystal_basis)
                if self.debug >= 2:
                    print ("Placing an atom at (%d+%f, %d+%f, %d+%f) ~ (%.3f, %.3f, %.3f)."
                           % (point[0], offset[0], point[1], offset[1],
                              point[2], offset[2], pos[0], pos[1], pos[2]))
                self.atoms[self.nput] = pos
                if self.element_basis is None:
                    self.elements[self.nput] = self.atomicnumber
                else:
                    self.elements[self.nput] = self.atomicnumber[i]
                self.nput += 1
                
    def find_directions(self, directions, miller):
        "Find missing directions and miller indices from the specified ones."
        directions = list(directions)
        miller = list(miller)
        # Now fill in missing directions and miller indices.  This is an
        # iterative process.
        change = 1
        while change:
            change = False
            missing = 0
            for i in (0,1,2):
                (j,k) = self.other[i]
                if directions[i] is None:
                    missing += 1
                    if miller[j] is not None and miller[k] is not None:
                        directions[i] = reduceindex(cross(miller[j],
                                                          miller[k]))
                        change = True
                        if self.debug >= 2:
                            print "Calculating directions[%d] from miller indices" % i
                if miller[i] is None:
                    missing += 1
                    if directions[j] is not None and directions[k] is not None:
                        miller[i] = reduceindex(cross(directions[j],
                                                      directions[k]))
                        change = True
                        if self.debug >= 2:
                            print "Calculating miller[%d] from directions" % i
        if missing:
            raise ValueError, "Specification of directions and miller indices is incomplete."
        # Make sure that everything is Numeric arrays
        self.directions = np.array(directions)
        self.miller = np.array(miller)
        # Check for left-handed coordinate system
        if determinant(self.directions) < 0:
            print "WARNING: Creating a left-handed coordinate system!"
            self.miller = -self.miller
            self.handedness = -1
        else:
            self.handedness = 1
        # Now check for consistency
        for i in (0,1,2):
            (j,k) = self.other[i]
            m = reduceindex(self.handedness *
                            cross(self.directions[j], self.directions[k]))
            if sum(not_equal(m, self.miller[i])):
                print "ERROR: Miller index %s is inconsisten with directions %d and %d" % (i,j,k)
                print "Miller indices:"
                print str(self.miller)
                print "Directions:"
                print str(self.directions)
                raise ValueError, "Inconsistent specification of miller indices and directions."

    def print_directions_and_miller(self, txt=""):
        "Print direction vectors and Miller indices."
        print "Direction vectors of unit cell%s:" % (txt,)
        for i in (0,1,2):
            print "   ", self.directions[i]
        print "Miller indices of surfaces%s:" % (txt,)
        for i in (0,1,2):
            print "   ", self.miller[i]

class MillerInfo:
    """Mixin class to provide information about Miller indices."""
    def miller_to_direction(self, miller):
	"""Returns the direction corresponding to a given Miller index."""
        return matrixmultiply(miller, self.millerbasis)

class Lattice(Atoms, MillerInfo):
    """List of atoms initially containing a regular lattice of atoms.

    A part from the usual list of atoms methods this list of atoms type
    also has a method, `miller_to_direction`, used to convert from Miller
    indices to directions in the coordinate system of the lattice.
    """
    pass

# Helper functions
def cross(a, b):
    """The cross product of two vectors."""
    return np.array((a[1]*b[2] - b[1]*a[2],
                     a[2]*b[0] - b[2]*a[0],
                     a[0]*b[1] - b[0]*a[1]))

def gcd(a,b):
    """Greatest Common Divisor of a and b."""
    while a != 0:
        a,b = b%a,a
    return b

def reduceindex(M):
    "Reduce Miller index to the lowest equivalent integers."
    oldM = M
    g = gcd(M[0], M[1])
    h = gcd(g, M[2])
    while h != 1:
        M = M/h
        g = gcd(M[0], M[1])
        h = gcd(g, M[2])
    if dot(oldM, M) > 0:
        return M
    else:
        return -M

