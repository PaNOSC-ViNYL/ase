# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""Definition of the Spacegroup class.

This module only depends on NumPy and the space group database.
"""

import os
import warnings

import numpy as np

__all__ = ['Spacegroup']


class SpacegroupError(Exception):
    """Base exception for the spacegroup module."""
    pass

class SpacegroupNotFoundError(SpacegroupError):
    """Raised when given space group cannot be found in data base."""
    pass

class SpacegroupValueError(SpacegroupError):
    """Raised when arguments have invalid value."""
    pass


class Spacegroup(object):
    """A space group class. 
    
    The instances of Spacegroup describes the symmetry operations for
    the given space group. 

    Example:
    
    >>> from spacegroup import Spacegroup
    >>> 
    >>> sg = Spacegroup(225)
    >>> print 'Space group', sg.no, sg.symbol
    Space group 225 F m -3 m
    >>> sg.scaled_primitive_cell
    array([[ 0. ,  0.5,  0.5],
           [ 0.5,  0. ,  0.5],
           [ 0.5,  0.5,  0. ]])
    >>> sites, kinds = sg.equivalent_sites([[0,0,0]])
    >>> sites
    array([[ 0. ,  0. ,  0. ],
           [ 0. ,  0.5,  0.5],
           [ 0.5,  0. ,  0.5],
           [ 0.5,  0.5,  0. ]])
    """
    no = property(
        lambda self: self._no,
        doc='Space group number in International Tables of Crystallography.')
    symbol = property(
        lambda self: self._symbol, 
        doc='Hermann-Mauguin (or international) symbol for the space group.')
    setting = property(
        lambda self: self._setting, 
        doc='Space group setting. Either one or two.')
    lattice = property(
        lambda self: self._symbol[0], 
        doc="""Lattice type: 
            P      primitive
            I      body centering, h+k+l=2n
            F      face centering, h,k,l all odd or even
            A,B,C  single face centering, k+l=2n, h+l=2n, h+k=2n
            R      rhombohedral centering, -h+k+l=3n (obverse); h-k+l=3n (reverse)
            """)
    centrosymmetric = property(
        lambda self: self._centrosymmetric, 
        doc='Whether a center of symmetry exists.')
    scaled_primitive_cell = property(
        lambda self: self._scaled_primitive_cell, 
        doc='Primitive cell in scaled coordinates as a matrix with the '
        'primitive vectors along the rows.')
    reciprocal_cell = property(
        lambda self: self._reciprocal_cell, 
        doc='Tree Miller indices that span all kinematically non-forbidden '
        'reflections as a matrix with the Miller indices along the rows.')
    nsubtrans = property(
        lambda self: len(self._subtrans), 
        doc='Number of cell-subtranslation vectors.')
    def _get_nsymop(self):
        if self.centrosymmetric:
            return 2 * len(self._rotations) * len(self._subtrans)
        else:
            return len(self._rotations) * len(self._subtrans)
    nsymop = property(_get_nsymop, doc='Total number of symmetry operations.')
    subtrans = property(
        lambda self: self._subtrans, 
        doc='Translations vectors belonging to cell-sub-translations.')
    rotations = property(
        lambda self: self._rotations, 
        doc='Symmetry rotation matrices. The invertions are not included '
        'for centrosymmetrical crystals.')
    translations = property(
        lambda self: self._translations, 
        doc='Symmetry translations. The invertions are not included '
        'for centrosymmetrical crystals.')

    def __init__(self, spacegroup, setting=1, datafile=None):
        """Returns a new Spacegroup instance. 

        Parameters:

        spacegroup : int | string | Spacegroup instance
            The space group number in International Tables of
            Crystallography or its Hermann-Mauguin symbol. E.g. 
            spacegroup=225 and spacegroup='F m -3 m' are equivalent.
        setting : 1 | 2
            Some space groups have more than one setting. `setting`
            determines Which of these should be used.
        datafile : None | string
            Path to database file. If `None`, the the default database 
            will be used.
        """
        if isinstance(spacegroup, Spacegroup):
            for k, v in spacegroup.__dict__.iteritems():
                setattr(self, k, v)
            return 
        if not datafile:
            datafile = get_datafile()
        f = open(datafile, 'r')
        try:
            _read_datafile(self, spacegroup, setting, f)
        finally:
            f.close()

    def __repr__(self):
        return 'Spacegroup(%d, setting=%d)'%(self.no, self.setting)
    
    def __str__(self):
        """Return a string representation of the space group data in
        the same format as found the database."""
        retval = []
        # no, symbol
        retval.append('%-3d   %s\n' % (self.no, self.symbol))
        # setting
        retval.append('  setting %d\n' % (self.setting))
        # centrosymmetric
        retval.append('  centrosymmetric %d\n' % (self.centrosymmetric))
        # primitive vectors
        retval.append('  primitive vectors\n')
        for i in range(3):
            retval.append('   ')
            for j in range(3):
                retval.append(' %13.10f' % (self.scaled_primitive_cell[i, j]))
            retval.append('\n')
        # primitive reciprocal vectors
        retval.append('  reciprocal vectors\n')
        for i in range(3):
            retval.append('   ')
            for j in range(3):
                retval.append(' %3d' % (self.reciprocal_cell[i, j]))
            retval.append('\n')
        # sublattice
        retval.append('  %d subtranslations\n' % self.nsubtrans)
        for i in range(self.nsubtrans):
            retval.append('   ')
            for j in range(3):
                retval.append(' %13.10f' % (self.subtrans[i, j]))
            retval.append('\n')
        # symmetry operations
        nrot = len(self.rotations)
        retval.append('  %d symmetry operations (rot+trans)\n' % nrot)
        for i in range(nrot):
            retval.append(' ')
            for j in range(3):
                retval.append(' ')
                for k in range(3):
                    retval.append(' %2d' % (self.rotations[i, j, k]))
                retval.append('  ')
            for j in range(3):
                retval.append(' %13.10f' % self.translations[i, j])
            retval.append('\n')
        retval.append('\n')
        return ''.join(retval)

    def __eq__(self, other):
        """Chech whether *self* and *other* refer to the same
        spacegroup number and setting."""
        if isinstance(other, int):
            return self.no == other
        elif isinstance(other, str):
            return self.symbol == other
        else:
            return self.no == other.no and self.setting == other.setting

    def __index__(self):
        return self.no

    def get_symop(self):
        """Returns all symmetry operations (including inversions and
        subtranslations) as a sequence of (rotation, translation)
        tuples."""
        symop = []
        parities = [1]
        if self.centrosymmetric:
            parities.append(-1)
        for parity in parities:
            for subtrans in self.subtrans:
                for rot,trans in zip(self.rotations, self.translations):
                    newtrans = np.mod(trans + subtrans, 1)
                    symop.append((parity*rot, newtrans))
        return symop

    def get_rotations(self):
        """Return all rotations, including inversions for
        centrosymmetric crystals."""
        if self.centrosymmetric:
            return np.vstack((self.rotations, -self.rotations))
        else:
            return self.rotations

    def equivalent_reflections(self, hkl):
        """Return all equivalent reflections to the list of Miller indices
        in hkl.

        Example:

        >>> from spacegroup import Spacegroup
        >>> sg = Spacegroup(225)  # fcc
        >>> sg.equivalent_reflections([[0,0,2]])
        array([[ 0,  0, -2],
               [ 0, -2,  0],
               [-2,  0,  0],
               [ 2,  0,  0],
               [ 0,  2,  0],
               [ 0,  0,  2]])
        """
        hkl = np.array(hkl, dtype='int')
        rot = self.get_rotations()
        n,nrot = len(hkl), len(rot)
        R = rot.transpose(0, 2, 1).reshape((3*nrot, 3)).T
        refl = np.dot(hkl, R).reshape((n*nrot, 3))
        ind = np.lexsort(refl.T)
        refl = refl[ind]
        diff = np.diff(refl, axis=0)
        mask = np.any(diff, axis=1)
        return np.vstack((refl[mask], refl[-1,:]))
            
    def equivalent_sites(self, scaled_positions, ondublicates='error', 
                         symprec=1e-3):
        """Returns the scaled positions and all their equivalent sites.

        Parameters:

        scaled_positions: list | array
            List of non-equivalent sites given in unit cell coordinates.
        ondublicates : 'keep' | 'replace' | 'warn' | 'error'
            Action if `scaled_positions` contain symmetry-equivalent
            positions:
            
            'keep'
               ignore additional symmetry-equivalent positions
            'replace'
                replace
            'warn'
                like 'keep', but issue an UserWarning
            'error'
                raises a SpacegroupValueError
                    
        symprec: float
            Minimum "distance" betweed two sites in scaled coordinates
            before they are counted as the same site.

        Returns:

        sites: array
            A NumPy array of equivalent sites.
        kinds: list
            A list of integer indices specifying which input site is 
            equivalent to the corresponding returned site.
        """
        kinds = []
        sites = []
        symprec2 = symprec**2
        for kind, pos in enumerate(scaled_positions):
            for rot, trans in self.get_symop():
                site = np.mod(np.dot(rot, pos) + trans, 1.)
                if not sites:
                    sites.append(site)
                    kinds.append(kind)
                    continue
                t = site - sites
                mask = np.sum(t*t, 1) < symprec2
                if np.any(mask):
                    ind = np.argwhere(mask)[0][0]
                    if kinds[ind] == kind:
                        pass
                    elif ondublicates == 'keep':
                        pass
                    elif ondublicates == 'replace':
                        kinds[ind] = kind
                    elif ondublicates == 'warn':
                        warnings.warn('scaled_positions %d and %d '
                                      'are equivalent'%(kinds[ind], kind))
                    elif ondublicates == 'error':
                        raise SpacegroupValueError(
                            'scaled_positions %d and %d are equivalent'%(
                                kinds[ind], kind))
                    else:
                        raise SpacegroupValueError(
                            'Argument "ondublicates" must be one of: '
                            '"keep", "replace", "warn" or "error".')
                else:
                    sites.append(site)
                    kinds.append(kind)
        return np.array(sites), kinds




def get_datafile():
    """Return default path to datafile."""
    return os.path.join(os.path.dirname(__file__), 'spacegroup.dat')




#-----------------------------------------------------------------
# Functions for parsing the database. They are moved outside the
# Spacegroup class in order to make it easier to later implement
# caching to avoid reading the database each time a new Spacegroup
# instance is created.
#-----------------------------------------------------------------

def _skip_to_blank(f, spacegroup, setting):
    """Read lines from f until a blank line is encountered."""
    while True:
        line = f.readline()
        if not line:
            raise SpacegroupNotFoundError(
                'invalid spacegroup %s, setting %i not found in data base' % 
                ( spacegroup, setting ) )
        if not line.strip():
            break


def _skip_to_nonblank(f, spacegroup, setting):
    """Read lines from f until a nonblank line not starting with a
    hash (#) is encountered and returns this and the next line."""
    while True:
        line1 = f.readline()
        if not line1:
            raise SpacegroupNotFoundError(
                'invalid spacegroup %s, setting %i not found in data base' %
                ( spacegroup, setting ) )
        line1.strip()
        if line1 and not line1.startswith('#'):
            line2 = f.readline()
            break
    return line1, line2


def _read_datafile_entry(spg, no, symbol, setting, f):
    """Read space group data from f to spg."""
    spg._no = no
    spg._symbol = symbol.strip()
    spg._setting = setting
    spg._centrosymmetric = bool(int(f.readline().split()[1]))
    # primitive vectors
    f.readline()
    spg._scaled_primitive_cell = np.array([map(float, f.readline().split()) 
                                       for i in range(3)], 
                                      dtype=np.float)
    # primitive reciprocal vectors
    f.readline()
    spg._reciprocal_cell = np.array([map(int, f.readline().split()) 
                                        for i in range(3)],
                                       dtype=np.int)
    # subtranslations
    spg._nsubtrans = int(f.readline().split()[0])
    spg._subtrans = np.array([map(float, f.readline().split()) 
                              for i in range(spg._nsubtrans)],
                             dtype=np.float)
    # symmetry operations
    nsym = int(f.readline().split()[0])
    symop = np.array([map(float, f.readline().split()) for i in range(nsym)],
                     dtype=np.float)
    spg._nsymop = nsym
    spg._rotations = np.array(symop[:,:9].reshape((nsym,3,3)), dtype=np.int)
    spg._translations = symop[:,9:]


def _read_datafile(spg, spacegroup, setting, f):
    if isinstance(spacegroup, int):
        pass
    elif isinstance(spacegroup, basestring):
        spacegroup = ' '.join(spacegroup.strip().split())
    else:
        raise SpacegroupValueError('`spacegroup` must be of type int or str')
    while True:
        line1, line2 = _skip_to_nonblank(f, spacegroup, setting)
        _no,_symbol = line1.strip().split(None, 1)
        _setting = int(line2.strip().split()[1])
        _no = int(_no)
        if ((isinstance(spacegroup, int) and _no == spacegroup) or
            (isinstance(spacegroup, basestring) and 
             _symbol == spacegroup)) and _setting == setting:
            _read_datafile_entry(spg, _no, _symbol, _setting, f)
            break
        else:
            _skip_to_blank(f, spacegroup, setting)
                


#-----------------------------------------------------------------
# Self test
if __name__ == '__main__':

    # Import spacegroup in order to ensure that __file__ is defined
    # such that the data base can be found.
    import spacegroup

    import doctest
    print 'doctest: ', doctest.testmod()
