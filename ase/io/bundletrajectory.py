"""bundletrajectory - a module for I/O from large MD simulations.

The BundleTrajectory class writes trajectory into a directory with the
following structure::

    filename.bundle (dir)
        metadata.pickle        Data about the file format, and about which
                               data is present.
        state.pickle           The number of frames
        F0 (dir)               Frame number 0
            small.pickle       Small data structures in a dictionary
                               (pbc, cell, ...)
            numbers.pickle     Atomic numbers
            positions.pickle   Positions
            momenta.pickle     Momenta
            ...
        F1 (dir)
"""

import ase.parallel 
from ase.parallel import paropen
import pupynere # Non-relative import ase.io.pupynere creates import loop!
import numpy as np
import os
import shutil
import time
import cPickle as pickle
import weakref
from copy import copy, deepcopy


class BundleTrajectory:
    """Reads and writes atoms into a .bundle directory.

    The BundleTrajectory is an alternative way of storing
    trajectories, intended for large-scale molecular dynamics
    simulations, where a single flat file becomes unwieldy.  Instead,
    the data is stored in directory, a 'bundle' (the name bundle is
    inspired from bundles in Mac OS, which are really just directories
    the user is supposed to think of as a single file-like unit).

    Parameters:

    filename:
        The name of the directory.  Preferably ending in .bundle.

    mode (optional):
        The file opening mode.  'r' means open for reading, 'w' for
        writing and 'a' for appending.  Default: 'r'.  If opening in
        write mode, and the filename already exists, the old file is
        renamed to .bak (any old .bak file is deleted), except if the
        existing file is empty.

    atoms (optional):
        The atoms that will be written.  Can only be specified in
        write or append mode.  If not specified, the atoms must be
        given as an argument to the .write() method instead.

    backup=True:
        Use backup=False to disable renaming of an existing file.
    """
    slavelog = True  # Log from all nodes
    def __init__(self, filename, mode='r', atoms=None, backup=True):
        self.state = 'constructing'
        self.filename = filename
        self.pre_observers = []   # Callback functions before write is performed
        self.post_observers = []  # Callback functions after write is performed
        self.master = ase.parallel.rank == 0
        self._set_defaults()
        self._set_backend()
        if mode == 'r':
            if atoms is not None:
                raise ValueError("You cannot specify atoms in read mode.")
            self._open_read()
        elif mode == 'w':
            self._open_write(atoms, backup)
        elif mode == 'a':
            self._open_append(atoms)

    def _set_defaults(self):
        "Set default values for internal parameters."
        self.version = 1
        self.subtype = 'normal'
        self.backend_name = 'pickle'
        self.datatypes = {'positions': True,
                          'numbers': 'once',
                          'tags': 'once',
                          'masses': 'once',
                          'momenta': True,
                          'forces': True,
                          'energy': True,
                          'energies': False,
                          'stress': False,
                          'magmoms': True
                          }

    def _set_backend(self, backend=None):
        """Set the backed doing the actual I/O.

        It should either be a backend object, or a string containing
        the name of the desired backend.
        """
        if backend is None:
            backend = self.backend_name
        if isinstance(backend, str):
            # A backend specified by name
            self.backend_name = backend
            if self.backend_name == 'pickle':
                self.backend = PickleBundleBackend(self.master)
            elif self.backend_name == 'netcdf':
                self.backend = NetCDFBundleBackend(self.master)
            else:
                raise NotImplementedError(
                    "This version of ASE cannot use BundleTrajectory with backend '%s'"
                    % self.backend_name)
        else:
            # A backend object
            self.backend = backend

    def write(self, atoms=None):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        # Check that we are in write mode
        if self.state == 'prewrite':
            self.state = 'write'
            assert self.nframes == 0
        elif self.state != 'write':
            raise RuntimeError('Cannot write in ' + self.state + ' mode.')

        if atoms is None:
            atoms = self.atoms

        if hasattr(atoms, 'interpolate'):
            # seems to be a NEB
            self.log('Beginning to write NEB data')
            neb = atoms
            assert not neb.parallel
            try:
                neb.get_energies_and_forces(all=True)
            except AttributeError:
                pass
            for image in neb.images:
                self.write(image)
            self.log('Done writing NEB data')
            return

        # OK, it is a real atoms object.  Write it.
        self._call_observers(self.pre_observers)
        self.log("Beginning to write frame "+str(self.nframes))
        framedir = os.path.join(self.filename, "F"+str(self.nframes))
        self._make_framedir(framedir)

        # Check which data should be written the first time:
        # Modify datatypes so any element of type 'once' becomes true
        # for the first frame but false for subsequent frames.
        datatypes = {}
        for k, v in self.datatypes.items():
            if v == 'once':
                v = (self.nframes == 0)
            datatypes[k] = v

        # Write 'small' data structures.  They are written jointly.
        smalldata = {'pbc': atoms.get_pbc(),
                     'cell': atoms.get_cell(),
                     'natoms': atoms.get_number_of_atoms(),
                     'constraints': atoms.constraints,
                     }
        if datatypes.get('energy'):
            try:
                smalldata['energy'] = atoms.get_potential_energy()
            except (RuntimeError, NotImplementedError):
                self.datatypes['energy'] = False
        if datatypes.get('stress'):
            try:
                smalldata['stress'] = atoms.get_stress()
            except NotImplementedError:
                self.datatypes['stress'] = False
        self.backend.write_small(framedir, smalldata)
        
        # Write the large arrays.
        if datatypes.get('positions'):
            self.backend.write(framedir, 'positions', atoms.get_positions())
        if datatypes.get('numbers'):
            self.backend.write(framedir, 'numbers', atoms.get_atomic_numbers())
        if datatypes.get('tags'):
            if atoms.has('tags'):
                self.backend.write(framedir, 'tags', atoms.get_tags())
            else:
                self.datatypes['tags'] = False
        if datatypes.get('masses'):
            if atoms.has('masses'):
                self.backend.write(framedir, 'masses', atoms.get_masses())
            else:
                self.datatypes['masses'] = False
        if datatypes.get('momenta'):
            if atoms.has('momenta'):
                self.backend.write(framedir, 'momenta', atoms.get_momenta())
            else:
                self.datatypes['momenta'] = False
        if datatypes.get('magmoms'):
            if atoms.has('magmoms'):
                self.backend.write(framedir, 'magmoms', atoms.get_magmoms())
            else:
                self.datatypes['magmoms'] = False
        if datatypes.get('forces'):
            try:
                x = atoms.get_forces()
            except (RuntimeError, NotImplementedError):
                self.datatypes['forces'] = False
            else:
                self.backend.write(framedir, 'forces', x) 
                del x
        if datatypes.get('energies'):
            try:
                x = atoms.get_potential_energies()
            except (RuntimeError, NotImplementedError):
                self.datatypes['energies'] = False
            else:
                self.backend.write(framedir, 'energies', x) 
                del x
        # Finally, write metadata if it is the first frame
        if self.nframes == 0:
            metadata = {'datatypes': self.datatypes}
            self._write_metadata(metadata)
        self._write_nframes(self.nframes + 1)
        self._call_observers(self.post_observers)
        self.log("Done writing frame "+str(self.nframes))
        self.nframes += 1

    def select_data(self, data, value):
        """Selects if a given data type should be written.

        Data can be written in every frame (specify True), in the
        first frame only (specify 'only') or not at all (specify
        False).  Not all data types support the 'only' keyword, if not
        supported it is interpreted as True.

        The following data types are supported, the letter in parenthesis
        indicates the default:

        positions (T), numbers (O), tags (O), masses (O), momenta (T),
        forces (T), energy (T), energies (F), stress (F), magmoms (T)

        If a given property is not present during the first write, it
        will be not be saved at all.
        """
        if value not in (True, False, 'once'):
            raise ValueError("Unknown write mode")
        if data not in self.datatypes:
            raise ValueError("Unsupported data type: " + data)
        self.datatypes[data] = value
        
    def close(self):
        "Closes the trajectory."
        self.state = 'closed'
        lf = getattr(self, 'logfile', None)
        if lf is not None:
            lf.close()
            del self.logfile
            
    def log(self, text):
        """Write to the log file in the bundle.

        Logging is only possible in write/append mode.

        This function is mainly for internal use, but can also be called by
        the user.
        """
        if not (self.master or self.slavelog):
            return
        text = time.asctime() + ': ' + text
        if hasattr(self, "logfile"):
            # Logging enabled
            if self.logfile is None:
                # Logfile not yet open
                try:
                    self.logdata.append(text)
                except AttributeError:
                    self.logdata = [text]
            else:
                self.logfile.write(text+'\n')
                self.logfile.flush()
        else:
            raise RuntimeError("Cannot write to log file in mode "+self.state)

    # __getitem__ is the main reading method.
    def __getitem__(self, n):
        "Read an atoms object from the BundleTrajectory."
        if self.state != 'read':
            raise IOError('Cannot read in %s mode' % (self.state,))
        if n < 0:
            n += self.nframes
        if n < 0 or n >= self.nframes:
            raise IndexError('Trajectory index %d out of range [0, %d['
                             % (n, self.nframes))
        framedir = os.path.join(self.filename, 'F'+str(n))
        framezero = os.path.join(self.filename, 'F0')
        smalldata = self.backend.read_small(framedir)
        data = {}
        data['pbc'] = smalldata['pbc']
        data['cell'] = smalldata['cell']
        data['constraint'] = smalldata['constraints']
        atoms = ase.Atoms(**data)
        natoms = smalldata['natoms']
        for name in ('positions', 'numbers', 'tags', 'masses',
                     'momenta', 'magmoms'):
            if self.datatypes.get(name):
                if self.datatypes[name] == 'once':
                    d = self.backend.read(framezero, name)
                else:
                    d = self.backend.read(framedir, name)
                assert len(d) == natoms
                atoms.arrays[name] = d
        # Create the atoms object
        if self.datatypes.get('energy'):
            if self.datatype.get('forces'):
                forces = self.backend.read(framedir, 'forces')
            else:
                forces = None
            calc = SinglePointCalculator(smalldata.get('energy'),
                                         forces,
                                         smalldata.get('stress'),
                                         magmoms, atoms)
            atoms.set_calculator(calc)
        return atoms

    def __len__(self):
        return self.nframes
    
    def _open_log(self):
        if not (self.master or self.slavelog):
            return
        if self.master:
            lfn = os.path.join(self.filename, "log.txt")
        else:
            lfn = os.path.join(self.filename, ("log-node%d.txt" %
                                               (ase.parallel.rank,)))
        self.logfile = open(lfn, "a")   # Append to log if it exists.
        if hasattr(self, 'logdata'):
            for text in self.logdata:
                self.logfile.write(text+'\n')
            self.logfile.flush()
            del self.logdata

    def _open_write(self, atoms, backup):
        "Open a bundle trajectory for writing."
        self.logfile = None # Enable delayed logging
        self.atoms = atoms
        if os.path.exists(self.filename):
            # The output directory already exists.
            if not self.is_bundle(self.filename):
                raise IOError("Filename '" + self.filename +
                              "' already exists, but is not a BundleTrajectory." +
                              "Cowardly refusing to remove it.")
            ase.parallel.barrier() # All must have time to see it exists.
            if self.is_empty_bundle(self.filename):
                ase.parallel.barrier()  # All must see it is empty
                self.log('Deleting old "%s" as it is empty' % (self.filename,)) 
                self.delete_bundle(self.filename)
            elif not backup:
                self.log('Deleting old "%s" as backup is turned off.' % (self.filename,))
                self.delete_bundle(self.filename)
            else:
                # Make a backup file
                bakname = self.filename + '.bak'
                if os.path.exists(bakname):
                    ase.parallel.barrier()  # All must see it exists
                    self.log('Deleting old backup file "%s"' % (bakname,))
                    self.delete_bundle(bakname)
                self.log('Renaming "%s" to "%s"' % (self.filename, bakname))
                self._rename_bundle(self.filename, bakname)
        # Ready to create a new bundle.
        self.log('Creating new "%s"' % (self.filename,))
        self._make_bundledir(self.filename)
        self.state = 'prewrite'
        self._write_metadata({})
        self._write_nframes(0)    # Mark new bundle as empty
        self._open_log()
        self.nframes = 0

    def _open_read(self):
        "Open a bundle trajectory for reading."
        if not os.path.exists(self.filename):
            raise IOError('File not found: ' + self.filename)
        if not self.is_bundle(self.filename):
            raise IOError('Not a BundleTrajectory: ' + self.filename)
        self.state = 'read'
        # Read the metadata
        metadata = self._read_metadata()
        self.metadata = metadata
        if metadata['version'] > self.version:
            raise NotImplementedError(
                "This version of ASE cannot read a BundleTrajectory version "
                + str(metadata['version']))
        if metadata['subtype'] != 'normal':
            raise NotImplementedError(
                "This version of ASE cannot read BundleTrajectory subtype "
                + metadata['subtype'])
        self._set_backend(metadata['backend'])
        self.nframes = self._read_nframes()
        if self.nframes == 0:
            raise IOError("Empty BundleTrajectory")
        self.datatypes = metadata['datatypes']
        
    def _write_nframes(self, n):
        "Write the number of frames in the bundle."
        assert self.state == 'write' or self.state == 'prewrite'
        f = paropen(os.path.join(self.filename, "frames"), "w")
        f.write(str(n)+'\n')
        f.close()

    def _read_nframes(self):
        "Read the number of frames."
        f = open(os.path.join(self.filename, 'frames'))
        n = int(f.read())
        return n

    def _write_metadata(self, metadata):
        """Write the metadata file of the bundle.

        Modifies the medadata dictionary!
        """
        # Add standard fields that must always be present.
        assert self.state == 'write' or self.state == 'prewrite'
        metadata['format'] = 'BundleTrajectory'
        metadata['version'] = self.version
        metadata['subtype'] = self.subtype
        metadata['backend'] = self.backend_name
        f = paropen(os.path.join(self.filename, "metadata"), "w")
        pickle.dump(metadata, f, -1)
        f.close()

    def _read_metadata(self):
        """Read the metadata."""
        assert self.state == 'read'
        f = open(os.path.join(self.filename, 'metadata'))
        metadata = pickle.load(f)
        f.close()
        return metadata

    @staticmethod
    def is_bundle(filename):
        """Check if a filename exists and is a BundleTrajectory."""
        if not os.path.isdir(filename):
            return False
        metaname = os.path.join(filename, 'metadata')
        if not os.path.isfile(metaname):
            return False
        f = open(metaname)
        mdata = pickle.load(f)
        f.close()
        try:
            return mdata['format'] == 'BundleTrajectory'
        except KeyError:
            return False

    @staticmethod
    def is_empty_bundle(filename):
        """Check if a filename is an empty bundle.  Assumes that it is a bundle."""
        f = open(os.path.join(filename, "frames"))
        nframes = int(f.read())
        return nframes == 0

    @classmethod
    def delete_bundle(cls, filename):
        "Deletes a bundle."
        if ase.parallel.rank == 0:
            # Only the master deletes
            if not cls.is_bundle(filename):
                raise IOError("Cannot remove '%s' as it is not a bundle trajectory."
                              % (filename,))
            if os.path.islink(filename):
                # A symbolic link to a bundle.  Only remove the link.
                os.remove(filename)
            else:
                # A real bundle
                shutil.rmtree(filename)
        else:
            # All other tasks wait for the directory to go away.
            while os.path.exists(filename):
                time.sleep(1)
        # The master may not proceed before all tasks have seen the
        # directory go away, as it might otherwise create a new bundle
        # with the same name, fooling the wait loop in _make_bundledir.
        ase.parallel.barrier()

    def _rename_bundle(self, oldname, newname):
        "Rename a bundle.  Used to create the .bak"
        if self.master:
            os.rename(oldname, newname)
        else:
            while os.path.exists(oldname):
                time.sleep(1)
        # The master may not proceed before all tasks have seen the
        # directory go away.
        ase.parallel.barrier()
            
    def _make_bundledir(self, filename):
        """Make the main bundle directory.

        Since all MPI tasks might write to it, all tasks must wait for
        the directory to appear.
        """
        self.log("Making directory "+filename)
        assert not os.path.isdir(filename)
        ase.parallel.barrier()
        if self.master:
            os.mkdir(filename)
        else:
            i = 0
            while not os.path.isdir(filename):
                time.sleep(1)
                i += 1
            if i > 10:
                self.log("Waiting %d seconds for %s to appear!"
                         % (i, filename))

    def _make_framedir(self, filename):
        """Make subdirectory for the frame.

        As only the master writes to it, no synchronization between
        MPI tasks is necessary.
        """
        if self.master:
            os.mkdir(filename)

    def pre_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called before writing begins.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError("Callback object must be callable.")
        self.pre_observers.append((function, interval, args, kwargs))

    def post_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called after writing ends.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError("Callback object must be callable.")
        self.post_observers.append((function, interval, args, kwargs))

    def _call_observers(self, obs):
        "Call pre/post write observers."
        for function, interval, args, kwargs in obs:
            if (self.nframes + 1) % interval == 0:
                function(*args, **kwargs)
    

class PickleBundleBackend:
    """Backend for writing BundleTrajectories stored as pickle files."""
    def __init__(self, master):
        self.master = master
        
    def write_small(self, framedir, smalldata):
        "Write small data to be written jointly."
        if self.master:
            f = open(os.path.join(framedir, "smalldata.pickle"), "w")
            pickle.dump(smalldata, f, -1)
            f.close()

    def write(self, framedir, name, data):
        "Write data to separate file."
        if self.master:
            fn = os.path.join(framedir, name + '.pickle')
            f = open(fn, "w")
            pickle.dump(data.shape, f, -1)
            pickle.dump(data, f, -1)
            f.close()

    def read_small(self, framedir):
        "Read small data."
        f = open(os.path.join(framedir, "smalldata.pickle"))
        data = pickle.load(f)
        f.close()
        return data

    def read(self, framedir, name):
        "Read data from separate file."
        fn = os.path.join(framedir, name + '.pickle')
        f = open(fn)
        shape = pickle.load(f)  # Discarded.
        data = pickle.load(f)
        f.close()
        return data

class NetCDFBundleBackend(PickleBundleBackend):
    dims = ('first', 'second', 'third', 'fourth', 'fifth')
    def write(self, framedir, name, data):
        if self.master:
            if data.dtype == np.int64:
                # Cannot write 64-bit integers!
                data = data.astype(np.int32)
            fn = os.path.join(framedir, name + '.nc')
            f = pupynere.netcdf_file(fn, "w")
            for i, n in enumerate(data.shape):
                f.createDimension(self.dims[i], n)
            var = f.createVariable(name, data.dtype,
                                   self.dims[:len(data.shape)])
            var[:] = data
            f.close()

    def read(self, framedir, name):
        fn = os.path.join(framedir, name + '.nc')
        f = pupynere.netcdf_file(fn)
        var = f.variables[name][:]
        f.close()
        return var
    
def read_bundletrajectory(filename, index=-1):
    """Reads one or more atoms objects from a BundleTrajectory.

    Arguments:

    filename: str
        The name of the bundle (really a directory!)
    index: int
        An integer specifying which frame to read, or an index object
        for reading multiple frames.  Default: -1 (reads the last
        frame).
    """
    traj = BundleTrajectory(filename, mode='r')
    if isinstance(index, int):
        return traj[index]
    else:
        # Here, we try to read only the configurations we need to read
        # and len(traj) should only be called if we need to as it will
        # read all configurations!

        # XXX there must be a simpler way?
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(traj)
            stop = index.stop or len(traj)
            if stop < 0:
                stop += len(traj)
        else:
            if index.start is None:
                start = len(traj) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(traj)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(traj)
                    
        return [traj[i] for i in range(start, stop, step)]

def write_bundletrajectory(filename, images):
    """Write image(s) to a BundleTrajectory.

    Write also energy, forces, and stress if they are already
    calculated.
    """

    traj = BundleTrajectory(filename, mode='w')

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    for atoms in images:
        # Avoid potentially expensive calculations:
        calc = atoms.get_calculator()
        if hasattr(calc, 'calculation_required'):
            for quantity in ('energy', 'forces', 'stress', 'magmoms'):
                traj.select_data(quantity,
                                 not calc.calculation_required(atoms, [quantity]))
        traj.write(atoms)
    traj.close()

        
if __name__ == '__main__':
    from ase.lattice.cubic import FaceCenteredCubic
    from ase.io import read, write
    atoms = FaceCenteredCubic(size=(5,5,5), symbol='Au')
    write('test.bundle', atoms)
    atoms2 = read('test.bundle')
    assert (atoms.get_positions() == atoms2.get_positions()).all()
    assert (atoms.get_atomic_numbers() == atoms2.get_atomic_numbers()).all()
    
    
