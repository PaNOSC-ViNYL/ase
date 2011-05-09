#!/usr//bin/env python

# lammps.py (2011/03/29)
# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable LAMMPS_COMMAND must be defined to point to the LAMMPS binary.
#
# Copyright (C) 2009 - 2011 Joerg Meyer, joerg.meyer@ch.tum.de
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or see <http://www.gnu.org/licenses/>.


import os 
import shutil
import shlex
import time
from subprocess import Popen, PIPE
from threading import Thread
from re import compile as re_compile, IGNORECASE
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
from ase import Atoms
from ase.parallel import paropen
from ase.units import GPa

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = "__end_of_ase_invoked_calculation__"

class LAMMPS:

    def __init__(self, label="lammps", tmp_dir=None, parameters={}, 
                 specorder=None, files=[], always_triclinic=False, 
                 keep_alive=True, keep_tmp_files=False):
        self.label = label
        self.tmp_dir = tmp_dir
        self.parameters = parameters
        self.specorder = specorder
        self.files = files
        self.always_triclinic = always_triclinic
        self.calls = 0
        self.forces = None
        self.keep_alive = keep_alive
        self.keep_tmp_files = keep_tmp_files
        self._lmp_handle = None        # To handle the lmp process

        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched aginst the log output. I.e. 
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        self._custom_thermo_args = ['step', 'temp', 'press', 'cpu', 
                                    'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                                    'ke', 'pe', 'etotal',
                                    'vol', 'lx', 'ly', 'lz']
        self._custom_thermo_mark = ' '.join([x.capitalize() for x in
                                             self._custom_thermo_args[0:3]])

        # Match something which can be converted to a float
        f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
        n = len(self._custom_thermo_args)
        # Create a re matching exactly N white space separated floatish things
        self._custom_thermo_re = re_compile(r'^\s*' + r'\s+'.join([f_re]*n) + r'\s*$',
                                            flags=IGNORECASE)
        # thermo_content contains data "writen by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each 
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be 
        # re-populated by the read_log method.
        self.thermo_content = []

        if self.tmp_dir is None:
            tmp_dir = self.tmp_dir = mkdtemp(prefix='LAMMPS')
        else:
            tmp_dir = self.tmp_dir = os.path.abspath(self.tmp_dir)
            if os.path.isdir(tmp_dir):
                if (os.listdir(tmp_dir) == []):
                    os.rmdir(tmp_dir)
                else:
                    os.rename(tmp_dir, 
                              tmp_dir+".bak-"+time.strftime("%Y%m%d-%H%M%S"))
            os.mkdir(tmp_dir, 0755)
        
        for f in files:
            shutil.copy(f, os.path.join(tmp_dir, f))

    def clean(self, force=False):

        self._lmp_end()

        if not self.keep_tmp_files:
            shutil.rmtree(self.tmp_dir)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.thermo_content[-1]['pe']

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        tc = self.thermo_content[-1]
        # 1 bar (used by lammps for metal units) = 1e-4 GPa
        return np.array([tc[i] for i in ('pxx','pyy','pzz',
                                         'pyz','pxz','pxy')])*(-1e-4*GPa)

    def update(self, atoms):
        # TODO: check if (re-)calculation is necessary
        self.calculate(atoms)

    def calculate(self, atoms):
        self.atoms = atoms.copy()
        self.run()

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running lammps process
        return self._lmp_handle and not isinstance(self._lmp_handle.poll(), int)

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            return self._lmp_handle.wait()
        
    def run(self):
        """Method which explicitely runs LAMMPS."""

        self.calls += 1

        # set LAMMPS command from environment variable
        if 'LAMMPS_COMMAND' in os.environ:
            lammps_cmd_line = shlex.split(os.environ['LAMMPS_COMMAND'])
            if len(lammps_cmd_line) == 0:
                self.clean()
                raise RuntimeError('The LAMMPS_COMMAND environment variable '
                                   'must not be empty')
            # want always an absolute path to LAMMPS binary when calling from self.dir                       
            lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])

        else:
	    self.clean()
            raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
        if 'LAMMPS_OPTIONS' in os.environ:
            lammps_options = shlex.split(os.environ['LAMMPS_OPTIONS'])
        else:
            lammps_options = shlex.split("-echo log -screen none")

        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)

        # setup file names for LAMMPS calculation
        label = "%s-%06d" % (self.label, self.calls)
        lammps_in = "in."+label
        lammps_log = label+".log"
        lammps_trj_fd = NamedTemporaryFile(prefix=label+".", suffix=".lammpstrj",
                                           dir=self.tmp_dir, 
                                           delete=(not self.keep_tmp_files))
        lammps_trj = lammps_trj_fd.name
 
        # see to it that LAMMPS is started
        if not self._lmp_alive():
            # Attempt to (re)start lammps
            self._lmp_handle = Popen(lammps_cmd_line+lammps_options+['-log', '/dev/stdout'], 
                                    stdin=PIPE, stdout=PIPE)
        lmp_handle = self._lmp_handle

        # Create thread reading lammps stdout (for reference, if requested,
        # also create lammps_log, although it is never used)
        if self.keep_tmp_files:
            lammps_log_fd = open(lammps_log, 'w')
            fd = special_tee(lmp_handle.stdout, lammps_log_fd)
        else:
            fd = lmp_handle.stdout
        thr_read_log = Thread(target=self.read_lammps_log, args=(fd,))
        thr_read_log.start()

        # write LAMMPS input (for reference, also create the file lammps_in, 
        # although it is never used)
        if self.keep_tmp_files:
            lammps_in_fd = open(lammps_in, 'w')
            fd = special_tee(lmp_handle.stdin, lammps_in_fd)
        else:
            fd = lmp_handle.stdin
        self.write_lammps_in(lammps_in=fd, lammps_trj=lammps_trj)

        if self.keep_tmp_files:
            lammps_in_fd.close()

        # Wait for log output to be read (i.e., for LAMMPS to finish)
        # and close the log file if there is one
        thr_read_log.join()
        if self.keep_tmp_files:
            lammps_log_fd.close()

        if not self.keep_alive:
            self._lmp_end()

        exitcode = lmp_handle.poll()
        if exitcode and exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError("LAMMPS exited in %s with exit code: %d." %\
                                   (cwd,exitcode))

        # read LAMMPS output files
        self.read_lammps_trj(lammps_trj=lammps_trj)
        lammps_trj_fd.close()

        # change back to previous working directory
        os.chdir(cwd)

    def write_lammps_data(self, lammps_data=None):
        """Method which writes a LAMMPS data file with atomic structure."""
        if (lammps_data == None):
            lammps_data = "data." + self.label
        write_lammps_data(lammps_data, self.atoms, self.specorder, 
                          force_skew=self.always_triclinic)

    def write_lammps_in(self, lammps_in=None, lammps_trj=None):
        """Method which writes a LAMMPS in file with run parameters and settings."""

        if (lammps_in == None):
            lammps_in = "in." + self.label
        if (lammps_trj == None):
            lammps_trj = self.label + ".lammpstrj"

        if isinstance(lammps_in, str):
            f = paropen(lammps_in, 'w')
            close_in_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_in
            close_in_file = False

        f.write("# (written by ASE) \n" +
                "clear\n" +
                "\n### variables \n" +
                ("variable dump_file string \"%s\" \n" % lammps_trj ))

        parameters = self.parameters
        pbc = self.atoms.get_pbc()
        f.write("\n### simulation box \n"
                "units metal \n")
        if ("boundary" in parameters):
            f.write("boundary %s \n" % parameters["boundary"])
        else:
            f.write("boundary %c %c %c \n" % tuple('sp'[x] for x in pbc))
        f.write("atom_modify sort 0 0.0 \n")
        for key in ("neighbor" ,"newton"):
            if key in parameters:
                f.write("%s %s \n" % (key, parameters[key]))
        f.write("\n")

        # Write the simulation box and the atoms
        p = prism(self.atoms.get_cell())
        f.write("lattice sc 1.0\n")
        xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism()
        if self.always_triclinic or p.is_skewed():
            f.write("region asecell prism 0.0 %15.10f 0.0 %15.10f 0.0 %15.10f " % \
                        (xhi, yhi, zhi))
            f.write("%15.10f %15.10f %15.10f side in units box\n" % \
                        (xy, xz, yz))
        else:
            f.write(("region asecell block 0.0 %15.10f 0.0 %15.10f 0.0 %15.10f "
                    "side in units box\n") % (xhi, yhi, zhi))
                    
        symbols = self.atoms.get_chemical_symbols()
        if self.specorder is None:
            # By default, atom types in alphabetic order
            species = sorted(list(set(symbols)))
        else:
            # By request, specific atom type ordering
            species = self.specorder

        n_atom_types = len(species)
        species_i = dict([(s,i+1) for i,s in enumerate(species)])

        f.write("create_box %i asecell\n" % n_atom_types)
        for s, r in zip(symbols, 
                        map(p.trans_pos_to_lammps, self.atoms.get_positions())):
            f.write("create_atoms %i single %15.10f %15.10f %15.10f units box\n" % \
                        ((species_i[s],)+tuple(r)))

        # Write interaction stuff
        f.write("\n### interactions \n")
        if ( ("pair_style" in parameters) and ("pair_coeff" in parameters) ):
            pair_style = parameters["pair_style"]
            f.write("pair_style %s \n" % pair_style)
            for pair_coeff in parameters["pair_coeff"]:
                f.write("pair_coeff %s \n" % pair_coeff)
            if "mass" in parameters:
                for mass in parameters["mass"]:
                    f.write("mass %s \n" % mass)
        else:
            # simple default parameters that should always make the LAMMPS 
            # calculation run
            f.write("pair_style lj/cut 2.5 \n" +
                    "pair_coeff * * 1 1 \n" +
                    "mass * 1.0 \n")

        f.write("\n### run\n" +
                "fix fix_nve all nve\n" +
                ("dump dump_all all custom 1 %s id type x y z vx vy vz fx fy fz\n" % lammps_trj) )
        f.write(("thermo_style custom %s\n" +
                "thermo_modify flush yes\n" +
                "thermo 1\n") % (' '.join(self._custom_thermo_args)))
        if "minimize" in parameters:
            f.write("minimize %s\n" % parameters["minimize"])
        if "run" in parameters:
            f.write("run %s\n" % parameters["run"])
        if not (("minimize" in parameters) or ("run" in parameters)):
            f.write("run 0\n")

        f.write('print "%s"\n' % CALCULATION_END_MARK)
        f.write('log /dev/stdout\n') # Force LAMMPS to flush log

        if close_in_file:
            f.close()

    def read_lammps_log(self, lammps_log=None, PotEng_first=False):
        """Method which reads a LAMMPS output log file."""

        if (lammps_log == None):
            lammps_log = self.label + ".log"

        if isinstance(lammps_log, str):
            f = paropen(lammps_log, 'w')
            close_log_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_log
            close_log_file = False

        thermo_content = []
        line = f.readline()
        while line and line.strip() != CALCULATION_END_MARK:
            # get thermo output
            if line.startswith(self._custom_thermo_mark):
                m = True
                while m:
                    line = f.readline()
                    m = self._custom_thermo_re.match(line)
                    if m:
                        # create a dictionary between each of the thermo_style args
                        # and it's corresponding value
                        thermo_content.append(dict(zip(self._custom_thermo_args, 
                                                       map(float, m.groups()))))
            else:
                line = f.readline()

        if close_log_file:
            f.close()

        self.thermo_content = thermo_content

    def read_lammps_trj(self, lammps_trj=None, set_atoms=False):
        """Method which reads a LAMMPS dump file."""
        if (lammps_trj == None):
            lammps_trj = self.label + ".lammpstrj"

        f = paropen(lammps_trj, 'r')
        while True:
            line = f.readline()

            if not line:
                break

            #TODO: extend to proper dealing with multiple steps in one trajectory file
            if "ITEM: TIMESTEP" in line:
                n_atoms = 0
                lo = [] ; hi = [] ; tilt = []
                id = [] ; type = []
                positions = [] ; velocities = [] ; forces = []

            if "ITEM: NUMBER OF ATOMS" in line:
                line = f.readline()
                n_atoms = int(line.split()[0])
            
            if "ITEM: BOX BOUNDS" in line:
                # save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
                tilt_items = line.split()[3:]
                for i in range(3):
                    line = f.readline()
                    fields = line.split()
                    lo.append(float(fields[0]))
                    hi.append(float(fields[1]))
                    if (len(fields) >= 3):
                        tilt.append(float(fields[2]))
            
            if "ITEM: ATOMS" in line:
                # (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
                # create corresponding index dictionary before iterating over atoms to (hopefully) speed up lookups...
                atom_attributes = {}
                for (i, x) in enumerate(line.split()[2:]):
                    atom_attributes[x] = i
                for n in range(n_atoms):
                    line = f.readline()
                    fields = line.split()
                    id.append( int(fields[atom_attributes['id']]) )
                    type.append( int(fields[atom_attributes['type']]) )
                    positions.append( [ float(fields[atom_attributes[x]]) for x in ['x', 'y', 'z'] ] )
                    velocities.append( [ float(fields[atom_attributes[x]]) for x in ['vx', 'vy', 'vz'] ] )
                    forces.append( [ float(fields[atom_attributes[x]]) for x in ['fx', 'fy', 'fz'] ] )
        f.close()

        # determine cell tilt (triclinic case!)
        if (len(tilt) >= 3):
            # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS" to assign tilt (vector) elements ...
            if (len(tilt_items) >= 3):
                xy = tilt[tilt_items.index('xy')]
                xz = tilt[tilt_items.index('xz')]
                yz = tilt[tilt_items.index('yz')]
            # ... otherwise assume default order in 3rd column (if the latter was present)
            else:
                xy = tilt[0]
                xz = tilt[1]
                yz = tilt[2]
        else:
            xy = xz = yz = 0
        xhilo = (hi[0] - lo[0]) - xy - xz
        yhilo = (hi[1] - lo[1]) - yz
        zhilo = (hi[2] - lo[2])
	
	# The simulation box bounds are included in each snapshot and if the box is triclinic (non-orthogonal), 
	# then the tilt factors are also printed; see the region prism command for a description of tilt factors. 
	# For triclinic boxes the box bounds themselves (first 2 quantities on each line) are a true "bounding box" 
	# around the simulation domain, which means they include the effect of any tilt.
	# [ http://lammps.sandia.gov/doc/dump.html , lammps-7Jul09 ]
	#
	# This *should* extract the lattice vectors that LAMMPS uses from the true "bounding box" printed in the dump file
	# It might fail in some cases (negative tilts?!) due to the MIN / MAX construction of these box corners:
	#
	#	void Domain::set_global_box() 
	#	[...]
	#	  if (triclinic) {
	#	    [...]
	#	    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
	#	    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
	#	    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
	#	    boxlo_bound[2] = boxlo[2];
	#
	#	    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
	#	    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
	#	    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
	#	    boxhi_bound[2] = boxhi[2];
	#	  }
	# [ lammps-7Jul09/src/domain.cpp ]
	#
        cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]

#       print n_atoms
#       print lo
#       print hi
#       print id
#       print type
#       print positions
#       print velocities
#       print forces

        # assume that LAMMPS does not reorder atoms internally
#        print np.array(positions)
#        print self.atoms.get_positions()
#        print np.array(positions) - self.atoms.get_positions()

        cell_atoms = np.array(cell)
        type_atoms = np.array(type)
#      velocities_atoms = np.array(velocities)
        positions_atoms = np.array(positions)
        forces_atoms = np.array(forces)

        if self.atoms:
            cell_atoms = self.atoms.get_cell()

            rotation_lammps2ase = np.dot(np.linalg.inv(np.array(cell)), cell_atoms)
#           print "rotation_lammps2ase:"
#           print rotation_lammps2ase
#           print np.transpose(rotation_lammps2ase)
#           print np.linalg.det(rotation_lammps2ase)                                            # check for orthogonality of matrix 'rotation'
#           print np.dot(rotation_lammps2ase, np.transpose(rotation_lammps2ase))                # check for orthogonality of matrix 'rotation'

            type_atoms = self.atoms.get_atomic_numbers()
            positions_atoms = np.array( [np.dot(np.array(r), rotation_lammps2ase) for r in positions] )
            velocities_atoms = np.array( [np.dot(np.array(v), rotation_lammps2ase) for v in velocities] )
            forces_atoms = np.array( [np.dot(np.array(f), rotation_lammps2ase) for f in forces] )

        if (set_atoms):
            # assume periodic boundary conditions here (like also below in write_lammps) <- TODO:?!
            self.atoms = Atoms(type_atoms, positions=positions_atoms, cell=cell_atoms, pbc=True)

        self.forces = forces_atoms



class special_tee:
    """A special purpose, with limited applicability, tee-like thing.

    A subset of stuff read from, or written to, orig_fd, 
    is also written to out_fd.
    It is used by the lammps calculator for creating file-logs of stuff read from,
    or written to, stdin and stdout, respectively.
    """
    def __init__(self, orig_fd, out_fd):
        self._orig_fd = orig_fd
        self._out_fd = out_fd
        self.name = orig_fd.name
    def write(self, data):
        self._orig_fd.write(data)
        self._out_fd.write(data)
    def read(self, *args, **kwargs):
        data = self._orig_fd.read(*args, **kwargs)
        self._out_fd.write(data)
        return data
    def readline(self, *args, **kwargs):
        data = self._orig_fd.readline(*args, **kwargs)
        self._out_fd.write(data)
        return data
    def readlines(self, *args, **kwargs):
        data = self._orig_fd.readlines(*args, **kwargs)
        self._out_fd.write(''.join(data))
        return data
    def flush(self):
        self._orig_fd.flush()
        self._out_fd.flush()


class prism:
    def __init__(self, cell):
        """Create a lammps-style triclinic prism object from a cell 
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]
        
        alpha = np.arccos(np.dot(b, c)/(bn*cn))
        beta  = np.arccos(np.dot(a, c)/(an*cn))
        gamma = np.arccos(np.dot(a, b)/(an*bn))
        
        xhi = an
        xyp = np.cos(gamma)*bn
        yhi = np.sin(gamma)*bn
        xzp = np.cos(beta)*cn
        yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    
        def fold_s(t, ref):
            return np.mod(0.5*ref+t, ref)-0.5*ref
        
        xy = fold_s(xyp, xhi)
        xz = fold_s(xzp, xhi)
        yz = fold_s(yzp, yhi)
        self.accuracy = np.finfo(xhi).eps
        self.prism = (xhi, yhi, zhi, xy, xz, yz)
        self.A = np.array(((xhi, 0,   0),
                           (xy,  yhi, 0),
                           (xz,  yz,  zhi)))
        self.R = np.dot(np.linalg.inv(cell), self.A)
        self.Ainv = np.linalg.inv(self.A)

    def dir2car(self, v):
        "Direct to cartesian coordinates"
        return np.dot(v, self.A)

    def car2dir(self, v):
        "Cartesian to direct coordinates"
        return np.dot(v, self.Ainv)

    def fold(self,v):
        "Fold a position into the (lammps) cell" 
        # This method is unfortunately ugly.
        # The "two-stage fold" is there to avoid atoms ending up at the far ends  
        # of the box, since those atoms would be ignored by lammps when using 
        # create_atoms (they would fail to fall within any sub-box bound; 
        # see lammps create_atoms.cpp). 
        x = self.dir2car(np.mod(self.car2dir(v), 1))
        return self.dir2car(np.mod(self.car2dir(x), 1-self.accuracy))

    def get_lammps_prism(self):
        return self.prism

    def trans_pos_to_lammps(self, position):
        "Rotate and fold a postion into the lammps cell"
        return self.fold(np.dot(self.R, position))

    def is_skewed(self):
        acc = self.accuracy
        xy, xz, yz = [np.abs(x) for x in self.prism[3:]]
        return (xy > acc) or (xz > acc) or (yz > acc)
        

def write_lammps_data(fileobj, atoms, specorder=[], force_skew=False):
    """Method which writes atomic structure data to a LAMMPS data file."""
    if isinstance(fileobj, file):
        f = fileobj
    elif isinstance(fileobj, str):
        f = paropen(fileobj, 'w')
    else :
        raise TypeError('fileobj is not of type file or str!')

    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration to a lammps data file!')
        atoms = atoms[0]

    f.write(f.name + " (written by ASE) \n\n")

    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    f.write("%d \t atoms \n" % n_atoms)

    if len(specorder):
        # To index elements in the LAMMPS data file (indices must
        # correspond to order in the potential file)
        species = specorder
    else :
        # This way it is assured that LAMMPS atom types are always
        # assigned predictively according to the alphabetic order 
        species = sorted(list(set(symbols)))
    n_atom_types = len(species)
    f.write("%d  atom types\n" % n_atom_types)

    p = prism(atoms.get_cell())
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism()

    f.write("%16.10f %16.10f  xlo xhi\n" % (0.0, xhi))
    f.write("%16.10f %16.10f  ylo yhi\n" % (0.0, yhi))
    f.write("%16.10f %16.10f  zlo zhi\n" % (0.0, zhi))
    
    if force_skew or p.is_skewed():
        f.write("%15.10f %15.10f %15.10f  xy xz yz\n" % (xy, xz, yz))
    f.write("\n\n")

    f.write("Atoms \n\n")
    for i, r in enumerate(map(p.trans_pos_to_lammps,
                              atoms.get_positions())):
        s = species.index(symbols[i]) + 1
        f.write("%6d %3d %16.10f %16.10f %16.10f\n" % ((i+1, s)+tuple(r)))

    f.close()


if __name__ == "__main__":

    pair_style = "eam"
    Pd_eam_file = "Pd_u3.eam"
    pair_coeff = [ "* * " + Pd_eam_file ]
    parameters = { "pair_style" : pair_style, "pair_coeff" : pair_coeff }
    files = [ Pd_eam_file ]
    calc = LAMMPS(parameters=parameters, files=files)
    from ase import Atoms
    a0 = 3.93
    b0 = a0 / 2.0
    if True:
        bulk = Atoms(['Pd']*4,
                     positions=[(0,0,0),(b0,b0,0),(b0,0,b0),(0,b0,b0)],
                     cell=[a0]*3,
                     pbc=True)
        # test get_forces
        print "forces for a = %f" % a0
        print calc.get_forces(bulk)
        # single points for various lattice constants
        bulk.set_calculator(calc)
        for n in range(-5,5,1):
            a = a0 * (1 + n/100.0)
            bulk.set_cell([a]*3)
            print "a : %f , total energy : %f" % (a, bulk.get_potential_energy())

    calc.clean()



