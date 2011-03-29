#!/usr//bin/env python

# lammps.py (2011/02/22)
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
import math
import decimal as dec
from subprocess import Popen
import numpy as np
from ase import Atoms
from ase.parallel import paropen


class LAMMPS:

    def __init__(self, label="lammps", dir="LAMMPS", parameters={}, files=[], always_triclinic=False):
        self.label = label
        self.dir = dir
        self.parameters = parameters
        self.files = files
        self.always_triclinic = always_triclinic
        self.calls = 0
        self.potential_energy = None
        self.forces = None

        if os.path.isdir(self.dir):
            if (os.listdir(self.dir) == []):
                os.rmdir(self.dir)
            else:
                os.rename(self.dir, self.dir + ".bak-" + time.strftime("%Y%m%d-%H%M%S"))
        os.mkdir(self.dir, 0755)
        
        for f in files:
            shutil.copy(f, os.path.join(self.dir, f))

    def clean(self, force=False):
        if os.path.isdir(self.dir):
            files_in_dir = os.listdir(self.dir)
            files_copied = self.files
            if (len(files_in_dir) == 0):
                os.rmdir(self.dir)
            elif (set(files_in_dir).issubset(set(files_copied)) or force):
                shutil.rmtree(self.dir)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.potential_energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        raise NotImplementedError

    def update(self, atoms):
        # TODO: check if (re-)calculation is necessary
        self.calculate(atoms)

    def calculate(self, atoms):
        self.atoms = atoms.copy()
        self.run()

    def run(self):
        """Method which explicitely runs LAMMPS."""

        self.calls += 1

        # set LAMMPS command from environment variable
        env = os.environ
        if 'LAMMPS_COMMAND' in env:
            lammps_cmd_line = shlex.split(env['LAMMPS_COMMAND'])
            if len(lammps_cmd_line) == 0:
                self.clean()
                raise RuntimeError('The LAMMPS_COMMAND environment variable '
                                   'must not be empty')
        else:
	    self.clean()
            raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
        if 'LAMMPS_OPTIONS' in env:
            lammps_options = shlex.split(env['LAMMPS_OPTIONS'])
        else:
            lammps_options = shlex.split("-echo log -screen none")

        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.dir)

        # setup file names for LAMMPS calculation
        label = "%s-%06d" % (self.label, self.calls)
        lammps_data = "data." + label
        lammps_in = "in." + label
        lammps_trj = label + ".lammpstrj"
        lammps_log = label + ".log"

        # write LAMMPS input files
        self.write_lammps_data(lammps_data=lammps_data)
        self.write_lammps_in(lammps_in=lammps_in, lammps_data=lammps_data, lammps_trj=lammps_trj, parameters=self.parameters)

        # run LAMMPS
        # TODO: check for successful completion (based on output files?!) of LAMMPS call...
        lmp_proc = Popen(lammps_cmd_line+lammps_options+\
                            ['-in', lammps_in, '-log', lammps_log])
        exitcode = lmp_proc.wait()
        if exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError("LAMMPS exited in %s with exit code: %d.  " % (cwd,exitcode))

        # read LAMMPS output files
        self.read_lammps_log(lammps_log=lammps_log)
        self.read_lammps_trj(lammps_trj=lammps_trj)

        # change back to previous working directory
        os.chdir(cwd)

    def write_lammps_data(self, lammps_data=None):
        """Method which writes a LAMMPS data file with atomic structure."""
        if (lammps_data == None):
            lammps_data = "data." + self.label
        write_lammps(lammps_data, self.atoms, force_skew=self.always_triclinic)

    def write_lammps_in(self, lammps_in=None, lammps_data=None, lammps_trj=None, parameters={}):
        """Method which writes a LAMMPS in file with run parameters and settings."""
        if (lammps_in == None):
            lammps_in = "in." + self.label
        if (lammps_data == None):
            lammps_data = "data." + self.label
        if (lammps_trj == None):
            lammps_trj = self.label + ".lammpstrj"

        f = paropen(lammps_in, 'w')
        f.write("# " + f.name + " (written by ASE) \n")
        f.write("\n")

        f.write("### variables \n")
        f.write("variable data_file index \"%s\" \n" % lammps_data)
        f.write("variable dump_file index \"%s\" \n" % lammps_trj)
        f.write("\n\n")

        pbc = self.atoms.get_pbc()
        f.write("### simulation box \n")
        f.write("units metal \n")
        if ("boundary" in parameters):
            f.write("boundary %s \n" % parameters["boundary"])
        else:
            f.write("boundary %c %c %c \n" % tuple('sp'[x] for x in pbc))
        f.write("atom_modify sort 0 0.0 \n")
        for key in ("neighbor" ,"newton"):
            if key in parameters:
                f.write("%s %s \n" % (key, parameters[key]))
        f.write("\n")
        f.write("read_data ${data_file} \n")
        f.write("\n\n")

        f.write("### interactions \n")
        if ( ("pair_style" in parameters) and ("pair_coeff" in parameters) ):
            pair_style = parameters["pair_style"]
            f.write("pair_style %s \n" % pair_style)
            for pair_coeff in parameters["pair_coeff"]:
                f.write("pair_coeff %s \n" % pair_coeff)
            if "mass" in parameters:
                for mass in parameters["mass"]:
                    f.write("mass %s \n" % mass)
        else:
            # simple default parameters that should always make the LAMMPS calculation run
            #    pair_style 	lj/cut 2.5
            #    pair_coeff 	* * 1 1
            #    mass           * 1.0        else:
            f.write("pair_style lj/cut 2.5 \n")
            f.write("pair_coeff * * 1 1 \n")
            f.write("mass * 1.0 \n")
        f.write("\n")

        f.write("### run \n")
        f.write("fix fix_nve all nve \n")
        f.write("\n")
        f.write("dump dump_all all custom 1 ${dump_file} id type x y z vx vy vz fx fy fz \n")
        f.write("\n")
        f.write("thermo_style custom step temp ke pe etotal cpu \n")
        f.write("thermo_modify format 1 %4d format 2 %9.3f format 3 %20.6f format 4 %20.6f format 5 %20.6f format 6 %9.3f \n")
        f.write("thermo 1 \n")
        f.write("\n")
        if ("minimize" in parameters):
            f.write("minimize %s \n" % parameters["minimize"])
        if ("run" in parameters):
            f.write("run %s \n" % parameters["run"])
        if not ( ("minimize" in parameters) or ("run" in parameters) ):
            f.write("run 0 \n")

        f.close()

    def read_lammps_log(self, lammps_log=None):
        """Method which reads a LAMMPS output log file."""
        if (lammps_log == None):
            lammps_log = self.label + ".log"

        epot = 0.0

        f = paropen(lammps_log, 'r')
        while True:
            line = f.readline()
            if not line:
                break
            # get potential energy of first step (if more are done)
            if "PotEng" in line:
                i = line.split().index("PotEng")
                line = f.readline()
                epot = float(line.split()[i])
        f.close()

#        print epot

        self.potential_energy = epot

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


# could perhaps go into io.lammps in the future...
#from ase.parallel import paropen
def write_lammps(fileobj, atoms, force_skew=False):
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

    # Uniqify 'symbols' list and sort to a new list 'species'.
    # Sorting is important in case of multi-component systems:
    # This way it is assured that LAMMPS atom types are always
    # assigned predictively according to the alphabetic order 
    # of the species present in the current system. 
    # (Hence e.g. the mapping in interaction statements for alloy
    #  potentials depending on the order of the elements in the
    #  external potential can also be safely adjusted accordingly.)
    species = sorted(list(set(symbols)))
    n_atom_types = len(species)
    f.write("%d \t atom types \n" % n_atom_types)


    ### cell -> bounding box for LAMMPS

    ## A) simple version
    # - For orthogonal cells which are properly aligned with the reference coordinate system only.
#    cell = atoms.get_cell()
#    lo = [0.0, 0.0, 0.0]
#    hi = sum(cell[0:3])
#    f.write("%15.8f %15.8f \t xlo xhi \n" % (lo[0], hi[0]))
#    f.write("%15.8f %15.8f \t ylo yhi \n" % (lo[1], hi[1]))
#    f.write("%15.8f %15.8f \t zlo zhi \n" % (lo[2], hi[2]))
#    f.write("\n\n")

    ## B) convert to (perhaps) triclinic cells in LAMMPS
    
    # B.1) calculate lengths of and angles between cell vectors
    #     (conventional crystallographic nomenclature)
    cell = atoms.get_cell()
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = np.arccos( np.vdot(cell[1], cell[2]) / (b * c) )
    beta = np.arccos( np.vdot(cell[0], cell[2]) / (a * c) )
    gamma = np.arccos( np.vdot(cell[0], cell[1]) / (a * b) )
#    print [a, b, c]
#    print map(np.degrees, [alpha, beta, gamma])

    # B.2) construct bounding box edges and skew vector of 
    #     corresponding (perhaps rotated!) LAMMPS cell
    # http://lammps.sandia.gov/doc/read_data.html :
    # 	a_LAMMPS = (xhi-xlo,0,0); b_LAMMPS = (xy,yhi-ylo,0); c_LAMMPS = (xz,yz,zhi-zlo)
    xlo = ylo = zlo = 0.0
    # this choice of origin simplifies things:
    # 	a_LAMMPS = (xhi,0,0); b_LAMMPS = (xy,yhi,0); c_LAMMPS = (xz,yz,zhi)
    xhi = a			# a_LAMMPS
    xy = np.cos(gamma) * b	# b_LAMMPS
    yhi = np.sin(gamma) * b
    xz = np.cos(beta) * c	# c_LAMMPS
    yz = ( b * c * np.cos(alpha) - xy * xz ) / yhi
    zhi = np.sqrt( c**2 - xz**2 - yz**2 )
    
    # B.3) obtain rotation of original cell with respect to LAMMPS cell
    #     to properly tranform the atoms contained within (see below!)
    # IMPORTANT: use same convention as in cell of ASE atoms object, i.e.
    #            cell vectors are ROW VECTORS in numpy array
    cell_lammps = np.array([[xhi-xlo,0,0],[xy,yhi-ylo,0],[xz,yz,zhi-zlo]])
#    a_lammps = np.linalg.norm(cell_lammps[0])
#    b_lammps = np.linalg.norm(cell_lammps[1])
#    c_lammps = np.linalg.norm(cell_lammps[2])
#    alpha_lammps = np.arccos( np.vdot(cell_lammps[1], cell_lammps[2]) / (b_lammps * c_lammps) )
#    beta_lammps = np.arccos( np.vdot(cell_lammps[0], cell_lammps[2]) / (a_lammps * c_lammps) )
#    gamma_lammps = np.arccos( np.vdot(cell_lammps[0], cell_lammps[1]) / (a_lammps * b_lammps) )
#    print [a_lammps, b_lammps, c_lammps]
#    print map(np.degrees, [alpha_lammps, beta_lammps, gamma_lammps])
    # IMPORTANT: need vector-(rotation-)matrix product (instead of matrix-vector product) here,
    #            cell vectors are ROW VECTORS (also see above)
    rotation = np.dot(np.linalg.inv(cell), cell_lammps)
#    print rotation
#    print np.transpose(rotation)
#    print np.linalg.det(rotation)				# check for orthogonality of matrix 'rotation'
#    print np.dot(rotation, np.transpose(rotation))		# check for orthogonality of matrix 'rotation'
#    print np.dot(rotation, cell[0])
#    print np.dot(rotation, cell[1])
#    print np.dot(rotation, cell[2])

    # B.4.1) write bounding box edges
    f.write("%15.8f %15.8f \t\t\t xlo xhi \n" % (xlo, xhi))
    f.write("%15.8f %15.8f \t\t\t ylo yhi \n" % (ylo, yhi))
    f.write("%15.8f %15.8f \t\t\t zlo zhi \n" % (zlo, zhi))
    
    # B.4.2) sanitize and write skew vector (if necessary)
    # The too simple version 
    #    f.write("%20.10f %20.10f %20.10f \t xy xz yz \n" % (xy, xz, yz))    
    # can make LAMMPS (easily) reject the definition of a triclinic box because
    #	a) skew vector components calculated above are outside the respective ranges 
    #	   which LAMMPS expects as described here
    #		http://lammps.sandia.gov/doc/read_data.html
    #	   and coded here
    #		domain.cpp -> Domain::set_initial_box()
    # 	b) rounding up in the skew vector output can make the result violate those 
    #      ranges for purely numeric reasons (perhaps even more frequent problem than a)
    #      and definitely more stupid!)
    # More elaborate solution addressing the issues described above:
    # -> a):
    # Fold back skew vector components calculated above to equivalent ones complying with
    # what LAMMPS expects (see references above) and describing another ("less skewed") 
    # unit cell for the same lattice:
    #		t -> t_folded in [-period/2, +period/2]
    def fold_skew_dec(t, period) :
        t_dec = dec.Decimal(repr(t))
        period_dec = dec.Decimal(repr(period))
	# dec.ROUND_HALF_DOWN leaves all t in [-period/2, +period/2] unchanged whereas
	# dec.ROUND_HALF_UP does the same for all t in (-period/2, +period/2) but
	# swaps the boundary values: 
	#	t=-period/2 -> t_folded=+period/2
	#	t=+period/2 -> t_folded=-period/2
        t_folded_dec = t_dec - period_dec * (t_dec/period_dec).to_integral_value(dec.ROUND_HALF_DOWN)
        return t_folded_dec
    skew_folded_dec = [ fold_skew_dec(xy, xhi-xlo), fold_skew_dec(xz, xhi-xlo), fold_skew_dec(yz, yhi-ylo) ]
    # -> b):
    # This should be kept in sync with the accuracy used for writing the bounding box edges above
    prec_dec = dec.Decimal('1E-8')
    # Make sure to always round down and hence avoid numerical problems with skew vector in LAMMPS.
    skew_folded_and_rounded_dec = [ d.quantize(prec_dec, dec.ROUND_DOWN) for d in skew_folded_dec ]
    # Do not write zero skew vector. 
    # (I.e. calculate orthogonal cell with LAMMPS in that case!)
    has_skew = any( [not d.is_zero() for d in skew_folded_and_rounded_dec] )
    if force_skew or has_skew:
        f.write( "%15s %15s %15s \t xy xz yz \n" % tuple( str(d) for d in skew_folded_and_rounded_dec ) )
    f.write("\n\n")


    ### atoms

    ## A) simple version
    # - For orthogonal cells which are properly aligned with the reference coordinate system only.
#    f.write("Atoms \n\n")
#    for i, (x, y, z) in enumerate(atoms.get_positions()):
#        s = species.index(symbols[i]) + 1
#        f.write("%6d %3d %22.15f %22.15f %22.15f \n" % (i+1, s, x, y, z))

    ## B) convert to (perhaps) triclinic cells in LAMMPS:
    # - Apply rotation determined above for this case.
    # - Lattice translation in case of a folded skew vector above should not be necessary,
    #   since the resulting (perhaps) new unit cell does describe the same lattice.
    #   Therefore the (rotated) atoms are still at correct lattice sites, only some of them
    #   probably are outside of that new unit cell, which LAMMPS can hopefully deal with.
    f.write("Atoms \n\n")
    for i, r in enumerate(atoms.get_positions()):
        s = species.index(symbols[i]) + 1
        [x, y, z] = np.dot(r, rotation)
#        print r, [x, y, z]
        f.write("%6d %3d %22.15f %22.15f %22.15f \n" % (i+1, s, x, y, z))

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

#    calc.clean(force=True)
    calc.clean()



