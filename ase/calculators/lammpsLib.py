"""ASE LAMMPS Calculator Library Verison"""

import numpy as np
from numpy.linalg import norm
from lammps import lammps
import ctypes

class LAMMPSLib:

	def __init__(self, parameters=[], atom_types=None, log_file = None):
		"""
		parameters is a sequence of strings that are full LAMMPS commands
		to setup things like pair_style, pair_coeff, pbc etc.

		atom_types is a dictionary of "atomic_symbol":lammps_atom_type pairs,
		e.g. "Cu":1 to bind copper to lammps atom type 1.
		Default assigns lammps atom types in alphabetic order by atomic symbol
		
		log_file is a string of the path to the desired LAMMPS log file
		"""
		self.parameters = parameters
		self.lmp = None
		self.atoms = None
		self.atom_types = atom_types
		self.coord_transform = None
		self.log_file = log_file
		self.f = None
		self.u = None

	def calculation_req(self, atoms):
		#returns 0 if no calculation needed
		#returns 1 if calculation needed
		#returns 2 if calculation and new lmp object are needed
		if  self.atoms == None:
			return 2
		if atoms.get_number_of_atoms() != self.atoms.get_number_of_atoms():
			return 2
		if atoms.get_cell().all() != self.atoms.get_cell().all():
			return 2
		if atoms!=self.atoms or self.f==None or self.u==None or atoms==None:
			return 1
		return 0

	def get_forces(self,atoms):
		calcReq = self.calculation_req(atoms)
		if calcReq == 1:
			self.atoms = atoms.copy()
			self.calculate()
		elif calcReq == 2:
			self.atoms = atoms.copy()
			self.make_new_lammps()
			self.calculate()
		return self.f	
	
	def get_potential_energy(self,atoms=None, force_consistent=False):
		calcReq = self.calculation_req(atoms)
		if calcReq == 1:
			self.atoms = atoms.copy()
			self.calculate()
		elif calcReq == 2:
			self.atoms = atoms.copy()
			self.make_new_lammps()
			self.calculate()
		return self.u
	
	def get_stress(self, atoms):
		raise NotImplementedError
	
	def set_atoms(self,atoms):
		pass

	def calculate(self):
		pos = self.atoms.get_positions()
		#If neccesary, transform the positions to new coordinate system
		if self.coord_transform != None:
			pos = self.coord_transform*np.matrix.transpose(pos)
			pos = np.matrix.transpose(pos)

		#Convert ase position matrix to lammps-style position array
		lmp_positions = list(pos.ravel())

		#Convert that lammps-style array into a C object
		lmp_c_positions =\
			(ctypes.c_double*len(lmp_positions))(*lmp_positions)
		self.lmp.put_coords(lmp_c_positions)

		#Run for 0 time to calculate
		self.lmp.command("run 0")
		
		#Extract the forces and energy
		f = np.zeros((len(self.atoms),3))
		f[:,0] = tuple(self.lmp.extract_variable("fx","all",1))
		f[:,1] = tuple(self.lmp.extract_variable("fy","all",1))
		f[:,2] = tuple(self.lmp.extract_variable("fz","all",1))
		u = float(self.lmp.extract_variable("pe",None,0))
		
		self.f = f
		self.u = u
	
	def is_upper_triangular(self,mat):
		"""test if 3x3 matrix is upper triangular"""

		def near0(x):
			"""Test if a float is within .00001 of 0"""
			return abs(x)<.00001

		return near0(mat[1,0]) and near0(mat[2,0]) and near0(mat[2,1])
			

	def convert_cell(self,ase_cell):
		"""
		Convert a parrallelepiped (forming right hand basis)
		to lower triangular matrix LAMMPS can accept. This 
		function transposes cell matrix so the bases are column vectors
		""" 
		cell = np.matrix.transpose(ase_cell)
		if not self.is_upper_triangular(cell):
			#rotate bases into triangular matrix
			tri_mat = np.zeros((3,3))
			A = cell[:,0]
			B = cell[:,1]
			C = cell[:,2]
			tri_mat[0,0] = norm(A)
			Ahat = A/norm(A)
			AxBhat = np.cross(A,B)/norm(np.cross(A,B))
			tri_mat[0,1] = np.dot(B,Ahat)
			tri_mat[1,1] = norm(np.cross(Ahat,B))
			tri_mat[0,2] = np.dot(C,Ahat)
			tri_mat[1,2] = np.dot(C,np.cross(AxBhat,Ahat))
			tri_mat[2,2] = norm(np.dot(C,AxBhat))
			
			#create and save the transformation for coordinates
			volume = np.linalg.det(ase_cell)
			trans = np.array([np.cross(B,C),np.cross(C,A),np.cross(A,B)])
			trans = trans/volume
			self.coord_transform = tri_mat*trans
			
			return tri_mat
		else:
			return cell
	def make_new_lammps(self):
		#Clear the old lammps object
		if(self.lmp != None):
			self.lmp.close()
			self.lmp = None
			self.coord_transform = None
		#create new lammps object
		if self.log_file == None:
			cmd_args = ["-echo", "log","-log","none","-screen","none"]
		else:
			cmd_args = ["-echo", "log","-log",self.log_file,"-screen","none"]
		self.lmp = lammps(cmd_args)

		###Initializing commands

		# Use metal units, Angstrom and eV
		self.lmp.command("units metal")
		self.lmp.command("atom_style atomic")
		self.lmp.command("atom_modify map array sort 0 0")

		#Initialize cell
		cell = self.convert_cell(self.atoms.get_cell())
		xhi = cell[0,0]
		yhi = cell[1,1]
		zhi = cell[2,2]
		xy = cell[0,1]
		xz = cell[0,2]
		yz = cell[1,2]
		cell_cmd = "region cell prism 0 {} 0 {} 0 {} {} {} {} units box"\
			.format(xhi,yhi,zhi,xy,xz,yz)
		self.lmp.command(cell_cmd)

		#The default atom_typeshas atom type in alphabetic order
		# by atomic symbol
		symbols = self.atoms.get_chemical_symbols()
		sorted_symbols = sorted(symbols)
		if self.atom_types == None:
			self.atom_types = {}
			for sym in sorted_symbols:
				if not sym in self.atom_types:
					self.atom_types[sym] =\
						 len(self.atom_types)+1

		#Initiaze box
		n_types = len(self.atom_types)
		types_command = "create_box {} cell".format(n_types)
		self.lmp.command(types_command)
		#Initialize the atoms with their types
		#positions do not matter here
		for sym in symbols:
			cmd = "create_atoms {} single {} {} {} units box"\
				.format(self.atom_types[sym],0.0,0.0,0.0)
			self.lmp.command(cmd)

		#Set masses, even though they don't matter
		self.lmp.command("mass * 1.0")

		#Read in the user parameters
		for cmd in self.parameters:
			self.lmp.command(cmd)


		#Define force&energy variables for extraction
		self.lmp.command("variable fx atom fx")
		self.lmp.command("variable fy atom fy")
		self.lmp.command("variable fz atom fz")
		self.lmp.command("variable pe equal pe")
