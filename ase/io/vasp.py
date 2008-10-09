"""
This module contains the class Atoms which extends ase Atoms
to handle VASP simulations. It will only work with ASE-3.0 and later
versions.
"""

from string import atoi,atof,split
import sys
import os
import re

def nlist(x):
    return range(len(x))

def read_vasp(filename='CONTCAR'):
    """
    Reads unitcell, atom positions (only Direct is implemented at the 
    moment) and constraints from the CONTCAR file and atomtypes 
    from the POTCAR file.
    """
 
    import os
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms, fix_scaled
    from ase.data import chemical_symbols

    dir, filename = os.path.split(filename)

    # Read from file
    file_in = ReadPOSCAR(filename, dir)
    unitcell = file_in.unit_cell()
    poscon = file_in.positions_constraints()
    cart=file_in.cartesian()
    # Read atomic species
    try:
        # Try to read chemical symbols from label in POSCAR/CONTCAR.
        atomtypes = file_in.atom_types()
    except KeyError:
        # If this failes, first try to find OUTCAR, then POTCAR.
        try:
            file_outcar = ReadOUTCAR(dir)
            atomtypes=file_outcar.atom_types()
        except IOError:
            file_potcar = ReadPOTCAR(dir)
            atomtypes=file_potcar.atom_types()
    # Set atomic positions
    atoms=Atoms()

    for ikind in nlist(poscon[0]):
        for nb in nlist(poscon[0][ikind][0:])[::3]:
            atoms.append(Atom(atomtypes[ikind], poscon[0][ikind][nb:nb+3]))
    # Set unitcell
    atoms.set_cell(unitcell, not(cart))
    # Set constraints
    constraints=[]
    indices=[]
    i=0
    for ikind in nlist(poscon[0]):
        for nb in nlist(poscon[0][ikind][0:])[::3]:
            if poscon[1][ikind][0:]:
                tmp = poscon[1][ikind][nb:nb+3]
                if any(tmp) and not all(tmp):
                    constraints.append(fix_scaled(atoms.get_cell(),i,tmp))
                elif all(tmp):
                    indices.append(i)
                i=i+1
    if indices:
        constraints.append(FixAtoms(indices=indices))
    if constraints:
        atoms.set_constraint(constraints)
    return atoms

def write_vasp(filename, atoms, label='', cart=1, dir='', sort=None):
        """Method to write VASP position (POSCAR/CONTCAR) files.

        Writes label, scalefactor, unitcell, # of various kinds of atoms,
        positions in scaled coordinates (Direct) and constraints on file
        POSCAR. Current directory is default and default label is: The 
        atomic species, e.g. 'C N H Cu'.
        """
        
        import numpy as np

        dir, filename = os.path.split(filename)
        
        # Read number of atoms of various kinds
        symbols = atoms.get_atomic_numbers()
        if sort:
            symbols.sort()
        atom_num=[[symbols[0],1]]
        for m in range(1,len(symbols)):
            if symbols[m]==symbols[m-1]:
                atom_num[-1][1]+=1
            else:
                atom_num.append([symbols[m], 1])
        atomlist=[]
        atomsort=[]
        for n in nlist(atom_num):
            atomlist.append(atom_num[n][1])
            atomsort.append(atom_num[n][0])
      
        # Create the label
        if label == '':
            label = element_string(atomsort)

        # Write unitcell in real coordinates and adapt to VASP convention 
        # for unit cell
        ucell=atoms.get_cell()
        
        # Write atom positions in scaled or cartesian coordinates
        if cart:
            coord = atoms.get_positions()
        else:
            coord = atoms.get_scaled_positions()
        con=[]
        constr=atoms.constraints
        if constr:
            con=np.zeros([len(atoms),3])
            for n_constr in nlist(constr):
                b=str(constr[n_constr])
                if b.startswith('FixScaled'):
                    print 1
                    con[constr[n_constr].a]=constr[n_constr].mask
                elif b.startswith('FixAtoms'):
                    con[constr[n_constr].index]=[1,1,1]
        # Write to file
        outfile = WritePOSCAR(filename, dir)
        outfile.label(label)
        outfile.unit_cell(ucell)
        outfile.atom_types(atomlist)
        outfile.coordinates_constraints(coord, con, cart)
        outfile.close()

class ReadPOSCAR:
    """Class to read data from CONTCAR and POSCAR files. 

    All read methods opens and closes POSCAR/CONTCAR. 
    """

    def __init__(self,filename='CONTCAR', dir='./'):

        self._file_= os.path.join(dir,filename)
        if not os.path.isfile(self._file_):
            print ' The file ' +self._file_+ ' does not exist. Please specify a valid filename.'
            sys.exit()
    
    def atom_types(self):
        from ase.data import chemical_symbols
        file=open(self._file_,'r')
        input = file.readline()
        str = split(input)
        atomtypes=[]
        for n in range(len(str)):
            atomtypes.append(str[n])
        for i in nlist(atomtypes):
           issymbol=False
           for symbol in chemical_symbols:
              if symbol==atomtypes[i]:
                 issymbol=True
           if not issymbol:
              raise KeyError
        return atomtypes

    def unit_cell(self):
        """Method that returns unitcell as a list of list of numbers"""
        file=open(self._file_,'r')
        input=file.readlines()
        file.close()
        scale= atof(input[1])
        unitcell=[]
        for i in range(2,5):
            unitcell.append([scale*float(x) for x in input[i].split()])
        return unitcell

    def cartesian(self):
        """Method the check if the coordinates are given in direct or carteisna
        coordinates
        """
        file=open(self._file_,'r')
        for i in range(6):
            skip=file.readline()
        str=file.readline()
        if str[0]=='S':
            str=file.readline()
        if str[0]=='D' or str[0]=='d':
            cartesian=False
        elif str[0]=='C' or str[0]=='c' or str[0]=='K' or str[0]=='k':
            cartesian=True
        return cartesian
        
    def positions_constraints(self):
        """Method that returns a list of list of atom positions and list of
           constraints. Atom positions is a list of the coordinates
           and the constraints is a string.
        """

        file=open(self._file_, 'r')
        input=file.readlines()
        str=split(input[5])
        natoms=[]
        for i in nlist(str):
            natoms.append(atoi(str[i]))
        if(input[6][0] == 'S'):
            nline=8
        else:
            nline=7
        pos=[]
        con=[]
        for ikind in nlist(natoms):
            posatom=[]
            conatom=[]
            for i in range(natoms[ikind]):
                str=split(input[nline])
                vec=[]
                for j in range(3):
                    posatom.append(atof(str[j]))
                if len(str) > 3:
                    for k in range(3):
                        conatom.append(-int(str[3+k]=='T')+1)
                nline+=1
            pos.append(posatom)
            con.append(conatom)
        file.close()
        return [pos,con]

class WritePOSCAR:
    """Class that generates POSCAR files.

    Directory can be specified, default is current directory.
    """

    def __init__(self,filename='POSCAR', dir='./'):
        self._file_ = open(os.path.join(dir,filename),'w')

    def close(self):
        """Method that closes POSCAR file"""
        self._file_.close()        

    def label(self,label):
        """Method that writes label on first line"""
        self._file_.write(label+'\n')

    def unit_cell(self,unitcell):
        """Method that writes scalefactor and unitcell with
          16 decimal places specifie in fix.
        """
        string = '    %.16f\n' % 1.0
        self._file_.write(string)
        for ivec in range(3):
            self._file_.write('    %.16f %.16f %.16f\n' %  tuple(unitcell[ivec]))

    def atom_types(self,atomlist):
        """Method that writes number of various atomtypes
           from 'atomlist'
        """
        str=''
        for ikind in range(len(atomlist)):
           str = str+'  %i' % atomlist[ikind]
        self._file_.write(str+'\n')

    def coordinates_constraints(self, coord, con, cart):
        """Method that writes the atom coordinates from the list 'coord' 
           with 16 decimal places and the constraints from the 
           list 'con'.
        """
        if con != []:
           self._file_.write('Selective dynamics\n')
        if cart :
            self._file_.write('Cartesian\n')
        else:
            self._file_.write('Direct\n')
        for iatom in nlist(coord):
           string='  %.16f %.16f %.16f' %  tuple(coord[iatom])
           if con != []:
               str_ = ''
               for c in range(3):
                   str_+=str(bool(-con[iatom][c]+1))[0]+' '
               string += '  ' + str_
           self._file_.write(string+'\n')

class ReadPOTCAR:
    """Class that read data from POTCAR file.

    Directory can be specified, default is current directory.
    """
    def __init__(self,dir='./'):
        self._file_ = os.path.join(dir, 'POTCAR')

    def atom_types(self):
        """Method that returns list of atomtypes."""
        file=open(self._file_,'r')
        lines=file.readlines()
        file.close()
        atomtypes=[]
        for line in lines:
            if re.search('TITEL',line):
                atomtypes.append(split(split(line.split()[3],'_')[0],'.')[0])
        return atomtypes

class ReadOUTCAR:
    """Class that read data from OUTCAR file. 

    Directory can be specified, default is current directory.
    """
    def __init__(self,dir='./'):
        
        self._file_ = os.path.join(dir, 'OUTCAR')

    def atom_types(self):
        """Method that returns list of atomtypes."""
        file=open(self._file_,'r')
        lines=file.readlines()
        file.close()
        atomtypes=[]
        for line in lines:
            if re.search('TITEL',line):
                atomtypes.append(split(split(line.split()[3],'_')[0],'.')[0])
        return atomtypes

def element_string(list):
    from ase.data import chemical_symbols
    elements = ''
    for n in range(len(list)):
        elements += chemical_symbols[list[n]]+' '
    return elements
