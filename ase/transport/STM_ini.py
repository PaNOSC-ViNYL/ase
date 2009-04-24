"""This program is to carry out STM simulation by using Green's functions"""

from gpaw.grid_descriptor import GridDescriptor
from ase.units import Bohr
from gpaw import *
from ase import *
import numpy as np
import os

class STM:
  """ The first step is to read the outputs of isolated tip and surface calculations """
  def __init__(self, surfacecalc, tipcalc, tipapex, surface_extention, tip_cutoff, vacuum,
               surface_z):
    self.scalc = surfacecalc  # calculator of surface calculation
    self.tcalc = tipcalc      # calculator of tip calculation
    self.tapex = tipapex      # 3D coordinate of tip-apex in tip unit cell
    self.sradius = surface_extention  # surface extension in pseduo surface cell (Angstrom)
    self.tcutoff = tip_cutoff         # tip cut off for pseduo tip cell (Angstrom)
    self.vacuum = vacuum              # tip-surface separation (Angstrom)
    self.spotential = surfacecalc.hamiltonian.vt_sG[0]  # get surface effective potential on surface grid
    self.tpotential = tipcalc.hamiltonian.vt_sG[0]      # get tip effective potential on tip grid
    self.satoms = surfacecalc.get_atoms() # get surface atoms
    self.tatoms = surfacecalc.get_atoms() # get tip atoms
    self.sgd = self.scalc.gd # grid points information of the surface unit cell
    self.tgd = self.tcalc.gd # grid points information of the tip unit cell
    self.surface_z = surface_z
    if self.tcutoff[0] > self.sradius:
      raise ValueError('The surface extension is too small comparing to the tip cutoff in x direction')
   
    if self.tcutoff[1] > self.sradius:
      raise ValueError('The surface extension is too small comparing to the tip cutoff in y direction')

    if self.tcutoff[1] > self.tatoms.get_cell()[1][1]:
      raise ValueError('The tip cutoff in y direction is too large for the original tip cell')

    if self.tcutoff[0] > self.tatoms.get_cell()[0][0]:
      raise ValueError('The tip cutoff in x direction is too large for the original tip cell')

    if self.tcutoff[2] > self.tatoms.get_cell()[2][2]:
      raise ValueError('The tip cutoff in z direction is too large for the original tip cell')

    if self.surface_z+self.tcutoff[2]+5.0 > self.satoms.get_cell()[2][2]: 
      raise ValueError('The original surface cell does not have enough vacuum to embed the tip')
    
    print os.system('date') 

    """ extend the surface cell so that the periodic boundary condition is included in a pseduo (large) surface cell """
    pseduo_surface_extension,pseduo_surface_gd = pseduo_surf(cell=self.satoms.get_cell(),h_c=self.sgd.h_c, 
                                                             N_c=self.sgd.N_c, radius=self.sradius)
    """ embed the effective potential of original surf cell to the pseduo surf cell """
    pseduo_surface_potential = pseduo_surf_potential(gd=pseduo_surface_gd,ext=pseduo_surface_extension,
                                                     potential = self.spotential)
    """ construct a pseduo tip cell """
    pseduo_tip_gd = pseduo_tgd(tgd=self.tgd, size=self.tcutoff)

    """ embed the effective potential of tip into the pseduo tip cell """
    pseduo_tip_potential = pseduo_tpotential(potential=self.tpotential, pseduo_tgd=pseduo_tip_gd,
                                             tapex=self.tapex, tcutoff=self.tcutoff)
    """ embed the functions of the original surf cell to the pseduo surf cell """
    pseduo_surface_wf = pseduo_surf_wannier_functions(calc=self.scalc, lgd=pseduo_surface_gd, 
                                                      ext=pseduo_surface_extension)

    """ embed the functions of the original tip cell to to the pseduo tip cell """
    pseduo_tip_wf = pseduo_tip_functions(calc=self.tcalc, pseduo_tgd=pseduo_tip_gd,
                             tapex=self.tapex, tcutoff=self.tcutoff)


#    print "start", os.system('date')
    """ calculate the second derivative of surface and tip functions """
    kinetic_surface, kinetic_tip = kinetic(sgrids=pseduo_surface_gd, sfunction=pseduo_surface_wf,
                                           tgrids=pseduo_tip_gd, tfunction=pseduo_tip_wf,
                                           scalc=self.scalc, tcalc=self.tcalc)

    """ set up STM and calculate the matrix element as the tip moves above the surface """

    setup_v, setup_phi = setup(distance=self.vacuum,
                                extension = pseduo_surface_extension,surface_z=self.surface_z, sgd=pseduo_surface_gd,
                                tcutoff=self.tcutoff)

    scan = array([0,0,0]) # Here you indicate the relative position of the tip apex to the (0,0,0) point of surface cell

    """ calculate the combined effective potential """
    V = couple_v(scan=scan,vtip=pseduo_tip_potential,vsurface=pseduo_surface_potential,
                 setup_v=setup_v) 

    print "start", os.system('date')
    """ calculate the combined kinetic energy """
    T = couple_t(scan=scan,ktip=kinetic_tip, ksurface=kinetic_surface,
                 setup=setup_v)
    print "finish",os.system('date')

    """ calculate the combined functions """
    Phi = couple_phi(scan=scan,swf=pseduo_surface_wf,twf=pseduo_tip_wf,setup=setup_phi)


    print "V",shape(V)
    print "T",shape(T)
    print "Phi",shape(Phi)

#    """ calculate the matrix element V_ij"""
#    H_V = V_element(V=V,Phi=Phi,gd=pseduo_surface_gd)
#    print "H_V",shape(H_V)
#  
#    """ calculate the matrix element T_ij"""
#    H_T = T_element(T=T,Phi=Phi,gd=pseduo_surface_gd)
#    print "H_T",shape(H_T)
#
#    """ calculate the overlap of the basis functions """
#    S = S_element(Phi,gd=pseduo_surface_gd)

##########################################################################################################################
def pseduo_surf(cell,h_c,N_c,radius):
  # cell: the original surface cell, unit is 'Angstrom'
  # h_c: the separation of grid points in the original surface cell, unit is 'Bohr'
  if radius > cell[1][1] or radius > cell[0][0]:
    raise ValueError('The surface extension you select to too large for the original surface cell')
  temp = np.round(radius/(h_c*Bohr)) # the number of grid points included in the radius in 3D
  extension = [int(temp[0]),int(temp[1]),0]   # convert data to integer type, and do not extend in z direction 
  extension_x = extension[0]*h_c[0]*Bohr  # extension in +/- x direction in Angstrom
  pseduo_x = cell[0,0]+2*extension_x # the new cell length in x direction in Angstrom
  extension_y = extension[1]*h_c[1]*Bohr  # extension in +/- y direction in Angstrom
  pseduo_y = cell[1,1]+2*extension_y # the new cell length in y direction in Angstrom
  pseduo_z = cell[2,2]  # the pseduo cell length in z direction in Angstrom
  pseduo_Nc = N_c + extension + extension # the No. of grid points in pseduo cell in 3D
  newcell = array([pseduo_x,pseduo_y,pseduo_z])/Bohr # the length of unit cell in x,y,z direction in Bohr
  """set up a pseduo surf cell"""
  pseduo_surf_gd = GridDescriptor(pseduo_Nc,newcell,pbc_c=False) # define the information of grid in pseduo cell
  print "Surface extention ok!"
  return extension,pseduo_surf_gd

###########################################################################################################################
def pseduo_surf_potential(gd,ext,potential):
  # gd: grid description of the pseduo surface cell
  # ext: the number of grid points extending in +/- x,y direction
  # potential: effective potential on original surface grid
  pseduo_spotential = empty((gd.N_c),dtype=float)
  ext_x = ext[0] # get the number of grid points shifted in x direction
  ext_y = ext[1] # get the number of grid points shifted in y direction
  ori_x = gd.N_c[0]-2*ext_x  # get the number of grid points in the original surface cell in x direction
  ori_y = gd.N_c[1]-2*ext_y  # get the number of grid points in the original surface cell in y direction
  """ write the calculated effective potential on the grids of the pseduo surface cell """
  pseduo_spotential[ext_x:(ext_x + ori_x), ext_y:(ext_y + ori_y), :] = potential
  pseduo_spotential[0:ext_x, ext_y:(ext_y+ori_y), :] = potential[(ori_x-ext_x):ori_x, :, :]
  pseduo_spotential[(ext_x+ori_x):(ext_x+ori_x+ext_x), ext_y:(ext_y+ori_y), :] = potential[0:ext_x, :, :]
  temp=pseduo_spotential.copy()
  pseduo_spotential[:,0:ext_y,:]=temp[:,ori_y:ori_y+ext_y,:]
  pseduo_spotential[:,ext_y+ori_y:,:]=temp[:,ext_y:ext_y+ext_y,:]
  spotential = pseduo_spotential[1:gd.N_c[0]-1,1:gd.N_c[1]-1,:]
  print "Surface effective potential embeding ok!" 
  return spotential
##########################################################################################################################
def pseduo_tgd(tgd, size):
  #tgd: the original tip grid description
  #size: the cutoff of tip (in Angstrom)
  h_c = array(tgd.h_c)   # get the separtion of grid points in surface cell in Bohr
  length = array(size)/Bohr   # get the pseduo tip cell box length in Bohr
  temp = around(length/h_c)+4   # get the number of grid points in 3D in the pseduo tip cell
  N_c = array(temp.copy(),dtype=integer)
  length = h_c*N_c  # the length of the pseduo tip cell in Bohr
  pseduo_tip_gd = GridDescriptor(N_c=N_c,cell_cv=length, pbc_c=False) # define the grids in large pseduo tip cell
  print "Pseduo tip gird ok!"
  return pseduo_tip_gd

##########################################################################################################################
def pseduo_tpotential(potential, pseduo_tgd, tapex, tcutoff):
  # potenital: effective pontential of STM tip in original grids
  # pseduo_tgd: grid description of the coarse grid
  # tapex: the apex of STM tip in origianl tip unit cell in Angstrom
  # tcutoff: the cutoff of tip in Angstrom
  
  tpotential = empty((pseduo_tgd.N_c),dtype=float) # set empty potential for all grid points
  tpotential*= 0.0
  """ the 3D coordinates of the starting grid points in the original tip grids is"""
  gd_start_c  = array(tapex)-array(tcutoff)*0.5
  if gd_start_c[0]<0. or gd_start_c[1]<0.:
    raise ValueError('The tip cutoff is too large for the original tip cell')
  """ the starting grid point is """
  temp = around(gd_start_c/pseduo_tgd.h_c)
  gd_start_g = array(temp.copy(),dtype=integer)

  sx = gd_start_g[0]
  sy = gd_start_g[1]
  sz = gd_start_g[2]
   
  nx = (pseduo_tgd.N_c-4+gd_start_g)[0]
  ny = (pseduo_tgd.N_c-4+gd_start_g)[1]
  nz = (pseduo_tgd.N_c-4+gd_start_g)[2]
  tpotential[2:-2,2:-2,2:-2] = potential[sx:nx,sy:ny,sz:nz]
  tip_potential = tpotential[1:pseduo_tgd.N_c[0]-1,1:pseduo_tgd.N_c[1]-1,1:pseduo_tgd.N_c[2]-1]
  print "Tip effective potential embeding ok!"

  return tip_potential
########################################################################################################################## 
def pseduo_surf_wannier_functions(calc, lgd, ext):
  #calc: the calculator of surface
  #lgd: the extended surface grid description
  #ext: the number of grid points in surface extension in 3D
  
  """ get surface functions from surface calculation """
  wfs = calc.wfs
  bfs = calc.wfs.basis_functions
  gd = calc.wfs.gd
  nao = wfs.setups.nao # get the total number of atomic orbitals
  calc.initialize_positions()
  C_MM = np.identity(nao)
  phi_MG = gd.zeros(nao)
  bfs.lcao_to_grid(C_MM, phi_MG, -1)

  """ get information of surface extension """
  ext_x = ext[0]  # get the number of grid points shifted in x direction
  ext_y = ext[1]  # get the number of grid points shifted in y direction
  ori_x = lgd.N_c[0]-2*ext_x  # get the number of grid points in the original surface cell in x direction
  ori_y = lgd.N_c[1]-2*ext_y  # get the number of grid points in the original surface cell in y direction
  """ embed each functions to the extend surface cell """   
  pseduo_swf = []
  for i in range(nao):
    temp_1 = empty((lgd.N_c),dtype=wfs.dtype)
    temp_1[ext_x:ext_x + ori_x, ext_y:ext_y+ori_y, :] = phi_MG[i] 
    temp_1[0:ext_x, ext_y:(ext_y+ori_y), :] = phi_MG[i][(ori_x-ext_x):ori_x, :, :]
    temp_1[(ext_x+ori_x):(ext_x+ori_x+ext_x), ext_y:(ext_y+ori_y), :] = phi_MG[i][0:ext_x, :, :]
    temp_2 = temp_1.copy()
    temp_1[:,0:ext_y,:]=temp_2[:,ori_y:ori_y+ext_y,:]
    temp_1[:,ext_y+ori_y:,:]=temp_2[:,ext_y:ext_y+ext_y,:]    
    pseduo_swf.append(temp_1.copy())
    
  print "Surface function embeding ok!"
  return pseduo_swf

########################################################################################################################## 
def pseduo_tip_functions(calc, pseduo_tgd, tapex, tcutoff):
  # calc: calculator of tip calculation
  # pseduo_tgd: grid description of large pseduo tip cell
  # tapex: tip coordinates in 3D in the original tip cell (in Angstrom)
  # tcutoff: tip cutoff for STM simulations (in Angstrom) 
  """ get tip functions from tip calculation """
  calc.initialize_positions()
  wfs = calc.wfs                   
  bfs = calc.wfs.basis_functions
  gd = calc.wfs.gd                   # get the grid description from original tip calculation
  nao = wfs.setups.nao # get the total number of orbitals
 
  atoms = calc.get_atoms()           # get the atoms from original tip calculation
  tipapex=array(tapex)             
  tip_cutoff=array(tcutoff)
  tip_from = tipapex - 0.5*tip_cutoff      # the lower bound of the pseduo tip volume in the the original tip cell
  tip_to = tipapex + 0.5*tip_cutoff        # the upper bound of the pseduo tip volume in the the original tip cell
  norbital = 0
  selection=[]
  global_index = []
  for atom in atoms:
    pos = atom.get_position()
    if tip_from[0] < pos[0] < tip_to[0] and tip_from[1] < pos[1] < tip_to[1] and tip_from[2] < pos[2] < tip_to[2]:
      selection.append(atom.index)              # select the atoms included in the pseduo tip 
      norbital += wfs.setups[atom.index].niAO   # add the number functions for each atom which is in the pseduo tip
      M1 = bfs.M_a[atom.index]
      M2 = wfs.setups[atom.index].niAO + M1
      for i in range(M1,M2):
        global_index.append(i)   # get the global index of the selected orbitals and put them in a list called global_index

  C_MM = np.identity(nao)
  phi_MG = gd.zeros(nao)
  bfs.lcao_to_grid(C_MM, phi_MG, -1)
  temp = np.round(tip_from/gd.h_c/Bohr)
  grid_from = array(temp.copy(),dtype=integer)
  nx = pseduo_tgd.N_c[0]
  ny = pseduo_tgd.N_c[1]
  nz = pseduo_tgd.N_c[2]
  pseduo_twf = []

  for j in global_index:
    temp_1 = empty((pseduo_tgd.N_c),dtype=wfs.dtype) # set empty values for all grid points
    temp_1 = phi_MG[j][grid_from[0]:grid_from[0]+nx,grid_from[1]:grid_from[1]+ny,grid_from[2]:grid_from[2]+nz]
    pseduo_twf.append(temp_1.copy())
  print "Tip function embeding ok!"
  return pseduo_twf


###############################################################################################################
def kinetic(sgrids,sfunction, tgrids, tfunction, scalc, tcalc):
  """ This function is to calculate the kinetic energy of a wavefunction on grid """
  # sgrids: the grid description of the pseduo surface cell  
  # fuction is the function on the extend surface cell
  # tgrids: the grid description of the pseduo tip cell
  # tfunction: the function on the pseduo tip cell
  # scalc and tcalc : calculators of original surface and tip calculations
  h_x = sgrids.h_c[0]
  h_y = sgrids.h_c[1]
  h_z = sgrids.h_c[2]
  ksurface = []
  for m in range(len(sfunction)):
    temp = empty((sgrids.N_c),dtype=sfunction[0].dtype)
    temp *= 0.0
    f = sfunction[m].copy()
    f_m = f.copy()
    f_p = f.copy()
    """Caclculate k_x"""
    f_m[1:sgrids.N_c[0],:,:] = f[0:sgrids.N_c[0]-1,:,:]
    f_p[0:sgrids.N_c[0]-1,:,:] = f[1:sgrids.N_c[0],:,:]
    k_x = (f_p+f_m-f-f)/h_x/h_x/2.0
    """Calculate k_y"""
    f=sfunction[m].copy()
    f_m[:,1:sgrids.N_c[1],:]=f[:,0:sgrids.N_c[1]-1,:]
    f_p[:,0:sgrids.N_c[1]-1,:]=f[:,1:sgrids.N_c[1],:]
    k_y = (f_p+f_m-f-f)/h_y/h_y/2.0
    """Calculate k_z"""
    f=sfunction[m].copy()
    f_m[:,:,1:] = f[:,:,0:sgrids.N_c[2]-1]
    f_m[:,:,0] = f[:,:,-1]
    f_p[:,:,0:sgrids.N_c[2]-1] = f[:,:,1:sgrids.N_c[2]]
    f_p[:,:,-1] = f[:,:,0]
    k_z = (f_p+f_m-f-f)/h_z/h_z/2.0
    """Calculate k """
    k = k_x + k_y + k_z   #Note that valid region is k[1:sgrids.N_c[0]-1,1:sgrids.N_c[1]-1,1:sgrids.N_c[2]-1]
    k_real = k[1:sgrids.N_c[0]-1,1:sgrids.N_c[1]-1,:]
    ksurface.append(k_real.copy())
  print "kinetic for surface functions done" 

  t_x = tgrids.h_c[0]
  t_y = tgrids.h_c[1]
  t_z = tgrids.h_c[2]
  ktip = []
  for n in range(len(tfunction)):
    temp = empty((tgrids.N_c),dtype=tfunction[0].dtype)
    temp*= 0.0
    f = tfunction[n].copy()
    f_m = f.copy()
    f_p = f.copy()
    """Caclculate k_x"""
    f_m[1:tgrids.N_c[0],:,:] = f[0:tgrids.N_c[0]-1,:,:]
    f_p[0:tgrids.N_c[0]-1,:,:] = f[1:tgrids.N_c[0],:,:]
    k_x = (f_p+f_m-f-f)/t_x/t_x/2.0
    """Calculate k_y"""
    f = tfunction[n].copy()
    f_m[:,1:tgrids.N_c[1],:]=f[:,0:tgrids.N_c[1]-1,:]
    f_p[:,0:tgrids.N_c[1]-1,:]=f[:,1:tgrids.N_c[1],:]
    k_y = (f_p+f_m-f-f)/t_y/t_y/2.0
    """Calculate k_z"""
    f=tfunction[n].copy()
    f_m[:,:,1:tgrids.N_c[2]] = f[:,:,0:tgrids.N_c[2]-1]
    f_p[:,:,0:tgrids.N_c[2]-1] = f[:,:,1:tgrids.N_c[2]]
    k_z = (f_p+f_m-f-f)/t_z/t_z/2.0
    """Calculate k """
    k = k_x + k_y + k_z   
    k_real = k[1:tgrids.N_c[0]-1,1:tgrids.N_c[1]-1,1:tgrids.N_c[2]-1]
    ktip.append(k_real.copy())
  print "kinetic for tip functions done"
  return ksurface, ktip
############################################################################################################
def setup(distance, extension, surface_z, sgd, tcutoff):
  # distance: tip-surface distance in Angstrom
  # extension: No. of grid points extened in orignal surface cell when we construct the pseduo surface cell
  # surface_z: the z coordinate of the top-most surface atom
  # sgd: grid description of the pseduo surface celli
  # tcutoff : size of the STM tip in Angstrom
  """Find out the grid index of (0,0,vacuum) point in the pseduo surface cell """
  gox = extension[0]-1
  goy = extension[0]-1 
  oz = surface_z + distance - (0.5 * tcutoff[2])
  temp = np.round(oz / Bohr /sgd.h_c[2])
  goz = int(temp)
  setup_v = array([gox,goy,goz])
  setup_phi = array([gox+1,goy+1,goz])
  return setup_v, setup_phi
############################################################################################################
def couple_v(scan, vtip, vsurface, setup_v):
  ox = scan[0]+setup_v[0]
  oy = scan[1]+setup_v[1]
  oz = scan[2]+setup_v[2]
  nx,ny,nz = shape(vtip)
  
  vcouple = empty((shape(vsurface)),dtype=vsurface.dtype) 
  vcouple*= 0.0 
  vcouple[ox:ox+nx,oy:oy+ny,oz:oz+nz] = vtip.copy()
  vcouple += vsurface.copy()

  return vcouple

##############################################################################################################
def couple_t(scan, ktip, ksurface, setup):
  ox = scan[0]+setup[0]
  oy = scan[1]+setup[1]
  oz = scan[2]+setup[2]
  nx = shape(ktip)[1]
  ny = shape(ktip)[2]
  nz = shape(ktip)[3]
  tcouple = ksurface[:]
  
  for i in range(shape(ktip)[0]):
    temp = empty((shape(ksurface)[1:4]),dtype=ksurface[0].dtype)
    temp *= 0.0
    temp[ox:ox+nx,oy:oy+ny,oz:oz+nz] = ktip[i].copy() 
    tcouple.append(temp.copy())
 
  return tcouple

###################################################################################
def couple_phi(scan, swf, twf, setup):
  ox = scan[0]+setup[0]
  oy = scan[1]+setup[1]
  oz = scan[2]+setup[2]
  nx = shape(twf)[1]-2
  ny = shape(twf)[2]-2
  nz = shape(twf)[3]-2
  
  phi_couple = []
  sx = shape(swf)[1]
  sy = shape(swf)[2]
  sz = shape(swf)[3]
  tx = shape(twf)[1]
  ty = shape(twf)[2]
  tz = shape(twf)[3]

  temp = empty((sx-2,sy-2,sz),dtype=swf[0].dtype)

  #add the the surface functions to phi_couple
  
  for i in range(shape(swf)[0]):
    temp = empty((sx-2,sy-2,sz),dtype=swf[0].dtype)
    temp = swf[i].copy()[1:sx-1,1:sy-1,:]
    phi_couple.append(temp.copy())     
  
  # add the tip functions to the phi_couple
  for j in range(shape(twf)[0]):
    temp = empty((sx-2,sy-2,sz),dtype=swf[0].dtype)
    temp *= 0.0
    temp[ox:ox+nx,oy:oy+ny,oz:oz+nz] = twf[j].copy()[1:tx-1,1:ty-1,1:tz-1]
    phi_couple.append(temp.copy())
  return phi_couple
###################################################################################
def V_element(V,Phi,gd):
  n = shape(Phi)[0]
  element_V = empty((n,n),dtype=Phi[0].dtype)
  print "#",shape(Phi[1].conj())
  print "#",shape(V)
  print "#",shape(Phi[1])

  ix,iy,iz = gd.h_c*Bohr
  print "ix,iy,iz",ix,iy,iz
  for i in range(n):
    for j in range(n):
      temp = Phi[i].conj()*(V*Phi[j])
      test = temp.sum()*ix*iy*iz
      element_V[i,j] = test
  return element_V
####################################################################################
def T_element(T,Phi,gd):
  n = shape(Phi)[0]
  element_T = empty((n,n),dtype=Phi[0].dtype)
  ix,iy,iz = gd.h_c*Bohr
  for i in range(n):
    for j in range(n):
      temp = Phi[i].conj()*T[j]
      test = temp.sum()*ix*iy*iz
      element_T[i,j] = test
  return element_T
#####################################################################################
def S_element(Phi,gd):
  n = shape(Phi)[0]
  element_S = empty((n,n),dtype=Phi[0].dtype)
  ix,iy,iz = gd.h_c*Bohr
  for i in range(n):
    for j in range(n):
      temp = Phi[i].conj()*Phi[j]
      test = temp.sum()*ix*iy*iz
      element_S[i,j] = test
  return element_S


