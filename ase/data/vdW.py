""" Van der Waals radii in [A] taken from
http://www.webelements.com/periodicity/van_der_waals_radius/
and the references given there:
1. A. Bondi, J. Phys. Chem., 1964, 68, 441.
2. L. Pauling, The Nature of the Chemical Bond, 
   Cornell University Press, USA, 1945.
3. J.E. Huheey, E.A. Keiter, and R.L. Keiter in Inorganic Chemistry : 
   Principles of Structure and Reactivity, 4th edition, HarperCollins, 
   New York, USA, 1993.W.W. Porterfield in Inorganic chemistry, 
   a unified approach, Addison Wesley Publishing Co., 
   Reading Massachusetts, USA, 1984.
4. A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data, 
   Macmillan, London, UK, 1992."""

import numpy as np

vdW_radii = np.array([
 np.nan, # X
 1.20, # H
 1.40, # He
 1.82, # Li
 np.nan, # Be
 np.nan, # B
 1.70, # C
 1.55, # N
 1.52, # O
 1.47, # F
 1.54, # Ne
 2.27, # Na
 1.73, # Mg
 np.nan, # Al
 2.10, # Si
 1.80, # P
 1.80, # S
 1.75, # Cl
 1.88, # Ar
 2.75, # K
 np.nan, # Ca
 np.nan, # Sc
 np.nan, # Ti
 np.nan, # V
 np.nan, # Cr
 np.nan, # Mn
 np.nan, # Fe
 np.nan, # Co
 1.63, # Ni
 1.40, # Cu
 1.39, # Zn
 1.87, # Ga
 np.nan, # Ge
 1.85, # As
 1.90, # Se
 1.85, # Br
 2.02, # Kr
 np.nan, # Rb
 np.nan, # Sr
 np.nan, # Y
 np.nan, # Zr
 np.nan, # Nb
 np.nan, # Mo
 np.nan, # Tc
 np.nan, # Ru
 np.nan, # Rh
 1.63, # Pd
 1.72, # Ag
 1.58, # Cd
 1.93, # In
 2.17, # Sn
 np.nan, # Sb
 2.06, # Te
 1.98, # I
 2.16, # Xe
 np.nan, # Cs
 np.nan, # Ba
 np.nan, # La
 np.nan, # Ce
 np.nan, # Pr
 np.nan, # Nd
 np.nan, # Pm
 np.nan, # Sm
 np.nan, # Eu
 np.nan, # Gd
 np.nan, # Tb
 np.nan, # Dy
 np.nan, # Ho
 np.nan, # Er
 np.nan, # Tm
 np.nan, # Yb
 np.nan, # Lu
 np.nan, # Hf
 np.nan, # Ta
 np.nan, # W
 np.nan, # Re
 np.nan, # Os
 np.nan, # Ir
 1.75, # Pt
 1.66, # Au
 1.55, # Hg
 1.96, # Tl
 2.02, # Pb
 np.nan, # Bi
 np.nan, # Po
 np.nan, # At
 np.nan, # Rn
 np.nan, # Fr
 np.nan, # Ra
 np.nan, # Ac
 np.nan, # Th
 np.nan, # Pa
 1.86, # U
 np.nan, # Np
 np.nan, # Pu
 np.nan, # Am
 np.nan, # Cm
 np.nan, # Bk
 np.nan, # Cf
 np.nan, # Es
 np.nan, # Fm
 np.nan, # Md
 np.nan, # No
 np.nan]) # Lr
