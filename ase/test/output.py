from ase import *
a = Atoms([Atom('H')], cell=(2, 3, 4))
write('H.eps', a)
write('H.nc', a)

          
