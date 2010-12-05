.. module:: phonons

Phonon calculations
-------------------

This is the first piece of documentation for the new phonon module.  You can
calculate the vibrational normal modes for periodic systems represented by an
:class:`~ase.atoms.Atoms` object in the harmonic approximation using the
:class:`~ase.phonons.Phonons` class.


Example
-------

Simple example showing how to calculate the phonon dispersion for bulk aluminum
using a 7x7x7 supercell within effective medium theory::

  from ase.structure import bulk
  from ase.calculators.emt import EMT
  from ase.dft.kpoints import ibz_points, get_bandpath
  from ase.phonons import Phonons
  
  # Setup crystal and EMT calculator
  atoms = bulk('Al', 'fcc', a=4.05)
  calc = EMT()
  
  # Phonon calculator
  N = 7
  ph = Phonons(atoms, calc, supercell=(N, N, N))
  # ph.run()
  
  # Read forces and assemble the dynamical matrix
  ph.read(acoustic=True)
  
  # High-symmetry points in the Brillouin zone
  points = ibz_points['fcc']
  G = points['Gamma']
  X = points['X']
  W = points['W']
  K = points['K']
  L = points['L']
  
  point_names = ['$\Gamma$', 'K', 'X', '$\Gamma$', 'L', 'X', 'W', 'L']    
  path = [G, K, X, G, L, X, W, L]
  
  path_kc, q, Q = get_bandpath(path, atoms.cell, 100)
  omega_kn = 1000 * ph.band_structure(path_kc)

  import pylab as plt
  for n in range(len(omega_kn[0])):
      omega_n = omega_kn[:, n]
      plt.plot(q, omega_n, 'k-', lw=2)

  plt.xticks(Q, point_names, fontsize=18)
  plt.yticks(fontsize=18)
  plt.xlim(q[0], q[-1])
  plt.ylabel("Frequency ($\mathrm{meV}$)", fontsize=22)
  plt.grid('on')
  plt.show()
  
  # Write modes for specific q-vector to trajectory files  
  ph.write_modes(K, branches=[0,1,2], repeat=(5, 5, 5))