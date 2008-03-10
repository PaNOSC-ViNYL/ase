"""ASE interface to SIESTA."""

import os

import numpy as npy

from ase.data import chemical_symbols


class Siesta:
    """Class for doing SIESTA calculations.

    """
    def __init__(self, label='siesta', xc='LDA', kpts=None, nbands=None,
                 width=None, meshcutoff=None, charge=None,
                 basis=None):
        """Construct SIESTA-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.fdf, label.txt, ...).
            Default is 'siesta'.
        xc: str
            Exchange-correlation functional.  Must be one of LDA, PBE,
            revPBE, RPBE.
        kpts: list of three int
            ...

        Examples
        ========
        Use default values:
        
        >>> h = Atoms('H', calculator=Siesta())
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()
        
        """
        self.label = label
        self.xc = xc
        self.kpts = kpts
        self.nbands = nbands
        self.width = width
        self.meshcutoff = meshcutoff
        self.charge = charge
        self.basis = basis
    
        self.etotal = None
        self.efree = None
        self.fdf = {}

    def update(self, atoms):
        if (self.etotal is None or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal + self.efree) / 2

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def initialize(self, atoms):
        pass
    
    def calculate(self, atoms):
        self.write_fdf(atoms)

        siesta = os.environ.get('SIESTA_SCRIPT', 'siesta_2.0')
        error = os.system('%s < %s.fdf > %s.txt' %
                          (siesta, self.label, self.label))
        if error != 0:
            raise RuntimeError('Siesta exited with error code: %d.' % error)
        
        self.read()

    def fdf(self, key, value):
        """Set FDF parameter."""
        self.fdf[key] = value
                
    def write_fdf(self, atoms):
        """Write input parameters to fdf-file."""
        self.positions = atoms.get_positions().copy()
        self.numbers = atoms.get_atomic_numbers().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()
        
        fdf = {
            'SystemLabel': self.label,
            'AtomicCoordinatesFormat': 'Ang',
            'LatticeConstant': 1.0,
            'NumberOfAtoms': len(atoms),
            'MeshCutoff': self.meshcutoff,
            'NetCharge': self.charge,
            'ElectronicTemperature': self.width,
            'NumberOfEigenStates': self.nbands,
            #'kgrid_cutoff': (Bohr, ''),??
            }
        
        if self.xc != 'LDA':
            fdf['xc.functional'] = 'GGA'
            fdf['xc.authors'] = self.xc

        try:
            magmoms = atoms.get_magnetic_moments()
        except KeyError:
            pass
        else:
            if magmoms is not None and magmoms.any():
                fdf['SpinPolarized'] = True
        
        species = {}
        n = 0
        for Z in self.numbers:
            if Z not in species:
                n += 1
                species[Z] = n
        fdf['Number_of_species'] = n

        #'user-basis': bool,
        #'PAO.BasisSize': ('standard',),
        #'PAO.BasisType': ('split',),
        #'PAO.EnergyShift': (Ry, 'eV'),
        #'PAO.SplitNorm': float,

        fdf.update(self.fdf)

        f = open(self.label + '.fdf', 'w')
        for key, value in fdf.items():
            if value is None:
                continue
            unit = keys_with_units.get(fdfify(key))
            if unit is None:
                f.write('%s %s\n' % (key, value))
            else:
                f.write('%s %f %s\n' % (key, value, unit))

        f.write('%block LatticeVectors\n')
        for v in self.cell:
            f.write('%.14f %.14f %.14f\n' %  tuple(v))
        f.write('%endblock\n')

        f.write('%block Chemical_Species_label\n')
        species = [(n, Z) for Z, n in species.items()]
        species.sort()
        for n, Z in species:
            f.write('%d %s %s\n' % (n, Z, chemical_symbols[Z]))
        f.write('%endblock\n')

        f.write('%block AtomicCoordinatesAndAtomicSpecies\n')
        for pos, Z in zip(self.positions, self.numbers):
            f.write('%.14f %.14f %.14f' %  tuple(pos))
            f.write(' %d\n' % Z)
        f.write('%endblock\n')
        f.close()
        
    def read(self):
        """Read results from SIESTA's text-output file."""
        lines = iter(open(self.label + '.txt', 'r').readlines())

        # Get the number of grid points used:
        for line in lines:
            if line.startswith('InitMesh: MESH ='):
                self.grid = [int(word) for word in line.split()[3:8:2]]
                break

        # Stress:
        for line in lines:
            if line.startswith('siesta: Stress tensor (total) (eV/Ang**3):'):
                self.stress = npy.empty((3, 3))
                for i in range(3):
                    self.stress[i] = [float(word)
                                      for word in lines.next().split()]
                break

        # Energy:
        for line in lines:
            if line.startswith('siesta: Etot    ='):
                self.etotal = float(line.split()[-1])
                self.efree = float(lines.next().split()[-1])
                break

        # Forces:
        for line in lines:
            if line.startswith('siesta: Atomic forces (eV/Ang):'):
                self.forces = npy.empty((len(self.numbers), 3))
                for force in self.forces:
                    force[:] = [float(word)
                                for word in lines.next().split()[-3:]]
                break
        

def fdfify(key):
    return key.lower().replace('_', '').replace('.', '').replace('-', '')


keys_with_units = {
    'paoenergyshift': 'eV',
    'zmunitslength': 'Bohr',
    'zmunitsangle': 'rad',
    'zmforcetollength': 'eV/Ang',
    'zmforcetolangle': 'eV/rad',
    'zmmaxdispllength': 'Ang',
    'zmmaxdisplangle': 'rad',
    'meshcutoff': 'eV',
    'dmenergytolerance': 'eV',
    'electronictemperature': 'eV',
    'oneta': 'eV',
    'onetaalpha': 'eV',
    'onetabeta': 'eV',
    'onrclwf': 'Ang',
    'onchemicalpotentialrc': 'Ang',
    'onchemicalpotentialtemperature': 'eV',
    'mdmaxcgdispl': 'Ang',
    'mdmaxforcetol': 'eV/Ang',
    'mdmaxstresstol': 'eV/Ang**3',
    'mdlengthtimestep': 'fs',
    #'mdinitialtemperature': (K, ,????
    #'mdtargettemperature': (K, '',????
    'mdtargetpressure': 'eV/Ang**3',
    'mdnosemass': 'eV*fs**2',
    'mdparrinellorahmanmass': 'eV*fs**2',
    'mdtaurelax': 'fs',
    'mdbulkmodulus': 'eV/Ang**3',
    'mdfcdispl': 'Ang',
    'warningminimumatomicdistance': 'Ang',
    'rcspatial': 'Ang',
    'kgridcutoff': 'Ang',
    'latticeconstant': 'Ang'}
