# this has been modified
def _write_atomic_coordinates(f, atoms):
    """Write atomic coordinates.
    
    Parameters:
    - f:     An open file object.
    - atoms: An atoms object.
    """
    #species, species_numbers = self.species(atoms)
    f.write('\n')
    f.write('#\n')
    f.write('GEOMETRY CARTESIAN ANGSTROMS\n')
    #f.write('%block AtomicCoordinatesAndAtomicSpecies\n')
    for i in range(len(atoms)):
        xyz = atoms.get_positions()[i]
        chem_symbol = atoms.get_chemical_symbols()[i]
        chem_symbol += str(i+1)
        
        # if tag is set to 1 then we have a ghost atom, set nuclear charge to 0
        if(atoms.get_tags()[i] == 1):
            nuc_charge = '0' 
        else:
            nuc_charge = str(atoms.get_atomic_numbers()[i])
            
        mass = atoms.get_masses()[i]
        #line = '{0:10s}  {1:.9f}  {2:.9f}  {3:.9f}  {4:5s}  {5:.9f}'.format(chem_symbol,
        #                                                                    xyz[0], xyz[1],xyz[2],
        #                                                                    nuc_charge, mass)

        line = '{0:10s}'.format(chem_symbol).rjust(16) + ' ' 
        line += '{0:.9f}'.format(xyz[0]).rjust(16) + ' ' 
        line += '{0:.9f}'.format(xyz[1]).rjust(16) + ' ' 
        line += '{0:.9f}'.format(xyz[2]).rjust(16) + ' ' 
        line += '{0:5s}'.format(nuc_charge).rjust(16) + ' ' 
        line += '{0:.9f}'.format(mass).rjust(16) + ' ' 

        f.write(line)
        f.write('\n')
    f.write('#\n')
