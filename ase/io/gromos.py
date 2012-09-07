""" write gromos96 geometry files 
(the exact file format is copied from the freely available 
gromacs package, http://www.gromacs.org
its procedure src/gmxlib/confio.c (write_g96_conf)
"""

from ase.parallel import paropen


def write_gromos(fileobj, images):
    """Write gromos geometry files (\*.g96).
    Writes:
    atom positions,
    and simulation cell (if present)
    """

    from ase import units

    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    natoms = images[-1].get_number_of_atoms()
    try:
        gromos_residuenames = images[-1].get_array('residuenames')
    except:
        gromos_residuenames = []
        for idum in range(natoms):
            gromos_residuenames.append('1DUM')
    try:
        gromos_atomtypes = images[-1].get_array('atomtypes')
    except:
        gromos_atomtypes = images[-1].get_chemical_symbols()
    #print "gromos_atomtypes", gromos_atomtypes

    pos = images[-1].get_positions()
    pos = pos / 10.0
    try:
        vel = images[-1].get_velocities()
        vel = vel * 1000.0 * units.fs / units.nm
    except:
        vel = pos
        vel = pos * 0.0

    fileobj.write('TITLE \n')
    fileobj.write('Gromos96 structure file written by ASE \n')
    fileobj.write('END \n')
    fileobj.write('POSITION \n')
    count = 1
    rescount = 0
    oldresname = ''
    for resname, atomtype, xyz in zip\
            (gromos_residuenames, gromos_atomtypes, pos):
        if resname != oldresname:
            oldresname = resname
            rescount = rescount + 1
        #print xyz
        okresname = resname.lstrip('0123456789 ')
        fileobj.write('%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n' % \
                          (rescount, okresname, atomtype, count, \
                               xyz[0], xyz[1], xyz[2]))
        count = count + 1
    fileobj.write('END \n')
    if images[-1].get_pbc().any():
        fileobj.write('BOX \n')
        mycell = images[-1].get_cell()
        fileobj.write('%15.9f%15.9f%15.9f' \
                          % (mycell[0, 0] * 0.1, \
                                 mycell[1, 1] * 0.1, \
                                 mycell[2, 2] * 0.1))
        fileobj.write('%15.9f%15.9f%15.9f' \
                          % (mycell[1, 0] * 0.1, \
                                 mycell[2, 0] * 0.1, \
                                 mycell[0, 1] * 0.1))
        fileobj.write('%15.9f%15.9f%15.9f\n' \
                          % (mycell[2, 1] * 0.1, \
                                 mycell[0, 2] * 0.1, \
                                 mycell[1, 2] * 0.1))
        fileobj.write('END \n')
    return
