def installed():
    import os
    from ase.test import NotAvailable
    lammpscmd = os.getenv('LAMMPS_COMMAND')
    if lammpscmd == None:
        raise NotAvailable('LAMMPS_COMMAND envronment variable not defined')
    return True
