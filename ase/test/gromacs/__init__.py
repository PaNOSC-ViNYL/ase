def installed():
    import os
    from ase.test import NotAvailable
    gromacs_mdrun = os.getenv('GMXCMD')
    #print 'gromacs_mdrun', gromacs_mdrun
    if gromacs_mdrun == None:
        raise NotAvailable(\
            'Gromacs_COMMAND ($GMXCMD) not defined (use: "setenv GMXCMD mdrun")')
    return True
