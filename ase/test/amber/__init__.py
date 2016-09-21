def installed():
    import os

    from ase.test import NotAvailable
    amberhome=os.getenv('AMBERHOME')
    if not os.path.isfile(amberhome+'/bin/sander'):
        raise NotAvailable(
            'No sander in $AMBERHOME/bin/sander (cannot run Amber)')
    return True
