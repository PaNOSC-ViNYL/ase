def installed():
    import os
    from ase.test import NotAvailable
    vase = os.getenv('ASE_VASP_COMMAND')
    vcmd = os.getenv('VASP_COMMAND')
    vscr = os.getenv('VASP_SCRIPT')
    if vcmd is None and vscr is None and vase is None:
        raise NotAvailable(('Neither ASE_VASP_COMMAND, VASP_COMMAND '
                            'nor VASP_SCRIPT defined'))
    return True
