def installed():
    import os
    import subprocess
    from ase.test import NotAvailable
    # check if CASTEP_COMMAND is set an environment variable
    if 'CASTEP_COMMAND' not in os.environ:
        print("WARNING: Environment variable CASTEP_COMMAND is not set")
        print("Will set CASTEP_COMMAND  = castep for the sake of this test")
        print("Please change it if this does not run castep in your environment")
        os.environ['CASTEP_COMMAND'] = 'castep'

    try:
        # not python 2.6-safe
        #subprocess.check_output([os.environ['CASTEP_COMMAND']])
        devnull = open(os.devnull)
        subprocess.Popen([os.environ['CASTEP_COMMAND']], stdout=devnull, stderr=devnull)
    except OSError:
        raise NotAvailable("""Could not find CASTEP. If you have it
                              installed make sure, you set the CASTEP_COMMAND
                              environment variable correctly""")
    return True
