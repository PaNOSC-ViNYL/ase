import numpy as np

def installed():
    try:
        import scipy
        return True
    except ImportError:
        raise NotAvailable()


# This test cross-checks our implementation of CODATA against the
# implementation that SCIPY brings with it

def test_units():
    installed()

    name_map = {'_c' : 'speed of light in vacuum',
                '_mu0' : 'mag. const.',
                '_Grav' : 'Newtonian constant of gravitation',
                '_hplanck' : 'Planck constant',
                '_e' : 'elementary charge',
                '_me' : 'electron mass',
                '_mp' : 'proton mass',
                '_Nav' : 'Avogadro constant',
                '_k' : 'Boltzmann constant',
                '_amu' : 'atomic mass unit-kilogram relationship'}

    from ase.units import CODATA
    import scipy.constants.codata
    for version in sorted(CODATA.keys()):
        print('Checking CODATA version "{}"'.format(version))

        try:
            scipy_CODATA = getattr(scipy.constants.codata, '_physical_constants_{}'.format(version))
        except AttributeError:
            print('\tNot available through scipy, skipping')
            continue

        for unit, scipyname in name_map.items():
            aseval = CODATA[version][unit]
            try:
                # 2002 in scipy contains too little data
                scipyval = scipy_CODATA[name_map[unit]][0]
                msg = 'Unit "{}" : '.format(name_map[unit])
                ok = True
                if np.isclose(aseval, scipyval):
                    msg += '[OK]'
                else:
                    msg += '[FALSE]'
                    ok = False
                print('\t'+msg)
                if not ok:
                    raise AssertionError

            except KeyError:
                continue
if __name__ == '__main__':
    test_units()
