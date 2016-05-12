from math import pi, sqrt

# this is the hard-coded CODATA values
# all other units are dynamically derived from these values upon import of the
# module
CODATA = {# the "original" CODATA version ase used ever since
          # Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA)
          # Add the constant pi used to define the mu0 and hbar here for reference as well
          '1986' : {'_c' : 299792458.,              # speed of light, m/s
                    '_mu0' : 4.e-7 * pi,            # permeability of vacuum
                    '_Grav' : 6.67259e-11,          # gravitational constant
                    '_hplanck' : 6.6260755e-34,     # Planck constant, J s
                    '_e' : 1.60217733e-19,          # elementary charge
                    '_me' : 9.1093897e-31,          # electron mass
                    '_mp' : 1.6726231e-27,          # proton mass
                    '_Nav' : 6.0221367e23,          # Avogadro number
                    '_k' : 1.380658e-23,            # Boltzmann constant, J/K
                    '_amu' : 1.6605402e-27,         # atomic mass unit, kg
                    '_pi' : pi                      # pi used for derived units
                    },

          # CODATA 1998 taken from
          # http://dx.doi.org/10.1103/RevModPhys.72.351
          '1998' : {'_c' : 299792458.,
                    '_mu0' : 4.0e-7*pi,
                    '_Grav' : 6.673e-11,
                    '_hplanck' : 6.62606876e-34,
                    '_e' : 1.602176462e-19,
                    '_me' : 9.10938188e-31,
                    '_mp' : 1.67262158e-27,
                    '_Nav' : 6.02214199e23,
                    '_k' : 1.3806503-23,
                    '_amu' : 1.66053873e-27,
                    '_pi' : pi
                    },

          # CODATA 2002 taken from
          # http://dx.doi.org/10.1103/RevModPhys.77.1
          '2002' : {'_c' : 299792458.,
                    '_mu0' : 4.0e-7*pi,
                    '_Grav' : 6.6742e-11,
                    '_hplanck' : 6.6260693e-34,
                    '_e' : 1.60217653e-19,
                    '_me' : 9.10938226e-31,
                    '_mp' : 1.67262171e-27,
                    '_Nav' : 6.0221415e23,
                    '_k' : 1.3806505-23,
                    '_amu' : 1.66053886e-27,
                    '_pi' : pi
                    },

          # CODATA 2006 taken from
          # http://dx.doi.org/10.1103/RevModPhys.80.633
          '2006' : {'_c' : 299792458.,
                    '_mu0' : 4.0e-7*pi,
                    '_Grav' : 6.67428e-11,
                    '_hplanck' : 6.62606896e-34,
                    '_e' : 1.602176487e-19,
                    '_me' : 9.10938215e-31,
                    '_mp' : 1.672621637e-27,
                    '_Nav' : 6.02214179e23,
                    '_k' : 1.3806504-23,
                    '_amu' : 1.660538782e-27,
                    '_pi' : pi
                    },

          # CODATA 2010 taken from
          # http://dx.doi.org/10.1103/RevModPhys.84.1527
          '2010' : {'_c' : 299792458.,
                    '_mu0' : 4.0e-7*pi,
                    '_Grav' : 6.67384e-11,
                    '_hplanck' : 6.62606957e-34,
                    '_e' : 1.602176565e-19,
                    '_me' : 9.10938291e-31,
                    '_mp' : 1.672621777e-27,
                    '_Nav' : 6.02214129e23,
                    '_k' : 1.3806488-23,
                    '_amu' : 1.660538921e-27,
                    '_pi' : pi
                    },

          # CODATA 2014 taken from
          # http://arxiv.org/pdf/1507.07956.pdf
          '2014' : {'_c' : 299792458.,
                    '_mu0' : 4.0e-7*pi,
                    '_Grav' : 6.67408e-11,
                    '_hplanck' : 6.626070040e-34,
                    '_e' : 1.6021766208e-19,
                    '_me' : 9.10938356e-31,
                    '_mp' : 1.672621898e-27,
                    '_Nav' : 6.022140857e23,
                    '_k' : 1.38064852e-23,
                    '_amu' : 1.660539040e-27,
                    '_pi' : pi
                    }
        }


# ok, this is a bit dirty and bad practice, any better suggestions?
def create_units(codata_version):
    """
    Function that creates a dictionary containing all units previously hard
    coded in ase.units depending on a certain CODATA version. Note that you can
    use the dictionary returned to update your local or global namespace.

    Parameters
    ----------
    codata_version : str
        The CODATA version to be used. Implemented are
            * '1986'
            * '1998'
            * '2002'
            * '2006'
            * '2010'
            * '2014'

    Returns
    -------
    units : dict
        Dictionary that contains all formerly hard coded variables from
        ase.units as key-value pairs.

    Raises
    ------
    NotImplementedError
        If the required CODATA version is not known.
    """

    import copy
    try:
        units = copy.deepcopy(CODATA[codata_version])
    except KeyError:
        raise NotImplementedError('CODATA version "{}" not implemented'.format(__codata_version__))

    # derived from the CODATA values
    units['_eps0'] = 1 / units['_mu0'] / units['_c']**2     # permittivity of vacuum
    units['_hbar'] = units['_hplanck'] / (2 * units['_pi'])  # Planck constant / 2pi, J s

    units['Ang'] = units['Angstrom'] = 1.0
    units['nm'] = 10.0
    units['Bohr'] = 4e10 * units['_pi'] * units['_eps0'] * units['_hbar']**2 / units['_me'] / units['_e']**2  # Bohr radius

    units['eV'] = 1.0
    units['Hartree'] = units['_me'] * units['_e']**3 / 16 / units['_pi']**2 / units['_eps0']**2 / units['_hbar']**2
    units['kJ'] = 1000.0 / units['_e']
    units['kcal'] = 4.184 * units['kJ']
    units['mol'] = units['_Nav']
    units['Rydberg'] = 0.5 * units['Hartree']
    units['Ry'] = units['Rydberg']
    units['Ha'] = units['Hartree']

    units['second'] = 1e10 * sqrt(units['_e'] / units['_amu'])
    units['fs'] = 1e-15 * units['second']

    units['kB'] = units['_k'] / units['_e']                 # Boltzmann constant, eV/K

    units['Pascal'] = (1 / units['_e']) / 1e30  # J/m^3
    units['GPa'] = 1e9 * units['Pascal']

    units['Debye'] = 1.0 / 1e11 / units['_e'] / units['_c']
    units['alpha'] = units['_e']**2 / (4 * units['_pi'] * units['_eps0']) / units['_hbar'] / units['_c']  # fine structure constant
    units['invcm'] = 100 * units['_c'] * units['_hplanck'] / units['_e']               # cm^-1 energy unit

    # Derived atomic units that have no assigned name:
    units['_aut'] = units['_hbar'] / (units['alpha']**2 * units['_me'] * units['_c']**2)      # atomic unit of time, s
    units['_auv'] = units['_e']**2 / units['_hbar'] / (4 * units['_pi'] * units['_eps0'])      # atomic unit of velocity, m/s
    units['_auf'] = units['alpha']**3 * units['_me']**2 * units['_c']**3 / units['_hbar']     # atomic unit of force, N
    units['_aup'] = units['alpha']**5 * units['_me']**4 * units['_c']**5 / units['_hbar']**3  # atomic unit of pressure, Pa

    units['AUT'] = units['second'] * units['_aut']

    # SI units
    units['m'] = 1e10 * units['Ang']    # metre
    units['kg'] = 1. / units['_amu']    # kilogram
    units['s'] = units['second']        # second
    units['A'] = 1.0 / units['_e'] / units['s']  # ampere
    # derived
    units['J'] = units['kJ'] / 1000  # Joule   = kg * m**2 / s**2
    units['C'] = 1.0 / units['_e']   # Coulomb = A * s

    return units


def __update_units(codata_version):
    globals().update(create_units(__codata_version__))


def __add_labelled_units():
    """
    Add units that carry the version number of the respective CODATA
    suggestions as as suffix, as suggested here:
        https://gitlab.com/ase/ase/issues/8
    """
    for codata_version in CODATA.keys():
        units = create_units(codata_version)
        globals().update({'{}_{}'.format(k, codata_version) : v for (k,v) in units.items()})
        globals().update({'units_{}'.format(codata_version) : units})

# the version we actually use
__codata_version__ = '2014'

# now update the module scope
__update_units(__codata_version__)

# add the labelled units
__add_labelled_units()
