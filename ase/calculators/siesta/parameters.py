from ase.calculators.calculator import LockedParameters


class PAOBasisBlock(LockedParameters):
    """
    Representing a block in PAO.Basis for one species.
    """
    def __init__(self, block):
        """
        Parameters:
            -block : String. A block defining the basis set of a single
                     specie using the format of a PAO.Basis block.
                     The initial label should be left out since it is
                     determined programatically.
                     Example1: 2 nodes 1.0
                               n=2 0 2 E 50.0 2.5
                               3.50 3.50
                               0.95 1.00
                               1 1 P 2
                               3.50
                     Example2: 1
                               0 2 S 0.2
                               5.00 0.00
                     See siesta manual for details.
        """
        assert isinstance(block, str)
        LockedParameters.__init__(self, block=block)

    def script(self, label):
        """
        Write the fdf script for the block.

        Parameters:
            -label : The label to insert in front of the block.
        """
        return label + ' ' + self['block']


class Specie(LockedParameters):
    """
    Parameters for specifying the behaviour for a single species in the
    calculation. If the tag argument is set to an integer then atoms with
    the specified element and tag will be a seperate species.

    Pseudopotential and basis set can be specified. Additionally the species
    can be set be a ghost species, meaning that they will not be considered
    atoms, but the corresponding basis set will be used.
    """
    def __init__(self,
                 symbol,
                 basis_set='DZP',
                 pseudopotential=None,
                 tag=None,
                 ghost=False,
                 ):
        kwargs = locals()
        kwargs.pop('self')
        LockedParameters.__init__(self, **kwargs)


class SiestaParameters(LockedParameters):
    """
    Parameter class which can write out its own fdf-script.
    """
    def write_fdf(self, f):
        for key, value in self.iteritems():
            key = self.prefix() + '.' + key
            f.write(format_fdf(key, value))


class SolutionMethod(SiestaParameters):
    """
    Collection of parameters related to a specific solution method.
    """
    def identitier(self):
        """
        The string which begins all fdf-keywords in the group.
        """
        raise NotImplementedError

    def write_fdf(self, f):
        """
        Write the SolutionMethod keyword to the fdf-script as well as all
        parameters in this group.
        """
        f.write(format_fdf('SolutionMethod', self.identifier()))
        SiestaParameters.write_fdf(self, f)


class Diag(SolutionMethod):
    """
    Parameters related to the diagonalization solution method.
    """
    def prefix(self):
        return 'Diag'

    def identifier(self):
        return 'diagon'

    def __init__(
            self,
            DivideAndConquer=False,
            AllInOne=False,
            NoExpert=False,
            PreRotate=False,
            Use2D=False,
            Memory=1.0,
            ParallelOverK=False, ):
        kwargs = locals()
        kwargs.pop('self')
        SolutionMethod.__init__(self, **kwargs)


class OrderN(SolutionMethod):
    """
    Parameters related to the OrderN solution method.
    """
    def prefix(self):
        return 'ON'

    def identifier(self):
        return 'ON'

    def __init__(
            self,
            functional='Kim',
            MaxNumIter=1000,
            etol=1e-8,
            eta='0.0 eV',
            eta_alpha='0.0 eV',
            eta_beta='0.0 eV',
            RcLWF='9.5 Bohr',
            ChemicalPotential=False,
            ChemicalPotentialUse=False,
            ChemicalPotentialRc='9.5 Bohr',
            ChemicalPotentialTemperature='0.05 Ry',
            ChemicalPotentialOrder=100,
            LowerMemory=False,
            UseSaveLWF=False,
            OccupationFunction='FD',
            OccupationMPOrder=1, ):

        kwargs = locals()
        kwargs.pop('self')
        SolutionMethod.__init__(self, **kwargs)


def format_fdf(key, value):
    """
    Write an fdf key-word value pair.

    Parameters:
        - key   : The fdf-key
        - value : The fdf value.
    """
    if isinstance(value, (list, tuple)) and len(value) == 0:
        return ''

    key = format_key(key)
    new_value = format_value(value)

    if isinstance(value, list):
        string = '%block ' + key + '\n' +\
            new_value + '\n' + \
            '%endblock ' + key + '\n'
    else:
        string = '%s  %s\n' % (key, new_value)

    return string


def format_value(value):
    """
    Format python values to fdf-format.

    Parameters:
        - value : The value to format.
    """
    if isinstance(value, tuple):
        sub_values = map(format_value, value)
        value = '\t'.join(sub_values)
    elif isinstance(value, list):
        sub_values = map(format_value, value)
        value = '\n'.join(sub_values)
    else:
        value = str(value)

    return value


def format_key(key):
    """ Fix the fdf-key replacing '_' with '.' """
    return key.replace('_', '.')
