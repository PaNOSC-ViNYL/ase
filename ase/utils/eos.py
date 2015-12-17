# -*- coding: utf-8 -*-


from ase.units import kJ
from ase.utils.sjeos import EquationOfStateSJEOS

try:
    # ase.utils.eosase2 requires scipy
    from ase.utils.eosase2 import EquationOfStateASE2

    class EquationOfState:
        """Fit equation of state for bulk systems.

        The following equation is used::

           sjeos (default)
               A third order inverse polynomial fit 10.1103/PhysRevB.67.026103

                               2      3        -1/3
           E(V) = c + c t + c t  + c t ,  t = V
                   0   1     2      3

           taylor
               A third order Taylor series expansion about the minimum volume

           murnaghan
               PRB 28, 5480 (1983)

           birch
               Intermetallic compounds: Principles and Practice,
               Vol I: Principles. pages 195-210

           birchmurnaghan
               PRB 70, 224107

           pouriertarantola
               PRB 70, 224107

           vinet
               PRB 70, 224107

           antonschmidt
               Intermetallics 11, 23-32 (2003)

           p3
               A third order polynomial fit

        Use::

           eos = EquationOfState(volumes, energies, eos='sjeos')
           v0, e0, B = eos.fit()
           eos.plot()

        """
        def __init__(self, volumes, energies, eos='sjeos'):
            if eos == 'sjeos':
                eosclass = EquationOfStateSJEOS
            else:
                eosclass = EquationOfStateASE2  # old ASE2 implementation
            self._impl = eosclass(volumes, energies, eos)

        def fit(self):
            return self._impl.fit()

        def plot(self, filename=None, show=None):
            return self._impl.plot(filename, show)

except ImportError:
    # ase.utils.sjeos requires only numpy
    EquationOfState = EquationOfStateSJEOS
    

def main():
    import optparse
    from ase.io import read
    parser = optparse.OptionParser(usage='python -m ase.io.eos [options] '
                                   'filename, ...',
                                   description='Calculate equation of state.')
    parser.add_option('-p', '--plot', action='store_true')
    opts, args = parser.parse_args()
    if not opts.plot:
        print('# filename                '
              'points     volume    energy  bulk modulus')
        print('#                         '
              '          [Ang^3]      [eV]         [GPa]')
    for name in args:
        if name == '-':
            # Special case - used by ase-gui:
            import pickle
            import sys
            v, e = pickle.load(sys.stdin)
        else:
            if '@' in name:
                index = None
            else:
                index = ':'
            images = read(name, index=index)
            v = [atoms.get_volume() for atoms in images]
            e = [atoms.get_potential_energy() for atoms in images]
        eos = EquationOfState(v, e)
        if opts.plot:
            eos.plot()
        else:
            try:
                v0, e0, B = eos.fit()
            except ValueError as ex:
                print('{0:30}{1:2}    {2}'.format(name, len(v), ex.message))
            else:
                print('{0:30}{1:2} {2:10.3f}{3:10.3f}{4:14.3f}'
                      .format(name, len(v), v0, e0, B / kJ * 1.0e24))
            
            
if __name__ == '__main__':
    main()
