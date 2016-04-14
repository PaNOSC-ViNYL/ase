import optparse
import sys

from ase.test import test


parser = optparse.OptionParser(
    usage='Usage: python -m ase.test [-c calc1,calc2,...]',
    description='Test ASE')

parser.add_option('-c', '--calculators',
                  help='Comma-separated list of calculators to test.')
parser.add_option('-v', '--verbosity', type=int, default=2, metavar='N',
                  help='Use 0, 1 or 2.')
parser.add_option('-g', '--test-also-gui', action='store_true',
                  help='Test also ase-gui.')

opts, args = parser.parse_args()

if opts.calculators:
    calculators = opts.calculators.split(',')
else:
    calculators = []

results = test(display=opts.test_also_gui,
               verbosity=opts.verbosity,
               calculators=calculators)
sys.exit(len(results.errors + results.failures))
