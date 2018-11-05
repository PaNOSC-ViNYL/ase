import sys

sys.path.insert(0, '../../doc/ase/')
import fix_symmetry_example

assert fix_symmetry_example.d_init6["number"] == 229
assert fix_symmetry_example.d_init8["number"] == 99
assert fix_symmetry_example.d_unsym["number"] == 1
assert fix_symmetry_example.d_sym["number"] == 229
