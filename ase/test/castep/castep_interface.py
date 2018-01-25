"""Simple shallow test of the CASTEP interface"""
import os
import shutil
import tempfile
import ase
import re
import ase.lattice.cubic
from ase.calculators.castep import (Castep, CastepParam,
                                    create_castep_keywords,
                                    import_castep_keywords)

tmp_dir = tempfile.mkdtemp()
cwd = os.getcwd()

c = Castep(directory=tmp_dir, label='test_label')
c.xc_functional = 'PBE'

lattice = ase.lattice.cubic.BodyCenteredCubic('Li')

print('For the sake of evaluating this test, warnings')
print('about auto-generating pseudo-potentials are')
print('normal behavior and can be safely ignored')

lattice.set_calculator(c)

create_castep_keywords(
    castep_command=os.environ['CASTEP_COMMAND'],
    path=tmp_dir,
    fetch_only=20)

param_fn = os.path.join(tmp_dir, 'myParam.param')
param = open(param_fn, 'w')
param.write('XC_FUNCTIONAL : PBE #comment\n')
param.write('XC_FUNCTIONAL : PBE #comment\n')
param.write('#comment\n')
param.write('CUT_OFF_ENERGY : 450.\n')
param.close()
c.merge_param(param_fn)

# check if the CastepOpt, CastepCell comparison mechanism works

castep_keywords = import_castep_keywords()
p1 = CastepParam(castep_keywords)
p2 = CastepParam(castep_keywords)

assert p1._options == p2._options

p1._options['xc_functional'].value = 'PBE'
p1.xc_functional = 'PBE'

assert p1._options != p2._options

assert c.calculation_required(lattice)

assert c.dryrun_ok()

c.prepare_input_files(lattice)


# detecting pseudopotentials tests

# typical filenames
files = ['Ag_00PBE.usp',
         'Ag_00.recpot',
         'Ag_C18_PBE_OTF.usp',
         'ag-optgga1.recpot',
         'Ag_OTF.usp',
         'ag_pbe_v1.4.uspp.F.UPF',
         'Ni_OTF.usp',
         'fe_pbe_v1.5.uspp.F.UPF',
         'Cu_01.recpot']

pp_path = os.path.join(tmp_dir, 'test_pp')
os.makedirs(pp_path)

for f in files:
    with open(os.path.join(pp_path, f), 'w') as _f:
        _f.write('DUMMY PP')


c = Castep(directory=tmp_dir, label='test_label_pspots', castep_pp_path=pp_path)
c._pedantic = True
atoms = ase.build.bulk('Ag')
atoms.set_calculator(c)

# I know, unittest would be nicer... maybe at a later point

# disabled, but may be useful still
# try:
    # # this should yield no files
    # atoms.calc.find_pspots(suffix='uspp')
    # raise AssertionError
# except RuntimeError as e:
    # #print(e)
    # pass

try:
    # this should yield non-unique files
    atoms.calc.find_pspots(suffix='recpot')
    raise AssertionError
except RuntimeError as e:
    #print(e)
    pass


# now let's see if we find all...
atoms.calc.find_pspots(pspot='00PBE', suffix='usp')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_00PBE.usp'

atoms.calc.find_pspots(pspot='00', suffix='recpot')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_00.recpot'

atoms.calc.find_pspots(pspot='C18_PBE_OTF', suffix='usp')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_C18_PBE_OTF.usp'

atoms.calc.find_pspots(pspot='optgga1', suffix='recpot')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'ag-optgga1.recpot'

atoms.calc.find_pspots(pspot='OTF', suffix='usp')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_OTF.usp'

atoms.calc.find_pspots(suffix='UPF')
assert atoms.calc.cell.species_pot.value.split()[-1] == 'ag_pbe_v1.4.uspp.F.UPF'


# testing regular workflow
c = Castep(directory=tmp_dir, label='test_label_pspots',
        castep_pp_path=pp_path, find_pspots=True)
c._build_missing_pspots = False
atoms = ase.build.bulk('Ag')
atoms.set_calculator(c)

# this should raise an error due to ambuiguity
try:
    c._fetch_pspots()
    raise AssertionError
except RuntimeError as e:
    #print(e)
    pass

for e in ['Ni', 'Fe', 'Cu']:
    atoms = ase.build.bulk(e)
    atoms.set_calculator(c)
    c._fetch_pspots()

# test writing to file
tmp_dir = os.path.join(tmp_dir, 'input_files')
c = Castep(directory=tmp_dir, label='test_label_pspots',
        find_pspots=True, castep_pp_path=pp_path)
c._label = 'test'
atoms = ase.build.bulk('Cu')
atoms.set_calculator(c)
c.prepare_input_files()

with open(os.path.join(tmp_dir, 'test.cell'), 'r') as f:
    assert re.search('Cu Cu_01\.recpot', ''.join(f.readlines())) is not None


os.chdir(cwd)
shutil.rmtree(tmp_dir)

