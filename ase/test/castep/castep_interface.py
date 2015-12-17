"""Simple shallow test of the CASTEP interface"""

from ase.test.castep import installed


import os
import shutil
import tempfile
import traceback

def test_castep_interface():
    # check if we have castep installed
    installed()

    # check if we can import everything
    ase_castep_dir = "ase"

    try:
        castep_calc = __import__(ase_castep_dir + ".calculators.castep", globals(), locals(), ["Castep", "CastepParam", "create_castep_keywords"])
        Castep = castep_calc.Castep
        CastepParam = castep_calc.CastepParam
        create_castep_keywords = castep_calc.create_castep_keywords

    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError('Castep calculator module could not be loaded')

    try:
        __import__(ase_castep_dir + ".io.castep")
    except Exception as e:
        raise AssertionError('Castep io module could not be loaded')


    tmp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()

    from ase.calculators.castep import Castep

    try:
        c = Castep(directory=tmp_dir, label='test_label')
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError('Could not instantiate castep calculator')


    try:
        c.xc_functional = 'PBE'
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError('Setting xc_functional  failed')

    import ase.lattice.cubic
    lattice = ase.lattice.cubic.BodyCenteredCubic('Li' )

    print('For the sake of evaluating this test, warnings')
    print('about auto-generating pseudo-potentials are')
    print('normal behavior and can be safely ignored')

    try:
        lattice.set_calculator(c)
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError('Setting the calculator %s failed' % c)


    try:
        create_castep_keywords(
            castep_command=os.environ['CASTEP_COMMAND'],
            path=tmp_dir,
            fetch_only=20)
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError("Cannot create castep_keywords, this usually means a  bug"\
        + " in the interface or the castep binary cannot be called")


    param_fn = os.path.join(tmp_dir, 'myParam.param')
    param = open(param_fn,'w')
    param.write('XC_FUNCTIONAL : PBE #comment\n')
    param.write('XC_FUNCTIONAL : PBE #comment\n')
    param.write('#comment\n')
    param.write('CUT_OFF_ENERGY : 450.\n')
    param.close()
    try:
        c.merge_param(param_fn)
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise AssertionError("Error in merge_param_filename, go figure")


    # check if the CastepOpt, CastepCell comparison mechanism works

    p1 = CastepParam()
    p2 = CastepParam()

    if not p1._options == p2._options:
        raise AssertionError("Print two newly created CastepParams are not the same")

    p1._options['xc_functional'].value = 'PBE'
    p1.xc_functional = 'PBE'

    if p1._options == p2._options:
        raise AssertionError("Changed one CastepParam, but the still look the same")

    if not c.calculation_required(lattice):
        raise AssertionError("Calculator does not fetch that a calculation is required")

    if not c.dryrun_ok():
        print(c._error)
        raise AssertionError("Dryrun_ok does not work, where it should")
    else:
        print("Dryrun is ok")

    c.prepare_input_files(lattice)

    os.chdir(cwd)
    shutil.rmtree(tmp_dir)


    print("Test finished without errors")

if __name__ == '__main__':
    test_castep_interface()
