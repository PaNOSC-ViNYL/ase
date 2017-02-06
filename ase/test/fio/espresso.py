"""Quantum ESPRESSO file parsers.

Implemented:
* Input file (pwi)

"""

import os

from ase import io


# This file is parsed correctly by pw.x, even though things are
# scattered all over the place with some namelist edge cases
pw_input_text = """
&CONTrol
   prefix           = 'surf_110_H2_md'
   calculation      = 'md'
   restart_mode     = 'from_scratch'
   pseudo_dir       = '.'
   outdir           = './surf_110_!H2_m=d_sc,ratch/'
   verbosity        = 'default'
   tprnfor          = .true.
   tstress          = .True.
!   disk_io          = 'low'
   wf_collect       = .false.
   max_seconds      = 82800
   forc_con!v_thr    = 1e-05
   etot_conv_thr    = 1e-06
   dt               = 41.3 , /

&SYSTEM ecutwfc     = 63,   ecutrho   = 577,  ibrav    = 0,
nat              = 8,   ntyp             = 2,  occupations      = 'smearing',
smearing         = 'marzari-vanderbilt',
degauss          = 0.01,   nspin            = 2,  !  nosym     = .true. ,
    starting_magnetization(2) = 0.32 /
&ELECTRONS
   electron_maxstep = 300
   mixing_beta      = 0.1
   conv_thr         = 1d-07
   mixing_mode      = 'local-TF'
   scf_must_converge = False
/
&IONS
   ion_dynamics     = 'verlet'
   ion_temperature  = 'rescaling'
   tolp             = 50.0
   tempw            = 500.0
/

ATOMIC_SPECIES
H 1.008 H.pbe-rrkjus_psl.0.1.UPF
Fe 55.845 Fe.pbe-spn-rrkjus_psl.0.2.1.UPF

K_POINTS automatic
2 2 2  1 1 1

CELL_PARAMETERS angstrom
5.6672000000000002 0.0000000000000000 0.0000000000000000
0.0000000000000000 8.0146311006808038 0.0000000000000000
0.0000000000000000 0.0000000000000000 27.0219466510212101

ATOMIC_POSITIONS angstrom
Fe 0.0000000000 0.0000000000 0.0000000000 0 0 0
Fe 1.4168000000 2.0036577752 -0.0000000000 0 0 0
Fe 0.0000000000 2.0036577752 2.0036577752 0 0 0
Fe 1.4168000000 0.0000000000 2.0036577752 0 0 0
Fe 0.0000000000 0.0000000000 4.0073155503
Fe 1.4168000000 2.0036577752 4.0073155503
H 0.0000000000 2.0036577752 6.0109733255
H 1.4168000000 0.0000000000 6.0109733255
"""


def test_pw_input():
    with open('pw_input.pwi', 'w') as pw_input_f:
        pw_input_f.write(pw_input_text)

    try:
        pw_input_atoms = io.read('pw_input.pwi', format='espresso-in')
        assert len(pw_input_atoms) == 8
    finally:
        os.unlink('pw_input.pwi')


if __name__ in ('__main__', '__builtin__'):
    test_pw_input()
