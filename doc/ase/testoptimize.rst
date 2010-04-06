.. _optimizer_tests:

===============
Optimizer tests
===============
This page shows benchmarks of optimizations done with our different optimizers.
Note that the iteration number (steps) is not the same as the number of force
evaluations. This is because some of the optimizers uses internal line searches
or similar.

N2Ru-surf
=========
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               55                55  14.024151 
LBFGS              55                55  14.024151 
FIRE               59                59  14.027805 
MDMin              19                19  14.039306 
SciPyFminCG        11                37  14.022259 
SciPyFminBFGS      26                28  14.022231 
BFGSLineSearch     13                13  14.027444 
=============== ===== ================= ========== ===============

N2Ru-N2
=======
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               88                88  20.658067 
LBFGS              88                88  20.658067 
FIRE              200               200  20.685025 Not converged in 200 steps
MDMin              49                49  32.053754 
SciPyFminCG        56               112  20.660197 
SciPyFminBFGS      69                71  20.659189 
BFGSLineSearch     34                35  20.657959 
=============== ===== ================= ========== ===============

Cu_bulk
=======
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               22                22  -0.031977 
LBFGS              22                22  -0.031977 
FIRE               58                58  -0.032252 
MDMin              12                12  -0.031269 
SciPyFminCG         8                23  -0.032246 
SciPyFminBFGS      21                22  -0.031997 
BFGSLineSearch      8                 8  -0.032193 
=============== ===== ================= ========== ===============

CO_Au111
========
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS              146               146  -0.289119 
LBFGS             125               125  -0.288257 
FIRE              112               112  -0.285567 
MDMin              21                21  15.077292 
SciPyFminCG        26                57  -0.287460 
SciPyFminBFGS      34                36  -0.285501 
BFGSLineSearch     36                37  -0.255474 
=============== ===== ================= ========== ===============

H2-emt
======
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               11                11   0.000000 
LBFGS               9                 9   4.417486 
FIRE               29                29   0.000019 
MDMin              31                31   0.000042 
SciPyFminCG         3                10   0.000000 
SciPyFminBFGS       4                 7   0.000002 
BFGSLineSearch      4                 5   0.000001 
=============== ===== ================= ========== ===============

H2-gpaw
=======
Calculator used: GPAW

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               10                10  -6.641238 
LBFGS              25                25  -0.579014 Not converged in 25 steps
FIRE               15                15  -6.641190 
MDMin              13                13        nan An exception occurred
SciPyFminCG         2                 3        nan An exception occurred
SciPyFminBFGS       2                 5  -6.641215 
BFGSLineSearch      4                 7  -6.641196 
=============== ===== ================= ========== ===============

nanoparticle
============
Calculator used: GPAW

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS              200               200 -16.487239 Not converged in 200 steps
LBFGS              20                20 -15.274295 
FIRE               83                83 -16.286412 
MDMin              45                45 -16.286837 
SciPyFminCG       187               755 -16.527323 
SciPyFminBFGS      13                22        nan An exception occurred
BFGSLineSearch     31                38 -16.526885 
=============== ===== ================= ========== ===============

neb-emt
=======
Calculator used: EMT

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               17                17   3.615078 
LBFGS              17                17   3.615078 
FIRE               35                35   3.616566 
MDMin              12                12   3.615979 
SciPyFminCG       101               774   4.158123 
SciPyFminBFGS       4                 9        nan An exception occurred
BFGSLineSearch     24               137   3.601026 
=============== ===== ================= ========== ===============

neb-gpaw
========
Calculator used: GPAW (lcao)

=============== ===== ================= ========== ===============
Optimizer       Steps Force evaluations Energy     Note           
=============== ===== ================= ========== ===============
BFGS               24                24 -46.365794 
LBFGS              24                24 -46.365794 
FIRE               47                47 -46.365964 
MDMin              10                10 -46.366941 
SciPyFminCG       101               297 -46.083070 
SciPyFminBFGS      13                84        nan An exception occurred
BFGSLineSearch     16                51 -46.366182 
=============== ===== ================= ========== ===============
