from ase.calculators.jacapo import Jacapo
from __init__ import *

Jacapo.script = 'dacapo_run'
Jacapo.calculate = queue_job
#Jacapo.update_from_database = update_from_database
