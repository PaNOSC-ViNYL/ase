# creates:  nitrogen.txt, ethane.txt, gold.txt
import io
import os
import sys


def output_to_string(pythonfile):
    """Returns the stdout of executing the code in pythonfile
    as a string."""
    if sys.version_info.major == 2:
        buffer = io.BytesIO()
    else:
        buffer = io.StringIO()
    sys.stdout = buffer
    exec(open(pythonfile).read())
    sys.stdout = sys.__stdout__
    return buffer.getvalue()

# Only save the parts relevant to thermochemistry
nitrogen = output_to_string('nitrogen.py')
nitrogen = nitrogen[nitrogen.find('Enthalpy'):]
with open('nitrogen.txt', 'w') as f:
    f.write(nitrogen)
ethane = output_to_string('ethane.py')
ethane = ethane[ethane.find('Internal'):]
with open('ethane.txt', 'w') as f:
    f.write(ethane)
gold = output_to_string('gold.py')
gold = gold[gold.find('Internal'):]
with open('gold.txt', 'w') as f:
    f.write(gold)

# Clean up.
vibfiles = [file for file in os.listdir(os.getcwd()) if
            file.startswith('vib.') or file.startswith('phonon.')]
for file in vibfiles:
    os.remove(file)
