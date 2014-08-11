# creates: mlab_options.txt
import sys
from ase.visualize.mlab import main
fd = open('mlab_options.txt', 'w')
try:
    sys.stdout, old = fd, sys.stdout
    main(['-h'])
    fd.close()
finally:
    sys.stdout = old

