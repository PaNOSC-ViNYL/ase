import urllib.request, urllib.error, urllib.parse
import urllib.request, urllib.parse, urllib.error

from ase.test import NotAvailable

dest = 'demo.ascii'
src = 'http://inac.cea.fr/L_Sim/V_Sim/files/' + dest

try:
    e = urllib.request.urlopen(src)
    urllib.request.urlretrieve(src, filename=dest)
except urllib.error.URLError:
    raise NotAvailable('Retrieval of ' + src + ' failed')

from ase.io import read

a = read(dest, format='v_sim')
