import urllib2
import urllib

from ase.test import NotAvailable

dest = 'demo.ascii'
src = 'http://inac.cea.fr/L_Sim/V_Sim/files/' + dest

try:
    e = urllib2.urlopen(src)
    urllib.urlretrieve(src, filename=dest)
except urllib2.URLError:
    raise NotAvailable('Retrieval of ' + src + ' failed')

from ase.io.v_sim import read_v_sim

a = read_v_sim(dest)
