try:
    from urllib.request import urlretrieve
    from urllib.error import URLError
except ImportError:
    from urllib import urlretrieve
    from urllib2 import URLError
from socket import error as SocketError

from ase.test import NotAvailable
from ase.io import read

raise NotAvailable

dest = 'demo.ascii'
src = 'http://inac.cea.fr/L_Sim/V_Sim/files/' + dest

try:
    urlretrieve(src, filename=dest)
except (IOError, URLError, SocketError):
    raise NotAvailable('Retrieval of ' + src + ' failed')

a = read(dest, format='v_sim')
