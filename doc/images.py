# creates: ag.png water_divide_surf.png
from urllib import urlretrieve
for name in ['ag.png', 'water_divide_surf.png']:
    urlretrieve('http://wiki.fysik.dtu.dk/ase-files/' + name, name)
