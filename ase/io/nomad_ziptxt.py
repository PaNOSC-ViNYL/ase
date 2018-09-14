import ase 
from ase.utils import basestring

def read_nomad_ziptxt(fd, index=':', only_atoms=False, skip_errors=False):
    images = []
    from ase.io.formats import string2index
    if isinstance(index, basestring):
        index = string2index(index)
    for bline in fd:
        line = bline.decode("utf-8")
        if line.startswith('#'):
            pass
        else:
            nmduri = line.split('/section_run')
            print('Requesting NOMAD archive at ' + nmduri[0])
            entry = ase.nomad.download(nmduri[0], only_atoms=only_atoms, skip_errors=skip_errors)
            nmd_entry_images = entry.toatoms()
            nmd_images = list(nmd_entry_images)
            if len(nmd_images)>0:
                print('Adding ' + str(len(nmd_images)) + ' structure(s) with ' + ','.join(
                    list(set([str(ni.get_chemical_formula('reduce')) for ni in nmd_images]))))
            else:
                print('No structures retrieved from this NOMAD archive!')
            images.extend(nmd_images)
    return list(images)[index]
