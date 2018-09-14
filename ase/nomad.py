import json
import numpy as np

import ase.units as units
from ase import Atoms
from ase.io import nomad_json
from ase.data import chemical_symbols
from ase.calculators.singlepoint import SinglePointCalculator


nomad_api_template = ('https://labdev-nomad.esc.rzg.mpg.de/'
                      'api/resolve/{hash}?format=recursiveJson')


def nmd2https(uri):
    assert uri.startswith('nmd://')
    return nomad_api_template.format(hash=uri[6:])


def nmd2dict(uri):
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen

    httpsuri = nmd2https(uri)
    response = urlopen(httpsuri)
    txt = response.read().decode('utf8')
    return json.loads(txt, object_hook=lambda dct: NomadEntry(dct))


def read(fd):
    dct = json.load(fd, object_hook=lambda dct: NomadEntry(dct))
    return dct


def download(uri):
    # Might want to look/return sections also
    dct = nmd2dict(uri)
    return NomadEntry(dct)


def section_method2metadata(method, methods, metainfo=None):
    # Collect all information starting from reference method
    if not metainfo:
        metainfo = {}
    xc_funcs = method.get('section_XC_functionals', [])
    if xc_funcs:
        xc_info = ','.join([
            xc_func['XC_functional_name'] for xc_func in xc_funcs])
        if 'nomad_XC_functionals' in metainfo:
            metainfo['nomad_XC_functionals'] = metainfo['nomad_XC_functionals'] + ',' + xc_info
        else:
            metainfo['nomad_XC_functionals'] = xc_info
    e_calc_method = method.get('electronic_structure_method', [])
    if e_calc_method:
        metainfo['nomad_electronic_structure_method'] = e_calc_method
    ref_methods = method.get('section_method_to_method_refs', [])
    if ref_methods:
        for ref_method in ref_methods:
            ref_id = ref_method.get('method_to_method_ref', [])
            if ref_id:
                metainfo.update(section_method2metadata(
                    methods[ref_id], methods, metainfo=metainfo))
    return metainfo


def add_nomad_metainfo(d, run, calc, system=[]):
    # More nomad metainfo can be add to key_value_pairs and 
    # key_value_pairs can also be stored at ASE db.
    info = {}
    info['nomad_metadata_type'] = run['type']
    info['nomad_run_gIndex'] = run['gIndex']
    if system:
        info['nomad_uri'] = system['uri']
        info['nomad_system_gIndex'] = system['gIndex']
    info['nomad_calculation_uri'] = d['uri']
    info['nomad_program_name'] = run['program_name']
    if 'program_version' in run:
        info['nomad_program_version'] = ' '.join(run['program_version'].split())
    if 'energy_total_T0' in calc:
        info['potential_energy'] = calc['energy_total_T0'] * units.J
    if 'energy_total' in calc:
        info['nomad_total_energy'] = calc['energy_total'] * units.J
        info['energy'] = calc['energy_total'] * units.J
    if 'energy_free' in calc:
        info['free_energy'] = calc['energy_free'] * units.J
    if 'single_configuration_calculation_converged' in calc:
        info['nomad_converged'] = calc['single_configuration_calculation_converged']
    # Checking the reference section_method for this calc, 
    # section_single_configuration_calculation
    ref_method = calc.get('single_configuration_to_calculation_method_ref') 
    methods = run.get('section_method', [])
    if methods:
        if ref_method is not None:
            info.update(section_method2metadata(
                methods[ref_method], 
                methods))
    # ?? In case there is no reference to section_method,
    # ?? can we assume section_method(s) is(are) nested in 
    # ?? section_single_configuration_calculation
    else:
        methods = calc.get('section_method', [])
        if methods:
            for method in methods:
                info.update(section_method2metadata(
                    method, 
                    methods))
    return info


def dict2images(d, only_atoms=False):
    # Check if server return with error or json file has error field.
    assert 'error' not in d, 'Request return with following error: ' + d['error']
    runs = d.get('section_run', [])
    assert 'section_run' in d, 'Missing section_run!'
    single_confs={}
    for run in runs:
        calculations = run.get('section_single_configuration_calculation', [])
        systems = run.get('section_system', [])
        if not only_atoms:
            assert 'section_system' in run, 'No section_system in section_run!'
        for nmd_calc in calculations:
            system_ref = nmd_calc.get('single_configuration_calculation_to_system_ref', -1)
            # if single calculation w/o system, the system ref is -1
            single_confs[run.get('gIndex'), system_ref] = nmd_calc
            nmd_system = []
            if systems and system_ref > -1:
                nmd_system = systems[system_ref]
            metainfo = add_nomad_metainfo(d, run, nmd_calc, nmd_system)
            if not only_atoms:
                calc = SinglePointCalculator(**metainfo)
            if not nmd_system: yield calc
            if 'atom_positions' not in nmd_system: yield calc
            atoms = section_system2atoms(nmd_system)
            if not only_atoms:
                calc.atoms = atoms.copy()
                yield calc
            else:
                info = atoms.info.get('key_value_pairs', {})
                info.update(metainfo)
                atoms.info['key_value_pairs'] = info
                yield atoms


def calcs2atoms(dct):
    for calc in list(dict2images(dct, 
                     only_atoms=dct.only_atoms)):
        if calc.atoms is not None:
            atm = calc.atoms.copy()
            atm.info['key_value_pairs'] = calc.results
            yield atm


class NomadEntry(dict):
    def __init__(self, dct, only_atoms=False):
        #assert dct['type'] == 'nomad_calculation_2_0'
        #assert dct['name'] == 'calculation_context'
        # We could implement NomadEntries that represent sections.
        dict.__init__(self, dct)
        self.only_atoms = only_atoms

    @property
    def hash(self):
        # The hash is a string, so not __hash__
        assert self['uri'].startswith('nmd://')
        return self['uri'][6:]

    def toatoms(self):
        if not self.only_atoms:
            return calcs2atoms(self)
        else:
            return list(dict2images(self, 
                        only_atoms=self.only_atoms))

    def iterimages(self):
        return dict2images(self, 
                only_atoms=self.only_atoms)


def section_system2atoms(section):
    assert section['name'] == 'section_system'
    numbers = section['atom_species']
    numbers = np.array(numbers, int)
    numbers[numbers < 0] = 0
    numbers[numbers > len(chemical_symbols)] = 0
    positions = section['atom_positions']['flatData']
    positions = np.array(positions).reshape(-1, 3) * units.m
    pbc = section.get('configuration_periodic_dimensions')
    cell = section.get('lattice_vectors')
    atoms = Atoms(numbers, positions=positions)
    if pbc is not None:
        assert len(pbc) == 1
        pbc = pbc[0]  # it's a list??
        pbc = pbc['flatData']
        assert len(pbc) == 3
        atoms.pbc = pbc

    # celldisp?
    if cell is not None:
        cell = cell['flatData']
        cell = np.array(cell).reshape(3, 3) * units.m
        atoms.cell = cell

    return atoms


def section_singleconfig2calc(metainfo):
    # Forces, total energy, ........
    # We should be able to extract e.g. a band structure as well.
    kwargs = {}
    if 'energy' in metainfo:
        kwargs['energy']=metainfo['energy']
    calc = SinglePointCalculator(**kwargs)
    return calc


def main(argv):
    uri = "nmd://N9Jqc1y-Bzf7sI1R9qhyyyoIosJDs/C74RJltyQeM9_WFuJYO49AR4gKuJ2"
    # TRY THESE ALSO:
    # nmd://N9GXsYFhz9mUMxVWxlvwnn7mCz3j1/CtiXMrvFRdyQV4fOxJkh8NwoCZR6Z
    # nmd://N9GXsYFhz9mUMxVWxlvwnn7mCz3j1/CHYLKLtXXU7w7VTzesEaWibL3_A7O 
    # nmd://NWApItBGtGUDsfMVlHKqrjUQ4rShT/C-1SH_T1kd13-U3MEB7Xz-_eToBHT
    if len(argv)>0:
        uri = argv[0]
    only_atoms = True if len(argv)>1 else False
    if uri.startswith('nmd://'):
        print(nmd2https(uri))
        entry = download(uri)
        nmd_images = entry.toatoms()
    else:
        with open(uri) as fd:
            nmd_images = nomad_json.read_nomad_json(fd, only_atoms=only_atoms)
    from ase.visualize import view
    nmd_atoms = []
    for image in nmd_images:
        if isinstance(image, SinglePointCalculator):
            print(image)
            if image.atoms:
                nmd_atoms.append(image.atoms)
        else:
            print(image)
            if image.info:
                print(image.info['key_value_pairs'])
                nmd_atoms.append(image)
    view(nmd_atoms)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
