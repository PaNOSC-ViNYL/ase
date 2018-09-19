import json
import numpy as np
import sys
import re
import requests

import ase.units as units
from ase import Atoms
from ase.atoms import symbols2numbers
from ase.io import nomad_json
from ase.io import nomad_ziptxt
from ase.data import chemical_symbols
from ase.calculators.singlepoint import SinglePointCalculator

if sys.version_info.major > 2:
    # For Python 3.0 and later
    from urllib.request import urlopen, Request
    from urllib.parse import quote, unquote_plus
else:
    # Fall back to Python 2's urllib2
    from urllib import quote, unquote_plus
    from urllib2 import urlopen, Request

nomad_api_url = 'https://labdev-nomad.esc.rzg.mpg.de'
nomad_query_url = 'https://analytics-toolkit.nomad-coe.eu'
nomad_api_template = (nomad_api_url + '/api/resolve/{hash}?format=recursiveJson')
nomad_nql_api_query_template = (nomad_api_url + '/dev/archive/nql-api/search?query={hash}')
# The next link for queries will be DEPRECATED from NOMAD!
nomad_api_query_template = (nomad_query_url + '/api/query/section_run?filter={hash}')

nomad_enc_url = 'https://encyclopedia.nomad-coe.eu/api/v1.0/materials'
nomad_enc_saml = 'https://encyclopedia.nomad-coe.eu/api/v1.0/saml/'
nomad_enc_calc_template = (nomad_enc_url + '/{}/calculations?pagination=off')
nomad_enc_sgrp_template = (nomad_enc_url + '/{}?property=space_group')
nomad_enc_cell_template = (nomad_enc_url + '/{}/cells')
nomad_enc_elmt_template = (nomad_enc_url + '/{}/elements')


def nmd2https(uri):
    assert uri.startswith('nmd://')
    return nomad_api_template.format(hash=uri[6:])


def nmd2dict(uri):
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen

    httpsuri = nmd2https(uri)
    try:
        response = urlopen(httpsuri)
        txt = response.read().decode('utf8')
        return json.loads(txt, object_hook=lambda dct: NomadEntry(dct))
    except Exception as exc:
        exc = sys.exc_info()[1]
        print('NOMAD Server ERROR: ' + str(exc))
        return dict()

def read(fd):
    dct = json.load(fd, object_hook=lambda dct: NomadEntry(dct))
    return dct


def download(uri, only_atoms=False, skip_errors=False):
    # Might want to look/return sections also
    dct = nmd2dict(uri)
    return NomadEntry(dct, only_atoms=only_atoms, skip_errors=skip_errors)


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
    if 'program_name' in run:
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
            try:
                info.update(section_method2metadata(
                    methods[ref_method], 
                    methods))
            except IndexError:
                pass
            return info
    # ?? In case there is no reference to section_method,
    # ?? can we assume section_method(s) is(are) nested in 
    # ?? section_single_configuration_calculation
    methods = calc.get('section_method', [])
    if methods:
        for method in methods:
            info.update(section_method2metadata(
                method, 
                methods))
    return info


def dict2images(d, only_atoms=False, skip_errors=False):
    extracted_systems = []
    if 'error' in d:
        if not skip_errors:
            # Check if server return with error or json file has error field.
            assert 'error' not in d, 'Request return with following error: ' + d['error']
    else:
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
                    # if the system is already written in the image list
                    # we can skip this system_ref link and just add calculation info
                    if system_ref not in extracted_systems:
                        try:
                            nmd_system = systems[system_ref]
                            if system_ref not in extracted_systems:
                                extracted_systems.append(system_ref)
                        except IndexError:
                            pass
                metainfo = add_nomad_metainfo(d, run, nmd_calc, nmd_system)
                calc = SinglePointCalculator(**metainfo)
                if not nmd_system: yield calc
                atom_pos_true = None
                try:
                    atom_pos_true = nmd_system['atom_positions']
                except (TypeError, KeyError):
                    yield calc
                if atom_pos_true is None:
                    yield calc
                else:
                    atoms = section_system2atoms(nmd_system)
                    if atoms is None:
                        if not only_atoms:
                            yield calc
                    else:
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
                     only_atoms=dct.only_atoms,
                     skip_errors=dct.skip_errors)):
        if calc.atoms is not None:
            atm = calc.atoms.copy()
            atm.info['key_value_pairs'] = calc.results
            yield atm


class NomadEntry(dict):
    def __init__(self, dct, only_atoms=False, skip_errors=False):
        #assert dct['type'] == 'nomad_calculation_2_0'
        #assert dct['name'] == 'calculation_context'
        # We could implement NomadEntries that represent sections.
        dict.__init__(self, dct)
        self.only_atoms = only_atoms
        self.skip_errors = skip_errors

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
                        only_atoms=self.only_atoms,
                        skip_errors=self.skip_errors))

    def iterimages(self):
        return dict2images(self, 
                only_atoms=self.only_atoms,
                skip_errors=self.skip_errors)


def section_system2atoms(section):
    #assert section['name'] == 'section_system'
    numbers = None
    if 'atom_species' in section:
        numbers = section['atom_species']
        numbers = np.array(numbers, int)
        numbers[numbers < 0] = 0
        numbers[numbers > len(chemical_symbols)] = 0
    if 'atom_positions' not in section:
        return None
    else:
        positions = section['atom_positions']['flatData']
        positions = np.array(positions).reshape(-1, 3) * units.m
        pbc = section.get('configuration_periodic_dimensions')
        cell = section.get('lattice_vectors')
        if numbers is not None:
            atoms = Atoms(numbers, positions=positions)
        else:
            atoms = Atoms(positions=positions)
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


def section_singleconfig2calc(dct, nmd_run, nmd_calc, nmd_system):
    # Forces, total energy, ........
    # We should be able to extract e.g. a band structure as well.
    kwargs = add_nomad_metainfo(dct, nmd_run, nmd_calc, nmd_system)
    calc = SinglePointCalculator(**kwargs)
    return calc


class NomadQuery(object):
    """
        NOMAD Query class. Requests archive info from NOMAD servers.

        Parameters:
        ===========
        atom_labels = String that includes atom element symbols seperated with commas.
        query = Raw query string for general usage
        nomad_interface = 'a' or 'archive' (Default NOMAD archive with NQL API)
                          'o' or 'old-archive' (Uses same NOMAD archive but with old API)
                          'e' or 'encyclopeadia' (Access to NOMAD encyclopeadia via its API)
        space_group = Integer. Supply with nomad_interface = e
        program_name = String to specify name of Ab-initio/MD program.
        exclusive = True to get archives only with specified atom symbols, set it to 
                    False if you would like to get all archive data that inlcudes these symbols.

        Returns:
        ========
        NomadQuery object.

        Methods:
        ========
        query() : To start a new query with given parameters.
        download() : Download archive data for nmd:// sequences that are
                     retrieved from query().
        save_db() : Save downloaded ASE.Atoms list to ASE.db file.

    """
    def __init__(self, atom_labels=None, query=None, 
                 nomad_interface='archive', nomad_token='',
                 space_group=None, program_name='', 
                 exclusive=True, number_of_results=None):
        self.response = None
        if nomad_interface.startswith('e'):
            self.nomad_interface = 'e'
        elif nomad_interface.startswith('o'):
            self.nomad_interface = 'o'
        else:
            self.nomad_interface = 'a'
        if nomad_token != '':
            self.auth = (nomad_token, '')
        else:
            if nomad_interface.startswith('e'):
                try:
                    response = requests.get(
                        nomad_enc_saml + '/user/',
                        verify=False
                        )
                    response = response.json()
                except Exception as exc:
                    response = self._handle_error(exc)
                if 'status' in response:
                    if 'Unauthenticated' in response['status']:
                        print("Your NOMAD Encyclopedia session is not authenticated."
                              " Type '" + nomad_enc_saml + "' in your browser.")
                if "token" in response:
                    nomad_token = response["token"]["data"]
                    self.auth = (nomad_token, '')
                else:
                    self.auth = ('', '')
            else:
                self.auth = ('', '')
        init_request = False
        if isinstance(atom_labels, (tuple, list, np.ndarray)):
            init_request = True
        elif isinstance(atom_labels, str):
            if ',' in atom_labels:
                init_request = True
            else:
                init_request = True
        else:
            if query is not None:
                if isinstance(query, str):
                    init_request = True
        if init_request:
            self.request(atom_labels=atom_labels, query=query, 
                nomad_interface=self.nomad_interface, nomad_token=self.auth[0],
                space_group=space_group, program_name=program_name, 
                exclusive=exclusive, number_of_results=number_of_results)

    def _reset_query(self):
        self.query = None
        self.number_of_results = None
        self.response = {}
        self.atom_labels = []

    def request(self, atom_labels=None, query=None, 
                nomad_interface='a', nomad_token='',
                space_group=None, program_name='', 
                exclusive=True, number_of_results=None):
        assert (atom_labels is not None or 
                query is not None), 'One of atom_labels or query should be given for NOMAD request.'

        self._reset_query()

        if isinstance(atom_labels, (tuple, list, np.ndarray)):
            self.atom_labels = atom_labels
        elif isinstance(atom_labels, str):
            if ',' in atom_labels:
                self.atom_labels = [c for c in atom_labels.split(',')]
            else:
                self.atom_labels = [atom_labels]

        if nomad_interface.startswith('e'):
            # encyclopedia
            self.nomad_interface = 'e'
        elif nomad_interface.startswith('o'):
            # old API for archive
            self.nomad_interface = 'o'
        else:
            # NQL API for archive
            self.nomad_interface = 'a'

        self.exclusive = exclusive

        if nomad_token != '':
            self.auth = (nomad_token, '')
        else:
            self.auth = ('', '')

        if isinstance(query, str):
            if len(query)>0:
                self.query = query 
            else:
                self.query = None

        if self.nomad_interface == 'o':
            if self.query is None:
                the_query = ''
                if self.atom_labels:
                    the_query = the_query + 'atom_symbols:' + ','.join(
                            [str(at) for at in self.atom_labels])
                if program_name != '':
                    the_query = the_query + ' AND ' + 'program_name:' + str(program_name)
                self.query = the_query
        elif self.nomad_interface == 'e':
            self.query = {
                    "search_by": {
                        "exclusive" : "0" if self.exclusive is False else "1",
                        "element": ','.join([at for at in self.atom_labels]),
                        "pagination": "off"
                        },
                    }
            if space_group:
                self.query["space_group"] = str(space_group)
        else:
            if self.query is None:
                the_query = ''
                if self.atom_labels:
                    the_query = the_query + 'all atom_species=' + ','.join([
                        str(num) for num in symbols2numbers(self.atom_labels)])
                if program_name != '':
                    the_query = the_query + ' and ' + 'program_name:' + str(program_name)
                if number_of_results is not None:
                    the_query = the_query + "&num_results="+str(int(number_of_results))
                self.query = the_query + '&sort_field=calculation_gid&format=nested'

        response = self._request()

        self.response = {}
        if self.nomad_interface == 'o':
            # Need LICENSE for this part ?
            # See accessing json file with URI list at 
            # https://analytics-toolkit.nomad-coe.eu/notebook-edit/data/shared/
            #     tutorialsNew/nomad-query/nomad-query.bkr
            download_path = response.get('path', '')
            nmd_uri_list = []
            if download_path:
                download_url = ''.join([
                    nomad_query_url + '/api/download/',
                    download_path.split('_')[-1],
                    '?file=', quote(download_path.encode('utf-8'))
                    ])
                download_json = download_url + '.json'
                print(download_json)
                json_file_request = self._request(url=download_json)
                if json_file_request['status'] == 'success':
                    regex = re.compile(r'(?<=/[a-zA-Z0-9\-_]{3}/)[^\.]+')
                    for uri_path in json_file_request['result']:
                        match = regex.search(uri_path)
                        if match:
                            # Substitute prefixes
                            groups = match.group(0).split('/')
                            groups[0] = 'N' + groups[0][1:]         # Normalized

                            if len(groups) == 2:
                                groups[1] = 'C' + groups[1][1:]     # Computed

                        nmd_uri_list.append('nmd://' + '/'.join(groups))
            self.response['data'] = response
            self.response['info'] = response['info']
            self.response['nmd_uri_list'] = nmd_uri_list
        elif self.nomad_interface == 'e':
            if 'result' in response and response['status'] != 'error':
                nomad_metarial_id = response["result"][0]['id']
                self.response['material_id'] = nomad_metarial_id
                try:
                    nomadmatdata = requests.get(
                            url=nomad_enc_calc_template.format(nomad_metarial_id),
                            auth=self.auth
                            )
                    nomadmatdata = nomadmatdata.json()
                    self.response['data'] = nomadmatdata['results']
                except Exception as exc:
                    nomadmatdata = self._handle_error(exc)
                    print(nomadmatdata)
                try:
                    nomadsgrpdata = requests.get(
                            url=nomad_enc_sgrp_template.format(nomad_metarial_id),
                            auth=self.auth
                            )
                    nomadsgrpdata = nomadsgrpdata.json()
                    self.response['space_group'] = nomadsgrpdata['space_group']
                except Exception as exc:
                    nomadsgrpdata = self._handle_error(exc)
                    print(nomadsgrpdata)
            else:
                print(response)
        else:
            nmd_uri_list = []
            data_list = response.get('result', [])
            for query_data in data_list:
                qdat = query_data["attributes"]["metadata"]
                archive_gid = qdat["archive_context"]['archive_gid'][0].replace('R','N',1)
                for cgid in qdat["calculation_context"]['calculation_gid']:
                    nmd_uri_list.append('nmd://' + str(archive_gid) + '/' + str(cgid))
            self.response['data'] = response
            self.response['info'] = response['info']
            self.response['nmd_uri_list'] = nmd_uri_list
            print('nmd_uri_list',len(nmd_uri_list))

    def _handle_error(self, exc):
        # Need LICENSE for this part ?
        # See handling of errors from NOMAD API at
        # https://analytics-toolkit.nomad-coe.eu/notebook-edit/data/shared/
        #     tutorialsNew/nomad-query/nomad-query.bkr
        error = {
                'status'  : 'error',
                }
        if self.nomad_interface == 'a':
            error['message'] = 'Unknown error for NOMAD Archive API.'
        elif self.nomad_interface == 'o':
            error['message'] = 'Unknown error for NOMAD Archive API.'
        else:
            error['message'] = 'Unknown error for NOMAD Encyclopedia API.'
        exc = sys.exc_info()[1]
                
        # Get error message
        message = exc
        if sys.version_info <= (2, 5) and hasattr(exc, 'message'):
            message = exc.message
        elif hasattr(exc, 'reason'):
            message = exc.reason
        error['message'] = str(message)
                
        # Fix error message
        if error['message'].endswith('timed out'):
            error['message'] = 'Connection timed out for NOMAD ' + \
                               'API Service. Service may currently unavailable.'
        return error

    def _handle_response(self, response, interface):
        # Need LICENSE for this part ?
        # See handling of responses from NOMAD APIs at
        # https://analytics-toolkit.nomad-coe.eu/notebook-edit/data/shared/
        #     tutorialsNew/nomad-query/nomad-query.bkr
        # and
        # https://encyclopedia.nomad-coe.eu/doc/
        error = {
                'status'  : 'error',
                }

        if interface == 'a':
            error['message'] = 'Unknown error for NOMAD Archive NQL API.'
        elif interface == 'o':
            error['message'] = 'Unknown error for NOMAD Archive API.'
        else:
            error['message'] = 'Unknown error for NOMAD Encyclopedia API.'

        # Extracts the response and its code
        if not isinstance(response, dict):
            if response.code != 200:
                response = error
            else:
                data = json.loads(response.read().decode('utf-8'))
                if 'error' in data:
                    response = error.copy()
                    response['message'] = 'Error:' + data['error'] + '.'
                elif 'errors' in data:
                    response = error.copy()
                    if isinstance(data['errors'], dict):
                        response['message'] = 'Error:' + '.'.join([
                            str(k)+str(v) for k,v in data['errors'].items()]) + '.'
                    elif isinstance(data['errors'], list):
                        response['message'] = 'Error:' + '.'.join(['.'.join([
                            str(k)+str(v) for k,v in errs.items()
                            ]) for errs in data['errors']]) + '.'
                    else:
                        response['message'] = 'Error:' + data['errors'] + '.'
                elif 'msg' in data:
                    response = error.copy()
                    response['message'] = response['message'] + ' ' + data['msg']
                elif 'message' in data:
                    response = error.copy()
                    response['message'] = response['message'] + ' ' + data['message']
                else:
                    # Everthing was ok
                    message = ''
                    status = 'success'
                        
                    if interface == 'o':
                        # Get status from backend
                        if data['status'] != 'success':
                            status = data['status']
                            message = data.get('message', '')
                            
                        # Check connection timeout
                        elif('timed_out' in data['result'] and
                                data['result']['timed_out']):
                            status = 'error'
                            message = 'Connection timed out.'
                        
                        # Construct response
                        response = {
                                'status': status,
                                'message': message,
                                'result': (data['result'] if status == 'success'
                                    else []),
                                'info': (data['result']['aggregations'] if 'aggregations' in data['result']
                                    else {}),
                                'path': data.get('path', '')
                                }
                    elif interface == 'e':
                        # Construct response
                        response = {
                                'status': status,
                                'message': message,
                                'result': (data['results'] if status == 'success'
                                    else []),
                                'path': data.get('login_url', '')
                                }
                    else:
                        # Construct response
                        response = {
                                'status': status,
                                'message': message,
                                'result': (data['data'] if status == 'success'
                                    else []),
                                'info': (data['meta'] if status == 'success'
                                    else ''),
                                'path': data.get('path', '')
                                }
        else:
            data = response
            if 'error' in data or 'msg' in data or 'message' in data:
                response = error.copy()
                response['login_url'] = data.get('login_url', '')
            if 'error' in data:
                response['message'] = 'Error:' + data['error'] + '.'
            elif 'errors' in data:
                    response = error.copy()
                    if isinstance(data['errors'], dict):
                        response['message'] = 'Error:' + '.'.join([
                            str(k)+str(v) for k,v in data['errors'].items()]) + '.'
                    elif isinstance(data['errors'], list):
                        response['message'] = 'Error:' + '.'.join(['.'.join([
                            str(k)+str(v) for k,v in errs.items()
                            ]) for errs in data['errors']]) + '.'
                    else:
                        response['message'] = 'Error:' + data['errors'] + '.'
            elif 'msg' in data:
                response['message'] = response['message'] + ' ' + data['msg']
            elif 'message' in data:
                response['message'] = response['message'] + ' ' + data['message']
            else:
                # Everthing was ok
                message = ''
                status = 'success'
                
                # Construct response
                if interface == 'e':
                    response = {
                            'status' : status,
                            'message': message,
                            'result' : data.get('results',[]),
                            'info'   : data.get('meta', ''),
                            'path'   : data.get('url', '')
                            }
                else:
                    response = {
                            'status' : status,
                            'message': message,
                            'result' : data.get('data', []),
                            'info'   : data.get('meta', ''),
                            'path'   : data.get('path', '')
                            }
        return response


    def _request(self, url=None, timeout=30):
        response = { 'status' : 'error' }
        # Resuest from NOMAD Archive API
        if url is not None:
            if self.nomad_interface == 'o':
                try:
                    response = urlopen(Request(url), timeout=timeout)
                except Exception as exc:
                    response = self._handle_error(exc)
                    #print(response)
                response = self._handle_response(response, interface='o')
            else:
                try:
                    response = requests.get(url)
                    response = response.json()
                except Exception as exc:
                    response = self._handle_error(exc)
                    print(response)
                response = self._handle_response(response, interface='a')
        else:
            if self.nomad_interface == 'o':
                url = nomad_api_query_template.format(hash=self.query)
                print(url)
                # Sends the request and catches the response
                try:
                    response = urlopen(Request(url), timeout=timeout)
                except Exception as exc:
                    response = self._handle_error(exc)
                    print(response)
            elif self.nomad_interface == 'e':
                if isinstance(self.query, dict):
                    req_json = self.query
                else:
                    req_json = {}
                # Sends the request and catches the response
                try:
                    response = requests.post(
                            url=nomad_enc_url,
                            json=req_json,
                            auth=self.auth
                            )
                    response = response.json()
                except Exception as exc:
                    response = self._handle_error(exc)
                    print(response)
            else:
                url = nomad_nql_api_query_template.format(hash=self.query)
                # Sends the request and catches the response
                try:
                    response = requests.get(url)
                    response = response.json()
                except Exception as exc:
                    response = self._handle_error(exc)
                    print(response)
            response = self._handle_response(response, interface=self.nomad_interface)
        return response


    def download(self, index=':', only_atoms=False,
                 skip_errors=False, no_dublicate_positions=False,
                 only_last_entries=False):
        images = []
        from ase.utils import basestring
        from ase.io.formats import string2index
        from math import pi
        from ase.spacegroup import crystal
        if isinstance(index, basestring):
            index = string2index(index)
        if(self.nomad_interface == 'a' or 
           self.nomad_interface == 'o'):
            for nmd_link in self.response["nmd_uri_list"]:
                ref_positions = None
                if 'section_run' in nmd_link:
                    nmd_uri = nmd_link.split('/section_run')
                    nmduri = nmd_uri[0]
                else:
                    nmduri = nmd_link
                print(nmduri)
                try:
                    entry = download(nmduri, only_atoms=only_atoms, skip_errors=skip_errors)
                    nmd_images_download = entry.toatoms()
                    nmd_images = list(nmd_images_download)
                except(TypeError, AttributeError, AssertionError, IndexError):
                    print(exc = sys.exc_info()[1])
                    nmd_images = []
                if self.atom_labels:
                    if self.exclusive:
                        nmd_images = [img for img in nmd_images if np.all(np.in1d(
                            img.get_chemical_symbols(), self.atom_labels
                            ))]
                if no_dublicate_positions:
                    nmd_new_images = []
                    if len(nmd_images)>0:
                        ref_positions = nmd_images[-1].positions
                        nmd_new_images.append(nmd_images[-1])
                        for img in reversed(nmd_images):
                            if np.linalg.norm(img.positions - ref_positions) > 0.:
                                nmd_new_images.append(img)
                                ref_positions = img.positions
                    nmd_images = nmd_new_images
                if only_last_entries:
                    formulas = list(set([','.join([str(i) for i in ni.numbers]) for ni in nmd_images]))
                    nmd_new_images = []
                    for formula in formulas:
                        for img in reversed(nmd_images):
                            if formula == ','.join([str(i) for i in img.numbers]):
                                nmd_new_images.append(img)
                                break
                    nmd_images = nmd_new_images
                #print(nmd_images)
                if len(nmd_images)>0:
                    print('Adding ' + str(len(nmd_images)) + ' structure(s) with ' + ','.join(
                        list(set([str(ni.get_chemical_formula('reduce')) for ni in nmd_images]))))
                else:
                    print('No structures is retrieved from this NOMAD archive!')
                images.extend(nmd_images)
        else:
            if 'data' in self.response:
                for result in self.response["data"]:
                    wyckoff_basis = []
                    nmd_elements = []
                    for item in result['wyckoff_groups_json']:
                        if item['variables']:
                            wyckoff_basis.append([
                                float(item['variables']['x']),
                                float(item['variables']['y']),
                                float(item['variables']['z'])])
                        if item['element']:
                            nmd_elements.append(str(item['element']))
                            lat_par = np.array([float(s) for s in result[
                                "lattice_parameters"][1:-1].split(",")])
                        nmd_cell = []
                        nmd_cell.extend(lat_par[0:3] * units.m)
                        nmd_cell.extend(lat_par[3:6] * 180.0/pi)
                        atoms = crystal(nmd_elements,
                                basis=wyckoff_basis,
                                spacegroup=int(self.response['space_group']),
                                cellpar=nmd_cell)
                        images.append(atoms)
        self.images = list(images)[index]

    def save_db(self, filename=None):
        self.db_file = filename if isinstance(
                filename, str) else None
        if self.db_file is None:
            self.db_file = 'nomad_asedb.db'
        if self.db_file and self.images:
            import ase.db
            with ase.db.connect(self.db_file) as atoms_db:
                for image in self.images:
                    if isinstance(image, Atoms):
                        atoms_db.write(image)
                        
    def __repr__(self):
        tokens = []
        
        if self.response:
            info = self.response.get('info', {})
            if info:
                if 'total_calculations' in info:
                    tokens.append("{0} calculations".format(
                        info['total_calculations']))
                if 'total_unique_geometries' in info:
                    tokens.append("{0} unique geometries".format(
                        info['total_unique_geometries']))
                if 'total_single_config_calculations' in info:
                    tokens.append("{0} configurations".format(
                        info['total_single_config_calculations']))
                if 'num_single_configurations' in info:
                    tokens.append("{0} configurations".format(
                        info['num_single_configurations']['value']))
            else:
                data = self.response.get('data', {})
                if 'status' in data:
                    if 'error' in data['status']:
                        if 'message' in data:
                            tokens.append("Status: {0}.".format(data.get('status', '')))
                            tokens.append("Message: {0}.".format(data.get('message', '')))

        return '{0}({1})'.format(self.__class__.__name__, ', '.join(tokens))


def main(argv):
    uri = "nmd://N9Jqc1y-Bzf7sI1R9qhyyyoIosJDs/C74RJltyQeM9_WFuJYO49AR4gKuJ2"
    # TRY THESE ALSO:
    # nmd://N9GXsYFhz9mUMxVWxlvwnn7mCz3j1/CtiXMrvFRdyQV4fOxJkh8NwoCZR6Z
    # nmd://N9GXsYFhz9mUMxVWxlvwnn7mCz3j1/CHYLKLtXXU7w7VTzesEaWibL3_A7O 
    # nmd://NWApItBGtGUDsfMVlHKqrjUQ4rShT/C-1SH_T1kd13-U3MEB7Xz-_eToBHT
    nmd_query = None
    if len(argv)>0:
        uri = argv[0]
    only_atoms = True if len(argv)>1 else False
    write_to_atomsdb = True if len(argv)>2 else False
    nmd_int = 'a'
    if len(argv)>3:
        nmd_int = 'e' if argv[3]=='e' else 'o'
    nmd_auth = argv[4] if len(argv)>4 else ''
    nmd_sgroup = int(argv[5]) if len(argv)>5 else None
    if uri.startswith('nmd://'):
        print(nmd2https(uri))
        entry = download(uri)
        nmd_images = entry.toatoms()
    else:
        if uri.endswith('.json'):
            with open(uri) as fd:
                nmd_images = nomad_json.read_nomad_json(fd, only_atoms=only_atoms)
        elif uri.endswith('.zip'):
            import zipfile
            zipfilelist = []
            ziptxts = []
            nmd_images = []
            with zipfile.ZipFile(uri) as zin:
                zipfilelist = zin.namelist()
            for zfile_name in zipfilelist:
                if zfile_name.startswith(uri.replace('.zip','')) and '.txt' in zfile_name:
                    ziptxts.append(zfile_name)
            with zipfile.ZipFile(uri) as zin:
                for zfile in ziptxts:
                    print('Found NMD txt file: ', str(zfile))
                    with zin.open(zfile, 'r') as fd:
                        nmd_txt_images = nomad_ziptxt.read_nomad_ziptxt(fd, 
                                only_atoms=only_atoms, skip_errors=True)
                        nmd_images.extend(nmd_txt_images)
        else:
            nmd_query = NomadQuery(atom_labels=uri, 
                    nomad_interface=nmd_int, 
                    nomad_token=nmd_auth,
                    space_group=nmd_sgroup)
            print(nmd_query)
            nmd_query.download(skip_errors=True, 
                    no_dublicate_positions=True,
                    only_last_entries=False)
    if write_to_atomsdb:
        if nmd_query:
            nmd_query.save_db(argv[2])
    else:
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
