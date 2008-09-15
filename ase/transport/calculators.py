from ase.transport.transmission import Transmission
from ase.transport.selfenergy import LeadSelfEnergy
from ase.transport.greenfunction import GreenFunction
from ase.transport.tools import subdiagonalize, cutcoupling, tri2full, dagger
import numpy as npy

class TransportCalculator:
    """Determine transport properties of device sandwiched between
    semi-infinite leads using nonequillibrium Green function methods.
    """

    def __init__(self, **kwargs):
        """args=(energies, h, h1, h2, s=None, s1=None, s2=None, align_bf=None)
        
        energies is the energy grid on which the transport properties should
        be determined.
        
        h1 (h2) is a matrix representation of the Hamiltonian of two
        principal layers of the left (right) lead, and the coupling between
        such layers.
        
        h is a matrix representation of the Hamiltonian of the scattering
        region. This must include at least on lead principal layer on each
        side. The coupling in (out) of the scattering region is by 
        default assumed to be identical to the coupling between left (right) 
        principal layers. However, these couplings can also be
        specified explicitly through hc1 and hc2. 
        
        s, s1, and s2 are the overlap matrices corresponding to h, h1, and
        h2. Default is the identity operator. sc1 and sc2 are the overlap
        matrices corresponding to the optional couplings hc1 and hc2.
        
        align_bf specifies the principal layer basis index used to
        align the fermi levels of the lead and scattering regions.
        """
        
        self.input_parameters = {'energies': None,
                                 'h': None,
                                 'h1': None,
                                 'h2': None,
                                 's': None,
                                 's1': None,
                                 's2': None,
                                 'hc1' : None,
                                 'hc2' : None,
                                 'sc1' : None,
                                 'sc2' :None,
                                 'align_bf': None,
                                 'eta1': 1e-3,
                                 'eta2': 1e-3,
                                 'eta': 1e-3,
                                 'logfile': None,
                                 'verbose': False}

        self.trans = Transmission()
        self.initialized =  False
        self.set(**kwargs)

    def set(self, **kwargs):
        p = self.input_parameters
        p.update(kwargs)

        self.pl1 = len(p['h1']) / 2
        self.pl2 = len(p['h2']) / 2
        assert p['h'] is not None
        
        if 'pdos' in p:
            self.trans.set(pdos=p['pdos'])
        
        self.initialize()

    def initialize(self):
        p = self.input_parameters
        self.verbose = p['verbose']
        if self.verbose:
            print "initializing calculator..."
        self.energies = p['energies']
        pl1 = self.pl1
        pl2 = self.pl2
        if p['s1'] == None:
            p['s1'] = npy.identity(len(p['h1']))
        if p['s2'] == None:
            p['s2'] = npy.identity(len(p['h2']))
        if p['s'] == None:
            p['s'] = npy.identity(len(p['h']))
            
        h1_ii = p['h1'][:pl1, :pl1]
        h1_ij = p['h1'][:pl1, pl1:2*pl1]
        s1_ii = p['s1'][:pl1, :pl1]
        s1_ij = p['s1'][:pl1, pl1:2*pl1]

        h2_ii = p['h2'][:pl2,:pl2]
        h2_ij = p['h2'][pl2:2*pl2,:pl2]
        s2_ii = p['s2'][:pl2,:pl2]
        s2_ij = p['s2'][pl2:2*pl2,:pl2]

        #h1_im = p['h'][:pl1,pl1:-pl2]
        #s1_im = p['s'][:pl1,pl1:-pl2]

        #h2_im = p['h'][-pl1:,pl1:-pl2]
        #s2_im = p['s'][-pl1:,pl1:-pl2]

        #h_mm = p['h'][pl1:-pl2,pl1:-pl2]
        #s_mm = p['s'][pl1:-pl2,pl1:-pl2]
        
        h_mm = p['h']
        s_mm = p['s']
        
        if p['hc1'] is None:
            nbf = len(h_mm)
            h1_im = npy.zeros((pl1, nbf), complex)
            s1_im = npy.zeros((pl1, nbf), complex)
            h1_im[:pl1, :pl1] = h1_ij
            s1_im[:pl1, :pl1] = s1_ij
        else:
            h1_im = p['hc1']
            if p['sc1'] is not None:
                s1_im = p['sc1']
            else:
                s1_im = npy.zeros(h1_im.shape, complex)

        if p['hc2'] is None:
            h2_im = npy.zeros((pl2, nbf), complex)
            s2_im = npy.zeros((pl2, nbf), complex)
            h2_im[-pl2:, -pl2:] = h2_ij
            s2_im[-pl2:, -pl2:] = s2_ij
        else:
            h2_im = p['hc2']
            if p['sc2'] is not None:
                s2_im[:] = p['sc2']
            else:
                s2_im = npy.zeros(h2_im.shape, complex)

        align_bf = p['align_bf']
        if align_bf != None:
            diff = (h_mm[align_bf, align_bf] - h1_ii[align_bf, align_bf]) \
                   / s_mm[align_bf, align_bf]
            if self.verbose:
                print "Alligning scattering H to left lead H. diff=", diff
            h_mm -= diff * s_mm

        self.h_pp = p['h']
        self.s_pp = p['s']
        #setup lead self-energies
        sigma1 = LeadSelfEnergy((h1_ii, s1_ii), 
                                (h1_ij, s1_ij),
                                (h1_im, s1_im),
                                p['eta1'])
        
        sigma2 = LeadSelfEnergy((h2_ii, s2_ii), 
                                (h2_ij, s2_ij),
                                (h2_im, s2_im),
                                p['eta2'])

        self.selfenergies = [sigma1, sigma2]
        #setup scattering green function
        self.gf = GreenFunction(selfenergies=self.selfenergies,
                                H=h_mm,
                                S=s_mm,
                                eta = p['eta'])

        #setup the basic transmission calculator 
        self.trans.set(energies=self.energies,
                       greenfunction=self.gf,
                       selfenergies=self.selfenergies,
                       logfile=p['logfile'],
                       verbose=p['verbose'])
        
        self.initialized = True

    def print_pl_convergence(self):
        p = self.input_parameters
        pl1 = self.pl1
        
        h_ii = self.selfenergies[0].h_ii
        s_ii = self.selfenergies[0].s_ii
        ha_ii = self.gf.H[:pl1, :pl1]
        sa_ii = self.gf.S[:pl1, :pl1]
        c1 = npy.abs(h_ii - ha_ii).max()
        c2 = npy.abs(s_ii - sa_ii).max()
        print "Conv (h,s)=%.2e, %2.e" % (c1, c2)

    def get_transmission(self):
        if not self.trans.uptodate:
            self.trans.update()
        return self.trans.T_e

    def get_dos(self):
        if not self.trans.uptodate:
            self.set(dos=True)
            self.trans.update()
        return self.trans.dos_e

    def get_eigenchannels(self,n=0):
        if n > self.trans.input_parameters['eigenchannels']:
            self.trans.set(eigenchannels=n)
            self.trans.initialize()
            self.trans.update()
        elif not self.trans.uptodate:
            self.trans.update()
        if n==0:
            n = self.trans.input_parameters['eigenchannels']
        return self.trans.eigenchannels_ne[:n]

    def get_pdos(self):
        if not self.trans.uptodate:
            self.trans.update()
        return self.trans.pdos_ne

    def subdiagonalize_bfs(self, bfs):
        bfs = npy.array(bfs)
        p = self.input_parameters
        h_pp = p['h']
        s_pp = p['s']
        ht_pp, st_pp, c_pp, e_p = subdiagonalize(h_pp, s_pp, bfs)
        c_pp = npy.take(c_pp, bfs, axis=0)
        c_pp = npy.take(c_pp, bfs, axis=1)
        return ht_pp, st_pp, e_p, c_pp

    def cutcoupling_bfs(self, bfs):
        bfs = npy.array(bfs)
        p = self.input_parameters
        h_pp = p['h'].copy()
        s_pp = p['s'].copy()
        cutcoupling(h_pp, s_pp, bfs)
        return h_pp, s_pp
        
    def get_left_channels(self, energy, n):
        return self.trans.get_left_channels(energy, n)
