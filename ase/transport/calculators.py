from ase.transport.transmission import Transmission
from ase.transport.selfenergy import LeadSelfEnergy
from ase.transport.greensfunction import GreensFunction
from ase.transport.tools import subdiagonalize, cutcoupling, tri2full, dagger
import numpy as npy

class TransportCalculator:

    def __init__(self, **kwargs):
        
        self.input_parameters = {'energies': None,
                                 'pl': None,
                                 'pl1' : None,
                                 'pl2' : None,
                                 'h1' : None,
                                 'h2' : None,
                                 's1' : None,
                                 's2' : None,
                                 'h' : None,
                                 's' : None,
                                 'align_bf' : None,
                                 'eta1' : 1.0e-3,
                                 'eta2' : 1.0e-3,
                                 'eta': 1.0e-3,
                                 'verbose' : False}

        self.trans = Transmission()
        self.initialized =  False
        self.set(**kwargs)
        #self.initialize()

    def set(self, **kwargs):
        p = self.input_parameters
        if 'pl' in kwargs: #using pl=pl1=pl2
            pl = kwargs['pl']
            p['pl1'] = pl
            p['pl2'] = pl
        p.update(kwargs)

        if 'pdos' in p:
            self.trans.set(pdos=p['pdos'])
        
        if p['h1'] != None and p['h2'] != None and p['h'] != None:
            self.initialize() #XXX more advanced log would be nicea

    def initialize(self):
        p = self.input_parameters
        self.verbose = p['verbose']
        if self.verbose:
            print "initializing calculator..."
        self.energies = p['energies']
        pl1 = p['pl1']
        pl2 = p['pl2']
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

        h1_im = p['h'][:pl1,pl1:-pl2]
        s1_im = p['s'][:pl1,pl1:-pl2]

        h2_im = p['h'][-pl1:,pl1:-pl2]
        s2_im = p['s'][-pl1:,pl1:-pl2]

        h_mm = p['h'][pl1:-pl2,pl1:-pl2]
        s_mm = p['s'][pl1:-pl2,pl1:-pl2]
        align_bf = p['align_bf']
        
        if align_bf != None:
            diff = (h_mm[align_bf,align_bf] - h1_ii[align_bf,align_bf]) / s_mm[align_bf,align_bf]
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
        #setup scattering green's function
        self.gf = GreensFunction(selfenergies=self.selfenergies,
                                 h_mm=h_mm,
                                 s_mm=s_mm,
                                 eta = p['eta'])
        self.gf.initialize()
        #setup the basic transmission calculator 
        self.trans.set(energies=self.energies,
                       greensfunction = self.gf,
                       selfenergies = self.selfenergies,
                       verbose=p['verbose'])
        
        self.initialized = True

    def print_pl_convergence(self):
        p = self.input_parameters
        pl1 = p['pl1']
        pl2 = p['pl2']
        for l in range(2):
            h_ii = self.selfenergies[l].h_ii
            s_ii = self.selfenergies[l].s_ii
            ha_ii = self.gf.h_mm[:pl1,:pl1]
            sa_ii = self.gf.s_mm[:pl1,:pl1]
            c1 = npy.abs(h_ii-ha_ii).max()
            c2 = npy.abs(s_ii-sa_ii).max()
            print "Conv %i: (h,s)=%.2e, %2.e" % (l,c1,c2)

    def get_transmission(self):
        if not self.trans.uptodate:
            self.trans.update()
        return self.trans.T_e

    def get_dos(self):
        if not self.trans.uptodate:
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
        bfs += p['pl1']
        h_pp = p['h']
        s_pp = p['s']
        ht_pp, st_pp, c_pp, e_p = subdiagonalize(h_pp, s_pp, bfs)
        c_pp = npy.take(c_pp,bfs,axis=0)
        c_pp = npy.take(c_pp,bfs,axis=1)
        return ht_pp, st_pp, e_p, c_pp

    def cutcoupling_bfs(self, bfs):
        bfs = npy.array(bfs)
        p = self.input_parameters
        bfs += p['pl1']
        h_pp = p['h'].copy()
        s_pp = p['s'].copy()
        cutcoupling(h_pp, s_pp, bfs)
        return h_pp, s_pp
        
    def get_left_channels(self,energy,n):
        return self.trans.get_left_channels(energy,n)
