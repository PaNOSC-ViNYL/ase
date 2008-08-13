# -*- coding: utf-8 -*-
# creates:  ener.png distance.png
import os
from ase import *
import matplotlib
matplotlib.use('Agg')
import pylab as plt


e_s = [0.01,0.1,0.2,0.3,0.4,0.5]
E = [-461.8360, -461.6443, -461.2166, -460.6679,
     -460.3595, -459.7719]
d = [1.1220, 1.1129, 1.1040, 1.0975, 1.0934,
     1.0877]
fig=plt.figure(figsize=(4.5, 3.6))
fig.subplots_adjust(left=.19, right=.97, top=.9, bottom=0.12)
plt.plot(e_s, E, 'o-')
plt.xlabel(u'Energy shift [eV]')
plt.ylabel(u'Energy [eV]')
plt.title('Total Energy vs Eshift')
plt.savefig('ener.png')

fig=plt.figure(figsize=(4.5, 3.6))
fig.subplots_adjust(left=.19, right=.97, top=.9, bottom=0.12)
plt.plot(e_s, d, 'o-')
plt.xlabel(u'O-H distance [Å]')
plt.ylabel(u'Energy [eV]')
limits = plt.axis('tight')
plt.title('O-H distance vs Eshift')
plt.savefig('distance.png')
