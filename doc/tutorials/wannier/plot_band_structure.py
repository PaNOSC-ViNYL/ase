import pylab as pl

k1, e1 = pl.load('KSbands.txt', unpack=True)
k2, e2 = pl.load('WANbands.txt', unpack=True)
k3, e3 = pl.load('NEWbands.txt', unpack=True)

fig = pl.figure(1, dpi=80, figsize=(4.2, 6))
fig.subplots_adjust(left=.16, right=.97, top=.97, bottom=.05)
pl.plot(k1, e1, 'ro', label='DFT', ms=9)
pl.plot(k2, e2, 'kx', label='Wannier', mew=2)
pl.plot(k3, e3, 'k.', label='Wannier')
pl.plot([-.5, .5], [1, 1], 'k:', label='_nolegend_')
pl.text(-.5, 1, 'fixedenergy', ha='left', va='bottom')
pl.axis('tight')
pl.xticks([-.5, -.25, 0, .25, .5],
          [ r'$X$', r'$\Delta$', r'$\Gamma$', r'$\Delta$', r'$X$'], size=16)
pl.ylabel(r'$E - E_F\  \rm{(eV)}$', size=16)
pl.legend()
pl.savefig('bands.png', dpi=80)
pl.show()
