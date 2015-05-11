.. module:: ase.dft.bandgap
   :synopsis: Band gap

========
Band gap
========

Example::

  from ase.dft.bandgap import get_band_gap
  calc = ...
  gap = get_band_gap(calc)

Here, gap is a list: [edir, ein, (sdir, kdir), (svin, kvin), (scin, kcin)], where edir and ein are the direct and indirect gaps respectively. sdir, kdir are the spin and k-point indices for the direct gap and (svin, kvin) and (scin, kcin) are the spin and k-point indices for the valence and conduction band respectively of the idirect gap.
