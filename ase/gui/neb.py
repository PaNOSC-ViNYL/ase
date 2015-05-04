# -*- coding: utf-8 -*-

from ase.neb import NEBtools


def NudgedElasticBand(images):
    nebtools = NEBtools(images)
    fig = nebtools.plot_band()
    fig.show()
