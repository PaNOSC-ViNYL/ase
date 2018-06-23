import numpy as np
from collections import OrderedDict


class PDOStype:
    def __init__(self, energy, weights,
                 **kwargs):
        """Basic PDOS derived type"""
        self.energy = energy
        self.weights = weights

        # Story any other information about this
        # PDOS that the user might want to use later.
        self.__dict__.update(**kwargs)

    def get_energy(self):
        return self.energy

    def get_weights(self):
        return self.weights


class PDOS:
    def __init__(self, dos, width=0.1, npts=401):

        self.width = width
        self.npts = npts

        # We maintain the order as well as we can
        self.pdos = OrderedDict()

        # This needs to be made more robust.
        # Just a quick implementation
        # Possible alternatives to consider: what if we're just
        # adding two lists (of lists) instead?
        # Possbility of telling that the energies are on the
        # same scale, so we only need to store it once?

        if isinstance(dos, dict):
            # assume it's a dict of dicts
            for key, value in dos.items():
                # The dict needs to have
                # "energy" and "weights"
                en = value.pop('energy')
                weights = value.pop('weights')
                self.pdos[key] = PDOStype(en, weights, **value)
        else:
            # Let's assume it's a zip(names, things)
            # and things is a dict
            for names, things in dos:
                en = things.pop('energy')
                weights = things.pop('weights')
                self.pdos[names] = PDOStype(en, weights, **things)

    def delta(self, x, x0):
        """Return a delta-function centered at 'x0'."""
        x1 = -((x - x0) / self.width)**2
        return np.exp(x1) / (np.sqrt(np.pi) * self.width)

    def smear(self, energy, weights):
        """Add Gaussian smearing, and map energies and weights on grid
        of length npts. Disabled for 0 width
        """

        if self.width == 0.0:
            return energy, weights
        else:
            dos = np.zeros(self.npts)
            energies = np.linspace(min(energy), max(energy), self.npts)

            for en, w in zip(energy, weights):
                dos += w * self.delta(energies, en)
            return energies, dos


class PDOSPlot:
    def __init__(self, pdos):
        self.pdos = pdos
        self.ax = None

    def plot(self, ax=None, emin=None, emax=None, ymin=None, ymax=None,
             filename=None, show=None, ylabel=None, colors=None,
             show_legend=True, loc=None, **plotkwargs):

        if self.ax is None:
            ax = self.prepare_plot(ax, emin, emax,
                                   ymin=ymin, ymax=ymax,
                                   ylabel=ylabel)

        for name, _pdos in self.pdos.pdos.items():
            energies = _pdos.get_energy()
            weights = _pdos.get_weights()
            en, w = self.pdos.smear(energies, weights)
            ax.plot(en, w, label=name)

        self.finish_plot(filename, show, show_legend, loc)

        return ax

    def prepare_plot(self, ax=None, emin=None, emax=None,
                     ymin=None, ymax=None, ylabel=None):
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.figure().add_subplot(111)

        ylabel = ylabel if ylabel is not None else 'DOS'
        ax.axis(xmin=emin, xmax=emax, ymin=ymin, ymax=ymax)
        ax.set_ylabel(ylabel)
        self.ax = ax
        return ax

    def finish_plot(self, filename, show, show_legend, loc):
        # Put in finalize_plot
        import matplotlib.pyplot as plt

        if show_legend:
            leg = plt.legend(loc='best')
            leg.get_frame().set_alpha(1)

        if filename:
            plt.savefig(filename)

        if show is None:
            show = not filename

        if show:
            plt.show()
