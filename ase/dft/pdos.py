import numpy as np
from collections import OrderedDict


class PDOStype:
    def __init__(self, weights,
                 energy=None, info=None):
        """Basic PDOS derived type"""
        self.energy = energy
        self.weights = weights

        # Store any other information about this
        # PDOS that the user might want to use later.
        # From atoms.info
        if info is None:
            self.info = {}
        else:
            self.info = dict(info)


class PDOS:
    def __init__(self, dos=None):

        # We maintain the order as well as we can
        self.pdos = OrderedDict()

        # This needs to be made more robust.
        # Just a quick implementation
        # Possible alternatives to consider: what if we're just
        # adding two lists (of lists) instead?
        # Possbility of telling that the energies are on the
        # same scale, so we only need to store it once?

        # Let's allow the user to set up an empty PDOS object,
        # and then add the energies and DOS afterwards
        if dos:
            if isinstance(dos, dict):
                # assume it's a dict of dicts
                for key, value in dos.items():
                    # The dict needs to have
                    # "energy" and "weights"
                    en = value.pop('energy')
                    weights = value.pop('weights')
                    self.pdos[key] = PDOStype(weights, energy=en, info=value)
            else:
                # Let's assume it's a zip(names, things)
                # and things is a dict
                for names, things in dos:
                    en = things.pop('energy')
                    weights = things.pop('weights')
                    self.pdos[names] = PDOStype(weights,
                                                energy=en, info=things)

    def add(self, name, weights, energy=None, **info):
        self.pdos[name] = PDOStype(weights, energy=energy, info=info)

    def __iter__(self):
        self._it = iter(self.pdos.items())
        return self

    def __next__(self):
        return next(self._it)

    def delta(self, x, x0, width):
        """Return a delta-function centered at 'x0'."""
        x1 = -((x - x0) / width)**2
        return np.exp(x1) / (np.sqrt(np.pi) * width)

    def smear(self, energy, weights, npts=401, width=0.1):
        """Add Gaussian smearing, and map energies and weights on grid
        of length npts. Disabled for 0 width
        """

        if width == 0.0:
            return energy, weights
        else:
            dos = np.zeros(npts)
            energies = np.linspace(min(energy), max(energy), npts)

            for en, w in zip(energy, weights):
                dos += w * self.delta(energies, en, width)
            return energies, dos

    def sample(self, npts=401, width=0.1, type='Gauss',
               window=None, grid=None, sampling={'type': 'raw'}):

        pdos_new = PDOS()

        if window is None:
            emin, emax = None, None
        else:
            emin, emax = window
        if emin is None:
            emin = -np.infty
        if emax is None:
            emax = np.infty

        for name, pd in self:
            energy = pd.energy

            idx = (emin <= energy) & (energy <= emax)
            energy = energy[idx]
            weights = pd.weights[idx]

            energy_lin, weights = self.smear(energy, weights,
                                             npts=npts, width=width)

            pdos_new.add(name, weights, energy=energy_lin, info=pd.info)

        return pdos_new

    def plot(self, *plotargs,
             # We need to grab the init keywords
             ax=None,
             emin=None, emax=None,
             ymin=None, ymax=None, ylabel=None,
             **plotkwargs):

        pdp = PDOSPlot(self, ax=None,
                       emin=None, emax=None,
                       ymin=None, ymax=None, ylabel=None)
        return pdp.plot(*plotargs, **plotkwargs)


class PDOSPlot:
    def __init__(self, pdos, ax=None,
                 emin=None, emax=None,
                 ymin=None, ymax=None, ylabel=None):
        self.pdos = pdos
        self.ax = ax
        if self.ax is None:
            self.ax = self.prepare_plot(ax, emin, emax,
                                        ymin=ymin, ymax=ymax,
                                        ylabel=ylabel)

    def plot(self, filename=None, show=None, colors=None,
             show_legend=True, loc=None, **plotkwargs):

        ax = self.ax

        for name, _pdos in self.pdos:
            ax.plot(_pdos.energy, _pdos.weights, label=name)

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
