import numpy as np
from collections import OrderedDict
import itertools


class PDOStype:
    def __init__(self, weights=None,
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
    def __init__(self, dos=None, sampling={'type': 'raw'}):

        # We maintain the order as well as we can
        self.pdos = OrderedDict()
        self.sampling = sampling

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
                for name, things in dos.items():
                    # The dict needs to have
                    # "energy" and "weights"
                    en = things.pop('energy', None)
                    weights = things.pop('weights', None)
                    self.add(name, weights=weights, energy=en, **things)
            else:
                # Let's assume it's a zip(names, things)
                # and things is a dict
                for names, things in dos:
                    en = things.pop('energy', None)
                    weights = things.pop('weights', None)
                    self.add(name, weights=weights, energy=en, **things)

    def add(self, name, weights=None, energy=None, **info):
        self.pdos[name] = PDOStype(weights=weights, energy=energy, info=info)

    def __iter__(self):
        self._it = iter(self.pdos.items())
        return self

    def __next__(self):
        return next(self._it)

    def delta(self, x, x0, width):
        """Return a delta-function centered at 'x0'."""
        x1 = -((x - x0) / width)**2
        return np.exp(x1) / (np.sqrt(np.pi) * width)

    def smear(self, energy, weights, npts=401, width=0.1, grid=None):
        """Add Gaussian smearing, and map energies and weights on grid
        of length npts. Disabled for 0 width
        """

        if width == 0.0:
            return energy, weights
        else:
            if grid is None:
                # Make new linear uniform grid
                dos = np.zeros(npts)
                energy_grid = np.linspace(min(energy), max(energy), npts)
            else:
                # Use the energy we specified as the grid
                dos = np.zeros(len(grid))
                energy_grid = grid

            for en, w in zip(energy, weights):
                dos += w * self.delta(energy_grid, en, width)
            return energy_grid, dos

    def sample(self, npts=401, width=0.1, type='Gauss',
               window=None, grid=None):

        # What exactly should 'sampling' do?
        if grid is not None:
            sampling = self.sampling  # Should this be something else?
        else:
            sampling = {'type': 'uniform'}

        pdos_new = PDOS(sampling=sampling)

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

            if energy is None:
                # Check if we can use the grid as the energies instead
                if grid is None:
                    msg = ('Either the PDOStype must contain'
                           ' energies or grid must be specified')
                    raise ValueError(msg)
                else:
                    # Should this raise a warning?
                    energy = grid

            energy_grid, weights = self.smear(energy, pd.weights,
                                              npts=npts, width=width,
                                              grid=grid)

            # Apply window
            idx = (emin <= energy_grid) & (energy_grid <= emax)
            energy_grid = energy_grid[idx]
            weights = weights[idx]

            pdos_new.add(name, weights=weights, energy=energy_grid, info=pd.info)

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

    @staticmethod
    def resample(doslist, names=None, grid=None):

        doslist = np.asarray(doslist)
        # Check correct dimensionality, either [...] or [[..], [..], ...]
        if doslist.ndim not in [1, 2]:
            msg = ('Incorrect number of dimensions for doslist.'
                   'Should be either 1 or 2, got {}'.format(doslist.ndim))
            raise ValueError(msg)

        if doslist.ndim == 1:
            # Add extra axis to preserve syntax
            doslist = doslist[np.newaxis]

        if grid is not None:
            grid = np.asarray(grid)
            # Check correct dimensionality, either [...] or [[..], [..], ...]
            if doslist.ndim not in [1, 2]:
                msg = ('Incorrect number of dimensions for energy grid.'
                       'Should be either 1 or 2, got {}'.format(grid.ndim))
                raise ValueError(msg)

            if grid.ndim == 1:
                # Use same grid for every dos
                grid = itertools.cycle(grid[np.newaxis])
        else:
            # Use None as energy grid
            grid = itertools.cycle([None])

        if names is None:
            names = range(len(doslist))
        else:
            # This could fail if names is just a string
            # but do we want to test every possibility?
            # Greater or equal, have at least enough names
            if len(names) >= len(doslist):
                msg = ('Not enough provided number of names.'
                       'Expected at least {}, got {}'.format(len(doslist),
                                                             len(names)))
                raise ValueError(msg)

        pdos = PDOS()
        for name, en, weights in zip(names, grid, doslist):
            pdos.add(name, weights=weights, energy=en)
        return pdos


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
