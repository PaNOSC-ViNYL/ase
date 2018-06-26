import numpy as np


class PDOS:
    def __init__(self, energy, weights, info=None, sampling={'type': 'raw'}):
        """
        Docstring here
        """
        self.energy = np.asarray(energy)
        self.weights = np.asarray(weights)
        self.sampling = sampling

        # Energy format: [e1, e2, ...]
        if self.energy.ndim > 1:
            msg = ('Incorrect Energy dimensionality. '
                   'Expected 1 got {}'.format(
                       self.energy.ndim))
            raise ValueError(msg)

        # Weights format: [[w1, w2, ...], [w1, w2, ..], ...]
        if self.weights.ndim != 2:
            msg = ('Incorrect weight dimensionality. '
                   'Expected 2, got {}'.format(
                       self.weights.ndim))
            raise ValueError(msg)

        # Check weight shape matches energy
        if self.weights.shape[1] != self.energy.shape[0]:
            msg = ('Weight dimensionality does not match energy.'
                   ' Expected {}, got {}'.format(self.energy.shape[0],
                                                 self.weights.shape[1]))
            raise ValueError(msg)

        # One entry for info for each weight
        if info is None:
            info = [None for _ in range(len(self.weights))]
        else:
            if len(info) != len(weights):
                msg = ('Incorrect number of entries in '
                       'info. Expected {}, got {}'.format(
                           len(self.weights), len(info)))
                raise ValueError(info)
        self.info = info

    def delta(self, x, x0, width):
        """Return a delta-function centered at 'x0'."""
        x1 = -((x - x0) / width)**2
        return np.exp(x1) / (np.sqrt(np.pi) * width)

    def smear(self, energy_grid, width=0.1):
        """Add Gaussian smearing, to all weights onto an energy grid.
        Disabled for width=0.0"""
        if width == 0.0:
            return self.weights

        en0 = self.energy[:, np.newaxis]  # Add axis to use NumPy broadcasting
        weights_grid = np.dot(self.weights,
                              self.delta(energy_grid, en0, width=width))

        return weights_grid

    def sample(self, grid, width=0.1, smearing='Gauss', gridtype='grid'):
        """Sample weights onto new specified grid"""

        npts = len(grid)
        sampling = {'width': width,
                    'smearing': smearing,
                    'npts': npts,
                    'type': gridtype}

        weights_grid = np.zeros((self.weights.shape[0], npts))

        weights_grid = self.smear(grid, width=width)

        pdos_new = PDOS(grid, weights_grid,
                        info=self.info, sampling=sampling)
        return pdos_new

    def sample_uniform(self, spacing=None, npts=None, width=0.1,
                       window=None, smearing='Gauss'):
        """Sample onto uniform grid"""

        if window is None:
            emin, emax = None, None
        else:
            emin, emax = window

        if emin is None:
            emin = self.energy.min()
        if emax is None:
            emax = self.energy.max()

        grid_uniform = PDOS._make_uniform_grid(emin, emax, spacing=spacing,
                                               npts=npts, width=width)

        return self.sample(grid_uniform, width=width,
                           smearing=smearing, gridtype='uniform')

    @staticmethod
    def resample(doslist, grid, width=0.1, smearing='Gauss',
                 gridtype='resample_grid'):
        """Take list of PDOS objects, and combine into 1, with same grid"""

        # Count the total number of weights
        n_weights = sum(len(dos.weights) for dos in doslist)

        npts = len(grid)

        weight_grid = np.zeros((n_weights, npts))
        info_new = []
        # Do sampling
        ii = 0
        for dos in doslist:
            pdos_sample = dos.sample(grid, width=width,
                                     smearing=smearing)
            info_new.extend(pdos_sample.info)
            for w_i in pdos_sample.weights:
                weight_grid[ii] = w_i
                ii += 1
        sampling = {'smearing': smearing,
                    'width': width,
                    'npts': npts,
                    'type': gridtype}
        return PDOS(energy=grid, weights=weight_grid, info=info_new,
                    sampling=sampling)

    @staticmethod
    def resample_uniform(doslist, window=None, spacing=None,
                         npts=None, width=0.1, smearing='Gauss'):
        """Resample list of PDOS objects onto uniform grid.
        Takes the lowest and highest energies as grid range, if
        no window is specified"""
        # Parse window
        if window is None:
            emin, emax = None, None
        else:
            emin, emax = window
        if emin is None:
            emin = -np.infty
        if emax is None:
            emax = np.infty
        dosen = [dos.energy for dos in doslist]

        # If needed, adjust emin and emax to be within
        # the range of sampled data
        emin = max(np.min(dosen), emin)
        emax = min(np.max(dosen), emax)

        grid_uniform = PDOS._make_uniform_grid(emin, emax, spacing=spacing,
                                               npts=npts, width=width)

        return PDOS.resample(doslist, grid_uniform, width=width,
                             smearing=smearing, gridtype='resample_uniform')

    @staticmethod
    def _make_uniform_grid(emin, emax, spacing=None, npts=None, width=0.1):
        if spacing and npts:
            msg = ('spacing and npts cannot both be defined'
                   ' at the same time.')
            raise ValueError(msg)
        if not spacing and not npts:
            # Default behavior
            spacing = 0.2 * width
        # Now either spacing or npts is defined
        if npts:
            grid_uniform = np.linspace(emin, emax, npts)
        else:
            grid_uniform = np.arange(emin, emax, spacing)
        return grid_uniform

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

        for ii, w_i in enumerate(self.pdos.weights):
            # We can add smater labeling later
            label = self.pdos.info[ii]
            ax.plot(self.pdos.energy, w_i, label=label, **plotkwargs)

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
