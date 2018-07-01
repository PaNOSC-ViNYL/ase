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
        if self.energy.ndim != 1:
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
            info = [{} for _ in self.weights]
        else:
            if len(info) != len(weights):
                msg = ('Incorrect number of entries in '
                       'info. Expected {}, got {}'.format(
                           len(self.weights), len(info)))
                raise ValueError(msg)
        self.info = np.asarray(info)  # Make info np array for slicing purposes

    def delta(self, x, x0, width, smearing='Gauss'):
        """Return a delta-function centered at 'x0'."""
        if smearing.lower() == 'gauss':
            x1 = -((x - x0) / width)**2
            return np.exp(x1) / (np.sqrt(np.pi) * width)
        else:
            msg = 'Requested smearing type not recognized. Got {}'.format(
                smearing)
            raise ValueError(msg)

    def smear(self, energy_grid, width=0.1, smearing='Gauss'):
        """Add Gaussian smearing, to all weights onto an energy grid.
        Disabled for width=0.0"""
        if width == 0.0:
            msg = 'Cannot add 0 width smearing'
            raise ValueError(msg)

        en0 = self.energy[:, np.newaxis]  # Add axis to use NumPy broadcasting
        weights_grid = np.dot(self.weights,
                              self.delta(energy_grid, en0, width,
                                         smearing=smearing))

        return weights_grid

    def sample(self, grid, width=0.1, smearing='Gauss', gridtype='general'):
        """Sample weights onto new specified grid"""

        npts = len(grid)
        sampling = {'width': width,
                    'smearing': smearing,
                    'npts': npts,
                    'type': gridtype}

        weights_grid = self.smear(grid, width=width, smearing=smearing)

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
        emin -= 5 * width
        emax += 5 * width

        grid_uniform = PDOS._make_uniform_grid(emin, emax, spacing=spacing,
                                               npts=npts, width=width)

        return self.sample(grid_uniform, width=width,
                           smearing=smearing, gridtype='uniform')

    @staticmethod
    def resample(doslist, grid, width=0.1, smearing='Gauss',
                 gridtype='general'):
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
        dosen = [dos.energy for dos in doslist]
        # Parse window
        if window is None:
            emin, emax = None, None
        else:
            emin, emax = window
        if emin is None:
            emin = np.min(dosen)
        if emax is None:
            emax = np.max(dosen)
        # Add a little extra to avoid stopping midpeak
        emin -= 5 * width
        emax += 5 * width

        grid_uniform = PDOS._make_uniform_grid(emin, emax, spacing=spacing,
                                               npts=npts, width=width)

        return PDOS.resample(doslist, grid_uniform, width=width,
                             smearing=smearing, gridtype='uniform')

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

    def plot(self,
             # We need to grab init keywords
             ax=None,
             emin=None, emax=None,
             ymin=None, ymax=None, ylabel=None,
             *plotargs, **plotkwargs):

        pdp = PDOSPlot(self, ax=None,
                       emin=None, emax=None,
                       ymin=None, ymax=None, ylabel=None)
        return pdp.plot(*plotargs, **plotkwargs)

    def sum(self):
        weights_sum = self.weights.sum(0)[np.newaxis]

        # Find shared (key, value) pairs
        # dict(set.intersection(*(set(d.items()) for d in info)))
        all_kv = []
        for d in self.info:
            kv_pairs = set()
            for key, value in d.items():
                try:
                    kv_pairs.add((key, value))
                except TypeError:
                    # Unhashable type, skip it
                    pass
            all_kv.append(kv_pairs)
        if all_kv:
            info_new = [dict(set.intersection(*all_kv))]
        else:
            # We didn't find any shared (key, value) pairs
            # This prevents set.intersection from blowing up
            info_new = None

        return PDOS(energy=self.energy, weights=weights_sum,
                    info=info_new, sampling=self.sampling)

    def pick(self, **kwargs):
        # Pick key/value pairs using logical AND
        # i.e., all conditions from kwargs must be met
        idx = [i for i, d in enumerate(self.info)
               if all(d.get(key) == value
                      for key, value in kwargs.items())]

        return self[idx]

    def split(self, key):
        # Find all unique instances of key in info
        unique = np.unique([info.get(key) for info in self.info
                            if info.get(key, None) is not None])

        pdos_lst = []
        for value in unique:
            # Use **{key: value} instead of key=value,
            # as key=value will litterally look up "key" in info.
            pdos_lst.append(self.pick(**{key: value}))
        return pdos_lst

    def __getitem__(self, i):
        if isinstance(i, int):
            n_weights = len(self.weights)
            if i < -n_weights or i >= n_weights:
                raise IndexError('Index out of range.')

        indices = np.arange(len(self.weights))[i]
        if len(indices.shape) == 0:
            indices = indices[np.newaxis]

        return PDOS(energy=self.energy,
                    weights=self.weights[indices],
                    info=self.info[indices],
                    sampling=self.sampling)


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
             labels=None, show_legend=True, loc='best', **plotkwargs):

        ax = self.ax

        for ii, w_i in enumerate(self.pdos.weights):
            # We can add smater labeling later
            kwargs = {}
            if colors is not None:
                kwargs['color'] = colors[ii]

            # We could possibly have some better label logic here
            if labels is not None:
                kwargs['label'] = labels[ii]
            else:
                kwargs['label'] = self.pdos.info[ii]
            kwargs.update(plotkwargs)
            ax.plot(self.pdos.energy, w_i,
                    **kwargs)

        self.finish_plot(filename, show, show_legend, loc)

        return ax

    def prepare_plot(self, ax=None, emin=None, emax=None,
                     ymin=None, ymax=None,
                     ylabel=None, xlabel=None):
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.figure().add_subplot(111)

        ylabel = ylabel if ylabel is not None else 'DOS'
        xlabel = xlabel if xlabel is not None else 'Energy [eV]'
        ax.axis(xmin=emin, xmax=emax, ymin=ymin, ymax=ymax)
        ax.set_ylabel(ylabel)
        self.ax = ax
        return ax

    def finish_plot(self, filename, show, show_legend, loc):
        import matplotlib.pyplot as plt

        if show_legend:
            leg = plt.legend(loc=loc)
            leg.get_frame().set_alpha(1)

        if filename:
            plt.savefig(filename)

        if show is None:
            show = not filename

        if show:
            plt.show()
