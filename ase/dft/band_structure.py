import numpy as np

from ase.dft.kpoints import labels_from_kpts
from ase.io.jsonio import encode, decode
from ase.parallel import paropen


def get_band_structure(atoms=None, calc=None):
    """Create band structure object from Atoms or calculator."""
    atoms = atoms if atoms is not None else calc.atoms
    calc = calc if calc is not None else atoms.calc

    kpts = calc.get_ibz_k_points()

    energies = []
    for s in range(calc.get_number_of_spins()):
        energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                         for k in range(len(kpts))])
    energies = np.array(energies)

    return BandStructure(cell=atoms.cell,
                         kpts=kpts,
                         energies=energies,
                         reference=calc.get_fermi_level())


class BandStructure:
    def __init__(self, *args, **kwargs):
        """Create band structure object from energies and k-points."""
        self.setvars(*args, **kwargs)

    def setvars(self, cell, kpts, energies, reference=0.0):
        assert cell.shape == (3, 3)
        self.cell = cell
        assert kpts.shape[1] == 3
        self.kpts = kpts
        self.energies = np.asarray(energies)
        self.reference = reference

    def get_labels(self):
        return labels_from_kpts(self.kpts, self.cell)

    def todict(self):
        return dict((key, getattr(self, key))
                    for key in
                    ['cell', 'kpts', 'energies', 'reference'])

    def write(self, filename):
        """Write to json file."""
        with paropen(filename, 'w') as f:
            f.write(encode(self))

    @staticmethod
    def read(filename):
        """Read from json file."""
        with open(filename, 'r') as f:
            dct = decode(f.read())
        return BandStructure(**dct)

    def plot(self, spin=None, emax=None, filename=None, ax=None, show=None,
             **plotkwargs):
        """Plot band-structure.

        spin: int or None
            Spin channel.  Default behaviour is to plot both spi up and down
            for spin-polarized calculations.
        emax: float
            Maximum energy above reference.
        filename: str
            Write image to a file.
        ax: Axes
            MatPlotLib Axes object.  Will be created if not supplied.
        show: bool
            Show the image.
        """

        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.gca()

        def pretty(kpt):
            if kpt == 'G':
                kpt = r'\Gamma'
            elif len(kpt) == 2:
                kpt = kpt[0] + '_' + kpt[1]
            return '$' + kpt + '$'

        if spin is None:
            e_skn = self.energies
        else:
            e_skn = self.energies[spin, None]

        emin = e_skn.min()
        if emax is not None:
            emax = emax + self.reference

        xcoords, label_xcoords, orig_labels = self.get_labels()

        labels = [pretty(name) for name in orig_labels]
        i = 1
        while i < len(labels):
            if label_xcoords[i - 1] == label_xcoords[i]:
                labels[i - 1] = labels[i - 1][:-1] + ',' + labels[i][1:]
                labels[i] = ''
            i += 1

        for spin, e_kn in enumerate(e_skn):
            color = 'br'[spin]
            for e_k in e_kn.T:
                kwargs = dict(color=color)
                kwargs.update(plotkwargs)
                ax.plot(xcoords, e_k, **kwargs)

        for x in label_xcoords[1:-1]:
            ax.axvline(x, color='0.5')

        ax.set_xticks(label_xcoords)
        ax.set_xticklabels(labels)
        ax.axis(xmin=0, xmax=xcoords[-1], ymin=emin, ymax=emax)
        ax.set_ylabel('eigenvalues [eV]')
        ax.axhline(self.reference, color='k')
        try:
            plt.tight_layout()
        except AttributeError:
            pass
        if filename:
            plt.savefig(filename)

        if show is None:
            show = not filename

        if show:
            plt.show()

        return ax
