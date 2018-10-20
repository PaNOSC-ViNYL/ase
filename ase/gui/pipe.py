from __future__ import print_function
import pickle
import sys


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    task, data = pickle.load(sys.stdin)
    if task == 'eos':
        from ase.eos import plot
        plot(*data)
    elif task == 'neb':
        from ase.neb import plot_band_from_fit
        plot_band_from_fit(*data)
    elif task == 'reciprocal':
        from ase.dft.bz import bz3d_plot
        bz3d_plot(**data)
    else:
        print('Invalid task {}'.format(task))
        sys.exit(17)

    # Magic string to tell GUI that things went okay:
    print('GUI:OK')

    plt.show()
