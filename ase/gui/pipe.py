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
    else:
        print('Invalid task {}'.format(task), file=sys.stderr)
        sys.exit(17)

    plt.show()
