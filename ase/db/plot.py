from ase.utils import basestring


def dct2plot(dct, filename=None, show=True):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    styles = ['k-', 'r-', 'g-', 'b-']
    for d in dct['data']:
        x = d['x']
        Y = d['y']
        if Y.ndim == 1:
            Y = Y[None]
        style = d.get('style')
        if not style:
            style = styles.pop()
        plt.plot(x, Y[0], style, label=d['label'])
        for y in Y[1:]:
            plt.plot(x, y, style)
        plt.legend()
    if isinstance(dct['xlabel'], basestring):
        plt.xlabel(dct['xlabel'])
    else:
        x, labels = dct['xlabel']
        plt.xticks(x, labels)
        plt.xlim(x[0], x[-1])
    plt.ylabel(dct['ylabel'])
    if 'ylim' in dct:
        plt.ylim(*dct['ylim'])
    plt.title(dct['title'])
    plt.tight_layout()
    if show:
        plt.show()
    if filename:
        plt.savefig(filename)
        plt.close(fig)
