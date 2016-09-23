import pickle
import subprocess
import sys

import ase.gui.ui as ui

graph_help_text = """\
Symbols:
<c>e</c>: total energy
<c>epot</c>: potential energy
<c>ekin</c>: kinetic energy
<c>fmax</c>: maximum force
<c>fave</c>: average force
<c>R[n,0-2]</c>: position of atom number <c>n</c>
<c>d(n<sub>1</sub>,n<sub>2</sub>)</c>: distance between two atoms <c>n<sub>1</sub></c> and <c>n<sub>2</sub></c>
<c>i</c>: current image number
<c>E[i]</c>: energy of image number <c>i</c>
<c>F[n,0-2]</c>: force on atom number <c>n</c>
<c>V[n,0-2]</c>: velocity of atom number <c>n</c>
<c>M[n]</c>: magnetic moment of atom number <c>n</c>
<c>A[0-2,0-2]</c>: unit-cell basis vectors
<c>s</c>: path length
<c>a(n1,n2,n3)</c>: angle between atoms <c>n<sub>1</sub></c>, <c>n<sub>2</sub></c> and <c>n<sub>3</sub></c>, centered on <c>n<sub>2</sub></c>
<c>dih(n1,n2,n3,n4)</c>: dihedral angle between <c>n<sub>1</sub></c>, <c>n<sub>2</sub></c>, <c>n<sub>3</sub></c> and <c>n<sub>4</sub></c>
<c>T</c>: temperature (K)\
"""


def help():
    win = ui.Window('Help for plot')
    win.add(ui.Text(graph_help_text))


class Graphs:
    def __init__(self, gui):
        win = ui.Window('Graphs')
        self.expr = ui.Entry('', 50, self.plot)
        win.add([self.expr, ui.Button('Help', help)])

        win.add([ui.Button('Plot', self.plot, 'xy'), ' x, y1, y2, ...'], 'w')
        win.add([ui.Button('Plot', self.plot, 'y'), ' y1, y2, ...'], 'w')

        win.add([ui.Button('Save', self.save),
                 ui.Button('Clear', self.clear)], 'w')

        self.gui = gui

    def plot(self, type=None, expr=None):
        if expr is None:
            expr = self.expr.value
        else:
            self.expr.value = expr

        data = self.gui.images.graph(expr)

        import matplotlib
        matplotlib.use('TkAgg')

        process = subprocess.Popen([sys.executable, '-m', 'ase.gui.graphs'],
                                   stdin=subprocess.PIPE)
        pickle.dump((data, self.gui.frame, expr, type), process.stdin)
        process.stdin.close()
        self.gui.graphs.append(process)

    python = plot

    def save(self, filename):
        chooser = ui.FileChooserDialog(
            _('Save data to file ... '), None, ui.FILE_CHOOSER_ACTION_SAVE,
            ('Cancel', ui.RESPONSE_CANCEL,
             'Save', ui.RESPONSE_OK))
        save = chooser.run()
        if save == ui.RESPONSE_OK:
            filename = chooser.get_filename()
            expr = self.expr.get_text()
            data = self.gui.images.graph(expr)
            expr = '# ' + expr
            fd = open(filename, 'w')
            fd.write("%s \n" % (expr))
            for s in range(len(data[0])):
                for i in range(len(data)):
                    val = data[i, s]
                    fd.write("%12.8e\t" % (val))
                fd.write("\n")
            fd.close()
        chooser.destroy()

    def clear(self):
        for process in self.gui.graphs:
            process.terminate()
        self.gui.graphs = []


def make_plot(data, i, expr, type):
    import matplotlib.pyplot as plt
    x = 4
    plt.figure(figsize=(x * 2.5**0.5, x))
    m = len(data)

    if type is None:
        if m == 1:
            type = 'y'
        else:
            type = 'xy'

    if type == 'y':
        for j in range(m):
            plt.plot(data[j])
            plt.plot([i], [data[j, i]], 'o')
    else:
        for j in range(1, m):
            plt.plot(data[0], data[j])
            plt.plot([data[0, i]], [data[j, i]], 'o')
    plt.title(expr)
    plt.show()


if __name__ == '__main__':
    fd = sys.stdin
    if sys.version_info[0] > 2:
        fd = fd.buffer
    make_plot(*pickle.load(fd))
