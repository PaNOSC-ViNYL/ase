#!/usr/bin/env python
from math import sqrt

import gtk
from gettext import gettext as _

from ase.gui.widgets import pack, help

graph_help_text = _("""\
Help for plot ...

Symbols:
<c>e</c>:\t\t\ttotal energy
<c>epot</c>:\t\tpotential energy
<c>ekin</c>:\t\tkinetic energy
<c>fmax</c>:\t\tmaximum force
<c>fave</c>:\t\taverage force
<c>R[n,0-2]</c>:\t\tposition of atom number <c>n</c>
<c>d(n<sub>1</sub>,n<sub>2</sub>)</c>:\t\tdistance between two atoms n1 and n2
<c>i</c>:\t\t\tcurrent image number
<c>E[i]</c>:\t\t\tenergy of image number <c>i</c>
<c>F[n,0-2]</c>:\t\tforce on atom number <c>n</c>
<c>V[n,0-2]</c>:\t\tvelocity of atom number <c>n</c>
<c>M[n]</c>:\t\tmagnetic moment of atom number <c>n</c>
<c>A[0-2,0-2]</c>:\tunit-cell basis vectors
<c>s</c>:\t\t\tpath length
<c>a(n1,n2,n3)</c>:\tangle between atoms n1, n2 and n3, centered on n2
<c>dih(n1,n2,n3,n4)</c>:\tdihedral angle between n1, n2, n3, and n4
<c>T</c>:\t\t\tTemperature (K)\
""")

class Graphs(gtk.Window):
    def __init__(self, gui):
        gtk.Window.__init__(self)
        #self.window.set_position(gtk.WIN_POS_CENTER)
        #self.window.connect("destroy", lambda w: gtk.main_quit())
        #self.window.connect('delete_event', self.exit)
        self.set_title('Graphs')
        vbox = gtk.VBox()
        self.expr = pack(vbox, [gtk.Entry(64),
                                help(graph_help_text)])[0]
        self.expr.connect('activate', self.plot)

        completion = gtk.EntryCompletion()
        self.liststore = gtk.ListStore(str)
        for s in ['fmax', 's, e-E[0]', 'i, d(0,1)']:
            self.liststore.append([s])
        completion.set_model(self.liststore)
        self.expr.set_completion(completion)
        completion.set_text_column(0)

        button = pack(vbox, [gtk.Button(_('Plot')),
                             gtk.Label(' x, y1, y2, ...')])[0]
        button.connect('clicked', self.plot, 'xy')
        button = pack(vbox, [gtk.Button(_('Plot')),
                             gtk.Label(' y1, y2, ...')])[0]
        button.connect('clicked', self.plot, 'y')
        save_button = gtk.Button(stock=gtk.STOCK_SAVE)
        save_button.connect('clicked',self.save)
        clear_button = gtk.Button(_('clear'))
        clear_button.connect('clicked', self.clear)
        pack(vbox, [save_button,clear_button])
        

        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def plot(self, button=None, type=None, expr=None):
        if expr is None:
            expr = self.expr.get_text()
        else:
            self.expr.set_text(expr)

        if expr not in [row[0] for row in self.liststore]:
            self.liststore.append([expr])

        data = self.gui.images.graph(expr)
        
        import matplotlib
        matplotlib.interactive(True)
        matplotlib.use('GTK')
        #matplotlib.use('GTK', warn=False)# Not avail. in 0.91 (it is in 0.98)
        import pylab
        pylab.ion()

        x = 2.5
        self.gui.graphs.append(pylab.figure(figsize=(x * 2.5**0.5, x)))
        i = self.gui.frame
        m = len(data)

        if type is None:
            if m == 1:
                type = 'y'
            else:
                type = 'xy'
                
        if type == 'y':
            for j in range(m):
                pylab.plot(data[j])
                pylab.plot([i], [data[j, i]], 'o')
        else:
            for j in range(1, m):
                pylab.plot(data[0], data[j])
                pylab.plot([data[0, i]], [data[j, i]], 'o')
        pylab.title(expr)
        #pylab.show()

    python = plot

    def save(self, filename):
        chooser = gtk.FileChooserDialog(
            _('Save data to file ... '), None, gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        save = chooser.run()
        if save == gtk.RESPONSE_OK:
            filename = chooser.get_filename()
            expr = self.expr.get_text()
            data = self.gui.images.graph(expr)
            expr = '# '+expr
            fd = open(filename,'w')
            fd.write("%s \n" % (expr))
            for s in range(len(data[0])):
                for i in range(len(data)):
                    val = data[i,s]
                    fd.write("%12.8e\t" % (val))
                fd.write("\n")
            fd.close()
        chooser.destroy()

    def clear(self, button):
        import pylab
        for graph in self.gui.graphs:
            pylab.close(graph)
        self.gui.graphs = []
        
