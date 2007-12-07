from ase.io.eps import EPS


class PNG(EPS):
    def write_header(self):
        from matplotlib.backends.backend_agg import RendererAgg
        from matplotlib.backend_bases import GraphicsContextBase
        from matplotlib.transforms import Value, identity_transform

        self.renderer = RendererAgg(self.w, self.h, Value(72))
        identity = identity_transform()
        def line(gc, x1, y1, x2, y2):
            self.renderer.draw_lines(gc,
                                     (round(x1), round(x2)),
                                     (round(y1), round(y2)), identity)
        self.line = line
        self.gc = GraphicsContextBase()
        self.gc.set_linewidth(2)

    def write_trailer(self):
        self.renderer._renderer.write_png(self.filename)


def write_png(filename, atoms, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    PNG(atoms, **parameters).write(filename)
