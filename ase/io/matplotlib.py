from ase.io.utils import generate_writer_variables, make_patch_list


class MATPLOTLIB:
    def __init__(self, atoms, ax,
                 rotation='', show_unit_cell=False, radii=None,
                 colors=None, scale=1, offset=(0, 0)):
        generate_writer_variables(
            self, atoms, rotation=rotation,
            show_unit_cell=show_unit_cell,
            radii=radii, bbox=None, colors=colors, scale=scale,
            extra_offset=offset)

        self.ax = ax
        self.figure = ax.figure
        self.ax.set_aspect('equal')

    def write(self):
        self.write_body()
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(0, self.h)

    def write_body(self):
        patch_list = make_patch_list(self)
        for patch in patch_list:
            self.ax.add_patch(patch)


def write_matplotlib(atoms, ax, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    MATPLOTLIB(atoms, ax, **parameters).write()
