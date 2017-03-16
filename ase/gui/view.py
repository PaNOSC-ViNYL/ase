from __future__ import division
from math import cos, sin, sqrt
from os.path import basename

import numpy as np

from ase.data.colors import jmol_colors
from ase.geometry import complete_cell
from ase.gui.repeat import Repeat
from ase.gui.rotate import Rotate
from ase.gui.render import Render
from ase.gui.colors import ColorWindow
from ase.utils import rotate


GREEN = '#DDFFDD'


class View:
    def __init__(self, rotations):
        self.colormode = 'jmol'  # The default colors
        self.nselected = 0
        self.labels = None
        self.axes = rotate(rotations)
        self.configured = False
        self.frame = None

        # XXX
        self.colormode = 'jmol'
        self.colors = {}

        for i, rgb in enumerate(jmol_colors):
            self.colors[i] = ('#{0:02X}{1:02X}{2:02X}'
                              .format(*(int(x * 255) for x in rgb)))

    @property
    def atoms(self):
        return self.images[self.frame]

    def set_coordinates(self, frame=None, focus=None):
        if frame is None:
            frame = self.frame
        self.make_box()
        self.bind(frame)
        atoms = self.images[frame]
        natoms = len(atoms)
        self.X = np.empty((natoms + len(self.B1) + len(self.bonds), 3))
        self.X_pos = self.X[:natoms]
        self.X_pos[:] = atoms.positions
        self.X_B1 = self.X[natoms:natoms + len(self.B1)]
        self.X_bonds = self.X[natoms + len(self.B1):]
        self.set_frame(frame, focus=focus, init=True)

    def set_frame(self, frame=None, focus=False, init=False):
        if frame is None:
            frame = self.frame
        atoms = self.images[frame]

        self.X[:len(atoms)] = atoms.positions

        # XXX this should raise an error; caller must provide valid numbers!
        if self.frame is not None and self.frame > len(self.images):
            self.frame = len(self.images) - 1

        if init or frame != self.frame:
            cell = atoms.cell
            nc = len(self.B1)
            nbonds = len(self.bonds)

            if init or (atoms.cell != self.atoms.cell).any():
                self.X_B1[:] = np.dot(self.B1, cell)
                self.B = np.empty((nc + nbonds, 3))
                self.B[:nc] = np.dot(self.B2, cell)

            if nbonds > 0:
                P = self.atoms.positions
                Af = self.images.repeat[:, np.newaxis] * cell
                a = P[self.bonds[:, 0]]
                b = P[self.bonds[:, 1]] + np.dot(self.bonds[:, 2:], Af) - a
                d = (b**2).sum(1)**0.5
                r = 0.65 * self.get_covalent_radii()
                x0 = (r[self.bonds[:, 0]] / d).reshape((-1, 1))
                x1 = (r[self.bonds[:, 1]] / d).reshape((-1, 1))
                self.X_bonds[:] = a + b * x0
                b *= 1.0 - x0 - x1
                b[self.bonds[:, 2:].any(1)] *= 0.5
                self.B[nc:] = self.X_bonds + b

            filenames = self.images.filenames
            filename = filenames[frame]
            if (self.frame is None or
                filename != filenames[self.frame] or
                filename is None):
                if filename is None:
                    filename = 'ase.gui'
            filename = basename(filename)
            self.window.title = filename

        self.frame = frame
        if focus:
            self.focus()
        else:
            self.draw()

    def make_box(self):
        if not self.window['toggle-show-unit-cell']:
            self.B1 = self.B2 = np.zeros((0, 3))
            return

        # This function uses the box of the first Atoms object!  How can this
        # be right??
        V = self.images[0].get_cell()
        nn = []
        for c in range(3):
            v = V[c]
            d = sqrt(np.dot(v, v))
            if d < 1e-12:
                n = 0
            else:
                n = max(2, int(d / 0.3))
            nn.append(n)
        self.B1 = np.zeros((2, 2, sum(nn), 3))
        self.B2 = np.zeros((2, 2, sum(nn), 3))
        n1 = 0
        for c, n in enumerate(nn):
            n2 = n1 + n
            h = 1.0 / (2 * n - 1)
            R = np.arange(n) * (2 * h)

            for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                self.B1[i, j, n1:n2, c] = R
                self.B1[i, j, n1:n2, (c + 1) % 3] = i
                self.B1[i, j, n1:n2, (c + 2) % 3] = j
            self.B2[:, :, n1:n2] = self.B1[:, :, n1:n2]
            self.B2[:, :, n1:n2, c] += h
            n1 = n2
        self.B1.shape = (-1, 3)
        self.B2.shape = (-1, 3)

    def bind(self, frame):
        if not self.window['toggle-show-bonds']:
            self.bonds = np.empty((0, 5), int)
            return

        from ase.neighborlist import NeighborList
        nl = NeighborList(self.get_covalent_radii() * 1.5,
                          skin=0, self_interaction=False)
        atomscopy = self.atoms.copy()
        atomscopy.cell *= self.images.repeat[:, np.newaxis]
        nl.update(atomscopy)
        nbonds = nl.nneighbors + nl.npbcneighbors

        bonds = np.empty((nbonds, 5), int)
        #self.coordination = np.zeros(len(self.atoms), dtype=int)
        if nbonds == 0:
            return

        n1 = 0
        for a in range(len(self.atoms)):
            indices, offsets = nl.get_neighbors(a)
            #self.coordination[a] += len(indices)
            #for a2 in indices:
            #    self.coordination[a2] += 1
            n2 = n1 + len(indices)
            bonds[n1:n2, 0] = a
            bonds[n1:n2, 1] = indices
            bonds[n1:n2, 2:] = offsets
            n1 = n2

        i = bonds[:n2, 2:].any(1)
        pbcbonds = bonds[:n2][i]
        bonds[n2:, 0] = pbcbonds[:, 1]
        bonds[n2:, 1] = pbcbonds[:, 0]
        bonds[n2:, 2:] = -pbcbonds[:, 2:]
        self.bonds = bonds

    def toggle_show_unit_cell(self, key=None):
        self.set_coordinates()

    def show_labels(self):
        index = self.window['show-labels']
        if index == 0:
            self.labels = None
        elif index == 1:
            self.labels = list(range(len(self.atoms)))
        elif index == 2:
            self.labels = list(self.images.get_magmoms(self.atoms))
        else:
            self.labels = self.atoms.get_chemical_symbols()

        self.draw()

    def toggle_show_axes(self, key=None):
        self.draw()

    def toggle_show_bonds(self, key=None):
        self.set_coordinates()

    def toggle_show_velocities(self, key=None):
        # XXX hard coded scale is ugly
        v = self.atoms.get_velocities()
        if v is not None:
            self.show_vectors(10 * v)
            self.draw()

    # transitional compat hack
    def get_forces(self):
        return self.images.get_forces(self.atoms)

    def toggle_show_forces(self, key=None):
        self.show_vectors(np.nan_to_num(self.get_forces()))
        self.draw()

    def hide_selected(self):
        self.images.visible[self.images.selected] = False
        self.draw()

    def show_selected(self):
        self.images.visible[self.images.selected] = True
        self.draw()

    def repeat_window(self, key=None):
        Repeat(self)

    def rotate_window(self):
        return Rotate(self)

    def colors_window(self, key=None):
        win = ColorWindow(self)
        self.register_vulnerable(win)
        return win

    def focus(self, x=None):
        cell = (self.window['toggle-show-unit-cell'] and
                self.images[0].cell.any())
        if (len(self.atoms) == 0 and not cell):
            self.scale = 1.0
            self.center = np.zeros(3)
            self.draw()
            return

        P = np.dot(self.X, self.axes)
        n = len(self.atoms)
        covalent_radii = self.get_covalent_radii()
        P[:n] -= covalent_radii[:, None]
        P1 = P.min(0)
        P[:n] += 2 * covalent_radii[:, None]
        P2 = P.max(0)
        self.center = np.dot(self.axes, (P1 + P2) / 2)
        S = 1.3 * (P2 - P1)
        w, h = self.window.size
        if S[0] * h < S[1] * w:
            self.scale = h / S[1]
        elif S[0] > 0.0001:
            self.scale = w / S[0]
        else:
            self.scale = 1.0
        self.draw()

    def reset_view(self, menuitem):
        self.axes = rotate('0.0x,0.0y,0.0z')
        self.set_coordinates()
        self.focus(self)

    def set_view(self, key):
        if key == 'Z':
            self.axes = rotate('0.0x,0.0y,0.0z')
        elif key == 'X':
            self.axes = rotate('-90.0x,-90.0y,0.0z')
        elif key == 'Y':
            self.axes = rotate('90.0x,0.0y,90.0z')
        elif key == 'Alt+Z':
            self.axes = rotate('180.0x,0.0y,90.0z')
        elif key == 'Alt+X':
            self.axes = rotate('0.0x,90.0y,0.0z')
        elif key == 'Alt+Y':
            self.axes = rotate('-90.0x,0.0y,0.0z')
        else:
            if key == '3':
                i, j = 0, 1
            elif key == '1':
                i, j = 1, 2
            elif key == '2':
                i, j = 2, 0
            elif key == 'Alt+3':
                i, j = 1, 0
            elif key == 'Alt+1':
                i, j = 2, 1
            elif key == 'Alt+2':
                i, j = 0, 2

            A = complete_cell(self.atoms.cell)
            x1 = A[i]
            x2 = A[j]

            norm = np.linalg.norm

            x1 = x1 / norm(x1)
            x2 = x2 - x1 * np.dot(x1, x2)
            x2 /= norm(x2)
            x3 = np.cross(x1, x2)

            self.axes = np.array([x1, x2, x3]).T

        self.set_coordinates()

    def get_colors(self, rgb=False):
        if rgb:
            return [tuple(int(_rgb[i:i + 2], 16) / 255 for i in range(1, 7, 2))
                    for _rgb in self.get_colors()]

        if self.colormode == 'jmol':
            return [self.colors[Z] for Z in self.atoms.numbers]

        scalars = self.get_color_scalars()
        colorscale, cmin, cmax = self.colormode_data
        N = len(colorscale)
        indices = np.clip(((scalars - cmin) / (cmax - cmin) * N +
                           0.5).astype(int),
                          0, N - 1)
        return [colorscale[i] for i in indices]

    def get_color_scalars(self, frame=None):
        if self.colormode == 'tag':
            return self.atoms.get_tags()
        if self.colormode == 'force':
            f = (self.get_forces()**2).sum(1)**0.5
            return f * self.images.get_dynamic(self.atoms)
        elif self.colormode == 'velocity':
            return (self.atoms.get_velocities()**2).sum(1)**0.5
        elif self.colormode == 'charge':
            return self.atoms.get_charges()
        elif self.colormode == 'magmom':
            return self.images.get_magmoms(self.atoms)

    def get_covalent_radii(self):
        return self.images.get_radii(self.atoms)

    def draw(self, status=True):
        self.window.clear()
        axes = self.scale * self.axes * (1, -1, 1)
        offset = np.dot(self.center, axes)
        offset[:2] -= 0.5 * self.window.size
        X = np.dot(self.X, axes) - offset
        n = len(self.atoms)
        # The indices enumerate drawable objects in z order:
        self.indices = X[:, 2].argsort()
        r = self.get_covalent_radii() * self.scale
        if self.window['toggle-show-bonds']:
            r *= 0.65
        P = self.P = X[:n, :2]
        A = (P - r[:, None]).round().astype(int)
        X1 = X[n:, :2].round().astype(int)
        X2 = (np.dot(self.B, axes) - offset).round().astype(int)
        disp = (np.dot(self.atoms.get_celldisp().reshape((3,)),
                       axes)).round().astype(int)
        d = (2 * r).round().astype(int)

        vectors = (self.window['toggle-show-velocities'] or
                   self.window['toggle-show-forces'])
        if vectors:
            V = np.dot(self.vectors[self.frame], axes) + X[:n]

        colors = self.get_colors()
        circle = self.window.circle
        line = self.window.line
        constrained = ~self.images.get_dynamic(self.atoms)

        selected = self.images.selected
        visible = self.images.visible
        ncell = len(self.B1)
        bond_linewidth = self.scale * 0.15
        for a in self.indices:
            if a < n:
                ra = d[a]
                if visible[a]:
                    # Draw the atoms
                    if self.moving and selected[a]:
                        circle(GREEN, False,
                               A[a, 0] - 4, A[a, 1] - 4,
                               A[a, 0] + ra + 4, A[a, 1] + ra + 4)

                    circle(colors[a], selected[a],
                           A[a, 0], A[a, 1], A[a, 0] + ra, A[a, 1] + ra)

                    # Draw labels on the atoms
                    if self.labels is not None:
                        self.window.text(A[a, 0], A[a, 1],
                                         str(self.labels[a]))

                    # Draw cross on constrained atoms
                    if constrained[a]:
                        R1 = int(0.14644 * ra)
                        R2 = int(0.85355 * ra)
                        line((A[a, 0] + R1, A[a, 1] + R1,
                              A[a, 0] + R2, A[a, 1] + R2))
                        line((A[a, 0] + R2, A[a, 1] + R1,
                              A[a, 0] + R1, A[a, 1] + R2))

                    # Draw velocities or forces
                    if vectors:
                        self.arrow((X[a, 0], X[a, 1], V[a, 0], V[a, 1]),
                                   width=2)
            else:
                # Draw unit cell and/or bonds:
                a -= n
                if a < ncell:
                    line((X1[a, 0] + disp[0], X1[a, 1] + disp[1],
                          X2[a, 0] + disp[0], X2[a, 1] + disp[1]))
                else:
                    line((X1[a, 0] + disp[0], X1[a, 1] + disp[1],
                          X2[a, 0] + disp[0], X2[a, 1] + disp[1]),
                         width=bond_linewidth)

        if self.window['toggle-show-axes']:
            self.draw_axes()

        if len(self.images) > 1:
            self.draw_frame_number()

        self.window.update()

        if status:
            self.status(self.atoms)

    def arrow(self, coords, width):
        line = self.window.line
        begin = np.array((coords[0], coords[1]))
        end = np.array((coords[2], coords[3]))
        line(coords, width)

        vec = end - begin
        length = np.sqrt((vec[:2]**2).sum())
        length = min(length, 0.3 * self.scale)

        angle = np.arctan2(end[1] - begin[1], end[0] - begin[0]) + np.pi
        x1 = (end[0] + length * np.cos(angle - 0.3)).round().astype(int)
        y1 = (end[1] + length * np.sin(angle - 0.3)).round().astype(int)
        x2 = (end[0] + length * np.cos(angle + 0.3)).round().astype(int)
        y2 = (end[1] + length * np.sin(angle + 0.3)).round().astype(int)
        line((x1, y1, end[0], end[1]), width)
        line((x2, y2, end[0], end[1]), width)

    def draw_axes(self):
        axes_length = 15

        rgb = ['red', 'green', 'blue']

        for i in self.axes[:, 2].argsort():
            a = 20
            b = self.window.size[1] - 20
            c = int(self.axes[i][0] * axes_length + a)
            d = int(-self.axes[i][1] * axes_length + b)
            self.window.line((a, b, c, d))
            self.window.text(c, d, 'XYZ'[i], color=rgb[i])

    def draw_frame_number(self):
        x, y = self.window.size
        self.window.text(x, y, '{0}/{1}'.format(self.frame,
                                                len(self.images)),
                         anchor='SE')

    def get_magmoms(self, init_magmom=False):
        try:
            if init_magmom:
                M = self.atoms.get_initial_magnetic_moments()
            else:
                M = self.atoms.get_magnetic_moments()
                if M.ndim == 2:
                    M = M[:, 2]  # XXX
        except (RuntimeError, AttributeError):
            M = self.atoms.get_initial_magnetic_moments()
        return M

    def release(self, event):
        if event.button in [4, 5]:
            self.scroll_event(event)
            return

        if event.button != 1:
            return

        selected = self.images.selected
        selected_ordered = self.images.selected_ordered

        if event.time < self.t0 + 200:  # 200 ms
            d = self.P - self.xy
            r = self.get_covalent_radii()
            hit = np.less((d**2).sum(1), (self.scale * r)**2)
            for a in self.indices[::-1]:
                if a < len(self.atoms) and hit[a]:
                    if event.modifier == 'ctrl':
                        selected[a] = not selected[a]
                        if selected[a]:
                            selected_ordered += [a]
                        elif len(selected_ordered) > 0:
                            if selected_ordered[-1] == a:
                                selected_ordered = selected_ordered[:-1]
                            else:
                                selected_ordered = []
                    else:
                        selected[:] = False
                        selected[a] = True
                        selected_ordered = [a]
                    break
            else:
                selected[:] = False
                selected_ordered = []
            self.draw()
        else:
            A = (event.x, event.y)
            C1 = np.minimum(A, self.xy)
            C2 = np.maximum(A, self.xy)
            hit = np.logical_and(self.P > C1, self.P < C2)
            indices = np.compress(hit.prod(1), np.arange(len(hit)))
            if event.modifier != 'ctrl':
                selected[:] = False
            selected[indices] = True
            if (len(indices) == 1 and
                indices[0] not in self.images.selected_ordered):
                selected_ordered += [indices[0]]
            elif len(indices) > 1:
                selected_ordered = []
            self.draw()

        # XXX check bounds
        indices = np.arange(len(self.atoms))[self.images.selected[:len(self.atoms)]]
        if len(indices) != len(selected_ordered):
            selected_ordered = []
        self.images.selected_ordered = selected_ordered

    def press(self, event):
        self.button = event.button
        self.xy = (event.x, event.y)
        self.t0 = event.time
        self.axes0 = self.axes
        self.center0 = self.center

    def move(self, event):
        x = event.x
        y = event.y
        x0, y0 = self.xy
        if self.button == 1:
            x0 = int(round(x0))
            y0 = int(round(y0))
            self.draw()
            self.window.canvas.create_rectangle((x, y, x0, y0))
            return

        if event.modifier == 'shift':
            self.center = (self.center0 -
                           np.dot(self.axes, (x - x0, y0 - y, 0)) / self.scale)
        else:
            # Snap mode: the a-b angle and t should multipla of 15 degrees ???
            a = x - x0
            b = y0 - y
            t = sqrt(a * a + b * b)
            if t > 0:
                a /= t
                b /= t
            else:
                a = 1.0
                b = 0.0
            c = cos(0.01 * t)
            s = -sin(0.01 * t)
            rotation = np.array([(c * a * a + b * b, (c - 1) * b * a, s * a),
                                 ((c - 1) * a * b, c * b * b + a * a, s * b),
                                 (-s * a, -s * b, c)])
            self.axes = np.dot(self.axes0, rotation)
            if len(self.atoms) > 0:
                com = self.X_pos.mean(0)
            else:
                com = self.atoms.cell.mean(0)
            self.center = com - np.dot(com - self.center0,
                                       np.dot(self.axes0, self.axes.T))
        self.draw(status=False)

    def render_window(self, action):
        Render(self)

    def show_vectors(self, vectors):
        self.vectors = vectors

    def resize(self, event):
        w, h = self.window.size
        self.scale *= (event.width * event.height / (w * h))**0.5
        self.window.size[:] = [event.width, event.height]
        self.draw()
