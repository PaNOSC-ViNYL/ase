import threading

import numpy as np

import ase.gui.ui as ui


class Movie:
    def __init__(self, gui):
        self.win = win = ui.Window('Movie', self.close)
        win.add('Image number:')
        self.frame_number = ui.Scale(gui.frame, 0,
                                     gui.images.nimages - 1,
                                     callback=self.new_frame)
        win.add(self.frame_number)

        win.add([ui.Button('First', self.click, -1, True),
                 ui.Button('Back', self.click, -1),
                 ui.Button('Forward', self.click, 1),
                 ui.Button('Last', self.click, 1, True)])

        play = ui.Button('Play', self.play)
        stop = ui.Button('Stop', self.stop)
        
        # TRANSLATORS: This function plays an animation forwards and backwards
        # alternatingly, e.g. for displaying vibrational movement
        self.rock = ui.CheckButton('Rock')

        win.add([play, stop, self.rock])

        if gui.images.nimages > 150:
            skipdefault = gui.images.nimages // 150
            tdefault = min(max(gui.images.nimages / (skipdefault * 5.0),
                               1.0), 30)
        else:
            skipdefault = 0
            tdefault = min(max(gui.images.nimages / 5.0, 1.0), 30)
        self.time = ui.SpinBox(tdefault, 1.0, 99, 0.1)
        self.skip = ui.SpinBox(skipdefault, 0, 99, 1)
        win.add([' Frame rate: ', self.time, ' Skip frames: ', self.skip])

        self.gui = gui
        self.direction = 1
        self.timer = None
        gui.register_vulnerable(self)

    def notify_atoms_changed(self):
        """Called by gui object when the atoms have changed."""
        self.close()

    def close(self):
        self.stop()
        self.win.close()

    def click(self, step, firstlast=False):
        if firstlast and step < 0:
            i = 0
        elif firstlast:
            i = self.gui.images.nimages - 1
        else:
            i = max(0, min(self.gui.images.nimages - 1, self.gui.frame + step))
        self.gui.set_frame(i)
        self.frame_number.value = i
        if firstlast:
            self.direction = np.sign(-step)
        else:
            self.direction = np.sign(step)

    def new_frame(self, value):
        self.gui.set_coordinates(value)

    def play(self):
        if self.timer is not None:
            self.timer.cancel()
        t = 1 / self.time.value
        self.timer = threading.Timer(t, self.step)
        self.timer.start()

    def stop(self):
        if self.timer is not None:
            self.timer.cancel()

    def step(self):
        i = self.gui.frame
        nimages = self.gui.images.nimages
        delta = int(self.skip.value + 1)

        if self.rock.value:
            if i <= self.skip.value:
                self.direction = 1
            elif i >= nimages - delta:
                self.direction = -1
            i += self.direction * delta
        else:
            i = (i + self.direction * delta + nimages) % nimages

        self.frame_number.value = i
        self.play()
