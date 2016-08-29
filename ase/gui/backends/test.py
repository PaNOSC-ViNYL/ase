import numpy as np


def name2str(name):
    return '-'.join(x.lower() for x in name.replace('_', '').split())


class MainWindow:
    def __init__(self, menu_description, config,
                 exit, scroll, scroll_event,
                 press, move, release):
        self.size = np.array([600, 600])

        for name, things in menu_description:
            for thing in things:
                if len(things) == 5 and isinstance(things[4], bool):
                    subname, key, text, callback, on = thing
                    # check
                elif len(thing) == 5:
                    subname, key, text, things, callback = thing
                    #submenu.add_radio(label=_(subname),
                    #                  command=callback)
                    #self.menu[name2str(subname)] = 0

    def update_status_line(self, text):
        pass

    def resize_event(self):
        self.scale *= sqrt(1.0 * self.width * self.height / (w * h))
        self.draw()
        self.configured = True

    def run(self):
        pass

    def __getitem__(self, name):
        return False

    def clear(self):
        pass

    def update(self):
        pass

    def circle(self, color, selected, *bbox):
        pass

    def line(self, bbox):
        pass

    def text(self, x, y, txt, anchor=None):
        pass
