try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

from gettext import gettext as _

import numpy as np


def name2str(name):
    return '-'.join(x.lower() for x in name.replace('_', '').split())


class MainWindow:
    def __init__(self, menu_description, config, a, b, c):
        self.size = np.array([600, 600])

        self.root = tk.Tk()

        #self.window.set_position(gtk.WIN_POS_CENTER)
        #self.window.connect('delete_event', self.exit)

        menu = tk.Menu(self.root)
        self.root.config(menu=menu)

        for name, things in menu_description:
            submenu = tk.Menu(menu)
            menu.add_cascade(label=_(name), menu=submenu)
            for thing in things:
                if thing == '---':
                    submenu.add_separator()
                elif len(things) == 5 and isinstance(things[4], bool):
                    subname, key, text, callback, on = thing
                    # check
                elif len(thing) == 5:
                    subname, key, text, things, callback = thing
                    #submenu.add_radio(label=_(subname),
                    #                  command=callback)
                    #self.menu[name2str(subname)] = 0
                elif isinstance(thing[3], list):
                    subname, key, text, subthings = thing
                    subsubmenu = tk.Menu(submenu)
                    submenu.add_cascade(label=_(subname), menu=subsubmenu)
                    for subsubname, key, text, callback in subthings:
                        subsubmenu.add_command(label=_(subsubname),
                                               command=callback,
                                               accelerator=key)
                else:
                    subname, key, text, callback = thing
                    submenu.add_command(label=_(subname),
                                        command=callback,
                                        accelerator=key)

        self.canvas = tk.Canvas(self.root,
                                     width=self.size[0],
                                     height=self.size[1],
                                     bg='white')
        self.canvas.pack(side=tk.TOP, fill=tk.X)# & tk.Y)

        status = tk.Label(self.root, text="asdgag", #bd=1,
                               #relief=tk.SUNKEN,
                               anchor=tk.W)
        status.pack(side=tk.BOTTOM, fill=tk.X)
        #self.window.connect('key-press-event', self.scroll)
        #self.window.connect('scroll_event', self.scroll_event)
        #self.drawing_area.connect('button_press_event', self.press)
        #self.drawing_area.connect('button_release_event', self.release)
        #self.drawing_area.connect('motion-notify-event', self.move)
        #self.drawing_area.connect('expose_event', self.expose_event)
        #self.drawing_area.connect('configure_event', self.configure_event)
        #self.drawing_area.set_events(gtk.gdk.BUTTON_PRESS_MASK |
        #                             gtk.gdk.BUTTON_RELEASE_MASK |
        #                             gtk.gdk.BUTTON_MOTION_MASK |
        #                             gtk.gdk.POINTER_MOTION_HINT_MASK)

        #    self.eventbox.set_tooltip_text(_('Tip for status box ...'))

        self.fg = config['gui_foreground_color']
        self.bg = config['gui_background_color']
        self.red = '#F30000'
        self.green = '#00D300'
        self.blue = '#0000D3'
        self.selected_color = '#004500'

    def update_status_line(self, text):
        pass

    def resize_event(self):
        self.scale *= sqrt(1.0 * self.width * self.height / (w * h))
        self.draw()
        self.configured = True

    def run(self):
        tk.mainloop()

    def __getitem__(self, name):
        return False

    def title(self, txt):
        self.root.title(txt)

    title = property(None, title)

    def clear(self):
        self.canvas.delete(tk.ALL)

    def update(self):
        self.canvas.update_idletasks()

    def circle(self, color, bbox):
        self.canvas.create_oval(*tuple(int(x) for x in bbox), fill=color)

    def line(self, bbox):
        self.canvas.create_line(*tuple(int(x) for x in bbox))

    def text(self, x, y, txt, anchor=tk.CENTER):
        anchor = {'SE': tk.SE}.get(anchor, anchor)
        self.canvas.create_text((x, y), text=txt, anchor=anchor)


class Window:
    def __init__(self, stuff):
        for line in stuff:
            box = Box()
            pack


class Button:
    pass

"""
Entry
Text
"""
