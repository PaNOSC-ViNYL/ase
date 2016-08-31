try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

from gettext import gettext

import numpy as np


def name2str(name):
    return '-'.join(x.lower() for x in name.replace('_', '').split())


def parselabel(label):
    label = gettext(label)
    i = label.find('_')
    if i >= 0:
        return i, label.replace('_', '')
    return None, label


class MainWindow:
    def __init__(self, menu_description, config,
                 exit, scroll, scroll_event,
                 press, move, release):
        self.size = np.array([600, 600])

        self.root = tk.Tk()

        # self.window.set_position(gtk.WIN_POS_CENTER)
        # self.window.connect('delete_event', self.exit)

        menu = tk.Menu(self.root)
        self.root.config(menu=menu)

        self.menu = {}

        for name, things in menu_description:
            submenu = tk.Menu(menu)
            underline, label = parselabel(name)
            menu.add_cascade(label=label, underline=underline, menu=submenu)
            for thing in things:
                if thing == '---':
                    submenu.add_separator()
                    continue
                subname, key, text, callback = thing[:4]
                underline, label = parselabel(subname)
                if len(thing) == 4:
                    submenu.add_command(label=label,
                                        underline=underline,
                                        command=callback,
                                        accelerator=key)
                    continue
                x = thing[4]
                if isinstance(x, bool):
                    on = x
                    var = tk.BooleanVar(value=on)
                    self.menu[name2str(subname)] = var
                    submenu.add_checkbutton(label=label,
                                            underline=underline,
                                            command=callback,
                                            accelerator=key,
                                            var=var)

                elif isinstance(x[0], str):
                    pass  # hmm = x
                    # submenu.add_radio(label=_(subname),
                    #                   command=callback)
                else:
                    subsubmenu = tk.Menu(submenu)
                    submenu.add_cascade(label=gettext(subname),
                                        menu=subsubmenu)
                    for subsubname, key, text, callback in x:
                        subsubmenu.add_command(label=gettext(subsubname),
                                               command=callback,
                                               accelerator=key)

        self.canvas = tk.Canvas(self.root,
                                width=self.size[0],
                                height=self.size[1],
                                bg='white')
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        status = tk.Label(self.root, text="asdgag",  # bd=1,
                          # relief=tk.SUNKEN,
                          anchor=tk.W)
        status.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas.bind('<ButtonPress>', Handler(press))
        self.canvas.bind('<B1-Motion>', Handler(move))
        self.canvas.bind('<B3-Motion>', Handler(move))
        self.canvas.bind('<ButtonRelease>', Handler(release))
        self.canvas.bind('<Control-ButtonRelease>',
                         Handler(release, 'ctrl'))
        self.root.bind('<Key>', Handler(scroll))
        self.root.bind('<Shift-Key>', Handler(scroll, 'shift'))
        self.root.bind('<Control-Key>', Handler(scroll, 'ctrl'))
        #self.canvas.bind('<B4>', Handler(scroll_event))
        #self.canvas.bind('<Shift-MouseWheel>', Handler(scroll_event, 'shift'))
        # self.root.bind('<Configure>', configure_event)
        # self.drawing_area.connect('expose_event', self.expose_event)

        #    self.eventbox.set_tooltip_text(_('Tip for status box ...'))

        self.fg = config['gui_foreground_color']
        self.bg = config['gui_background_color']

    def update_status_line(self, text):
        pass

    def resize_event(self):
        # self.scale *= sqrt(1.0 * self.width * self.height / (w * h))
        self.draw()
        self.configured = True

    def run(self):
        tk.mainloop()

    def __getitem__(self, name):
        return self.menu[name].get()

    def title(self, txt):
        self.root.title(txt)

    title = property(None, title)

    def clear(self):
        self.canvas.delete(tk.ALL)

    def update(self):
        self.canvas.update_idletasks()

    def circle(self, color, selected, *bbox):
        if selected:
            outline = '#004500'
            width = 3
        else:
            outline = 'black'
            width = 1
        self.canvas.create_oval(*tuple(int(x) for x in bbox), fill=color,
                                outline=outline, width=width)

    def line(self, bbox, width=1):
        self.canvas.create_line(*tuple(int(x) for x in bbox), width=width)

    def text(self, x, y, txt, anchor=tk.CENTER, color='black'):
        anchor = {'SE': tk.SE}.get(anchor, anchor)
        self.canvas.create_text((x, y), text=txt, anchor=anchor, fill=color)


class Handler:
    def __init__(self, callback, modifier=None):
        self.callback = callback
        self.modifier = modifier

    def __call__(self, event):
        self.callback(Event(event, self.modifier))


class Event:
    def __init__(self, event, modifier):
        print(event.delta, event.keysym)
        self.x = event.x
        self.y = event.y
        self.time = event.time
        self.button = event.num
        self.modifier = modifier
        self.key = event.keysym.lower()

        if 1:
            print(self.x,
                  self.y,
                  self.time,
                  self.button,
                  self.modifier)


class Window:
    def __init__(self, stuff):
        for line in stuff:
            pass  # box = Box()
            # pack


class Button:
    pass

"""
Entry
Text
"""
