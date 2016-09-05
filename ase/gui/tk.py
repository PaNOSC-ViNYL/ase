try:
    import tkinter as tk
    from tkinter.messagebox import askokcancel
except ImportError:
    import Tkinter as tk
    from tkMessageBox import askokcancel

from functools import partial
from gettext import gettext

import numpy as np


font = ('Helvetica', 10)
        

def ask_question(title, question):
    return askokcancel(gettext(title), gettext(question))
    

def parse(name, key):
    label = gettext(name)
    name = name.replace('_', '').replace('.', '').strip()
    id = '-'.join(x.lower() for x in name.split())
    underline = -1
    if key[:4] == 'Ctrl':
        key = key[-1].lower()
        underline = label.lower().find(key)
        keyname = '<Control-{0}>'.format(key)
    else:
        keyname = {'Home': '<Home>',
                   'End': '<End>',
                   'Page-Up': '<Prior>',
                   'Page-Down': '<Next>',
                   'Backspace': '<BackSpace>'}.get(key, key.lower())
    return id, label, underline, keyname


class MainWindow:
    def __init__(self, menu_description, config,
                 exit, scroll, scroll_event,
                 press, move, release):
        self.size = np.array([450, 450])

        self.root = tk.Tk()
        # self.root.tk.call('tk', 'scaling', 3.0)
        # self.root.tk.call('tk', 'scaling', '-displayof', '.', 7)

        self.root.protocol('WM_DELETE_WINDOW', exit)

        menu = tk.Menu(self.root, font=font)
        self.root.config(menu=menu)

        self.menu = {}
        self.callbacks = {}

        for name, things in menu_description:
            submenu = tk.Menu(menu, font=font)
            label = gettext(name)
            menu.add_cascade(label=label, menu=submenu)
            for thing in things:
                if thing == '---':
                    submenu.add_separator()
                    continue
                subname, key, text, callback = thing[:4]
                id, label, underline, keyname = parse(subname, key)
                self.callbacks[id] = callback
                if len(thing) == 4:
                    submenu.add_command(label=label,
                                        underline=underline,
                                        command=callback,
                                        accelerator=key)
                    if key:
                        print(keyname, callback)
                        self.root.bind(keyname, callback)
                    continue
                    
                x = thing[4]
                if isinstance(x, bool):
                    on = x
                    var = tk.BooleanVar(value=on)
                    self.menu[id] = var
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
                        id, label, underline, keyname = parse(subsubname, key)
                        subsubmenu.add_command(label=label,
                                               underline=underline,
                                               command=callback,
                                               accelerator=key)
                        if key:
                            self.root.bind(keyname, callback)

        self.canvas = tk.Canvas(self.root,
                                width=self.size[0],
                                height=self.size[1],
                                bg='white')
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.status = tk.Label(self.root, text='',  # bd=1,
                               # relief=tk.SUNKEN,
                               anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas.bind('<ButtonPress>', bind(press))
        self.canvas.bind('<B1-Motion>', bind(move))
        self.canvas.bind('<B3-Motion>', bind(move))
        self.canvas.bind('<ButtonRelease>', bind(release))
        self.canvas.bind('<Control-ButtonRelease>',
                         bind(release, 'ctrl'))
        self.root.bind('<Key>', bind(scroll))
        self.root.bind('<Shift-Key>', bind(scroll, 'shift'))
        self.root.bind('<Control-Key>', bind(scroll, 'ctrl'))
        # self.canvas.bind('<B4>', bind(scroll_event))
        # self.canvas.bind('<Shift-MouseWheel>', bind(scroll_event, 'shift'))
        # self.root.bind('<Configure>', configure_event)
        # self.drawing_area.connect('expose_event', self.expose_event)

        #    self.eventbox.set_tooltip_text(_('Tip for status box ...'))

        self.fg = config['gui_foreground_color']
        self.bg = config['gui_background_color']

    def close(self):
        self.root.destroy()

    def update_status_line(self, text):
        self.status.config(text=text)

    def resize_event(self):
        # self.scale *= sqrt(1.0 * self.width * self.height / (w * h))
        self.draw()
        self.configured = True

    def run(self, click):
        if click:
            self.root.after_idle(self.click, click)
        tk.mainloop()

    def click(self, name):
        self.callbacks[name]()

    def __getitem__(self, name):
        return self.menu[name].get()

    def __setitem__(self, name, value):
        return self.menu[name].set(value)

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
        self.canvas.create_text((x, y), text=txt, anchor=anchor, fill=color,
                                font=font)


def bind(callback, modifier=None):
    def handle(event):
        event.button = event.num
        event.key = event.keysym.lower()
        event.modifier = modifier
        callback(event)
    return handle


class Window:
    def __init__(self, title, close=None):
        self.top = tk.Toplevel()
        self.top.title(gettext(title))
        if close:
            self.top.protocol('WM_DELETE_WINDOW', close)
        else:
            self.top.protocol('WM_DELETE_WINDOW', self.close)

    def close(self):
        self.top.destroy()

    def add(self, stuff, anchor='center'):
        if isinstance(stuff, str):
            stuff = Label(stuff)
        if isinstance(stuff, list):
            frame = tk.Frame(self.top)
            for thing in stuff:
                if isinstance(thing, str):
                    thing = Label(thing)
                thing.pack(frame, 'left')
            frame.pack(anchor=anchor)
        else:
            stuff.pack(self.top, anchor=anchor)


class Widget(object):
    def pack(self, parent, side='top', anchor='center'):
        widget = self.create(parent)
        widget.pack(side=side, anchor=anchor)
        widget['font'] = font

    def create(self, parent):
        self.widget = self.creator(parent)
        return self.widget


class Label(Widget):
    def __init__(self, text):
        self.text = gettext(text)

    def create(self, parent):
        return tk.Label(parent, text=self.text)


class Button(Widget):
    def __init__(self, text, on_press, *args, **kwargs):
        self.text = gettext(text)
        self.on_press = partial(on_press, *args, **kwargs)

    def create(self, parent):
        return tk.Button(parent, text=self.text, command=self.on_press)


class CheckButton(Widget):
    def __init__(self, text, value=False):
        self.text = gettext(text)
        self.var = tk.BooleanVar(value=value)

    def create(self, parent):
        self.check = tk.Checkbutton(parent, text=self.text, var=self.var)
        return self.check

    @property
    def value(self):
        return self.var.get()


class SpinBox(Widget):
    def __init__(self, value, start, end, step, on_change=None):
        self.creator = partial(tk.Spinbox,
                               from_=start,
                               to=end,
                               increment=step,
                               command=on_change,
                               width=6)
        self.initial = str(value)

    def create(self, parent):
        self.spin = self.creator(parent)
        self.value = self.initial
        return self.spin

    @property
    def value(self):
        return float(self.spin.get().replace(',', '.'))

    @value.setter
    def value(self, x):
        self.spin.delete(0, 'end')
        self.spin.insert(0, x)


class Entry(Widget):
    def __init__(self, value='', width=20, callback=None):
        self.creator = partial(tk.Entry,
                               width=width)
        self.callback = lambda event: callback()
        self.initial = value

    def create(self, parent):
        self.entry = self.creator(parent)
        self.value = self.initial
        self.entry.bind('<Return>', self.callback)
        return self.entry

    @property
    def value(self):
        return self.entry.get()

    @value.setter
    def value(self, x):
        self.entry.delete(0, 'end')
        self.entry.insert(0, x)


class Scale(Widget):
    def __init__(self, value, start, end, callback):
        def command(val):
            callback(int(val))

        self.creator = partial(tk.Scale,
                               from_=start,
                               to=end,
                               orient='horizontal',
                               command=command)
        self.initial = value

    def create(self, parent):
        self.scale = self.creator(parent)
        self.value = self.initial
        return self.scale

    @property
    def value(self):
        return self.scale.get()

    @value.setter
    def value(self, x):
        self.scale.set(x)

        
class RadioButtons(list):
    def __init__(self, labels, values=None, toggle=None):
        self.var = tk.IntVar
        
        def callback():
            toggle(self.value)
            
        self.values = values or list(range(len(labels)))
        list.__init__(self, [RadioButton(label, i, self.var, callback)
                             for i, label in enumerate(labels)])
        
    @property
    def value(self):
        return self.values[self.var.get()]

    
class RadioButton(Widget):
    def __init__(self, label, i, var, callback):
        self.creator = partial(tk.Radiobutton,
                               text=gettext(label),
                               var=var,
                               value=i,
                               callback=callback)
