try:
    # Python 3
    import tkinter as tk
    import tkinter.ttk as ttk
    from tkinter.messagebox import askokcancel as ask_question
    from tkinter.messagebox import showerror as oops
except ImportError:
    # Python 2
    import Tkinter as tk
    import ttk
    from tkMessageBox import askokcancel as ask_question

import re
from functools import partial
from gettext import gettext as _

import numpy as np

from ase.utils import basestring

__all__ = [
    'oops', 'ask_question', 'MainWindow',
    'ASEGUIWindow', 'Button', 'CheckButton', 'ComboBox', 'Entry', 'Label',
    'Window', 'MenuItem', 'RadioButton', 'RadioButtons', 'Rows', 'Scale',
    'SpinBox', 'Text']

font = ('Helvetica', 10)


def about(name, version, webpage):
    text = [name,
            '',
            _('Version') + ': ' + version,
            _('Web-page') + ': ' + webpage]
    win = Window(_('About'))
    win.add(Text('\n'.join(text)))


def helpbutton(text):
    return Button(_('Help'), helpwindow, text)


def helpwindow(text):
    win = Window(_('Help'))
    win.add(Text(text))


class BaseWindow:
    def __init__(self, title, close=None):
        self.title = title
        if close:
            self.win.protocol('WM_DELETE_WINDOW', close)
        else:
            self.win.protocol('WM_DELETE_WINDOW', self.close)

        self.things = []

    def close(self):
        self.win.destroy()

    def title(self, txt):
        self.win.title(txt)

    title = property(None, title)

    def add(self, stuff, anchor='w'):  # 'center'):
        if isinstance(stuff, basestring):
            stuff = Label(stuff)
        elif isinstance(stuff, list):
            stuff = Row(stuff)
        stuff.pack(self.win, anchor=anchor)
        self.things.append(stuff)


class Window(BaseWindow):
    def __init__(self, title, close=None):
        self.win = tk.Toplevel()
        BaseWindow.__init__(self, title, close)


class Widget(object):
    def pack(self, parent, side='top', anchor='center'):
        widget = self.create(parent)
        widget.pack(side=side, anchor=anchor)
        if not isinstance(self, (Rows, RadioButtons)):
            pass  # widget['font'] = font

    def grid(self, parent):
        widget = self.create(parent)
        widget.grid()
        # widget['font'] = font

    def create(self, parent):
        self.widget = self.creator(parent)
        return self.widget

    @property
    def active(self):
        return self.widget['state'] == 'normal'

    @active.setter
    def active(self, value):
        self.widget['state'] = ['disabled', 'normal'][bool(value)]


class Row(Widget):
    def __init__(self, things):
        self.things = things

    def create(self, parent):
        self.widget = tk.Frame(parent)
        for thing in self.things:
            if isinstance(thing, basestring):
                thing = Label(thing)
            thing.pack(self.widget, 'left')
        return self.widget

    def __getitem__(self, i):
        return self.things[i]


class Label(Widget):
    def __init__(self, text='', color=None):
        self.creator = partial(tk.Label, text=text, fg=color)

    def text(self, new):
        self.widget.config(text=new)

    text = property(None, text)


class Text(Widget):
    def __init__(self, text):
        self.creator = partial(tk.Text, height=text.count('\n') + 1)
        s = re.split('<(.*?)>', text)
        self.text = [(s[0], ())]
        i = 1
        tags = []
        while i < len(s):
            tag = s[i]
            if tag[0] != '/':
                tags.append(tag)
            else:
                tags.pop()
            self.text.append((s[i + 1], tuple(tags)))
            i += 2

    def create(self, parent):
        widget = Widget.create(self, parent)
        widget.tag_configure('sub', offset=-6)
        widget.tag_configure('sup', offset=6)
        widget.tag_configure('c', foreground='blue')
        for text, tags in self.text:
            widget.insert('insert', text, tags)
        widget.configure(state='disabled', background=parent['bg'])
        return widget


class Button(Widget):
    def __init__(self, text, callback, *args, **kwargs):
        self.callback = partial(callback, *args, **kwargs)
        self.creator = partial(tk.Button,
                               text=text,
                               command=self.callback)


class CheckButton(Widget):
    def __init__(self, text, value=False):
        self.text = text
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
        self.widget = self.creator(parent)
        self.value = self.initial
        return self.widget

    @property
    def value(self):
        x = self.widget.get().replace(',', '.')
        if '.' in x:
            return float(x)
        return int(x)

    @value.setter
    def value(self, x):
        self.widget.delete(0, 'end')
        self.widget.insert(0, x)


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


class RadioButtons(Widget):
    def __init__(self, labels, values=None, toggle=None, vertical=True):
        self.var = tk.IntVar()

        if toggle:
            def callback():
                toggle(self.value)
        else:
            callback = None

        self.values = values or list(range(len(labels)))
        self.buttons = [RadioButton(label, i, self.var, callback)
                        for i, label in enumerate(labels)]

    def create(self, parrent):
        frame = tk.Frame(parrent)
        for button in self.buttons:
            button.create(frame).pack(side='top')  # 'left')
        return frame

    @property
    def value(self):
        return self.values[self.var.get()]

    @value.setter
    def value(self, value):
        self.var.set(self.values.index(value))

    def __getitem__(self, value):
        return self.buttons[self.values.index(value)]


class RadioButton(Widget):
    def __init__(self, label, i, var, callback):
        self.creator = partial(tk.Radiobutton,
                               text=label,
                               var=var,
                               value=i,
                               command=callback)


class ComboBox(Widget):
    def __init__(self, labels, values=None, callback=None):
        self.var = tk.StringVar()
        self.values = dict(zip(labels, values or list(range(len(labels)))))
        self.callback = callback
        self.creator = partial(ttk.Combobox,
                               textvariable=self.var, values=labels)

    def create(self, parrent):
        widget = Widget.create(self, parrent)
        widget.current(0)
        if self.callback:
            def callback(event):
                self.callback(self.value)
            widget.bind('<<ComboboxSelected>>', callback)
        return widget

    @property
    def value(self):
        return self.values[self.var.get()]


class Rows(Widget):
    def __init__(self, rows=None):
        self.rows_to_be_added = rows or []
        self.creator = tk.Frame
        self.rows = []

    def create(self, parent):
        widget = Widget.create(self, parent)
        for row in self.rows_to_be_added:
            self.add(row)
        self.rows_to_be_added = []
        return widget

    def add(self, row):
        if isinstance(row, basestring):
            row = Label(row)
        elif isinstance(row, list):
            row = Row(row)
        row.grid(self.widget)
        self.rows.append(row)

    def clear(self):
        while self.rows:
            del self[0]

    def __getitem__(self, i):
        return self.rows[i]

    def __delitem__(self, i):
        widget = self.rows.pop(i).widget
        widget.grid_remove()
        widget.destroy()

    def __len__(self):
        return len(self.rows)


class MenuItem:
    def __init__(self, label, callback=None, key=None,
                 value=None, choices=None, submenu=None):
        self.underline = label.find('_')
        self.label = label.replace('_', '')

        if key:
            if key[:4] == 'Ctrl':
                self.keyname = '<Control-{0}>'.format(key[-1].lower())
            else:
                self.keyname = {
                    'Home': '<Home>',
                    'End': '<End>',
                    'Page-Up': '<Prior>',
                    'Page-Down': '<Next>',
                    'Backspace': '<BackSpace>'}.get(key, key.lower())

        if key:
            def callback2(event=None):
                callback(key)

            callback2.__name__ = callback.__name__
            self.callback = callback2
        else:
            self.callback = callback

        self.key = key
        self.value = value
        self.choices = choices
        self.submenu = submenu

    def addto(self, menu, window, stuff=None):
        if self.label == '---':
            menu.add_separator()
        elif self.value is not None:
            var = tk.BooleanVar(value=self.value)
            stuff[self.callback.__name__.replace('_', '-')] = var
            menu.add_checkbutton(label=self.label,
                                 underline=self.underline,
                                 command=self.callback,
                                 accelerator=self.key,
                                 var=var)
        elif self.choices:
            submenu = tk.Menu(menu)
            menu.add_cascade(label=self.label, menu=submenu)
            var = tk.IntVar()
            var.set(0)
            stuff[self.callback.__name__.replace('_', '-')] = var
            for i, choice in enumerate(self.choices):
                submenu.add_radiobutton(label=choice.replace('_', ''),
                                        underline=choice.find('_'),
                                        command=self.callback,
                                        value=i,
                                        var=var)
        elif self.submenu:
            submenu = tk.Menu(menu)
            menu.add_cascade(label=self.label,
                             menu=submenu)
            for thing in self.submenu:
                thing.addto(submenu, window)
        else:
            menu.add_command(label=self.label,
                             underline=self.underline,
                             command=self.callback,
                             accelerator=self.key)
        if self.key:
            window.bind(self.keyname, self.callback)


class MainWindow(BaseWindow):
    def __init__(self, title, close=None, menu=[]):
        self.win = tk.Tk()
        BaseWindow.__init__(self, title, close)

        # self.win.tk.call('tk', 'scaling', 3.0)
        # self.win.tk.call('tk', 'scaling', '-displayof', '.', 7)

        self.menu = {}

        if menu:
            self.create_menu(menu)

    def create_menu(self, menu_description):
        menu = tk.Menu(self.win, font=font)
        self.win.config(menu=menu)

        for label, things in menu_description:
            submenu = tk.Menu(menu, font=font)
            menu.add_cascade(label=label.replace('_', ''),
                             underline=label.find('_'),
                             menu=submenu)
            for thing in things:
                thing.addto(submenu, self.win, self.menu)

    def resize_event(self):
        # self.scale *= sqrt(1.0 * self.width * self.height / (w * h))
        self.draw()
        self.configured = True

    def run(self):
        tk.mainloop()

    def test(self, test, close_after_test=False):
        def callback():
            try:
                next(test)
            except StopIteration:
                if close_after_test:
                    self.close()
            else:
                self.win.after_idle(callback)

        self.win.after_idle(test)#callback)
        self.run()

    def __getitem__(self, name):
        return self.menu[name].get()

    def __setitem__(self, name, value):
        return self.menu[name].set(value)


def bind(callback, modifier=None):
    def handle(event):
        event.button = event.num
        event.key = event.keysym.lower()
        event.modifier = modifier
        callback(event)
    return handle


class ASEGUIWindow(MainWindow):
    def __init__(self, close, menu, config,
                 scroll, scroll_event,
                 press, move, release):
        MainWindow.__init__(self, 'ASE-GUI', close, menu)

        self.size = np.array([450, 450])

        self.canvas = tk.Canvas(self.win,
                                width=self.size[0],
                                height=self.size[1],
                                bg='white')
        self.canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.status = tk.Label(self.win, text='',  # bd=1,
                               # relief=tk.SUNKEN,
                               anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas.bind('<ButtonPress>', bind(press))
        self.canvas.bind('<B1-Motion>', bind(move))
        self.canvas.bind('<B3-Motion>', bind(move))
        self.canvas.bind('<ButtonRelease>', bind(release))
        self.canvas.bind('<Control-ButtonRelease>', bind(release, 'ctrl'))
        self.canvas.bind('<Shift-ButtonRelease>', bind(release, 'shift'))
        self.win.bind('<Key>', bind(scroll))
        self.win.bind('<Shift-Key>', bind(scroll, 'shift'))
        self.win.bind('<Control-Key>', bind(scroll, 'ctrl'))
        # self.canvas.bind('<B4>', bind(scroll_event))
        # self.canvas.bind('<Shift-MouseWheel>', bind(scroll_event, 'shift'))
        # self.win.bind('<Configure>', configure_event)
        # self.drawing_area.connect('expose_event', self.expose_event)

        #    self.eventbox.set_tooltip_text(_('Tip for status box ...'))

        self.fg = config['gui_foreground_color']
        self.bg = config['gui_background_color']

    def update_status_line(self, text):
        self.status.config(text=text)

    def run(self):
        MainWindow.run(self)

    def click(self, name):
        self.callbacks[name]()

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
