try:
    import tkinter
except ImportError:
    import Tkinter as tkinter
    
from gettext import gettext as _

import numpy as np

from ase.gui.defaults import read_defaults


class MainWindow:
    def __init__(self, menu_description, a, b, c):
        self.size = np.array([600, 600])
        
        root = tkinter.Tk()

        #self.window.set_position(gtk.WIN_POS_CENTER)
        #self.window.connect('delete_event', self.exit)

        menu = tkinter.Menu(root)
        root.config(menu=menu)

        for name, things in menu_description:
            submenu = tkinter.Menu(menu)
            menu.add_cascade(label=_(name), menu=submenu)
            for thing in things:
                if thing == '---':
                    submenu.add_separator()
                else:
                    subname, key, text, callback = thing
                    submenu.add_command(label=_(name), command=callback)
        
        self.canvas = tkinter.Canvas(root,
                                     width=self.size[0],
                                     height=self.size[1],
                                     bg='white')
        self.canvas.pack(side=tkinter.TOP, fill=tkinter.X)# & tkinter.Y)
        
        status = tkinter.Label(root, text="asdgag", bd=1,
                               relief=tkinter.SUNKEN, anchor=tkinter.W)
        status.pack(side=tkinter.BOTTOM, fill=tkinter.X)
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

        config = read_defaults()
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
        tkinter.mainloop()
        
    def __getitem__(self, name):
        return False

    def title(self, txt):
        pass
        
    title = property(None, title)
    
    def clear(self):
        self.canvas.delete(tkinter.ALL)

    def update(self):
        self.canvas.update_idletasks()
        
    def circle(self, color, bbox):
        self.canvas.create_oval(*bbox, fill=color)

    def line(self):
        self.canvas.create_line(0, 0, 200, 100)
        
        
class Button:
    pass
