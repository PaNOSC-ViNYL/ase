from __future__ import print_function
import sys
import ase.gui.ui as ui
import re




class Help:
    __instance = None

    def __new__(cls, *args, **kwargs):
        # Make this a singleton.
        if Help.__instance is None:
            Help.__instance = ui.Window.__new__(cls, *args, **kwargs)
        return Help.__instance

    def __init__(self, text):
        # Now, __init__ may be called multiple times!
        if not hasattr(self, '_initialized'):
            self.initialize(text)
        else:
            self.set_text(text)
        self.present()  # Show the window.
        
    def initialize(self, text):
        ui.Window.__init__(self)
        self.set_title(_("Help"))
        self._initialized = True
        vbox = ui.VBox()
        self.add(vbox)
        self.label = pack(vbox, ui.Label())
        self.label.set_line_wrap(True)
        self.set_text(text)
        close = ui.Button(_('Close'))
        pack(vbox, [close])
        close.connect('clicked', self.destroy)
        self.connect("delete-event", self.destroy)
        self.show_all()

    def set_text(self, text):
        # Count line length without all the markup tags
        plaintext = ''.join(re.split('<[^>]+>', text))
        linelen = max([len(x) for x in plaintext.split('\n')])
        text = text.replace('<c>', '<span foreground="blue">')
        text = text.replace('</c>', '</span>')
        self.label.set_width_chars(linelen)
        self.label.set_line_wrap(False)
        self.label.set_markup(text)

    def destroy(self, *args):
        self.hide()
        return True  # Prevents destruction of the window.
    
        
def help(text):
    button = ui.Button(_('Help'))
    button.connect('clicked', lambda widget, text=text: Help(text))
    return button


class Window:
    def __init__(self, gui):
        self.gui = gui
        ui.Window.__init__(self)
        self.set_title(_('Constraints'))
        vbox = ui.VBox()
        b = pack(vbox, [ui.Button(_('Constrain')),
                        ui.Label(_(' selected atoms'))])[0]
        b.connect('clicked', self.selected)
        b = pack(vbox, [ui.Button(_('Constrain')),
                        ui.Label(_(' immobile atoms:'))])[0]
        b.connect('clicked', self.immobile)
        b = pack(vbox, ui.Button(_('Clear constraint')))
        b.connect('clicked', self.clear)
        close = pack(vbox, ui.Button(_('Close')))
        close.connect('clicked', lambda widget: self.destroy())
        self.add(vbox)
        vbox.show()
        self.show()

        
def pack(vbox, widgets, end=False, bottom=False, expand=False, padding=0):
    if not isinstance(widgets, list):
        widgets.show()
        if bottom:
            vbox.pack_end(widgets, expand, expand, padding)
        else:
            vbox.pack_start(widgets, expand, expand, padding)
        return widgets
    hbox = ui.HBox(0, 0)
    hbox.show()
    if bottom:
        vbox.pack_end(hbox, expand, expand, padding)
    else:
        vbox.pack_start(hbox, expand, expand, padding)
    for widget in widgets:
        if type(widget) is ui.Entry:  # isinstance does not work here
            widget.set_size_request(widget.get_max_length() * 9, 24)
        widget.show()
        if end and widget is widgets[-1]:
            hbox.pack_end(widget, expand, expand, padding)
        else:
            hbox.pack_start(widget, expand, expand, padding)
    return widgets

    
class cancel_apply_ok:
    "Widget with Cancel, Apply and OK buttons.  The arguments are callbacks."
    def __init__(self, cancel, apply, ok):
        ui.HButtonBox.__init__(self)
        cancel_but = ui.Button('Cancel')
        cancel_but.connect('clicked', cancel)
        apply_but = ui.Button('Apply')
        apply_but.connect('clicked', apply)
        ok_but = ui.Button('OK')
        ok_but.connect('clicked', ok)
        for w in (cancel_but, apply_but, ok_but):
            self.pack_start(w, 0, 0)
            w.show()
        # self.show_all()
       
        
def oops(message, message2=None):
    dialog = ui.MessageDialog(flags=ui.DIALOG_MODAL,
                               type=ui.MESSAGE_WARNING,
                               buttons=ui.BUTTONS_CLOSE,
                               message_format=message)
    try:
        dialog.format_secondary_text(message2)
    except AttributeError:
        print(message, file=sys.stderr)
        print(message2, file=sys.stderr)
    dialog.connect('response', lambda x, y, dialog=dialog: dialog.destroy())
    dialog.show()

    
class AseGuiCancelException(Exception):
    pass
