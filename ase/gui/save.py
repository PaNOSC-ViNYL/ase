"""Dialog for saving one or more configurations."""
import os

import numpy as np
import gtk
from gettext import gettext as _

from ase.io.formats import (write, parse_filename, get_ioformat, string2index,
                            filetype)


def save_dialog(gui):
    dialog = gtk.FileChooserDialog(_('Save ...'),
                                   None,
                                   gtk.FILE_CHOOSER_ACTION_SAVE,
                                   (gtk.STOCK_SAVE, gtk.RESPONSE_OK,
                                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))
    dialog.set_current_name('')
    dialog.set_current_folder(os.getcwd())
    text = _('Append name with "@n" in order to write image number "n" '
             'instead of the current image.\n'
             'Append "@start:stop" or "@start:stop:step" if you want to write '
             'a range of images.\n'
             'You can leave out "start" and "stop" so that "name@:" will '
             'give you all images.\n'
             'Negative numbers count from the last image '
             '("name@-1": last image, "name@-2:": last two).')
    dialog.set_extra_widget(gtk.Label(text))
    response = dialog.run()
    if response == gtk.RESPONSE_OK:
        filename = dialog.get_filename()
        dialog.destroy()
    else:
        dialog.destroy()
        return
        
    filename, index = parse_filename(filename)
    if index is None:
        index = slice(gui.frame, gui.frame + 1)
    if isinstance(index, str):
        index = string2index(index)
    format = filetype(filename, read=False)
    io = get_ioformat(format)
    
    extra = {}
    remove_hidden = False
    if format in ['png', 'eps', 'pov']:
        bbox = np.empty(4)
        size = np.array([gui.width, gui.height]) / gui.scale
        bbox[0:2] = np.dot(gui.center, gui.axes[:, :2]) - size / 2
        bbox[2:] = bbox[:2] + size
        extra['rotation'] = gui.axes
        extra['show_unit_cell'] = gui.ui.get_widget(
            '/MenuBar/ViewMenu/ShowUnitCell').get_active()
        extra['bbox'] = bbox
        extra['colors'] = gui.get_colors(rgb=True)[gui.images.visible]
        remove_hidden = True

    images = [gui.images.get_atoms(i, remove_hidden=remove_hidden)
              for i in range(*index.indices(gui.images.nimages))]
    
    if len(images) > 1 and io.single:
        # We want to write multiple images, but the file format does not
        # support it.  The solution is to write multiple files, inserting
        # a number in the file name before the suffix.
        j = filename.rfind('.')
        filename = filename[:j] + '{0:05d}' + filename[j:]
        for i, atoms in enumerate(images):
            write(filename.format(i), atoms, **extra)
    else:
        write(filename, images, **extra)
