# encoding: utf-8

"Module for displaying information about the system."

import ase.gui.ui as ui
from ase.gui.widgets import pack

singleimage = _("Single image loaded.")
multiimage = _("Image %d loaded (0 - %d).")
ucconst = _("Unit cell is fixed.")
ucvaries = _("Unit cell varies.")

format = _("""\
%s

Number of atoms: %d.

Unit cell:
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f
  %8.3f  %8.3f  %8.3f

%s
%s
""")

class QuickInfo:
    def __init__(self, gui):
        ui.Window.__init__(self)
        self.set_title(_("Quick Info"))
        vbox = ui.VBox()
        images = gui.images
        if images.natoms < 1:
            txt = _("No atoms loaded.")
        else:
            (nimg, natoms, three) = images.P.shape
            assert three == 3
            img = gui.frame
            uc = images.A[img]
            if nimg > 1:
                equal = True
                for i in range(nimg):
                    equal = equal and (uc == images.A[i]).all()
                if equal:
                    uctxt = ucconst
                else:
                    uctxt = ucvaries
            else:
                uctxt = ""
            if nimg == 1:
                imgtxt = singleimage
            else:
                imgtxt = multiimage % (img, nimg - 1)
            
            periodic = [[_('no'), _('yes')][periodic]
                        for periodic in images.pbc]
            
            # TRANSLATORS: This has the form Periodic: no, no, yes
            pbcstring = _('Periodic: %s, %s, %s') % tuple(periodic)
            txt = format % ((imgtxt, natoms) + tuple(uc.flat) +
                            (pbcstring,) + (uctxt,))
        label = ui.Label(txt)
        pack(vbox, [label])
        but = ui.Button('Close')
        but.connect('clicked', self.close)
        pack(vbox, [but], end=True)
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

    def close(self, *args):
        self.destroy()

    
