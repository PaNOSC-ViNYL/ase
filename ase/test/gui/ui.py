import ase.gui.ui as ui


def window():

    def hello(event=None):
        print('hello', event)

    menu = [('Hi', [ui.MenuItem('_Hello', hello, 'Ctrl+H')]),
            ('Hell_o', [ui.MenuItem('ABC', hello, choices='ABC')])]
    win = ui.MainWindow('Test', menu=menu)

    win.add(ui.Label('Hello'))
    win.add(ui.Button('Hello', hello))

    r = ui.Rows([ui.Label(x * 7) for x in 'abcd'])
    win.add(r)
    r.add('11111\n2222\n333\n44\n5')

    def abc(x):
        print(x, r.rows)

    cb = ui.ComboBox(['Aa', 'Bb', 'Cc'], callback=abc)
    win.add(cb)

    rb = ui.RadioButtons(['A', 'B', 'C'], 'ABC', abc)
    win.add(rb)

    b = ui.CheckButton('Hello')

    def hi():
        print(b.value, rb.value, cb.value)
        del r[2]
        r.add('-------------')

    win.add([b, ui.Button('Hi', hi)])

    return win


def run():
    win = window()
    t = test(win)
    win.test(t)


def test(w):
    w.things[1].callback()
    yield
    w.things[1].callback()


run()
