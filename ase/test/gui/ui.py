def test():
    import ase.gui.ui as ui

    #ui.select_backend('test')

    def hello(event=None):
        print('hello')

    menu = [('Hi', [('Hello', 'Ctrl+H', '', hello)])]
    win = ui.MainWindow('Test', menu=menu)
    win.add(ui.Label('Hello'))
    win.add(ui.Button('Hello', hello))

    def abc(x):
        print(x)

    c = ui.ComboBox(['Aa', 'Bb', 'Cc'], abc)
    win.add(c)

    r = ui.RadioButtons(['A', 'B', 'C'], 'ABC', abc)
    win.add(r)

    b = ui.CheckButton('Hello')

    def hi():
        print(b.value, r.value, c.value)

    win.add([b, ui.Button('Hi', hi)])

    win.run()


test()
