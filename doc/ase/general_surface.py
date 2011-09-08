# creates: general_surface.pdf
import os

for i in range(2):
    error = os.system('pdflatex -interaction=nonstopmode general_surface ' +
                      '> /dev/null')
    assert error == 0, 'pdflatex failed'
