#!/usr/bin/python
import os
import sys
import time
import glob
import trace
import tempfile

def build():
    if os.system('svn export ' +
                 'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
        raise RuntimeError('Checkout of ASE failed!')
    os.chdir('ase')
    if os.system('python setup.py install --home=.') != 0:
        raise RuntimeError('Installation failed!')
    sys.path.insert(0, 'lib/python')
    from ase.test import test

    # Run test-suite:
    results = test(verbosity=2, dir='ase/test')
    if len(results.failures) > 0 or len(results.errors) > 0:
        raise RuntimeError('Testsuite failed!')

    # Generate tar-file:
    assert os.system('python setup.py sdist') == 0

    if os.system('epydoc --docformat restructuredtext --parse-only ' +
                 '--name ASE ' +
                 '--url http://wiki.fysik.dtu.dk/ase ' +
                 '--show-imports --no-frames -v ase >& epydoc.out') != 0:
        raise RuntimeError('Epydoc failed!')

    if ' Warning:' in open('epydoc.out').read():
        raise RuntimeError('Warning(s) from epydoc!')

    os.chdir('doc')
    os.mkdir('_build')
    if os.system('sphinx-build . _build') != 0:
        raise RuntimeError('Sphinx failed!')
    os.system('cd _build; cp _static/searchtools.js .')

    if 1:
        if os.system('sphinx-build -b latex . _build') != 0:
            raise RuntimeError('Sphinx failed!')
        os.chdir('_build')
        #os.system('cd ../..; ln -s doc/_static')
        if os.system('make ase-manual.pdf') != 0:
            raise RuntimeError('pdflatex failed!')
    else:
        os.chdir('_build')

    assert os.system('mv ../../html epydoc;' +
                     'mv ../../dist/python-ase-3.0.0.tar.gz .') == 0
    
build()
