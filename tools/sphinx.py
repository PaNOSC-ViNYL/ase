#!/usr/bin/python
import os
import sys
import time
import glob
import trace
import tempfile

def build():
    if os.system('svn export ' +
                 'https://svn.fysik.dtu.dk/projects/ase3000/trunk ase3') != 0:
        raise RuntimeError('Checkout of ASE failed!')
    os.chdir('ase3')
    if os.system('python setup.py install --home=.') != 0:
        raise RuntimeError('Installation failed!')

    # Run test-suite:
    os.mkdir('test')
    os.chdir('test')
    from ase.test import test
    results = test()
    if len(results.failures) > 0 or len(results.errors) > 0:
        raise RuntimeError('Testsuite failed!')

    # Generate tar-file:
    os.chdir('..')
    assert os.system('python setup.py sdist') == 0

    if os.system('epydoc --docformat restructuredtext --parse-only ' +
                 '--name ASE ' +
                 '--url http://wiki.fysik.dtu.dk/stuff/ase/sphinx ' +
                 '--show-imports --no-frames -v ase >& epydoc.out') != 0:
        raise RuntimeError('Epydoc failed!')

    if ' Warning:' in open('epydoc.out').read():
        raise RuntimeError('Warning(s) from epydoc!')

    os.chdir('doc')
    os.mkdir('.static')
    os.mkdir('.build')
    if os.system('sphinx-build . .build') != 0:
        raise RuntimeError('Sphinx failed!')

    if os.system('sphinx-build -b latex . .build') != 0:
        raise RuntimeError('Sphinx failed!')
    os.chdir('.build')
    if os.system('make ase.pdf') != 0:
        raise RuntimeError('pdflatex failed!')
    
build()
