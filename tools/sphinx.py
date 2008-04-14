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
    os.mkdir('test')
    os.chdir('test')
    results = test(verbosity=2)
    if len(results.failures) > 0 or len(results.errors) > 0:
        raise RuntimeError('Testsuite failed!')

    # Generate tar-file:
    os.chdir('..')
    assert os.system('python setup.py sdist') == 0

    if os.system('epydoc --docformat restructuredtext --parse-only ' +
                 '--name ASE ' +
                 '--url http://web2.fysik.dtu.dk/ase ' +
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

    assert os.system('mv ../../html epydoc;' +
                     'mv ../../dist/python-ase-3.0.0.tar.gz .') == 0
    
build()
