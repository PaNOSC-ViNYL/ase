from __future__ import print_function
import glob
import optparse
import os
import shutil
import subprocess
import sys
sys.path.insert(0, 'ase')

from ase.version import version


def svn_update(name='ase'):
    os.chdir(name)
    output = subprocess.check_output('svn update', shell=True)
    os.chdir('..')
    lastline = output.splitlines()[-1]
    return not lastline.startswith('At revision')

        
def build(force_build, name='ase', env='', version=version):
    changes = svn_update()
    if not (force_build or changes):
        return
        
    os.chdir(name)

    # Clean up:
    shutil.rmtree('doc')
    subprocess.check_call('svn update', shell=True)

    # Create development snapshot tar-file and install:
    try:
        shutil.rmtree('dist')
    except OSError:
        pass
    subprocess.check_call('python setup.py sdist install --home=..',
                          shell=True)

    # Build web-page:
    os.chdir('doc')
    subprocess.check_call(env + ' ' +
                          'PYTHONPATH=../../lib/python:$PYTHONPATH '
                          'PATH=../../bin:$PATH '
                          'make html', shell=True)
           
    # Use https for mathjax:
    subprocess.check_call(
        'find build -name "*.html" | '
        'xargs sed -i "s|http://cdn.mathjax.org|https://cdn.mathjax.org|"',
        shell=True)
        
    # Set correct version of snapshot tar-file:
    subprocess.check_call(
        'find build/html -name download.html | '
        'xargs sed -i s/snapshot.tar/{}.tar/g'.format(version),
        shell=True)
    
    tar = glob.glob('../dist/*.tar.gz')[0]
    os.rename(tar, 'build/html/' + tar.split('/')[-1])
    
    os.chdir('..')
    output = subprocess.check_output(
        'epydoc --docformat restructuredtext --parse-only '
        '--name {0} --url http://wiki.fysik.dtu.dk/{1} '
        '--show-imports --no-frames -v {1}'.format(name.upper(), name),
        shell=True)
    
    # Check for warnings:
    for line in output.splitlines():
        if line.startswith('|'):
            print(line, file=sys.stderr)

    os.rename('html', 'doc/build/html/epydoc')
    
    os.chdir('doc/build')
    os.rename('html', name)
    subprocess.check_call('tar -czf {0}.tar.gz {0}'.format(name),
                          shell=True)
    os.rename('{}.tar.gz'.format(name), '../../../{}.tar.gz'.format(name))
    os.chdir('../../..')
    shutil.rmtree('lib')
    shutil.rmtree('bin')
    shutil.rmtree('share')


def main(build=build):
    if os.path.isfile('build-web-page.lock'):
        print('Locked', file=sys.stderr)
        return
    try:
        open('build-web-page.lock', 'w').close()
            
        parser = optparse.OptionParser(usage='Usage: %prog [-f]',
                                       description='Build web-page')
        parser.add_option('-f', '--force-build', action='store_true',
                          help='Force build instead of building only when '
                          'there are changes to the docs or code.')
        opts, args = parser.parse_args()
        assert len(args) == 0
        build(opts.force_build)
    finally:
        os.remove('build-web-page.lock')

        
if __name__ == '__main__':
    main()
