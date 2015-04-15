from __future__ import print_function
import glob
import optparse
import os
import shutil
import subprocess
import sys


def svn_update(name='ase'):
    os.chdir(name)
    output = subprocess.check_output('svn update', shell=True)
    os.chdir('..')
    lastline = output.splitlines()[-1]
    return not lastline.startswith('At revision')

        
def build(force_build, name='ase', env=''):
    if not force_build:
        return
        
    home = os.getcwd()
    
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
    subprocess.check_call(env + ' PYTHONPATH='
                          '{0}/lib/python:{0}/lib64/python:$PYTHONPATH '
                          'PATH={0}/bin:$PATH '.format(home) +
                          'make html', shell=True)
           
    # Use https for mathjax:
    subprocess.check_call(
        'find build -name "*.html" | '
        'xargs sed -i "s|http://cdn.mathjax.org|https://cdn.mathjax.org|"',
        shell=True)
        
    tar = glob.glob('../dist/*.tar.gz')[0].split('/')[-1]
    os.rename('../dist/' + tar, 'build/html/' + tar)
    
    # Set correct version of snapshot tar-file:
    subprocess.check_call(
        'find build/html -name download.html | '
        'xargs sed -i s/snapshot.tar.gz/{}/g'.format(tar),
        shell=True)
    
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
    dir = name + '-web-page'
    os.rename('html', dir)
    subprocess.check_call('tar -czf {0}.tar.gz {0}'.format(dir),
                          shell=True)
    os.rename('{}.tar.gz'.format(dir), '../../../{}.tar.gz'.format(dir))
    os.chdir('../../..')
    try:
        shutil.rmtree('lib64')
    except OSError:
        pass
    shutil.rmtree('lib')
    shutil.rmtree('bin')
    shutil.rmtree('share')


def main(build=build):
    """Build web-page if there are changes in the source.
    
    The optional build function is used by GPAW to build its web-page.
    """
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
        changes = svn_update('ase')
        build(opts.force_build or changes)
    finally:
        os.remove('build-web-page.lock')

        
if __name__ == '__main__':
    main()
