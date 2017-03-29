"""
cd
python3 -m venv web-page
cd web-page
. bin/activate
pip install sphinx-rtd-theme
pip install Sphinx
pip install matplotlib scipy
git clone git@gitlab.com:ase/ase
cd ase
pip install -e .

cd ~/web-page
. bin/activate
cd ase
"""

from __future__ import print_function
import os
import subprocess
import sys
import time


def git_pull(name='ase'):
    os.chdir(name)
    try:
        sleep = 1  # 1 minute
        while True:
            try:
                output = subprocess.check_output(
                    'GIT_HTTP_LOW_SPEED_LIMIT=1000 '
                    'GIT_HTTP_LOW_SPEED_TIME=20 '  # make sure we get a timeout
                    'git pull 2>> pull.err', shell=True)
            except subprocess.CalledProcessError:
                if sleep > 16:
                    raise
                time.sleep(sleep * 60)
                sleep *= 2
            else:
                break
    finally:
        os.chdir('..')
    lastline = output.splitlines()[-1]
    return not lastline.startswith('Already up-to-date')


cmds = """\
touch build-web-page.lock
git clean -fdx
git checkout web-page
git pull
cd doc; sphinx-build -b html -d build/doctrees . build/html
mv build/html web-page
git clean -fdx doc
git checkout master
git pull
cd doc; sphinx-build -b html -d build/doctrees . build/html
mv build/html web-page/dev
python setup.py sdist
cp dist/ase-*.tar.gz web-page/
cp dist/ase-*.tar.gz web-page/dev/
tarfile=`(cd dist; ls ase-*.tar.gz)`
find html -name install.html | xargs sed -i s/snapshot.tar.gz/$tarfile/g
tar -cf web-page.tar.gz web-page"""


def build():
    if os.path.isfile('build-web-page.lock'):
        print('Locked', file=sys.stderr)
        return
    try:
        for cmd in cmds.splitlines():
            subprocess.check_call(cmd, shell=True)
    finally:
        os.remove('build-web-page.lock')


if __name__ == '__main__':
    build()
