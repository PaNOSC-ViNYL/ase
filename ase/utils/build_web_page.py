"""Build ASE's web-page.

Initial setup::

    cd ~
    python3 -m venv web-page
    cd web-page
    . bin/activate
    pip install sphinx-rtd-theme
    pip install Sphinx
    pip install matplotlib scipy
    git clone git@gitlab.com:ase/ase
    cd ase
    pip install -e .

Crontab::

    build="python -m ase.utils.build_web_page"
    10 * * * * cd ~/web-page; . bin/activate; cd ase; $build > ../ase.log

"""

import os
import subprocess
import sys

from ase import __version__


cmds = """\
touch ../ase-web-page.lock
git clean -fdx
git checkout web-page
git pull
cd doc; sphinx-build -b html -d build/doctrees . build/html
mv doc/build/html web-page
git clean -fdx doc
git checkout master
git pull
cd doc; sphinx-build -b html -d build/doctrees . build/html
mv doc/build/html web-page/dev
python setup.py sdist
cp dist/ase-*.tar.gz web-page/
cp dist/ase-*.tar.gz web-page/dev/
find web-page -name install.html | xargs sed -i s/snapshot.tar.gz/{}/g
tar -cf web-page.tar.gz web-page""".format('ase-' + __version__ + '.tar.gz')


def build():
    if os.path.isfile('../ase-web-page.lock'):
        print('Locked', file=sys.stderr)
        return
    try:
        for cmd in cmds.splitlines():
            subprocess.check_call(cmd, shell=True)
    finally:
        os.remove('../ase-web-page.lock')


if __name__ == '__main__':
    build()
