if [ -f build-ase-webpage.lock ]
then
    >&2 echo "Locked"
    exit
fi
touch build-ase-webpage.lock

# Update ASE:
cd ase
svn up
cd ..

# Name of tar-file containing the whole web-page:
html=html.tar.gz

changes=`find ase -newer $html | grep -v .svn`
echo Changes:
echo $changes

if [ -n "$changes" ]
then
    # Clean up:
    rm -rf ase/doc/
    
    # Create development snapshot tar-file:
    cd ase
    svn up
    rm -r dist
    python setup.py sdist > sdist.out
    
    ase=`pwd`
    cd doc
    export PYTHONPATH=$ase:$PYTHONPATH
    export PATH=$PATH:$ase/tools
    
    # Clean up:
    python -m ase.utils.sphinx clean
    
    make html
           
    # Use https for mathjax:
    cd build
    find -name "*.html" | xargs \
        sed -i "s%http://cdn.mathjax.org%https://cdn.mathjax.org%"
        
    # Set correct version of snapshot tar-file:
    version=`python -c "from ase.version import version; print(version)"`
    sed -i s/snapshot.tar/$version.tar/g html/download.html
    
    mv ../../dist/python-ase-*.tar.gz html
    
    cd ../..
    epydoc --docformat restructuredtext --parse-only \
           --name ASE --url http://wiki.fysik.dtu.dk/ase \
           --show-imports --no-frames -v ase > epydoc.out
    # Check for warnings:
    >&2 grep '^|' epydoc.out

    mv html doc/build/html/epydoc
    
    cd doc/build
    tar -czf $html html
    cp $html ../../..
    cd ../../..
fi

rm build-ase-webpage.lock
