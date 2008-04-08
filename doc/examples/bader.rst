Bader Analysis
--------------

Henkelman et al. have implemented a fast and robust algorithm for calculating the electronic charges on individual atoms in molecules or crystals with the Bader scheme (G. Henkelman, A. Arnaldsson, and H. Jónsson, `A fast and robust algorithm for Bader decomposition of charge density`_, Comput. Mater. Sci. (in press, 2005)). In that method the electron density is analyzed and a so-called zero-flux surfaces are used to divide a system into atoms. A more detailed description about the method can be found at `Richard Bader's`_ homepage at McMaster University and more details about the algorithm (how to use, about the output files, discussion forum and an article about the algorithm) can be found at `Graeme Henkelman's`_ homepage at the University of Texas at Austin. A presentation about the Bader method and this algorithm can be found `here`_.

.. _A fast and robust algorithm for Bader decomposition of charge density: http://theory.cm.utexas.edu/bader/bader.pdf
.. _Richard Bader's: http://www.chemistry.mcmaster.ca/aim/aim_0.html 
.. _Graeme Henkelman's: http://theory.cm.utexas.edu/bader/
.. _here: http://www.hi.is/~egillsk/stuff/annad/egillsk_bader_pres_150206.pdf

This algorithm suites very well for large solid state physical systems as well as large biomolecular systems. The computational time depends only on the size of the 3D grid used in the calculations and typically it takes less than a minute to do the analysis. To get more accurate results in the charge analysis it is recommended to use a finer grid than the default values, however, usually it is enough to calculate the electronic structure with a finer grid after the optimization has been done with the default values. 

To be able to use this algorithm one has to create cube files from the NetCDFiles of e.g. a dacapo calculation. These cube files contain the size of the unitcell, coordinates of the atoms, and the electronic charge density on a grid. This can be done with a simple python script:

>>> from Dacapo import Dacapo
>>> from ASE.IO.Cube import *
>>> atoms = Dacapo.read_atoms('filename.nc')
>>> calc = atoms.get_calculator()
>>> dens = calc.get_density_array()
>>> density =  dens * (0.529177)**3     # changes from A^-3 units to bohr^-3 units
>>> WriteCube(atoms, density, 'filename.cube')

After the cube files have been created, one can run the bader code on the cube files with four possible options (1, 2, 3 or 4) where the choose of options is a choose of output files from the analysis as explained on the code's homepage. The bader code has been installed on all camp's desktop and will also be installed on niflheim in near future. Examples of using the program in the command line are::

    bader filename.cube 4

where no Bader volumes are written at all. It is recommended to use this option first to see if the analysis is correct because writing the Bader volumes take a longer time. If one uses e.g. ::

    bader filename.cube 3

then all the Bader volumes around each Bader maxima are written in the same cube file. If option 1 is used then all Bader volumes are written in separate files for each Bader region.

Few things should be mentioned for the people that wants to use this analysis tool; 

- The unitcell used in the calculations of the electronic structure has to be rectangular.

- You can only create the cube files when your calculations are finished, because it is not until then that the charge density is written on the grid. 

- There can be a problem with O-H bonds (and possibly C-H bonds too) when the electron density is coming from pseudopotential code like dacapo. 


For questions about the Bader algorithm in general should be directed to the code's `discussion forum`_ but questions of how to use the program with cube files written from NetCDFiles in e.g. Dacapo should be directed to egillsk@fysik.dtu.dk

.. _discussion forum: http://theory.cm.utexas.edu/forum/ 
