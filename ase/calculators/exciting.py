"""
========
exciting
========

.. image:: ../../_static/exciting.png

Introduction
============

``exciting`` is a full-potential *all-electron*
density-functional-theory (:term:`DFT`) package based on the
linearized augmented planewave (:term:`LAPW`) method. It can be
applied to all kinds of materials, irrespective of the atomic species
involved, and also allows for the investigation of the core
region. The website is http://exciting-code.org/

The module depends on lxml  http://codespeak.net/lxml/


Exciting Calculator Class
=========================

.. autoclass:: ase.calculators.exciting.Exciting

There are two ways to construct the exciting calculator.

Constructor Keywords
--------------------

One is by giving parameters of the ground state in the
constructor. The possible attributes can be found at
http://exciting-code.org/input-reference#groundstate.

.. class:: Exciting(bin='excitingser', kpts=(4, 4, 4), xctype='GGArevPBE')

XSLT Template
-------------
The other way is to use an XSLT template

.. class:: Exciting( template='template.xsl'  , bin='excitingser')

The template may look like:

.. highlight:: xml

::

    <?xml version="1.0" encoding="UTF-8" ?>
    <xsl:stylesheet version="1.0"
      xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
      <xsl:output method="xml" />
      <xsl:template match="/">
        <xsl:comment>
          created from template
        </xsl:comment>
    
    
    <!-- ############# -->
        <input>
          <title></title>
          <structure speciespath="./">
            <xsl:copy-of select="/input/structure/crystal" />
            <xsl:copy-of select="/input/structure/species" />
          </structure>
          <groundstate ngridk="4  4  4"  vkloff="0.5  0.5  0.5" tforce="true" />
        </input>
    <!-- ############# -->
    
    
    
      </xsl:template>
    </xsl:stylesheet>
    
.. highlight:: python    

Current Limitations
===================

The calculator supports only total energy and forces no stress strain
implemented in exciting yet.  However it is planned to be implemented
soon. http://exciting-code.org/current-developments.

The keywords on the constructor are converted into attributes 
of the ``groundstate`` element. No check is done there and not all 
options are accessible in that way. By using the template one can 
exceed this limitation.
"""
import os
import numpy as np
from ase.units import Bohr, Hartree


class Exciting:
    """exciting calculator object."""
   
    def __init__(self, dir='.', template=None, speciespath=None,
                 bin='excitingser', kpts=(1, 1, 1), **kwargs):
        """Exciting calculator object constructor
        
        Parameters
        ----------
        dir : string
            directory in which to excecute exciting
        template :string
            Path to XSLT templat if it schould be used
            default: none
        bin :string
            Path or executable name of exciting 
            default: ``excitingser`` 
        kpts:integer list length 3
            Number of kpoints
        **kwargs : dictionary like
            list of key value pairs to be converted into groundstate attributes
        
        """
        self.dir = dir
        self.energy = None
        self.template = template
        if speciespath is None:
            self.speciespath = os.environ.get('EXCITING_SPECIES_PATH', './')
        self.converged = False
        self.excitingbinary = bin
        self.groundstate_attributes = kwargs
        if not 'ngridk' in kwargs.keys():
            self.groundstate_attributes['ngridk'] = ' '.join(map(str, kpts))

    def update(self, atoms):
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.write(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        self.initialize(atoms)
        syscall = ('cd %(dir)s; %(bin)s;' %
                   {'dir': self.dir, 'bin': self.excitingbinary})
        print syscall
        assert os.system(syscall ) == 0
        self.read()

    def write(self, atoms):
        from lxml import etree as ET
        from ase.io.exciting import  atoms2etree
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        root = atoms2etree(atoms)
        root.find('structure').attrib['speciespath'] = self.speciespath
        groundstate = ET.SubElement(root, 'groundstate', tforce='true')
        for key, value in self.groundstate_attributes.items():
            groundstate.attrib[key] = str(value)
        if self.template:
            xslf = open(self.template, 'r')
            xslt_doc = ET.parse(xslf)
            transform = ET.XSLT(xslt_doc)
            result = transform(root)
            fd = open('%s/input.xml' % self.dir, 'w')
            fd.write(ET.tostring(result, method='xml', pretty_print=True,
                                 xml_declaration=True, encoding='UTF-8'))
            fd.close()
        else:
            fd = open('%s/input.xml' % self.dir, 'w')
            fd.write(ET.tostring(root, method='xml', pretty_print=True,
                                 xml_declaration=True, encoding='UTF-8'))
            fd.close()
        
    def read(self):
        """ reads Total energy and forces from info.xml
        """
        from lxml import etree as ET
        INFO_file = '%s/info.xml' % self.dir
   
        try:
            fd = open(INFO_file)
        except IOError:
            print "file doesn't exist"
        info = ET.parse(fd)
        self.energy = float(info.xpath('//@totalEnergy')[-1]) * Hartree
        forces = []
        forcesnodes = info.xpath(
            '//structure[last()]/species/atom/forces/totalforce/@*')
        for force in forcesnodes:
            forces.append(np.array(float(force)))
        self.forces = np.reshape(forces, (-1, 3))
        
        if str(info.xpath('//groundstate/@status')[0]) == 'finished':
            self.converged = True
        else:
            raise RuntimeError('calculation did not finish correctly')
  
        # Stress
        self.stress = np.empty((3, 3))
