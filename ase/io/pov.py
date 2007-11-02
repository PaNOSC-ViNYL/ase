from ase.atoms import Atoms
from ase.data import cpk_colors, covalent_radii, symbols

def write_py(fileobj, atoms):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')

    fileobj.write("""\
//EXAMPLE OF SPHERE

//Files with predefined colors and textures
#include "colors.inc"
#include "glass.inc"
#include "golds.inc"
#include "metals.inc"
#include "stones.inc"
#include "woods.inc"

//Place the camera
camera {
  sky <0,0,1>           //Don't change this
  direction <-1,0,0>    //Don't change this  
  right <-4/3,0,0>      //Don't change this
  location <30,10,1.5> //Camera location
  look_at <0,0,0>     //Where camera is pointing
  angle 15      //Angle of the view--increase to see more, decrease to see less
}

//Ambient light to "brighten up" darker pictures
global_settings { ambient_light White }

//Place a light--you can have more than one!
light_source {
  <10,-10,20>   //Change this if you want to put the light at a different point
  color White*2         //Multiplying by 2 doubles the brightness
}

//Set a background color
background { color White }

//Create a "floor"
//plane {
//  <0,0,1>, 0            //This represents the plane 0x+0y+z=0
//  texture { T_Silver_3A }       //The texture comes from the file "metals.inc"
//}

//Sphere with specified center point and radius
//The texture comes from the file "stones.inc""""
                  )
    numbers = atoms.get_atomic_numbers()
    done = {}
    for Z in numbers:
        if Z not in done:
            symbol = symbols[Z]
            fileobj.write('#define C_%s ' % symbol +
                          'color rgb <%.2f, %.2f, %.2f>\n' %
                          tuple(cpk_colors[Z]))
            fileobj.write('#define R_%s %.2f\n' % (symbol, covalent_radii[Z]))
            done[Z] = True

    for Z, p in zip(numbers, atoms.get_positions()):
        symbol = symbols[Z]
        fileobj.write('sphere{<%.2f, %.2f, %.2f>, ' % tuple(p) +
                      'R_%s C_%s}\n' % (symbol, symbol))
