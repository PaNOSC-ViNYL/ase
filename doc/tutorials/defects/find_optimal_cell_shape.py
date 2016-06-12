#!/usr/bin/python -u
# PE, 2016/06/11

import numpy as np
from ase.io import read

# Handle the command line
import argparse
parser = argparse.ArgumentParser(description='Find most compact supercell with a given number of atoms.')
parser.add_argument('target_size', type=int, help='Size in multiples of primitive cell.')
parser.add_argument('primfile', help='File with primitive unit cell in ASE readable format.')
defval = 'sc'
parser.add_argument('target_shape', choices=['fcc','sc'], help='type of supercell to create, choices: fcc, sc [default: %s]' % defval)
defval = -6
parser.add_argument('-l', '--lower', dest='lower_limit', default=defval, type=int,
                    help='Lower limit of range of elements of the transformation matrix [default: %d]' % defval)
defval = 6
parser.add_argument('-u', '--upper', dest='upper_limit', default=defval, type=int,
                    help='Upper limit of range of elements of the transformation matrix [default: %d]' % defval)
args = parser.parse_args()

print "target size in number of unit cells: %d" % args.target_size
print "primitive cell read from file: %s" % args.primfile

#---------------------------

# Read primitive unit cell.
conf_prim = read(args.primfile)
cell_prim = conf_prim.get_cell()
print "primitive cell metric:"
print cell_prim

# Inline C code that does the heavy work.
# It iterates over all possible matrices and finds the best match.
code = """
#include <math.h>

// The code below will search over a large number of matrices
// that are generated from a set of integer numbers between
// imin and imax.
int imin = range[0];
int imax = range[1];
long int ntotal = pow(imax - imin + 1, 9);
printf("searching over %ld matrices with coefficients in the range from %d to %d\\n", ntotal, imin, imax);

// For the sake of readability the primitive cell is copied
// to a set of double variables.
double h11 = cell_prim[0];   // [3 * 0 + 0];
double h12 = cell_prim[1];   // [3 * 0 + 1];
double h13 = cell_prim[2];   // [3 * 0 + 2];
double h21 = cell_prim[3];   // [3 * 1 + 0];
double h22 = cell_prim[4];   // [3 * 1 + 1];
double h23 = cell_prim[5];   // [3 * 1 + 2];
double h31 = cell_prim[6];   // [3 * 2 + 0];
double h32 = cell_prim[7];   // [3 * 2 + 1];
double h33 = cell_prim[8];   // [3 * 2 + 2];

double det_P;  // will store determinant of P matrix
double m;      // auxiliary variable
double current_score; // l2-norm of (Ph - h_cub)

for (int xx=imin ; xx<=imax ; xx++) 
for (int xy=imin ; xy<=imax ; xy++) 
for (int xz=imin ; xz<=imax ; xz++) 
for (int yx=imin ; yx<=imax ; yx++) 
for (int yy=imin ; yy<=imax ; yy++) 
for (int yz=imin ; yz<=imax ; yz++) 
for (int zx=imin ; zx<=imax ; zx++) 
for (int zy=imin ; zy<=imax ; zy++) 
for (int zz=imin ; zz<=imax ; zz++) {

  // compute determinant
  det_P = 0;
  det_P += xx*yy*zz;
  det_P += xy*yz*zx;
  det_P += xz*yx*zy;
  det_P -= xx*yz*zy;
  det_P -= xy*yx*zz;
  det_P -= xz*yy*zx;

  if (det_P == target_size[0]) {

    // compute second moment
    current_score = 0.0;
    current_score += pow(xx * h11 + xy * h21 + xz * h31 - target_h[0], 2); // 1-1
    current_score += pow(xx * h12 + xy * h22 + xz * h32 - target_h[1], 2);
    current_score += pow(xx * h13 + xy * h23 + xz * h33 - target_h[2], 2);
    current_score += pow(yx * h11 + yy * h21 + yz * h31 - target_h[3], 2);
    current_score += pow(yx * h12 + yy * h22 + yz * h32 - target_h[4], 2); // 2-2
    current_score += pow(yx * h13 + yy * h23 + yz * h33 - target_h[5], 2);
    current_score += pow(zx * h11 + zy * h21 + zz * h31 - target_h[6], 2);
    current_score += pow(zx * h12 + zy * h22 + zz * h32 - target_h[7], 2);
    current_score += pow(zx * h13 + zy * h23 + zz * h33 - target_h[8], 2); // 3-3
    
    if (current_score < best_score[0]) {
      best_score[0] = current_score;
      optimal_P[0] = xx;   // [3 * 0 + 0];
      optimal_P[1] = xy;   // [3 * 0 + 1];
      optimal_P[2] = xz;   // [3 * 0 + 2];
      optimal_P[3] = yx;   // [3 * 1 + 0];
      optimal_P[4] = yy;   // [3 * 1 + 1];
      optimal_P[5] = yz;   // [3 * 1 + 2];
      optimal_P[6] = zx;   // [3 * 2 + 0];
      optimal_P[7] = zy;   // [3 * 2 + 1];
      optimal_P[8] = zz;   // [3 * 2 + 2];
    }
  }
}
"""

# Set up target metric
if args.target_shape in ['sc', 'simple-cubic']:
    p = np.eye(3)
elif args.target_shape in ['fcc']:
    p = 0.5 * np.array([[0,1,1],
                        [1,0,1],
                        [1,1,0]], dtype=float)

L = (1.0 / np.linalg.det(p) * args.target_size * conf_prim.get_volume()) ** (1./3)
target_h = L * p
print 'effective target dimension: %g' % L
print 'target metric:'
print target_h

# Initialize run ...
best_score = np.array([100000.0])
optimal_P = np.zeros((3,3))
range = np.array([args.lower_limit,args.upper_limit], dtype=int)
target_size = np.array(args.target_size, dtype=int)
# ... and execute the inline C code.
from scipy import weave
weave.inline(code, ['range', 'cell_prim', 'best_score', 'optimal_P', 'target_h', 'target_size'])

# Finalize.
print 'smallest |P h_p - h_target|_2 value: %f' % best_score
print 'optimal P:'
print optimal_P
print 'cell:'
print np.round(np.dot(optimal_P, cell_prim), 4)
print 'determinant of optimal P: ', np.linalg.det(optimal_P)
