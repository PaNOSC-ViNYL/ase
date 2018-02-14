#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL ASE_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#define DOUBLEP(a) ((double*)PyArray_DATA(a))
#define INTP(a) ((int*)PyArray_DATA(a))

#define VELOCITY_VERLET_ADJUST_POSITION_TOL 1e-13
#define VELOCITY_VERLET_ADJUST_VELOCITY_TOL 1e-13
#define COULOMB_CONSTANT 14.399645351950548
// #define FF_DEBUG_M 1

#ifdef __clang__
#define INLINE static
#else
#define INLINE inline
#endif

// Index notation
// v is coordinate (x,y or z)
// x is distance index (1-2, 2-3, or 3-1)
// n is water index
// i is atom index (O, H1 or H2)

INLINE void vec3_sub(double* target_v, double* a_v, double* b_v)
{
  target_v[0] = a_v[0] - b_v[0];
  target_v[1] = a_v[1] - b_v[1];
  target_v[2] = a_v[2] - b_v[2];
}

INLINE void vec3_zero(double* target_v)
{
  target_v[0] = 0.0;
  target_v[1] = 0.0;
  target_v[2] = 0.0;
}

INLINE void vec3_submin(double* target_v, double* a_v, double* b_v, unsigned char* pbc, double* celldiag)
{
  for (unsigned int v=0; v<3; v++)
  if (pbc[v]) {
    double s = a_v[v] - b_v[v];
    s -= round(s / celldiag[v])*celldiag[v];
    target_v[v] = s;
  } else {
    target_v[v] = a_v[v] - b_v[v];
  }
}

INLINE void vec9_diffs(double* target_xv, double* R_iv)
{
   // R1-R2
   target_xv[0] = R_iv[0*3 + 0] - R_iv[1*3 + 0]; // dOH1_x
   target_xv[1] = R_iv[0*3 + 1] - R_iv[1*3 + 1]; // dOH1_y
   target_xv[2] = R_iv[0*3 + 2] - R_iv[1*3 + 2]; // dOH1_z
   // R2-R3
   target_xv[3] = R_iv[1*3 + 0] - R_iv[2*3 + 0]; // dOH2_x
   target_xv[4] = R_iv[1*3 + 1] - R_iv[2*3 + 1]; // dOH2_y
   target_xv[5] = R_iv[1*3 + 2] - R_iv[2*3 + 2]; // dOH2_z
   // R3-R1
   target_xv[6] = R_iv[2*3 + 0] - R_iv[0*3 + 0]; // dH1H2_x
   target_xv[7] = R_iv[2*3 + 1] - R_iv[0*3 + 1]; // dH1H2_y
   target_xv[8] = R_iv[2*3 + 2] - R_iv[0*3 + 2]; // dH1H2_z
}

INLINE void vec9_massdiffs(double* target_xv, double* im_i, double* P_iv)
{
   // R1-R2
   target_xv[0] = im_i[0] * P_iv[0*3 + 0] - im_i[1] * P_iv[1*3 + 0]; // dOH1_x
   target_xv[1] = im_i[0] * P_iv[0*3 + 1] - im_i[1] * P_iv[1*3 + 1]; // dOH1_y
   target_xv[2] = im_i[0] * P_iv[0*3 + 2] - im_i[1] * P_iv[1*3 + 2]; // dOH1_z
   // R2-R3
   target_xv[3] = im_i[1] * P_iv[1*3 + 0] - im_i[2] * P_iv[2*3 + 0]; // dOH2_x
   target_xv[4] = im_i[1] * P_iv[1*3 + 1] - im_i[2] * P_iv[2*3 + 1]; // dOH2_y
   target_xv[5] = im_i[1] * P_iv[1*3 + 2] - im_i[2] * P_iv[2*3 + 2]; // dOH2_z
   // R3-R1
   target_xv[6] = im_i[2] * P_iv[2*3 + 0] - im_i[0] * P_iv[0*3 + 0]; // dH1H2_x
   target_xv[7] = im_i[2] * P_iv[2*3 + 1] - im_i[0] * P_iv[0*3 + 1]; // dH1H2_y
   target_xv[8] = im_i[2] * P_iv[2*3 + 2] - im_i[0] * P_iv[0*3 + 2]; // dH1H2_z
}

INLINE void vec9_dot(double* lambda_x, double* newd_xv, double* d_xv)
{
   lambda_x[0] = newd_xv[0*3 + 0] * d_xv[0*3 + 0] +
                 newd_xv[0*3 + 1] * d_xv[0*3 + 1] +
                 newd_xv[0*3 + 2] * d_xv[0*3 + 2];
   lambda_x[1] = newd_xv[1*3 + 0] * d_xv[1*3 + 0] +
                 newd_xv[1*3 + 1] * d_xv[1*3 + 1] +
                 newd_xv[1*3 + 2] * d_xv[1*3 + 2];
   lambda_x[2] = newd_xv[2*3 + 0] * d_xv[2*3 + 0] +
                 newd_xv[2*3 + 1] * d_xv[2*3 + 1] +
                 newd_xv[2*3 + 2] * d_xv[2*3 + 2];
}

INLINE double vec3_dot(double* newd_xv, double* d_xv)
{
   return newd_xv[0] * d_xv[0] +
          newd_xv[1] * d_xv[1] +
          newd_xv[2] * d_xv[2];
}

INLINE void vec3_imul(double* target_x, double* a_x)
{
  target_x[0]*= a_x[0];
  target_x[1]*= a_x[1];
  target_x[2]*= a_x[2];
}

INLINE void vec3_div(double* target_x, double* a_x, double* b_x)
{
  target_x[0]= a_x[0] / b_x[0];
  target_x[1]= a_x[1] / b_x[1];
  target_x[2]= a_x[2] / b_x[2];
}

INLINE void vec3_axpy(double* target_v, double coeff, double* a_v)
{
  target_v[0] += coeff * a_v[0];
  target_v[1] += coeff * a_v[1];
  target_v[2] += coeff * a_v[2];
}

INLINE double sqr(double a)
{
  return a*a;
}

INLINE void vec9_sqrsum(double* target_x, double* R_xv)
{
  target_x[0] = sqr(R_xv[0*3 + 0]) + sqr(R_xv[0*3 + 1]) + sqr(R_xv[0*3 + 2]);
  target_x[1] = sqr(R_xv[1*3 + 0]) + sqr(R_xv[1*3 + 1]) + sqr(R_xv[1*3 + 2]);
  target_x[2] = sqr(R_xv[2*3 + 0]) + sqr(R_xv[2*3 + 1]) + sqr(R_xv[2*3 + 2]);
}

INLINE double vec3_sqrsum(double* R_v)
{
  return sqr(R_v[0]) + sqr(R_v[1]) + sqr(R_v[2]);
}

void vec9_print(char* title, double* target_iv)
{
  printf("%s:\n", title);
  printf("%20.15f %20.15f %20.15f\n", target_iv[0], target_iv[1], target_iv[2]);
  printf("%20.15f %20.15f %20.15f\n", target_iv[3], target_iv[4], target_iv[5]);
  printf("%20.15f %20.15f %20.15f\n", target_iv[6], target_iv[7], target_iv[8]);
}

void vec3_print(char* title, double* target_v)
{
  printf("%s:\n", title);
  printf("%20.15f %20.15f %20.15f\n", target_v[0], target_v[1], target_v[2]);
}

INLINE double coulomb(double Za, double Zb, double* Ra_v, double* Rb_v, double* Fa_v, double* Fb_v, unsigned char* pbc, double* celldiag)
{
  double d_v[3];
  vec3_submin(d_v, Ra_v, Rb_v, pbc, celldiag);
  double r2 = vec3_sqrsum(d_v);
  double r = sqrt(r2);
  double E = COULOMB_CONSTANT * Za * Zb / r;
  double str = E / r2;
  vec3_axpy(Fa_v, str, d_v);
  vec3_axpy(Fb_v, -str, d_v);
  return E;
}

INLINE double coulomb_cutoff(double Za, double Zb, double t, double* Ra_v, double* Rb_v, double* Fa_v, double* Fb_v, unsigned char* pbc, double* celldiag)
{
  double d_v[3];
  vec3_submin(d_v, Ra_v, Rb_v, pbc, celldiag);
  double r2 = vec3_sqrsum(d_v);
  double r = sqrt(r2);
  double E = COULOMB_CONSTANT * Za * Zb / r;
  double str = E / r2 * t;
  vec3_axpy(Fa_v, str, d_v);
  vec3_axpy(Fb_v, -str, d_v);
  return E; // The energy is not correct energy, it is multiplied later with t
}

INLINE double LJ(double A, double B, double* Ra_v, double* Rb_v, double* Fa_v, double* Fb_v, unsigned char* pbc, double* celldiag)
{
  double d_v[3];
  vec3_submin(d_v, Ra_v, Rb_v, pbc, celldiag);
  double r2 = vec3_sqrsum(d_v);
  double r4 = r2*r2;
  double r6 = r4*r2;
  double r8 = r4*r4;
  double r12 = r8*r4;
  double r14 = r12*r2;
  double str = 12*A/r14+6*B/r8;
  vec3_axpy(Fa_v, str, d_v);
  vec3_axpy(Fb_v, -str, d_v);
  return A / r12 + B/r6;
}

INLINE double pair_interaction_cutoff(double A, double B, double cutoff, double width, double* Z_i, double* Ra_iv, double* Rb_iv, double* Fa_iv, double* Fb_iv, unsigned char* pbc, double* celldiag)
{
  double d_v[3];
  vec3_submin(d_v, Ra_iv + 0*3, Rb_iv + 0*3, pbc, celldiag);
  double r2 = vec3_sqrsum(d_v);
  double t;
  double dtdr;
  double r;
  if (r2 > (cutoff*cutoff)) { return 0.0; }
  if (r2 < (cutoff-width)*(cutoff-width)) {
    t = 1;
    dtdr = 0;
    r = 1; // XXX Dummy value, this is not used because dtdr = 0
  } else {
    r = sqrt(r2);
    double y = (r - cutoff + width) / width;
    t = 1.0-y*y * (3.0 - 2.0 * y);
    dtdr = -6.0 / width * y * (1.0 - y);
  }

  double r4 = r2*r2;
  double r6 = r4*r2;
  double r8 = r4*r4;
  double r12 = r8*r4;
  double r14 = r12*r2;
  double E = (A / r12 + B/r6);
  double str = (12*A/r14+6*B/r8)* t - E*dtdr / r;
  E *= t;
  vec3_axpy(Fa_iv + 0*3, str, d_v);
  vec3_axpy(Fb_iv + 0*3, -str, d_v);

  for (unsigned int i1 = 0; i1 < 3; i1++)
    for (unsigned int i2 = 0; i2 < 3; i2++) {
       double Ept = coulomb_cutoff(Z_i[i1], Z_i[i2], t, Ra_iv+i1*3, Rb_iv+i2*3, Fa_iv+i1*3, Fb_iv+i2*3, pbc, celldiag);
       E+= Ept * t;
       double str = Ept / r * dtdr;
       vec3_axpy(Fa_iv + 0*3, -str, d_v); // Note that these are added to O, even for H-H interaction.
       vec3_axpy(Fb_iv + 0*3, str, d_v);
    }

  return E;

}

PyObject* adjust_positions(PyObject *self, PyObject *args)
{
  PyArrayObject* arraylen_x = 0;     // Input: the 3 constraint lengths
  PyArrayObject* arraymass_i = 0;    // Input: the 3 masses
  PyArrayObject* arrayR_niv = 0;     // Output: Adjusted positions will be written here.
  PyArrayObject* arraynewR_niv = 0;  // Input: gives the positions to be adjusted.

  if (!PyArg_ParseTuple(args, "OOOO", &arraylen_x, &arraymass_i, 
                        &arrayR_niv, &arraynewR_niv))
    {
      return NULL;
    }

  unsigned int NA = PyArray_DIM(arrayR_niv, 0);

  if (NA % 3 != 0)
  {
    PyErr_SetString(PyExc_TypeError, "Number of atoms not divisible with 3.");
    return NULL;
  }

  unsigned int NW = NA / 3;

  if (!(PyArray_NDIM(arraymass_i) == 1 &&
        PyArray_DIM(arraymass_i,0) == 3))
  {
    PyErr_SetString(PyExc_TypeError, "mass_i should be array with length 3.");
    return NULL;
  }

  if (!(PyArray_NDIM(arraylen_x) == 1 &&
        PyArray_DIM(arraylen_x,0) == 3))
  {
    PyErr_SetString(PyExc_TypeError, "len_x should be array with length 3.");
    return NULL;
  }

  double* len_x = DOUBLEP(arraylen_x);
  double* mass_i = DOUBLEP(arraymass_i);

  double len2_x[3];
  len2_x[0] = sqr(len_x[0]);
  len2_x[1] = sqr(len_x[1]);
  len2_x[2] = sqr(len_x[2]);

  double mu_x[3]; // Reduced masses of pairs
  mu_x[0] = 1.0 / ( (1.0/mass_i[0]) + (1.0/mass_i[1]) );
  mu_x[1] = 1.0 / ( (1.0/mass_i[1]) + (1.0/mass_i[2]) );
  mu_x[2] = 1.0 / ( (1.0/mass_i[2]) + (1.0/mass_i[0]) );

  double invm_i[3]; // Inverse masses of atoms
  invm_i[0] = 0.5 / mass_i[0]; // Note: /2 is embedded here
  invm_i[1] = 0.5 / mass_i[1];
  invm_i[2] = 0.5 / mass_i[2];

  double* R_niv = DOUBLEP(arrayR_niv);
  double* newR_niv = DOUBLEP(arraynewR_niv);

  for (int n=0; n<NW; n++) // For each water molecule
  {
      #ifdef FF_DEBUG
      vec9_print("R_iv", newR_niv+n*9);
      #endif
      double d_xv[9];
      // 1. Calculate (R1-R2), (R2-R3) and (R3-R1)
      vec9_diffs(d_xv, R_niv+n*9);
      #ifdef FF_DEBUG
      vec9_print("d_xv", d_xv);
      #endif

      // 2. Calculate (R1-R2)^2, (R2-R3)^2, (R3-R1)^2
      // double R_x[3];
      // vec9_sqrsum(R_x, d_xv);

      int iteration = 0;
      while (1) {
        #ifdef FF_DEBUG
        printf("Iteration: %d\n", iteration);
        vec9_print("newR_iv", newR_niv+n*9);
        #endif
        double newd_xv[9];
        // 1. Calculate (R1'-R2'), (R2'-R3') and (R3'-R1')
        vec9_diffs(newd_xv, newR_niv+n*9);
        #ifdef FF_DEBUG
        vec9_print("newd_xv", newd_xv);
        #endif

        // 2. Calculate (R1'-R2')^2, (R2'-R3')^2, (R3'-R1')^2
        double newR_x[3];
        vec9_sqrsum(newR_x, newd_xv);
        #ifdef FF_DEBUG
        vec3_print("newR_x", newR_x);
        #endif

        // 3. Calculate f12 = (R1'-R2')^2-D12^2, f23 and f31
        double f_x[3];
        vec3_sub(f_x, newR_x, len2_x);
        #ifdef FF_DEBUG
        vec3_print("f_x", f_x);
        #endif

        if (iteration++ > 1000) {
            printf("Warning: Adjust positions did not converge.\n");
            break;
        }
        if ((fabs(f_x[0]) < VELOCITY_VERLET_ADJUST_POSITION_TOL) &&
            (fabs(f_x[1]) < VELOCITY_VERLET_ADJUST_POSITION_TOL) &&
            (fabs(f_x[2]) < VELOCITY_VERLET_ADJUST_POSITION_TOL)) {
            break;
        }

        double lambda_x[3];
        double denom_x[3];
        // Calculate lambdas
        vec9_dot(denom_x, newd_xv, d_xv); // (R1-R2) . (R1'-R2') etc.
        vec3_div(lambda_x, f_x, denom_x); // f_12 / ( (R1-R2) . (R1'-R2') etc.
        vec3_imul(lambda_x, mu_x); // lambda /= (m_1^-1 + m_2^-1) etc.
        #ifdef FF_DEBUG
        vec3_print("lambda_x", lambda_x);
        #endif

        // Update newR's
        // newR_1 += -lambda12 * (R1-R2)
        vec3_axpy(newR_niv + n*9 + 0*3, -lambda_x[0] * invm_i[0], d_xv + 0*3);
        // newR_1 += lambda31 * (R3-R1)
        vec3_axpy(newR_niv + n*9 + 0*3, lambda_x[2] * invm_i[0], d_xv + 2*3);

        // newR_2 += +lambda12 * (R1-R2)
        vec3_axpy(newR_niv + n*9 + 1*3, +lambda_x[0] * invm_i[1], d_xv + 0*3);
        // newR_2 += -lambda23 * (R2-R3)
        vec3_axpy(newR_niv + n*9 + 1*3, -lambda_x[1] * invm_i[1], d_xv + 1*3);

        // newR_3 += +lambda23 * (R2-R3)
        vec3_axpy(newR_niv + n*9 + 2*3, +lambda_x[1] * invm_i[2], d_xv + 1*3);
        // newR_3 += -lambda31 * (R3-R1)
        vec3_axpy(newR_niv + n*9 + 2*3, -lambda_x[2] * invm_i[2], d_xv + 2*3);
      }
  }
  Py_RETURN_NONE;
}

PyObject* adjust_momenta(PyObject *self, PyObject *args)
{
  PyArrayObject* arraymass_i = 0;    // Input: the 3 masses
  PyArrayObject* arrayR_niv = 0;     //
  PyArrayObject* arraynewP_niv = 0;  // 

  if (!PyArg_ParseTuple(args, "OOO", &arraymass_i, &arrayR_niv,  
                        &arraynewP_niv))
    {
      return NULL;
    }

  unsigned int NA = PyArray_DIM(arrayR_niv, 0);

  if (NA % 3 != 0)
  {
    PyErr_SetString(PyExc_TypeError, "Number of atoms not divisible with 3.");
    return NULL;
  }

  unsigned int NW = NA / 3;

  if (!(PyArray_NDIM(arraymass_i) == 1 &&
        PyArray_DIM(arraymass_i,0) == 3))
  {
    PyErr_SetString(PyExc_TypeError, "mass_i should be array with length 3.");
    return NULL;
  }

  double* mass_i = DOUBLEP(arraymass_i);

  double invm_i[3]; // Inverse masses of atoms
  invm_i[0] = 1.0 / mass_i[0];
  invm_i[1] = 1.0 / mass_i[1];
  invm_i[2] = 1.0 / mass_i[2];

  double mu_x[3]; // Reduced masses of pairs
  mu_x[0] = 1.0 / ( (1.0/mass_i[0]) + (1.0/mass_i[1]) );
  mu_x[1] = 1.0 / ( (1.0/mass_i[1]) + (1.0/mass_i[2]) );
  mu_x[2] = 1.0 / ( (1.0/mass_i[2]) + (1.0/mass_i[0]) );

  double* R_niv = DOUBLEP(arrayR_niv);
  double* newP_niv = DOUBLEP(arraynewP_niv);

  for (int n=0; n<NW; n++) // For each water molecule
  {
      #ifdef FF_DEBUG_M
      vec9_print("R_iv", newP_niv+n*9);
      #endif
      double d_xv[9];
      // 1. Calculate (R1-R2), (R2-R3) and (R3-R1)
      vec9_diffs(d_xv, R_niv+n*9);
      #ifdef FF_DEBUG_M
      vec9_print("d_xv", d_xv);
      #endif

      int iteration = 0;
      while (1) {
        #ifdef FF_DEBUG_M
        printf("Iteration: %d\n", iteration);
        vec9_print("newP_iv", newP_niv+n*9);
        #endif
        double newd_xv[9];
        // 1. Calculate (m_1^-1 P1'-m_2^-1 P2'), etc.
        vec9_massdiffs(newd_xv, invm_i, newP_niv+n*9);

        #ifdef FF_DEBUG_M
        vec9_print("(m_1^-1 P1'-m_2^-1 P2')", newd_xv);
        #endif

        // 3. Calculate g12 = (R1-R2) . (m_1^-1 P1'-m_2^-1 P2') etc.
        double g_x[3];
        vec9_dot(g_x, newd_xv, d_xv);

        #ifdef FF_DEBUG_M
        vec3_print("g_x", g_x);
        #endif

        if (iteration++ > 1000) {
            printf("Warning: Adjust velocities did not converge.\n");
            break;
        }

        if ((fabs(g_x[0]) < VELOCITY_VERLET_ADJUST_VELOCITY_TOL) &&
            (fabs(g_x[1]) < VELOCITY_VERLET_ADJUST_VELOCITY_TOL) &&
            (fabs(g_x[2]) < VELOCITY_VERLET_ADJUST_VELOCITY_TOL)) {
            break;
        }

        double lambda_x[3];
        double denom_x[3];

        vec9_dot(denom_x, d_xv, d_xv); // (R1-R2)^2 etc.
        vec3_div(lambda_x, g_x, denom_x); // g_12 / (R1-R2)^2
        vec3_imul(lambda_x, mu_x); // lambda /= (m_1^-1 + m_2^-1) etc.
        
        #ifdef FF_DEBUG_M
        vec3_print("lambda_x", lambda_x);
        #endif

        // Update newR's
        // newR_1 += -lambda12 * (R1-R2)
        vec3_axpy(newP_niv + n*9 + 0*3, -lambda_x[0], d_xv + 0*3);
        // newR_1 += lambda31 * (R3-R1)
        vec3_axpy(newP_niv + n*9 + 0*3, lambda_x[2], d_xv + 2*3);

        // newR_2 += +lambda12 * (R1-R2)
        vec3_axpy(newP_niv + n*9 + 1*3, +lambda_x[0], d_xv + 0*3);
        // newR_2 += -lambda23 * (R2-R3)
        vec3_axpy(newP_niv + n*9 + 1*3, -lambda_x[1], d_xv + 1*3);

        // newR_3 += +lambda23 * (R2-R3)
        vec3_axpy(newP_niv + n*9 + 2*3, +lambda_x[1], d_xv + 1*3);
        // newR_3 += -lambda31 * (R3-R1)
        vec3_axpy(newP_niv + n*9 + 2*3, -lambda_x[2], d_xv + 2*3);
      }
  }
  Py_RETURN_NONE;
}

PyObject* calculate_forces_H2O(PyObject *self, PyObject *args)
{
  double A;
  double B;
  PyArrayObject* arraypbc = 0;
  PyArrayObject* arrayZ_i = 0;
  PyArrayObject* arrayR_niv = 0;
  PyArrayObject* arrayF_niv = 0;
  PyArrayObject* arraycell = 0;

  double cutoff=0;
  double width=0;
 
  if (!PyArg_ParseTuple(args, "OOddddOOO", &arraypbc, &arraycell, &A, &B, &cutoff, &width, &arrayZ_i, &arrayR_niv, &arrayF_niv))
    {
      return NULL;
    }

  if (!(PyArray_NDIM(arraypbc) == 1 &&
        PyArray_DIM(arraypbc,0) == 3)) {
    PyErr_SetString(PyExc_TypeError, "pbc should be array with length 3.");
    return NULL;
  }

  if (!(PyArray_NDIM(arrayZ_i) == 1 &&
        PyArray_DIM(arrayZ_i,0) == 3)) {
    PyErr_SetString(PyExc_TypeError, "Z_i should be array with length 3.");
    return NULL;
  }

  if (!(PyArray_NDIM(arraycell) == 2 && 
       (PyArray_DIM(arraycell,0) == 3) &&
       (PyArray_DIM(arraycell,1) == 3))) {
    PyErr_SetString(PyExc_TypeError, "Cell should be array with size 3x3.");
    return NULL;
  }

  double* cell_vc = DOUBLEP(arraycell);
  //printf("Cell\n");
  for (unsigned int v1=0; v1<3; v1++) 
    for (unsigned int v2=0; v2<3; v2++) {
      //printf("%20.15f \n", cell_vc[v1+v2*3]);
      if (v1 != v2) 
        if (fabs(cell_vc[v1+v2*3]) > 1e-10) {
          PyErr_SetString(PyExc_TypeError, "Cell array should be diagonal.");
          return NULL;
        }
    }

  double celldiag[3];
  celldiag[0] = cell_vc[0];
  celldiag[1] = cell_vc[1+1*3];
  celldiag[2] = cell_vc[2+2*3];

  unsigned int NA = PyArray_DIM(arrayR_niv, 0);

  if (NA % 3 != 0)
  {
    PyErr_SetString(PyExc_TypeError, "Number of atoms not divisible with 3.");
    return NULL;
  }

  char *pbc = (char*) PyArray_DATA(arraypbc);
  //printf("Boundary conditions %u %u %u \n", pbc[0], pbc[1], pbc[2]);

  unsigned int NW = NA / 3;

  double* Z_i = DOUBLEP(arrayZ_i);
  double* R_niv = DOUBLEP(arrayR_niv);
  double* F_niv = DOUBLEP(arrayF_niv);

  double E = 0.0;
  if (cutoff <= 0.0) { // No cutoff 
    // Double loop of atoms
    for (unsigned int n1 = 0; n1 < NW; n1++)
      for (unsigned int n2 = n1+1; n2 < NW; n2++) {
        E += LJ(A, B, R_niv+n1*9, R_niv+n2*9, F_niv+n1*9, F_niv+n2*9, pbc, celldiag);
        for (unsigned int i1 = 0; i1 < 3; i1++)
           for (unsigned int i2 = 0; i2 < 3; i2++)
             E+= coulomb(Z_i[i1], Z_i[i2], R_niv+n1*9+i1*3, R_niv+n2*9+i2*3, F_niv+n1*9+i1*3, F_niv+n2*9+i2*3, pbc, celldiag);
        }
  } else { // With cutoff
    // Double loop of atoms
    for (unsigned int n1 = 0; n1 < NW; n1++)
      for (unsigned int n2 = n1+1; n2 < NW; n2++) {
        E += pair_interaction_cutoff(A, B, cutoff, width, Z_i, R_niv+n1*9, R_niv+n2*9, F_niv+n1*9, F_niv+n2*9, pbc, celldiag);
      }
 }

  return Py_BuildValue("d", E);
}


PyObject* adjust_positions_general(PyObject *self, PyObject *args)
{
  PyArrayObject* arraypairs_px = 0;
  PyArrayObject* arraymasses_a = 0; 
  PyArrayObject* arraybondlengths_p = 0; 
  PyArrayObject* arrayR_av = 0; 
  PyArrayObject* arraynewR_av = 0; 

  if (!PyArg_ParseTuple(args, "OOOOO", &arraypairs_px, &arraymasses_a, &arraybondlengths_p, &arrayR_av, &arraynewR_av))
    {
      return NULL;
    }

  if (!((PyArray_NDIM(arraypairs_px) == 2) &&
        (PyArray_DIM(arraypairs_px,1) == 2)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of pairs should be (N, 2).");
    return NULL;
  }

  if (!((PyArray_NDIM(arrayR_av) == 2) &&
        (PyArray_DIM(arrayR_av,1) == 3)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of positions should be (NA, 3).");
    return NULL;
  }

  if (!((PyArray_NDIM(arraynewR_av) == 2) &&
        (PyArray_DIM(arraynewR_av,1) == 3)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of new positions should be (NA, 3).");
    return NULL;
  }

  unsigned int NA = PyArray_DIM(arrayR_av, 0); // Number of atoms
  unsigned int NP = PyArray_DIM(arraypairs_px, 0); // Number of pairs

  if (!((PyArray_NDIM(arraymasses_a) == 1) &&
        (PyArray_DIM(arraymasses_a,0) == NA)))
  {
    PyErr_SetString(PyExc_TypeError, "masses should be array with length NA.");
    return NULL;
  }

  int* pairs_px = INTP(arraypairs_px);
  double* masses_a = DOUBLEP(arraymasses_a);
  double* bondlengths_p = DOUBLEP(arraybondlengths_p);
  double* R_av = DOUBLEP(arrayR_av);
  double* newR_av = DOUBLEP(arraynewR_av);

  unsigned int iteration = 0;

  while (1) {
    int converged = 1;
    for (unsigned int p=0; p<NP; p++) { 
      int a1 = pairs_px[2*p];
      int a2 = pairs_px[2*p+1];
      double cd = bondlengths_p[p];

      double d0_v[3]; // Omitting mic here for now
      vec3_sub(d0_v, R_av + a1*3, R_av + a2*3);
      double d1_v[3];
      vec3_sub(d1_v, newR_av + a1*3, newR_av + a2*3);
      
      double m = 1.0 / (1.0 / masses_a[a1] + 1.0 / masses_a[a2]);
      double temp2 = vec3_sqrsum(d1_v);
      double temp3 = vec3_sqrsum(d0_v);
      double x = 0.5 * (cd*cd - vec3_sqrsum(d1_v)) / vec3_dot(d0_v,d1_v);
      if (fabs(x) > VELOCITY_VERLET_ADJUST_POSITION_TOL) { converged = 0; }
      vec3_axpy(newR_av + a1*3, x * m / masses_a[a1], d0_v);
      vec3_axpy(newR_av + a2*3, -x * m / masses_a[a2], d0_v);
    }
    if (converged) break;
    if (iteration++ > 1000) break;
   }

  Py_RETURN_NONE;
}

PyObject* adjust_momenta_general(PyObject *self, PyObject *args)
{
  PyArrayObject* arraypairs_px = 0;
  PyArrayObject* arraymasses_a = 0; 
  PyArrayObject* arraybondlengths_p = 0; 
  PyArrayObject* arrayR_av = 0; 
  PyArrayObject* arraynewP_av = 0; 

  if (!PyArg_ParseTuple(args, "OOOOO", &arraypairs_px, &arraymasses_a, &arraybondlengths_p, &arrayR_av, &arraynewP_av))
    {
      return NULL;
    }

  if (!((PyArray_NDIM(arraypairs_px) == 2) &&
        (PyArray_DIM(arraypairs_px,1) == 2)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of pairs should be (N, 2).");
    return NULL;
  }

  if (!((PyArray_NDIM(arrayR_av) == 2) &&
        (PyArray_DIM(arrayR_av,1) == 3)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of positions should be (NA, 3).");
    return NULL;
  }

  if (!((PyArray_NDIM(arraynewP_av) == 2) &&
        (PyArray_DIM(arraynewP_av,1) == 3)))
  {
    PyErr_SetString(PyExc_TypeError, "Dimension of new momentum should be (NA, 3).");
    return NULL;
  }

  unsigned int NA = PyArray_DIM(arrayR_av, 0); // Number of atoms
  unsigned int NP = PyArray_DIM(arraypairs_px, 0); // Number of pairs

  if (!((PyArray_NDIM(arraymasses_a) == 1) &&
        (PyArray_DIM(arraymasses_a,0) == NA)))
  {
    PyErr_SetString(PyExc_TypeError, "masses should be array with length NA.");
    return NULL;
  }

  int* pairs_px = INTP(arraypairs_px);
  double* masses_a = DOUBLEP(arraymasses_a);
  double* bondlengths_p = DOUBLEP(arraybondlengths_p);
  double* R_av = DOUBLEP(arrayR_av);
  double* newP_av = DOUBLEP(arraynewP_av);

  unsigned int iteration = 0;

  while (1) {
    int converged = 1;
    for (unsigned int p=0; p<NP; p++) {
      int a1 = pairs_px[2*p];
      int a2 = pairs_px[2*p+1];
      double cd = bondlengths_p[p];
      double d0_v[3]; // Omitting mic here for now
      vec3_sub(d0_v, R_av + a1*3, R_av + a2*3);
      double dv_v[3];
      vec3_zero(dv_v);
      vec3_axpy(dv_v, 1.0 / masses_a[a1], newP_av + 3*a1);       
      vec3_axpy(dv_v, -1.0 / masses_a[a2], newP_av + 3*a2);       
      double m = 1.0 / (1.0 / masses_a[a1] + 1.0 / masses_a[a2]);
      double x = -vec3_dot(d0_v, dv_v) / (cd*cd);

      if (fabs(x) > VELOCITY_VERLET_ADJUST_VELOCITY_TOL) { converged = 0; }
  
      vec3_axpy(newP_av + a1*3, x * m, d0_v);
      vec3_axpy(newP_av + a2*3, -x * m, d0_v);
    }
    if (converged) break;
    if (iteration++ > 1000) break;
   }

  Py_RETURN_NONE;
}
