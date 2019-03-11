/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(cosine/vdwl13, AngleCosineVdwl13)

#else

#ifndef LMP_ANGLE_COSINE_VDWL13_H
#define LMP_ANGLE_COSINE_VDWL13_H

#include <stdio.h>
#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineVdwl13 : public Angle {
 public:
  AngleCosineVdwl13(class LAMMPS *);
  virtual ~AngleCosineVdwl13();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *k,*th0;
  double *flj;
  double *eps,*d0;
  int *multiplicity;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
