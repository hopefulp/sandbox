/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hendrik Heenen (Technical University of Munich)
                        (hendrik.heenen at mytum.com)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "angle_cosine_vdwl13.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineVdwl13::AngleCosineVdwl13(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCosineVdwl13::~AngleCosineVdwl13()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(multiplicity);
    memory->destroy(th0);
    memory->destroy(flj);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineVdwl13::compute(int eflag, int vflag)
{
  int i,i1,i2,i3,n,type,itype,jtype;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double tk;

  // extra lj variables
  double delx3,dely3,delz3,rsq3,fpair,evdwl;

  eangle = evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;
  
  // insure pair->ev_tally() will use 1-3 virial contribution

  if (vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int newton_bond = force->newton_bond;
  int *atomtype = atom->type;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // c = cosine of angle

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;
    
    // force & energy

    // explicit lj-contribution

    itype = atomtype[i1];
    jtype = atomtype[i3];

    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

    evdwl = force->pair->single(i1,i3,itype,jtype,rsq3,0.0,flj[type],fpair);

    // add forces of additional LJ interaction
    if (newton_pair || i1 < nlocal) {
      f[i1][0] += delx3*fpair;
      f[i1][1] += dely3*fpair;
      f[i1][2] += delz3*fpair;
    }
    if (newton_pair || i3 < nlocal) {
      f[i3][0] -= delx3*fpair;
      f[i3][1] -= dely3*fpair;
      f[i3][2] -= delz3*fpair;
    }

    //update pair energy and velocities

    if (evflag) force->pair->ev_tally(i1,i3,nlocal,newton_pair,
                                      evdwl,0.0,fpair,delx3,dely3,delz3); 
    
    tk = multiplicity[type]*acos(c)-th0[type];

    if (eflag) eangle = k[type]*(1.0+cos(tk));

    a = k[type]*multiplicity[type]*sin(tk)*s;

    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineVdwl13::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(multiplicity,n+1,"angle:multiplicity");
  memory->create(th0,n+1,"angle:th0");
  memory->create(flj,n+1,"angle:flj");
  
  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleCosineVdwl13::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  double c_one = force->numeric(FLERR,arg[1]);
  int n_one = force->inumeric(FLERR,arg[2]);
  int th0_one = force->numeric(FLERR,arg[3]);
  int flj_one = force->numeric(FLERR,arg[4]);
  if (n_one <= 0) error->all(FLERR,"Incorrect args for angle coefficients");
  if (flj_one < 0.0 || flj_one > 1.0) 
    error->all(FLERR,"Incorrect args for angle coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = c_one;
    multiplicity[i] = n_one;
    // transform offset angle to radians
    th0[i] = th0_one/180.0 * MY_PI;
    flj[i] = flj_one;

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   check special_bond settings are valid and initialize vdwl parameters
------------------------------------------------------------------------- */

void AngleCosineVdwl13::init_style()
{
  if (force->pair == NULL)
    error->all(FLERR,"No pair style is defined for angle style");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support angle style");
  
  // special bond multi-count warning
  if (force->special_lj[2] != 0.0)
    error->warning(FLERR,"LJ 1-3 interactions are active and may "
                           " be counted twice");
}

/* ---------------------------------------------------------------------- */

double AngleCosineVdwl13::equilibrium_angle(int i)
{
  return MY_PI;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCosineVdwl13::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&multiplicity[1],sizeof(int),atom->nangletypes,fp);
  fwrite(&th0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&flj[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCosineVdwl13::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&multiplicity[1],sizeof(int),atom->nangletypes,fp);
    fread(&th0[1],sizeof(double),atom->nangletypes,fp);
    fread(&flj[1],sizeof(double),atom->nangletypes,fp);
  }

  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&multiplicity[1],atom->nangletypes,MPI_INT,0,world);
  MPI_Bcast(&th0[1],atom->nangletypes,MPI_INT,0,world);
  MPI_Bcast(&flj[1],atom->nangletypes,MPI_INT,0,world);
  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCosineVdwl13::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++) {
    fprintf(fp,"%d %g %d %d %d\n",i,k[i],multiplicity[i],th0[i],flj[i]);
  }
}

/* ---------------------------------------------------------------------- */

double AngleCosineVdwl13::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double tk = multiplicity[type]*acos(c)-th0[type];

  return k[type]*(1.0+cos(tk));
}
