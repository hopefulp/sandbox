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
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_buck6d_coul_dsf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "math_special.h"

#include <iostream>
#include <iomanip>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBuck6dCoulDSF::PairBuck6dCoulDSF(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairBuck6dCoulDSF::~PairBuck6dCoulDSF()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);

      memory->destroy(cut_lj);
      memory->destroy(cut_ljsq);
      memory->destroy(epsilon);
      memory->destroy(d0);
      memory->destroy(d02);
      memory->destroy(d06);
      memory->destroy(d014);
      memory->destroy(alpha_ij);
      memory->destroy(f_shift_ij);
      memory->destroy(e_shift_ij);
      memory->destroy(buck6d1);
      memory->destroy(buck6d2);
      memory->destroy(buck6d3);
      memory->destroy(buck6d4);
      memory->destroy(offset);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBuck6dCoulDSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,r14inv,rexp,forcecoul,forcebuck6d,factor_coul,factor_lj;
  double prefactor,erfcc,erfcd,t,arg;
  double term1,term2,term3,term4,term5;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;

        if (rsq < cut_ljsq[itype][jtype]) {
          r = sqrt(rsq);
          r6inv = r2inv*r2inv*r2inv;
          r14inv = r6inv*r6inv*r2inv;
          rexp = exp(-r*buck6d2[itype][jtype]);
          term1 = buck6d3[itype][jtype]*r6inv;
          term2 = buck6d4[itype][jtype]*r14inv;
          term3 = term2*term2;
          term4 = 1.0/(1.0 + term2);
          term5 = 1.0/(1.0 + 2.0*term2 + term3);
          forcebuck6d = buck6d1[itype][jtype]*buck6d2[itype][jtype]*r*rexp; 
          forcebuck6d -= term1*(6.0*term4 - term5*14.0*term2);
        } else forcebuck6d = 0.0;

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j]/r; // dlpoly constant offset*1.0000001960334732;


        //below parametrization for: erfcc = erf(alpha_ij[itype][jtype]*r);
          arg = alpha_ij[itype][jtype]*r;
          erfcd = MathSpecial::expmsq(arg);
          erfcc = 1 - (MathSpecial::my_erfcx(arg) * erfcd);
          
          forcecoul = prefactor * ((erfcc/r) - (2.0/MY_PIS*alpha_ij[itype][jtype]*erfcd) + 
                                                r*f_shift_ij[itype][jtype]) * r;
        
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        fpair = (forcecoul + factor_lj*forcebuck6d) * r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = buck6d1[itype][jtype]*rexp - term1*term4 -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;

          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r*e_shift_ij[itype][jtype] - 
                                 rsq*f_shift_ij[itype][jtype]);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(d02,n+1,n+1,"pair:d02");
  memory->create(d06,n+1,n+1,"pair:d06");
  memory->create(d014,n+1,n+1,"pair:d014");
  memory->create(alpha_ij,n+1,n+1,"pair:alpha_ij");
  memory->create(f_shift_ij,n+1,n+1,"pair:f_shift_ij");
  memory->create(e_shift_ij,n+1,n+1,"pair:e_shift_ij");
  memory->create(buck6d1,n+1,n+1,"pair:buck6d1");
  memory->create(buck6d2,n+1,n+1,"pair:buck6d2");
  memory->create(buck6d3,n+1,n+1,"pair:buck6d3");
  memory->create(buck6d4,n+1,n+1,"pair:buck6d4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::coeff(int narg, char **arg)
{
//if (narg < 5 || narg > 6)
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

//double epsilon_one = force->numeric(FLERR,arg[2]);
//double d0_one = force->numeric(FLERR,arg[3]);
//double alpha_one = force->numeric(FLERR,arg[4]);
  double buck6d1_one = force->numeric(FLERR,arg[2]);
  double buck6d2_one = force->numeric(FLERR,arg[3]);
  double buck6d3_one = force->numeric(FLERR,arg[4]);
  double buck6d4_one = force->numeric(FLERR,arg[5]);
  double alpha_one = force->numeric(FLERR,arg[6]);

  double cut_lj_one = cut_lj_global;
  if (narg == 6) cut_lj_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
    //epsilon[i][j] = epsilon_one;
    //d0[i][j] = d0_one;
      buck6d1[i][j] = buck6d1_one;
      buck6d2[i][j] = buck6d2_one;
      buck6d3[i][j] = buck6d3_one;
      buck6d4[i][j] = buck6d4_one;
      alpha_ij[i][j] = alpha_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style buck6d/coul/dsf requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuck6dCoulDSF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  
//d02[i][j] = d0[i][j]*d0[i][j];
//d06[i][j] = d02[i][j]*d02[i][j]*d02[i][j];
//d014[i][j] = d06[i][j]*d06[i][j]*d02[i][j];
  
//buck6d1[i][j] = epsilon[i][j]*1.84e5;
//buck6d2[i][j] = 12.0/d0[i][j];
//buck6d3[i][j] = epsilon[i][j]*2.25*d06[i][j];
//buck6d4[i][j] = 6.0*3.725290298461914e-9*d014[i][j];

  // check if this is correct
  if (offset_flag) {
    double term1 = buck6d3[i][j]/pow(cut_lj[i][j],6.0);
    double term4 = 1.0/(1.0 + (buck6d4[i][j]/pow(cut_lj[i][j],14.0)));
    double rexp = exp(-cut_lj[i][j]*buck6d2[i][j]);
    offset[i][j] = buck6d1[i][j]*rexp - term1*term4;
  } else offset[i][j] = 0.0;
  
  double erfcd_c = exp(-alpha_ij[i][j]*alpha_ij[i][j]*cut_coul*cut_coul);
  double erfcc_c = erf(alpha_ij[i][j]*cut_coul);
  //f_shift = -(erfcc_c/cut_coulsq + 2.0/MY_PIS*alpha_ij[itype][jtype]*erfcd_c/cut_coul); //NOTE old f_shift from erfc //keep until PW is without problem
  f_shift_ij[i][j] = -erfcc_c/cut_coulsq + 2.0/MY_PIS*alpha_ij[i][j]*erfcd_c/cut_coul;
  e_shift_ij[i][j] = erfcc_c/cut_coul - f_shift_ij[i][j]*cut_coul;
  
  cut_ljsq[j][i] = cut_ljsq[i][j];
  alpha_ij[j][i] = alpha_ij[i][j];
  f_shift_ij[j][i] = f_shift_ij[i][j];
  e_shift_ij[j][i] = e_shift_ij[i][j];
  buck6d1[j][i] = buck6d1[i][j];
  buck6d2[j][i] = buck6d2[i][j];
  buck6d3[j][i] = buck6d3[i][j];
  buck6d4[j][i] = buck6d4[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  /*
  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
               sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
               sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }
  */
  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha_ij[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
	    }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&d0[i][j],sizeof(double),1,fp);
          fread(&alpha_ij[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha_ij[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,
            epsilon[i][i],d0[i][i],alpha_ij[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulDSF::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,
              epsilon[i][j],d0[i][j],alpha_ij[i][j],cut_lj[i][j]);
}
/* ---------------------------------------------------------------------- */
//write when works
/*
double PairBuck6dCoulDSF::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,erfcc,erfcd,prefactor;
  double forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  } else forcelj = 0.0;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = factor_coul * force->qqrd2e * atom->q[i]*atom->q[j]/r;
    erfcc = erfc(alpha*r);
    erfcd = exp(-alpha*alpha*r*r);
    forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +
      r*f_shift) * r;
  } else forcecoul = 0.0;

  fforce = (forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
      offset[itype][jtype];
    eng += factor_lj*philj;
  }

  if (rsq < cut_coulsq) {
    phicoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
    eng += phicoul;
  }

  return eng;
}
*/
/* ---------------------------------------------------------------------- */

void *PairBuck6dCoulDSF::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}
