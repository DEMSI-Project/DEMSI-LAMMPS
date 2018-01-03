/* ----------------------------------------------------------------------
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Dan S. Bolintineanu (SNL), Adrian K. Turner (LANL)
   ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_gran_hopkins.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-10

/* ---------------------------------------------------------------------- */

PairGranHopkins::PairGranHopkins(LAMMPS *lmp) :
      PairGranHookeHistory(lmp, 12)
{
  history_ndim = 12;
}

/* ---------------------------------------------------------------------- */
PairGranHopkins::~PairGranHopkins()
{
}
/* ---------------------------------------------------------------------- */

void PairGranHopkins::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  int itype,jtype;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body
  // Not yet applicable for DEMSI, but may be something to look into

  if (fix_rigid && neighbor->ago == 0){
      int tmp;
      int *body = (int *) fix_rigid->extract("body",tmp);
      double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
      if (atom->nmax > nmax) {
	  memory->destroy(mass_rigid);
	  nmax = atom->nmax;
	  memory->create(mass_rigid,nmax,"pair:mass_rigid");
      }
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++)
	if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
	else mass_rigid[i] = 0.0;
      comm->forward_comm_pair(this);
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  // loop over my elements
  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = atom->type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      allhistory = firsthistory[i];

      // loop over neighbors of each element
      for (jj = 0; jj < jnum; jj++) {
	  j = jlist[jj];
	  jtype = atom->type[j];
	  j &= NEIGHMASK;

	  history = &allhistory[history_ndim*jj];
	  //'history' now points to the ii-jj array that stores 
	  //all the history associated with pair ii-jj
	  //For bonded pairs:
	  // history[0-7]: x,y components of s1i, s2i, s1j, s2j
	  // history[8,9]: chi1, chi2
	  // history[10]: bond length
	  // history[11]: bond thickness, h
	  //For unbonded pairs:
	  // history[0,1]: normal force at previous time step, x and y components
	  // history[2,3]: accumulated tangential displacement at contact, x and y
	  // history[4] : contact thickness h

	  if (history[8] >= history[9]){ // Un-bonded, chi1 >= chi2
	      compute_nonbonded(history, firsttouch[i], i, j);
	  }
	  else { //Bonded	    
	      compute_bonded(history, i, j);
	  }
      }
  }
}

void PairGranHopkins::compute_nonbonded(double *history, int* touch, int i, int j){

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delx, dely;
  double r, rsq, rinv, radsum, radmin;
  double nx, ny, vrx, vry, vnnr, vnx, vny;
  double wrz, vtrx, vtry, vtx, vty, vrel;
  double delta, delta_dot;

  double fnx, fny;
  double sig_c, hmin;
  double hprime, kp, kr, ke, L, hstar;
  double num, denom, fnmag_plastic, fnmag_elastic, fnmag;

  double ndisp, dispmag, scalefac;
  double ftx, fty, ftmag, ftcrit;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  //Hard-coded, ask Adrian what this should actually be??
  hstar = 0.3;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  rsq = delx*delx + dely*dely;
  radsum = radius[i] + radius[j];

  if (rsq >= radsum*radsum){
      touch[j] = 0;
      for (int i = 0; i < 4; i++)
	history[i] = 0;
  }
  else{
      touch [j] = 1;
      r = sqrt(rsq);
      rinv = 1.0/r;
      nx = delx/r;
      ny = dely/r;

      radmin = MIN(radius[i], radius[j]);
      L = 2*(1+radmin)*(r/abs(radius[i] - radius[j]));

      // relative translational velocity
      vrx = v[i][0] - v[j][0];
      vry = v[i][1] - v[j][1];

      delta = radsum - r;
      delta_dot = -vnnr;

      // Compute plastic normal force
      hprime = history[4];
      ke = Emod/L*(1/(1/atom->mean_thickness[i] + 1/atom->mean_thickness[j]));
      hmin = MIN(atom->min_thickness[i], atom->min_thickness[j]);
      if (hprime < hstar){
	  kr = 26126*hprime;
	  kp = 928*hprime*hprime;
      }
      else{
	  kr = kp = hprime*sig_c;
      }

      num = history[0]/(kp*update->dt) + delta_dot*L + delta*L*ke/damp_normal + ke*kr/(damp_normal*kp);
      denom = 1/(kp*update->dt) + 1/damp_normal*(1+ke/kp);
      fnmag_plastic = num/denom;

      // Elastic normal force
      fnmag_elastic = ke*delta*L + damp_normal*delta_dot*L;

      fnmag = MIN(fnmag_elastic, fnmag_plastic);

      fnx = fnmag*nx;
      fny = fnmag*ny;
      f[i][0] += fnx;
      f[i][1] += fny;
      history[0] = fnx;
      history[1] = fny;

      // Compute tangential force
      // normal component of relative translational velocity
      vnnr = vrx*nx + vry*ny;
      vnx = nx*vnnr;
      vny = ny*vnnr;

      // subtract to compute tangential component of relative translational velocity
      vtrx = vrx - vnx;
      vtry = vry - vny;

      // total relative tangential velocities at contact
      wrz = radius[i]*omega[i][2] + radius[j]*omega[j][2];
      vtx = vtx + ny*wrz;
      vty = vty - nx*wrz;
      vrel = vtx*vtx + vty*vty;
      vrel = sqrt(vrel);

      // update tangential displacement, rotate if needed
      if (historyupdate){
	  ndisp = nx*history[2] + ny*history[3];
	  dispmag =sqrt( history[2]*history[2] + history[3]*history[3]);
	  denom = dispmag - ndisp;
	  if (ndisp > EPSILON && denom != 0){
	      scalefac = dispmag/denom;
	      history[2] -= ndisp*nx;
	      history[3] -= ndisp*ny;
	      history[2] *= scalefac;
	      history[3] *= scalefac;
	  }
	  history[2] += vtx*update->dt;
	  history[3] += vty*update->dt;
      }

      // total tangential force
      ftx = - (kt*history[2] + damp_tangential*vtx);
      fty = - (kt*history[3] + damp_tangential*vty);

      ftmag = sqrt(ftx*ftx + fty*fty);
      ftcrit = friction_tangential*fabs(fnmag);
      if (ftmag > ftcrit){
	  if (dispmag != 0){
	      ftx *= ftcrit/ftmag;
	      fty *= ftcrit/ftmag;
	      history[2] = -(ftcrit + damp_tangential*vtx)/kt;
	      history[3] = -(ftcrit + damp_tangential*vty)/kt;
	  }
	  else ftx = fty = 0;
      }

      //Apply tangential forces
      f[i][0] += ftx;
      f[i][1] += fty;

      // torque induced by tangential force (don't apply to particle j, since we require newton_off for now)
      torque[i][0] += -radius[i]*(nx*fty - ny*ftx);

  }

}

void PairGranHopkins::compute_bonded(double *history, int i, int j){
  double s1x, s1y, s2x, s2y, mx, my, mmag, mex, mey, bex, bey, rx, ry;
  double An, Bn, Cn, Dn, Bt, Ct, Dt;
  double Fnmag, Ftmag, Nn, Nt;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;

  double chi1 = history[8];
  double chi2 = history[9];
  double chidiff, chidiff2, chidiff3;

  double sig_n1, sig_s1;
  double chi_c, chi_t, chi_s1, chi_s2;
  double nprefac, sprefac;
  double damp_prefac, area_bond, fdampx, fdampy, torquedamp;
  double hmin;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  if (historyupdate){
      // Update bond end points based on particle translations
      history[0] += dt*v[i][0];
      history[1] += dt*v[i][1];
      history[2] += dt*v[i][0];
      history[3] += dt*v[i][1];

      history[4] += dt*v[j][0];
      history[5] += dt*v[j][1];
      history[6] += dt*v[j][0];
      history[7] += dt*v[j][1];

      // Update bond end points based on particle rotations
      history[0] += -dt*omega[i][2]*(history[1]-x[i][1]);
      history[1] +=  dt*omega[i][2]*(history[0]-x[i][0]);
      history[2] += -dt*omega[i][2]*(history[3]-x[i][1]);
      history[3] +=  dt*omega[i][2]*(history[2]-x[i][0]);

      history[4] += -dt*omega[j][2]*(history[5]-x[j][1]);
      history[5] +=  dt*omega[j][2]*(history[4]-x[j][0]);
      history[6] += -dt*omega[j][2]*(history[7]-x[j][1]);
      history[7] +=  dt*omega[j][2]*(history[6]-x[j][0]);
  }

  //Compute s_1, s_2, m, m_e, b_e
  s1x = history[4] - history[0];
  s1y = history[5] - history[1];
  s2x = history[6] - history[2];
  s2y = history[7] - history[3];

  mx = history[2] + 0.5*s2x - history[0] - 0.5*s1x;
  my = history[3] + 0.5*s2y - history[1] - 0.5*s1y;
  mmag = sqrt(mx*mx + my*my);
  mex = mx/mmag;
  mey = my/mmag;

  bex = -mey;
  bey = mex;

  rx = history[0] + 0.5*s1x - x[i][0];
  ry = history[1] + 0.5*s1y - x[i][1];

  //Compute forces and torques
  Dn = s1x*bex + s1y*bey;
  Cn = s2x*bex + s2y*bey - Dn;
  Bn = rx*bey - ry*bex;
  An = mx*bey - my*bex;

  Dt = s1x*mex + s1y*mey;
  Ct = s2x*mex + s2y*mey - Dt;
  Bt = rx*mey - ry*mex;

  chidiff = chi2-chi1;
  chidiff2 = 0.5*(chi2*chi2 - chi1*chi1);
  chidiff3 = MathConst::THIRD*(chi2*chi2*chi2 - chi1*chi1*chi1);
  nprefac = history[11]*history[10]*kn;
  sprefac = history[11]*history[10]*kt;

  Fnmag = nprefac*(Dn*chidiff + Cn*chidiff2);
  Ftmag = sprefac*(Dt*chidiff + Ct*chidiff2);
  Nn = nprefac*(An*Dn*chidiff3 + (Bn*Cn+An*Dn)*chidiff2 + Bn*Dn*chidiff);
  Nt = sprefac*(Bt*Ct*chidiff2 + Bt*Dt*chidiff);

  //Damping force
  area_bond = history[10]*history[11]*chidiff;
  damp_prefac = damp_bonded*area_bond;
  fdampx = -damp_prefac*(v[i][0] - v[j][0]);
  fdampy = -damp_prefac*(v[i][1] - v[j][1]);

  //Damping torque
  torquedamp = -damp_bonded*area_bond*area_bond*(omega[i][2]-omega[j][2]);

  //Update forces and torque
  f[i][0] += Fnmag*bex + Ftmag*mex + fdampx;
  f[i][1] += Fnmag*bey + Ftmag*mey + fdampy;

  torque[i][2] += Nn + Nt + torquedamp;

  //Update chi1, chi2
  hmin = MIN(atom->min_thickness[i], atom->min_thickness[j]);
  if (historyupdate) update_chi(Dn, Cn, Dt, Ct, hmin, history[8], history[9]);
}


void PairGranHopkins::update_chi(double Dn, double Cn, double Dt, double Ct, double hmin, double &chi1, double &chi2){
  double sig_n1 = Dn + Cn*chi1;
  double sig_s1 = Dt + Ct*chi1;
  double sig_n2 = Dn + Cn*chi2;
  double sig_s2 = Dt + Ct*chi2;

  double sig_c = sig_c0*pow(hmin,(2.0/3.0));
  double sig_t = sig_t0*sig_c;

  //Check for purely tensile/compressive failure
  if (sig_n1 > sig_t || sig_n2 > sig_t || sig_n1 < sig_c || sig_n2 < sig_c){
      if (sig_n1 > sig_t){
	  if (Cn != 0) chi1 = (sig_t - Dn)/Cn;
	  else {
	      chi1 = chi2;
	      return;
	  }
      }
      if (sig_n2 > sig_t) chi2 = (sig_t - Dn)/Cn; //if Cn==0, function would've returned above

      if (sig_n1 < sig_c){
	  if (Cn != 0) chi1 = (sig_c - Dn)/Cn;
	  else{
	      chi1 = chi2;
	      return;
	  }
      }
      if (sig_n2 < sig_c) chi2 = (sig_c - Dn)/Cn; //if Cn==0, function would've returned above
  }

  //Check for 'cohesion' shear failure at chi1
  if (sig_s1 > -tanphi*sig_n1 + tanphi*sig_t){
      if (Ct + tanphi*Cn != 0) chi1 = (tanphi*(sig_t-Dn)-Dt)/(Ct + tanphi*Cn);
      else if (Cn != 0 && Ct != 0){
	  //No need to treat case where Cn = Ct = 0, since sig_s = 0 for that case,
	  // and it would've been picked up above
	  //For the rare case of Ct = -tanphi*Cn, yield criterion is independent of chi,
	  // therefore bond breaks.
	  chi1 = chi2;
	  return;
      }
  }
  if (sig_s1 < tanphi*sig_n1 - tanphi*sig_t){
      if (tanphi*Cn - Ct != 0) chi1 = (Dt-tanphi*(Dn-sig_t))/(tanphi*Cn - Ct);
      else if (Cn != 0 && Ct != 0){
	  chi1 = chi2;
	  return;
      }
  }

  //Check for 'cohesion' shear failure at chi2
  if (sig_s2 > -tanphi*sig_n2 + tanphi*sig_t){
      if (Ct + tanphi*Cn != 0) chi2 = (tanphi*(sig_t-Dn)-Dt)/(Ct + tanphi*Cn);
      else if (Cn != 0 && Ct != 0){
	  //No need to treat case where Cn = Ct = 0, since sig_s = 0 for that case,
	  // and it would've been picked up above
	  //For the rare case of Ct = -tanphi*Cn, yield criterion is independent of chi,
	  // therefore bond breaks.
	  chi2 = chi1;
	  return;
      }
  }
  if (sig_s2 < tanphi*sig_n2 - tanphi*sig_t){
      if (tanphi*Cn - Ct != 0) chi1 = (Dt-tanphi*(Dn-sig_t))/(tanphi*Cn - Ct);
      else if (Cn != 0 && Ct != 0){
	  chi2 = chi1;
	  return;
      }
  }  
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairGranHopkins::settings(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Illegal pair_style command");

  Emod = force->numeric(FLERR, arg[0]);
  poiss = force->numeric(FLERR, arg[1]);
  sig_c0 = force->numeric(FLERR, arg[2]);
  sig_t0 = force->numeric(FLERR, arg[3]);
  phi = force->numeric(FLERR, arg[4]);
  damp_bonded = force->numeric(FLERR, arg[5]);
  friction_tangential = force->numeric(FLERR, arg[6]);
  damp_normal = force->numeric(FLERR, arg[7]);
  damp_tangential = force->numeric(FLERR, arg[8]);

  tanphi = tan(phi*MathConst::MY_PI/180.0);
}

/* ---------------------------------------------------------------------- */

double PairGranHopkins::single(int i, int j, int itype, int jtype,
			       double rsq,
			       double factor_coul, double factor_lj,
			       double &fforce)
{
  return 0.0;
}
