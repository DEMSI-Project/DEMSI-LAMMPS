/* ----------------------------------------------------------------------
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include "pair_gran_bpm.h"
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

#define EPSILON 1e-10
/* ---------------------------------------------------------------------- */

PairGranBPM::PairGranBPM(LAMMPS *lmp) : PairGranHookeHistory(lmp, 6) {
  nondefault_history_transfer = 1;
  beyond_contact = 1;
}

void PairGranBPM::compute(int eflag, int vflag) {
  if (eflag || vflag) 
    ev_setup(eflag,vflag);
  else 
    evflag = vflag_fdotr = 0;
  
  // loop over my elements
  for (int ii = 0; ii < list->inum; ii++) {
    int i = list->ilist[ii];
    int* jlist = list->firstneigh[i];
    double* allhistory = fix_history->firstvalue[i];
    int* firsttouch = fix_history->firstflag[i];

    // loop over neighbors of each element
    for (int jj = 0; jj < list->numneigh[i]; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK; // TODO: What does this do?

      double* history = &allhistory[size_history*jj];
      // 'history' now points to the ii-jj array that stores
      // all the history associated with pair ii-jj
      // For bonded pairs:
      // history[0]: Bonded flag (off or on)
      // history[1]: Fn_bar (normal bond force)
      // history[2]: Fs_bar_x (shear force in the x-dir)
      // history[3]: Fs_bar_y (shear force in the y-dir)
      // history[4]: Ms_bar_i (home bond moment)
      // history[5]: Ms_bar_j (neighbor bond moment)
      // history[6]: Bond thickness

      // For unbonded pairs:
      // history[0]: Bonded flag (off or on)
      // history[1,2]: normal force at previous time step, x and y components
      // history[3,4]: accumulated tangential displacement at contact, x and y
      // history[5] : delta_0 (initial overlap at bond breakage)

      std::vector<double> bondInfoIJ;
      if (history[0]) {
        compute_bonded(history, &firsttouch[jj], i, j, bondInfoIJ);
      } else { 
        compute_nonbonded(history, &firsttouch[jj], i, j, bondInfoIJ);
      }

      if (write_bonds) {
        // bondGlobalIDs
        std::pair<tagint, tagint> particle_pair{atom->tag[i], atom->tag[j]};
        bondGlobalIDs.push_back(particle_pair);

        // bondContactHistory
        std::vector<double> bondContactHistoryIJ;
        for (int m=0; m<size_history; m++)
          bondContactHistoryIJ.push_back(history[m]);
        bondContactHistory.push_back(bondContactHistoryIJ);

        // bondInfo
        bondInfo.push_back(bondInfoIJ);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

void PairGranBPM::compute_nonbonded(double* history, int* touch, int i, int j, std::vector<double>& bondInfoIJ) {
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  int historyupdate = 1;
  if (update->setupflag || write_bonds) 
    historyupdate = 0;

  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];

  double rsq = delx*delx + dely*dely;
  double radsum = radius[i] + radius[j];

  double fx, fy, fnx, fny, ftx, fty, torque_i, torque_j;
  if (rsq >= radsum*radsum) {
    if (historyupdate) {
      *touch = 0;
      for (int k = 0; k < size_history; k++) {
        history[k] = 0;
      }
    }
    delx = 0.; dely = 0.; 
    fx = 0.; fy = 0.;
    fnx = 0.; fny = 0.; 
    ftx = 0.; fty = 0.;
    torque_i = 0.; torque_j = 0.;

  } else {
    if (!*touch && historyupdate) {
      *touch = 1;
      for (int k = 0; k < size_history; k++) {
        history[k] = 0;
      }
    }
    double r = std::sqrt(rsq);
    double nx = delx/r;
    double ny = dely/r;

    double vrx = v[i][0] - v[j][0];
    double vry = v[i][1] - v[j][1];

    double delta = radsum - r - history[5];
    if (delta <= 0.) {
      if (historyupdate) {
        *touch = 0;
        for (int k = 0; k < size_history; k++) {
          history[k] = 0;
        }
      }
      delx = 0.; dely = 0.; 
      fx = 0.; fy = 0.;

      if (write_bonds) {
        bondInfoIJ.push_back(fx);
        bondInfoIJ.push_back(fy);
        bondInfoIJ.push_back(0.);
        bondInfoIJ.push_back(0.);
        bondInfoIJ.push_back(0.);
        bondInfoIJ.push_back(0.);
        bondInfoIJ.push_back(0.);
        bondInfoIJ.push_back(0.);
      }
      
      if (evflag) 
        ev_tally_xyz(i, j, atom->nlocal, force->newton_pair, 0.0, 0.0, fx, fy, 0., delx, dely, 0.);

      return;
    }

    // subtract to compute tangential component of relative translational velocity
    double vnnr = vrx*nx + vry*ny;
    double vtrx = vrx - nx*vnnr;
    double vtry = vry - ny*vnnr;

    // Add in the rotational component
    // TODO: Calculate the relative angular velocity at the contact point
    double wrz = radius[i]*omega[i][2] + radius[j]*omega[j][2];
    double vtx = vtrx + ny*wrz;
    double vty = vtry - nx*wrz;

    // TODO: Is mean thickness from the thickness distribution?
    double k_n_i = 2. * atom->mean_thickness[i] * E_c;
    double k_n_j = 2. * atom->mean_thickness[j] * E_c;
    kn = k_n_i * k_n_j / (k_n_i + k_n_j);

    // Elastic normal force
    // TODO: Create a critical damping?
    double fnmag = kn*delta;
    fnx = (fnmag - damp_normal*vnnr)*nx;
    fny = (fnmag - damp_normal*vnnr)*ny;

    // update tangential displacement, rotate if needed
    double disp_tx = history[3];
    double disp_ty = history[4];
    rotate_tangential_vector(disp_tx, disp_ty, nx, ny);
    disp_tx += vtx*update->dt;
    disp_ty += vty*update->dt;

    double dispmag = std::sqrt(disp_tx*disp_tx + disp_ty*disp_ty);

    // total tangential force
    double k_s_i = k_n_i/kn_ks_ratio;
    double k_s_j = k_n_j/kn_ks_ratio;
    kt = k_s_i * k_s_j / (k_s_i + k_s_j);

    ftx = -kt*disp_tx;
    fty = -kt*disp_ty;

    double ftmag = std::sqrt(ftx*ftx + fty*fty);
    double ftcrit = mu*std::fabs(fnmag);
    if (ftmag > ftcrit) {
      if (dispmag != 0.) {
        ftx *= ftcrit/ftmag;
        fty *= ftcrit/ftmag;
        // TODO: Do we also modify the tangential displacment?
        disp_tx = -ftx/kt;
        disp_ty = -fty/kt;
      } else { 
        ftx = 0.; fty = 0.;
      }
    }

    // Apply forces
    fx = fnx + ftx;
    fy = fny + fty;

    f[i][0] += fx;
    f[i][1] += fy;

    // Torque induced by tangential force
    torque_i = -radius[i]*(nx*fty - ny*ftx);
    torque[i][2] += torque_i;

    torque_j = -radius[j]*(nx*fty - ny*ftx);
    if (force->newton_pair || j < atom->nlocal) {
      f[j][0] -= fx;
      f[j][1] -= fy;
      torque[j][2] += torque_j;
    }

    if (historyupdate) {
      *touch = 1;
      history[1] = fnx; 
      history[2] = fny;
      history[3] = disp_tx;
      history[4] = disp_ty;
    }
  }

  if (write_bonds) {
    bondInfoIJ.push_back(fx);
    bondInfoIJ.push_back(fy);
    bondInfoIJ.push_back(fnx);
    bondInfoIJ.push_back(fny);
    bondInfoIJ.push_back(ftx);
    bondInfoIJ.push_back(fty);
    bondInfoIJ.push_back(torque_i);
    bondInfoIJ.push_back(torque_j);
  }

  if (evflag) 
    ev_tally_xyz(i, j, atom->nlocal, force->newton_pair, 0.0, 0.0, fx, fy, 0., delx, dely, 0.);
}

void PairGranBPM::compute_bonded(double *history, int* touch, int i, int j, std::vector<double>& bondInfoIJ) {
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;

  int historyupdate = 1;
  if (update->setupflag || write_bonds) 
    historyupdate = 0;

  // Get previous history
  double Fn_bar_old = history[1];
  double Fs_bar_x_old = history[2];
  double Fs_bar_y_old = history[3];
  double Ms_bar_old_i = history[4];
  double Ms_bar_old_j = history[5];

  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double r = std::sqrt(delx*delx + dely*dely);
  if (std::abs(r) < 1e-12)
    error->all(FLERR,"Particles have overlapped and centroids are coincident.");

  double nx = delx/r;
  double ny = dely/r;

  double vrx = v[i][0] - v[j][0];
  double vry = v[i][1] - v[j][1];

  double vr_dot_n = vrx*nx + vry*ny;

  // Add in the rotational component
  // TODO: Calculate the relative angular velocity at the contact point
  double wrz = radius[i]*omega[i][2] + radius[j]*omega[j][2];
  double vtx = vrx - nx*vr_dot_n + ny*wrz;
  double vty = vry - ny*vr_dot_n - nx*wrz;

  double r_bar = lambda_bar * std::min(radius[i], radius[j]);
  double area_bond = 2.*r_bar*history[6];

  kn = E_c_bar / (radius[i] + radius[j]);
  kt = kn/kn_ks_ratio_bar;
  double moi = 2./3. * std::pow(r_bar, 3) * history[6];

  // Update Fn_bar magnitude
  double Fn_bar = Fn_bar_old - kn*area_bond*vr_dot_n*update->dt;

  // Update x,y component of Fs_bar
  rotate_tangential_vector(Fs_bar_x_old, Fs_bar_y_old, nx, ny);
  double Fs_bar_x = Fs_bar_x_old - kt*area_bond*vtx*update->dt;
  double Fs_bar_y = Fs_bar_y_old - kt*area_bond*vty*update->dt;

  double delta_theta_s = (omega[j][2] - omega[i][2])*update->dt;
  double delta_Ms_bar = kn*moi*delta_theta_s;

  double Ms_bar_i = Ms_bar_old_i + delta_Ms_bar;
  double Ms_bar_j = Ms_bar_old_j - delta_Ms_bar;

  double sigma_bar_max = -Fn_bar/area_bond + std::fabs(Ms_bar_i)*r_bar/moi;
  double tau_bar_max = std::sqrt(Fs_bar_x*Fs_bar_x + Fs_bar_y*Fs_bar_y)/area_bond;

  double tau_c_bar;
  if (std::strcmp(failure_model, "mohrCoulomb") == 0)
    tau_c_bar = tau_c_bar_0*(1. + (Fn_bar/area_bond)/sigma_c_bar);
  else if (std::strcmp(failure_model, "constant") == 0)
    tau_c_bar = tau_c_bar_0;
  else
    error->all(FLERR,"Unknown shearBreakingStressType");

  bool failure = false;
  if (sigma_bar_max >= sigma_c_bar || tau_bar_max >= tau_c_bar) {
    failure = true;
    // TODO: Do we set all forces and tourque to zero now?
    Fn_bar = 0.; Fs_bar_x = 0.; Fs_bar_y = 0.; Ms_bar_i = 0.; Ms_bar_j = 0.;
  }

  // Local bond damping
  double Fn_bar_damp, Fs_bar_damp_x, Fs_bar_damp_y, Ms_bar_damp; 
  if (!failure) {
    // TODO: Can I use atom->rmass?
    // TODO: Is this the actual mass? Ice area from voronoi diagram?
    double m1 = 900.*M_PI*radius[i]*radius[i]*atom->mean_thickness[i];
    double m2 = 900.*M_PI*radius[j]*radius[j]*atom->mean_thickness[j];
    double m_eq = m1*m2/(m1 + m2);

    Fn_bar_damp = vr_dot_n*normal_damping_ratio*2.*std::sqrt(kn*area_bond*m_eq);

    Fs_bar_damp_x = vtx*shear_damping_ratio*2.*std::sqrt(kt*area_bond*m_eq);
    Fs_bar_damp_y = vty*shear_damping_ratio*2.*std::sqrt(kt*area_bond*m_eq);

    double moi_eq = atom->momentOfInertia[i]*atom->momentOfInertia[j]/(atom->momentOfInertia[i] + atom->momentOfInertia[j]);
    Ms_bar_damp = (omega[j][2] - omega[i][2])*moment_damping_ratio*2.*std::sqrt(kn*moi*moi_eq);
  } else {
    Fn_bar_damp = 0.;
    Fs_bar_damp_x = 0.;
    Fs_bar_damp_y = 0.;
    Ms_bar_damp = 0.;
  }

  double fnx = (Fn_bar - Fn_bar_damp)*nx;
  double fny = (Fn_bar - Fn_bar_damp)*ny;
  double ftx = Fs_bar_x - Fs_bar_damp_x; 
  double fty = Fs_bar_y - Fs_bar_damp_y; 
  double tor_i = Ms_bar_i + Ms_bar_damp; // TODO: check the sign on this
  double tor_j = Ms_bar_j - Ms_bar_damp; // TODO: check the sign on this

  double fx = fnx + ftx;
  double fy = fny + fty;

  f[i][0] += fx;
  f[i][1] += fy;
  torque[i][2] += tor_i;

  if (force->newton_pair || j < atom->nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;

    // TODO: check the sign on this
    torque[j][2] += tor_j;
  }

  if (write_bonds) {
    bondInfoIJ.push_back(fx);
    bondInfoIJ.push_back(fy);
    bondInfoIJ.push_back(fnx);
    bondInfoIJ.push_back(fny);
    bondInfoIJ.push_back(ftx);
    bondInfoIJ.push_back(fty);
    bondInfoIJ.push_back(tor_i);
    bondInfoIJ.push_back(tor_j);
  }

  if (historyupdate) {
    *touch = 1;
    history[1] = Fn_bar;
    history[2] = Fs_bar_x;
    history[3] = Fs_bar_y;
    history[4] = Ms_bar_i;
    history[5] = Ms_bar_j;

    if (failure) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double rij = sqrt(dx*dx + dy*dy);
      double delta_0 = atom->radius[i] + atom->radius[j] - rij;

      if (delta_0 < 0) 
        delta_0 = 0.;

      for (int k = 0; k < size_history; ++k)
        history[k] = 0;

      history[5] = delta_0;
    }
  }

  if (evflag) 
    ev_tally_xyz(i,j,atom->nlocal, force->newton_pair, 0.0, 0.0, fx, fy, 0, x[i][0]-x[j][0], x[i][1]-x[j][1], 0);
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairGranBPM::settings(int narg, char **arg) {
  if (narg != 11) error->all(FLERR,"Illegal pair_style command");

  double youngs = utils::numeric(FLERR,arg[0],false,lmp); // elasticModulus
  double nu = utils::numeric(FLERR,arg[1],false,lmp); // poissonRatio

  // Grain
  E_c = youngs;
  mu = utils::numeric(FLERR,arg[2],false,lmp); // frictionCoefficient
  damp_normal = utils::numeric(FLERR,arg[3],false,lmp); // nonbondedNormalDamping
  kn_ks_ratio = 2.*(1. + nu);

  // Bond
  E_c_bar = youngs;
  lambda_bar = utils::numeric(FLERR,arg[4],false,lmp); // lambdaBar
  sigma_c_bar = utils::numeric(FLERR,arg[5],false,lmp); // tensileBreakingStress
  tau_c_bar_0 = utils::numeric(FLERR,arg[6],false,lmp); // shearBreakingStress
  normal_damping_ratio = utils::numeric(FLERR,arg[7],false,lmp); // normalDampingRatio
  shear_damping_ratio = utils::numeric(FLERR,arg[8],false,lmp); // shearDampingRatio
  moment_damping_ratio = utils::numeric(FLERR,arg[9],false,lmp); // momentDampingRatio
  std::strcpy(failure_model,arg[10]); // failureModel ("constant" or "mohrCoulomb")
  kn_ks_ratio_bar = 2.*(1. + nu);
}

/* ---------------------------------------------------------------------- */
double PairGranBPM::init_one(int i, int j) {
  double cutoff;
  cutoff = PairGranHookeHistory::init_one(i, j);
  cutoff += maxrad_dynamic[i]*0.1; // This could be an input parameter?
  return cutoff;
}

/* ---------------------------------------------------------------------- */

double PairGranBPM::single(int i, int j, int itype, int jtype, double rsq, double factor_coul, double factor_lj, double &fforce) {
  return 0.0;
}

/* ---------------------------------------------------------------------- */
void PairGranBPM::transfer_history(double* sourcevalues, double* targetvalues) {
  // history[0]: Bonded flag (off or on)
  // history[1]: Fn_bar (normal bond force)
  // history[2]: Fs_bar_x (shear force in the x-dir)
  // history[3]: Fs_bar_y (shear force in the y-dir)
  // history[4]: Ms_bar_i (previous home bond moment)
  // history[5]: Ms_bar_j (previous neighbor bond moment)
  // history[6]: Bond thickness

  // For unbonded pairs:
  // history[0]: Bonded flag (off or on)
  // history[1,2]: normal force at previous time step, x and y components
  // history[3,4]: accumulated tangential displacement at contact, x and y
  // history[5] : delta_0 (initial overlap at bond breakage)

  if (sourcevalues[0]) {
    // Transfers for bonded particles
    targetvalues[0] = sourcevalues[0];
    targetvalues[1] = sourcevalues[1];
    targetvalues[2] = -sourcevalues[2];
    targetvalues[3] = -sourcevalues[3];
    targetvalues[4] = sourcevalues[5];
    targetvalues[5] = sourcevalues[4];
    targetvalues[6] = sourcevalues[6];
  } else {
    // Transfers for non-bonded particles
    targetvalues[0] = sourcevalues[0];
    for (int i = 1; i < 5; i++) {
      targetvalues[i] = -sourcevalues[i];
    }
    targetvalues[5] = sourcevalues[5];
  }
}

void PairGranBPM::rotate_tangential_vector(double& x, double& y, const double nx, const double ny) {
  // Tangential displacement in normal direction
  double ndisp = nx*x + ny*y;

  if (std::fabs(ndisp) > EPSILON) {
    double dispmag = std::sqrt(x*x + y*y);
    // Remove component along normal direction
    x -= ndisp*nx;
    y -= ndisp*ny;

    // Rescale to preserve magnitude
    double prjmag = std::sqrt(x*x + y*y);
    double scalefac = 0.0;
    if (prjmag > 0.0) {
      scalefac = dispmag/prjmag;
    }
    x *= scalefac;
    y *= scalefac;
  }
}

void PairGranBPM::get_bond_info() {
  bondGlobalIDs.clear();
  bondContactHistory.clear();
  bondInfo.clear();
  
  write_bonds = 1;
  compute(0,0);
  write_bonds = 0;
}
