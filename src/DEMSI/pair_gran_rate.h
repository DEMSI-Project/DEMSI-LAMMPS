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

#ifdef PAIR_CLASS

PairStyle(gran/rate,PairGranRate)

#else

#ifndef LMP_PAIR_GRAN_RATE_H
#define LMP_PAIR_GRAN_RATE_H

#include "pair_gran_hooke_history.h"
#include <cmath>

namespace LAMMPS_NS {

class PairGranRate : public PairGranHookeHistory {
public:
  PairGranRate(class LAMMPS *);
  virtual ~PairGranRate();
  virtual void compute(int, int);
  void settings(int, char **);
  double single(int, int, int, int, double, double, double, double &);
  virtual void transfer_history(double*, double*);
  double init_one(int, int);
protected:
  void compute_TV_Transport(int);
  void compute_rate_explicit(int);
  void compute_rate_implicit(int);
  void load_new_forces(int);

  void allocate();

  int history_ndim;
  int icount;

  double  bulkModulus;
  double shearModulus;
  double       shearYieldStress;
  double compressiveYieldStress;
  double  tensileFractureStress;

  int timeIntegrationFlag;
  double bulkCFL;
  bool UseTVTransport;
  bool UseGradO2;
  double contactMargin;

  double ultimateTensileStrain;
  double elastic_wavespeed;

  const double OverLap = 1.050075135808664;   // R^\prime/R = \sqrt(2\sqrt(3) \slash \pi)
  const double LapOver = 0.9523128068639574;  // R/R^\prime = \sqrt(\pi \slash 2\sqrt(3))
  const double rho0 = 900.;  // kg/m^3
  const double MY_PI3 = 1.0471975511965976; // pi/3
  const double PIby2 = asin(1.);  // (1/2) pi
  const double PI3by2 = 3.*PIby2; // (3/2) pi
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

 */
