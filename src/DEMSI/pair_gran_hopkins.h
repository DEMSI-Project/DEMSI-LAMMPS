/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gran/hopkins,PairGranHopkins)

#else

#ifndef LMP_PAIR_GRAN_HOPKINS_H
#define LMP_PAIR_GRAN_HOPKINS_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHopkins : public PairGranHookeHistory {
public:
  PairGranHopkins(class LAMMPS *);
  virtual ~PairGranHopkins();
  virtual void compute(int, int);
  void settings(int, char **);
  double single(int, int, int, int, double, double, double, double &);
  virtual void transfer_history(double*, double*);
  double init_one(int, int);
protected:
  void compute_bonded(double*, int*, int, int);
  void compute_nonbonded(double*, int*, int, int);
  void update_chi(double, double, double, double, double, double, double, double&, double&);

  int history_ndim;
  double Emod;
  double poiss;
  double sig_c0;
  double sig_t0;
  double phi;
  double damp_bonded;
  double damp_tangential;
  double damp_normal;
  double tanphi;
  double friction_tangential;
  double hprime_0;
  double plasticFrictionCoeff;
  double plasticHardeningCoeff;
  double exponentialIceStrengthCoeff;
  double Gmod;
  char sig_c0_type[256];
  char sig_t0_type[256];
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
