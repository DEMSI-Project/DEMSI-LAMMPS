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

PairStyle(gran/bpm,PairGranBPM)

#else

#ifndef LMP_PAIR_GRAN_BPM_H
#define LMP_PAIR_GRAN_BPM_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranBPM : public PairGranHookeHistory {
 public:

  PairGranBPM(class LAMMPS *);
  virtual ~PairGranBPM() {}

  virtual void compute(int, int) override;
  virtual void settings(int, char **) override;
  virtual double single(int, int, int, int, double, double, double, double &) override;
  virtual void transfer_history(double*, double*) override;
  virtual double init_one(int, int) override;
  void get_bond_info();

  std::vector<std::pair<tagint,tagint>> bondGlobalIDs;
  std::vector<std::vector<double>> bondContactHistory;
  std::vector<std::vector<double>> bondInfo;

 private:
  void compute_bonded(double* history, int* firsttouch, const int& i, const int& j, std::vector<double>& bondInfoIJ);
  void compute_nonbonded(double* history, int* firsttouch, const int& i, const int& j, std::vector<double>& bondInfoIJ);
  void rotate_tangential_vector(double& x, double& y, const double nx, const double ny);

  // history update flag
  int write_bonds = 0;

  // non-bonded
  double E_c;
  double kn_ks_ratio;
  double mu;
  double damp_normal;

  // bonded
  double lambda_bar;
  double E_c_bar;
  double normal_damping_ratio;
  double shear_damping_ratio;
  double moment_damping_ratio;
  double sigma_c_bar;
  double tau_c_bar_0;
  double kn_ks_ratio_bar;
  double bond_thickness;
  char failure_model[256];
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
