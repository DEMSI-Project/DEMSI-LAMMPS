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

PairStyle(gran/hooke/thickness/kk,PairGranHookeThicknessKokkos<LMPDeviceType>)
PairStyle(gran/hooke/thickness/kk/device,PairGranHookeThicknessKokkos<LMPDeviceType>)
PairStyle(gran/hooke/thickness/kk/host,PairGranHookeThicknessKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_THICKNESS_KOKKOS_H
#define LMP_PAIR_GRAN_HOOKE_THICKNESS_KOKKOS_H

#include "pair_gran_hooke_thickness.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixNeighHistoryKokkos;
  
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairGranHookeThicknessCompute {};

template <class DeviceType>
class PairGranHookeThicknessKokkos : public PairGranHookeThickness {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairGranHookeThicknessKokkos(class LAMMPS *);
  virtual ~PairGranHookeThicknessKokkos();
  virtual void compute(int, int);
  void init_style();

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeThicknessCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeThicknessCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int) const;

  template<int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
		    F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
		    X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;
  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz_atom(EV_FLOAT &ev, int i, int j,
			 F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
			 X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;
    
 protected:
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_v_array_randomread v;
  typename AT::t_v_array_randomread omega;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_float_1d_randomread rmass;
  typename AT::t_float_1d_randomread radius;
  typename AT::t_float_1d_randomread mean_thickness;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;
  typename AT::t_tagint_1d tag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int newton_pair;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;

  friend void pair_virial_fdotr_compute<PairGranHookeThicknessKokkos>(PairGranHookeThicknessKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair granular requires atom attributes radius, rmass

The atom style defined does not have these attributes.

E: Pair granular requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Could not find pair fix neigh history ID

UNDOCUMENTED

U: Pair granular with shear history requires newton pair off

This is a current restriction of the implementation of pair
granular styles with history.

U: Could not find pair fix ID

A fix is created internally by the pair style to store shear
history information.  You cannot delete it.

*/
