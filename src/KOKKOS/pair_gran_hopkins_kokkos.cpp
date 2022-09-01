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

//For bonded pairs:
// history[0,1]: chi1, chi2
// history[2-9]: x,y components of s1i, s2i, s1j, s2j
// history[10]: bond length
// history[11]: bond thickness, h

//For unbonded pairs:
// history[0,1]: chi1, chi2
// history[2,3]: accumulated tangential displacement at contact, x and y
// history[4]  : delta_0: initial overlap at bond break
// history[5]  : previousForceConv
// history[6]  : ridgeSlipConv
// history[7]  : ridgeSlipUsedConv
// history[8]  :
// history[9]  :
// history[10] :

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include "pair_gran_hopkins_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "update.h"
#include "force.h"
#include "fix_neigh_history_kokkos.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-10

#define DEBUGID_1 35
#define DEBUGID_2 55
#define DEBUG_TIMESTEP 32696
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHopkinsKokkos<DeviceType>::PairGranHopkinsKokkos(LAMMPS *lmp) : PairGranHopkins(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MASK_MASK | 
                  ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK | THICKNESS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK | THICKNESS_MASK;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
PairGranHopkinsKokkos<DeviceType>::~PairGranHopkinsKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairGranHopkinsKokkos<DeviceType>::init_style()
{
  if (history && fix_history == nullptr) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",12);
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY_HH";
    fixarg[1] = (char *) "all";
    if (execution_space == Device)
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/DEVICE";
    else
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/HOST";
    fixarg[3] = dnumstr;
    modify->replace_fix("NEIGH_HISTORY_HH_DUMMY",4,fixarg,1);
    delete [] fixarg;
    int ifix = modify->find_fix("NEIGH_HISTORY_HH");
    fix_history = (FixNeighHistory *) modify->fix[ifix];
    fix_history->pair = this;
    fix_historyKK = (FixNeighHistoryKokkos<DeviceType> *)fix_history;
  }

  // PairGranHopkins uses PairGranHookeHistory::init_style
  PairGranHookeHistory::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with gran/hopkins/kk");
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairGranHopkinsKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  // reallocate per-atom arrays if necessary
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  min_thickness = atomKK->k_min_thickness.view<DeviceType>();
  mean_thickness = atomKK->k_mean_thickness.view<DeviceType>();
  iceConcentration = atomKK->k_iceConcentration.view<DeviceType>();
  ridgingIceThickness = atomKK->k_ridgingIceThickness.view<DeviceType>();
  ridgingIceThicknessWeight = atomKK->k_ridgingIceThicknessWeight.view<DeviceType>();
  netToGrossClosingRatio = atomKK->k_netToGrossClosingRatio.view<DeviceType>();
  changeEffectiveElementArea = atomKK->k_changeEffectiveElementArea.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  d_firsttouch = fix_historyKK->d_firstflag;
  d_firsthistory = fix_historyKK->d_firstvalue;

  EV_FLOAT ev;

  if (strcmp(sig_c0_type,"constant") == 0) {
    strcmp_sig_c0_type_constant = true;
    strcmp_sig_c0_type_KovacsSodhi = false;
  } else if (strcmp(sig_c0_type,"KovacsSodhi") == 0) {
    strcmp_sig_c0_type_constant = false;
    strcmp_sig_c0_type_KovacsSodhi = true;
  } else {
    error->all(FLERR,"Unknown sig_c0_type");
  }
  if (strcmp(sig_t0_type,"constant") == 0) {
    strcmp_sig_t0_type_constant = true;
    strcmp_sig_t0_type_multiply_sig_c0 = false;
  } else if (strcmp(sig_t0_type,"multiply_sig_c0") == 0) {
    strcmp_sig_t0_type_constant = false;
    strcmp_sig_t0_type_multiply_sig_c0 = true;
  } else {
    error->all(FLERR,"Unknown sig_t0_type");
  }
  update_dt = update->dt;


 if (lmp->kokkos->neighflag == HALF) {
    if (force->newton_pair) {
       if (historyupdate) {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,1,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,1,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,1,0>>(0,inum),*this, ev);
          }
       } else {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,0,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,0,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,0,0>>(0,inum),*this, ev);
          }
       }
    } else {
       if (historyupdate) {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,1,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,1,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,1,0>>(0,inum),*this, ev);
          }
       } else {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,0,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,0,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,0,0>>(0,inum),*this, ev);
          }
       }
    }
 }
 else if (lmp->kokkos->neighflag == HALFTHREAD) {
    if (force->newton_pair) {
      if (historyupdate) {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,1,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,1,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,1,0>>(0,inum),*this, ev);
          }
       } else {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,0,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,0,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,0,0>>(0,inum),*this, ev);
          }
       }
    } else {
       if (historyupdate) {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,1,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,1,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,1,0>>(0,inum),*this, ev);
          }
       } else {
          if (vflag_atom) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,0,2>>(0,inum),*this, ev);
          } else  if (vflag_global) {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,0,1>>(0,inum),*this, ev);
          } else {
             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,0,0>>(0,inum),*this, ev);
          }
       }
    }
 }

   if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::operator()(TagPairGranHopkinsCompute<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE,EVFLAG>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
//  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
//  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;
    Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
    Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  int i = d_ilist[ii];
  int itype = type[i];
  int jnum = d_numneigh[i];

  F_FLOAT fx_sum = 0.0;
  F_FLOAT fy_sum = 0.0;
  F_FLOAT torque_i_sum = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    int jtype = type[j];

    F_FLOAT fx = 0.0;
    F_FLOAT fy = 0.0;
    F_FLOAT torque_i = 0.0;
    F_FLOAT torque_j = 0.0;

    F_FLOAT fnx,fny;
    F_FLOAT ftx,fty;

    // calculate the force for a single bond
    compute_single_bond<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE>(i,j,jj,
							     fx,fy,
							     fnx,fny,
							     ftx,fty,
							     torque_i,torque_j,
							     true);

    fx_sum += fx;
    fy_sum += fy;
    torque_i_sum += torque_i;

    if (NEWTON_PAIR || j < nlocal){
       a_f(j,0) -= fx;
       a_f(j,1) -= fy;
       a_torque(j,2) += torque_j;
    }

    F_FLOAT delx = x(i,0) - x(j,0);
    F_FLOAT dely = x(i,1) - x(j,1);
    if (EVFLAG == 2)
      ev_tally_xyz_atom<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, fx_sum, fy_sum, 0, delx, dely, 0);
    if (EVFLAG == 1)
      ev_tally_xyz<NEWTON_PAIR>(ev, i, j, fx_sum, fy_sum, 0, delx, dely, 0);
  }

  a_f(i,0) += fx_sum;
  a_f(i,1) += fy_sum;

  // torque induced by tangential force
  a_torque(i,2) += torque_i_sum;

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::demsi_single_bond(int i,
						    int j,
						    int jj,
						    F_FLOAT &fx,
						    F_FLOAT &fy,
						    F_FLOAT &fnx,
						    F_FLOAT &fny,
						    F_FLOAT &ftx,
						    F_FLOAT &fty,
						    F_FLOAT &torque_i,
						    F_FLOAT &torque_j) {

  if (lmp->kokkos->neighflag == HALF) {
    if (force->newton_pair) {
      compute_single_bond<HALF,1,0>(i,j,jj,
				    fx,fy,
				    fnx,fny,
				    ftx,fty,
				    torque_i,torque_j,
				    false);
    } else {
      compute_single_bond<HALF,0,0>(i,j,jj,
				    fx,fy,
				    fnx,fny,
				    ftx,fty,
				    torque_i,torque_j,
				    false);
    }
  } else if (lmp->kokkos->neighflag == HALFTHREAD) {
    if (force->newton_pair) {
      compute_single_bond<HALFTHREAD,1,0>(i,j,jj,
					  fx,fy,
					  fnx,fny,
					  ftx,fty,
					  torque_i,torque_j,
					  false);
    } else {
      compute_single_bond<HALFTHREAD,0,0>(i,j,jj,
					  fx,fy,
					  fnx,fny,
					  ftx,fty,
					  torque_i,torque_j,
					  false);
    }
  }

}

//-----------------------------------------------------------------------------

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::compute_single_bond(int i,
							    int j,
							    int jj,
							    F_FLOAT &fx,
							    F_FLOAT &fy,
							    F_FLOAT &fnx,
							    F_FLOAT &fny,
							    F_FLOAT &ftx,
							    F_FLOAT &fty,
							    F_FLOAT &torque_i,
							    F_FLOAT &torque_j,
							    bool modifyState) const {

  F_FLOAT chi1 = d_firsthistory(i,size_history*jj);
  F_FLOAT chi2 = d_firsthistory(i,size_history*jj+1);
  if (chi1 >= chi2) {
    // Un-bonded, chi1 >= chi2
    compute_nonbonded_kokkos<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE>(i,j,jj,
								  fx,fy,
								  fnx,fny,
								  ftx,fty,
								  torque_i,torque_j,
								  modifyState);
  } else {
    // Bonded
    compute_bonded_kokkos<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE>(i,j,jj,
							       fx,fy,
							       fnx,fny,
							       ftx,fty,
							       torque_i,torque_j,
							       modifyState);
  }

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::elastic_stiffness(F_FLOAT meanIceThickness1,
							  F_FLOAT meanIceThickness2,
							  F_FLOAT radius1,
							  F_FLOAT radius2,
							  F_FLOAT mass1,
							  F_FLOAT mass2,
							  F_FLOAT bondLength,
							  F_FLOAT &elasticStiffness,
							  F_FLOAT &elasticDamping) const {

  F_FLOAT stiffness1 = (Emod * meanIceThickness1) / (2.0 * radius1);
  F_FLOAT stiffness2 = (Emod * meanIceThickness2) / (2.0 * radius2);
  elasticStiffness = MIN(stiffness1,stiffness2);

  elasticDamping = std::sqrt(2.0 * elasticStiffness * bondLength * MIN(mass1,mass2)) / bondLength;
  elasticDamping *= damp_normal;

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::plastic_parameters(F_FLOAT particleRadius,
							   F_FLOAT plasticFrictionCoeff,
							   F_FLOAT plasticHardeningCoeff,
							   F_FLOAT exponentialIceStrengthCoeff,
							   F_FLOAT iceConcentration1,
							   F_FLOAT iceConcentration2,
							   F_FLOAT ridgingIceThickness1,
							   F_FLOAT ridgingIceThickness2,
							   F_FLOAT ridgingIceThicknessWeight1,
							   F_FLOAT ridgingIceThicknessWeight2,
							   F_FLOAT radius1,
							   F_FLOAT radius2,
							   F_FLOAT &plasticFriction,
							   F_FLOAT &plasticHardeningStiffness) const {

  F_FLOAT ridgingThickness =
    (ridgingIceThickness1       + ridgingIceThickness2      ) /
    (ridgingIceThicknessWeight1 + ridgingIceThicknessWeight2);

  F_FLOAT resolutionScaling = 10000.0 / (2.0*particleRadius);

  F_FLOAT iceConcentration = MIN(iceConcentration1, iceConcentration2);
  F_FLOAT iceConcFactor = std::exp(-exponentialIceStrengthCoeff * (1.0 - iceConcentration));
  plasticFriction           = plasticFrictionCoeff  * ridgingThickness * iceConcFactor;
  plasticHardeningStiffness = plasticHardeningCoeff * ridgingThickness*ridgingThickness * resolutionScaling;

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::elastic_plastic_model(F_FLOAT bondLength,
							      F_FLOAT previousForce,
							      F_FLOAT overlap,
							      F_FLOAT convergence,
							      F_FLOAT elasticStiffness,
							      F_FLOAT elasticDamping,
							      F_FLOAT plasticFriction,
							      F_FLOAT plasticHardeningStiffness,
							      F_FLOAT &ridgingForce,
							      F_FLOAT &elasticOverlap,
							      F_FLOAT &plasticOverlap,
							      F_FLOAT &elasticConvergence,
							      F_FLOAT &plasticConvergence)  const{

  // pressure ridging
  F_FLOAT A = 1.0 / (plasticHardeningStiffness * update_dt);
  F_FLOAT B = 1.0;
  F_FLOAT C = (elasticStiffness * bondLength) / elasticDamping;
  F_FLOAT D = (elasticStiffness * plasticFriction * bondLength) / (elasticDamping * plasticHardeningStiffness);
  F_FLOAT denominator = (A + (1.0 + (elasticStiffness / plasticHardeningStiffness)) / elasticDamping);

  ridgingForce = (A * previousForce + B * convergence * bondLength + C * overlap + D) / denominator;

  plasticOverlap = (ridgingForce - plasticFriction * bondLength) / (plasticHardeningStiffness * bondLength);
  elasticOverlap = overlap - plasticOverlap;

  elasticConvergence = (ridgingForce - elasticStiffness * bondLength * elasticOverlap) / (elasticDamping * bondLength);
  plasticConvergence = convergence - elasticConvergence;

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::elastic_model(F_FLOAT bondLength,
						      F_FLOAT elasticOverlap,
						      F_FLOAT elasticConvergence,
						      F_FLOAT elasticStiffness,
						      F_FLOAT elasticDamping,
						      F_FLOAT &elasticForce) const {

  elasticForce = (elasticStiffness * elasticOverlap + elasticDamping * elasticConvergence) * bondLength;

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::geometry_change(bool modifyOtherElement,
							F_FLOAT bondLength,
							F_FLOAT netToGrossClosingRatio1,
							F_FLOAT netToGrossClosingRatio2,
							F_FLOAT ridgeSlip,
							F_FLOAT &ridgeSlipUsed,
							F_FLOAT &changeEffectiveElementArea1,
							F_FLOAT &changeEffectiveElementArea2) const {

  if (ridgeSlip > ridgeSlipUsed) {

    F_FLOAT weight1 = netToGrossClosingRatio1 / (netToGrossClosingRatio1 + netToGrossClosingRatio2);
    F_FLOAT weight2 = netToGrossClosingRatio2 / (netToGrossClosingRatio1 + netToGrossClosingRatio2);

    F_FLOAT areaDecrease = bondLength * (ridgeSlip - ridgeSlipUsed);

    changeEffectiveElementArea1 = changeEffectiveElementArea1 - weight1 * areaDecrease;
    if (modifyOtherElement) {
      changeEffectiveElementArea2 = changeEffectiveElementArea2 - weight2 * areaDecrease;
    } // modifyOtherElement

    ridgeSlipUsed = ridgeSlip;

  }

}

//-----------------------------------------------------------------------------

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::hopkins_ridging_model(bool modifyOtherElement,
							      bool updateGeometry,
							      F_FLOAT overlap,
							      F_FLOAT convergence,
							      F_FLOAT iceConcentration1,
							      F_FLOAT iceConcentration2,
							      F_FLOAT meanIceThickness1,
							      F_FLOAT meanIceThickness2,
							      F_FLOAT radius1,
							      F_FLOAT radius2,
							      F_FLOAT mass1,
							      F_FLOAT mass2,
							      F_FLOAT ridgingIceThickness1,
							      F_FLOAT ridgingIceThickness2,
							      F_FLOAT ridgingIceThicknessWeight1,
							      F_FLOAT ridgingIceThicknessWeight2,
							      F_FLOAT netToGrossClosingRatio1,
							      F_FLOAT netToGrossClosingRatio2,
							      F_FLOAT &changeEffectiveElementArea1,
							      F_FLOAT &changeEffectiveElementArea2,
							      F_FLOAT particleRadius,
							      F_FLOAT plasticFrictionCoeff,
							      F_FLOAT plasticHardeningCoeff,
							      F_FLOAT exponentialIceStrengthCoeff,
							      F_FLOAT bondLength,
							      F_FLOAT &ridgeSlip,
							      F_FLOAT &ridgeSlipUsed,
							      F_FLOAT &previousForce,
							      F_FLOAT &contactForce) const {

  // overlaps
  F_FLOAT totalOverlap = overlap;
  F_FLOAT elasticOverlap = overlap - ridgeSlip;
  F_FLOAT plasticOverlap = ridgeSlip;

  // only have contact force if overlap after plastic slip taken into account
  if (elasticOverlap > 0.0) {

    // elastic stiffness
    F_FLOAT elasticStiffness;
    F_FLOAT elasticDamping;
    elastic_stiffness(meanIceThickness1,
		      meanIceThickness2,
		      radius1,
		      radius2,
		      mass1,
		      mass2,
		      bondLength,
		      elasticStiffness,
		      elasticDamping);

    // plastic parameters
    F_FLOAT plasticFriction;
    F_FLOAT plasticHardeningStiffness;
    plastic_parameters(particleRadius,
		       plasticFrictionCoeff,
		       plasticHardeningCoeff,
		       exponentialIceStrengthCoeff,
		       iceConcentration1,
		       iceConcentration2,
		       ridgingIceThickness1,
		       ridgingIceThickness2,
		       ridgingIceThicknessWeight1,
		       ridgingIceThicknessWeight2,
		       radius1,
		       radius2,
		       plasticFriction,
		       plasticHardeningStiffness);

    // ridging force
    F_FLOAT ridgingForce;
    F_FLOAT elasticOverlapRidging;
    F_FLOAT plasticOverlapRidging;
    F_FLOAT elasticConvergenceRidging;
    F_FLOAT plasticConvergenceRidging;
    elastic_plastic_model(bondLength,
			  previousForce,
			  totalOverlap,
			  convergence,
			  elasticStiffness,
			  elasticDamping,
			  plasticFriction,
			  plasticHardeningStiffness,
			  ridgingForce,
			  elasticOverlapRidging,
			  plasticOverlapRidging,
			  elasticConvergenceRidging,
			  plasticConvergenceRidging);

    // previous force for ridging
    previousForce = ridgingForce;

    // elastic force
    F_FLOAT elasticForce;
    elastic_model(bondLength,
		  elasticOverlap,
		  convergence,
		  elasticStiffness,
		  elasticDamping,
		  elasticForce);

    // use elastic/plastic force if elastic/plastic force is less than elastic
    // force and ridging is positive and converging
    if (ridgingForce < elasticForce and
	plasticOverlapRidging > 0.0 and
	plasticConvergenceRidging > 0.0) {

      // use ridging force
      ridgeSlip = plasticOverlapRidging;

      contactForce = ridgingForce;

    } else {

      // use elastic force
      contactForce = elasticForce;

    }

  } else {

    // no element overlap so no interaction
    contactForce = 0.0;
    previousForce = 0.0;

  }

  // change the element geometry
  if (updateGeometry) {
    geometry_change(modifyOtherElement,
		    bondLength,
		    netToGrossClosingRatio1,
		    netToGrossClosingRatio2,
		    ridgeSlip,
		    ridgeSlipUsed,
		    changeEffectiveElementArea1,
		    changeEffectiveElementArea2);
  } // updateGeometry

}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::compute_nonbonded_kokkos(int i,
								 int j,
								 int jj,
								 F_FLOAT &fx,
								 F_FLOAT &fy,
								 F_FLOAT &fnx,
								 F_FLOAT &fny,
								 F_FLOAT &ftx,
								 F_FLOAT &fty,
								 F_FLOAT &torque_i,
								 F_FLOAT &torque_j,
								 bool modifyState) const {

   F_FLOAT r, rinv, nx, ny, radmin;
   F_FLOAT vnnr, vnx, vny;
   F_FLOAT wrz, vtrx, vtry, vtx, vty, vrel;
   F_FLOAT delta, delta_dot;

   F_FLOAT sig_c = 0, hmin;
   F_FLOAT hprime, kp, kr, ke, L, kt0;
   F_FLOAT num, denom, fnmag_plastic, fnmag_elastic, fnmag;
   F_FLOAT ncrossF;

   F_FLOAT ndisp, dispmag, scalefac;
   F_FLOAT ftmag, ftcrit;
   F_FLOAT var1, var2, disp_tx, disp_ty, prjmag;

   X_FLOAT delx = x(i,0) - x(j,0);
   X_FLOAT dely = x(i,1) - x(j,1);
   X_FLOAT rsq = delx*delx + dely*dely;
   F_FLOAT radsum = radius[i] + radius[j];

   if (rsq >= radsum*radsum){
     if (modifyState) {
       d_firsttouch(i,jj) = 0;
       for (int k = 0; k < size_history; k++) {
	 d_firsthistory(i,size_history*jj+k) = 0;
       }
     }
     fx = fy = 0;
     fnx = fny = 0;
     ftx = fty = 0;
     torque_i = torque_j = 0;
   }
   else{
     if (modifyState) {
       if (!d_firsttouch(i,jj)){ //If this is first contact
	 d_firsttouch(i,jj) = 1;
	 for (int k = 2; k < size_history; k++) {
	   d_firsthistory(i,size_history*jj+k) = 0;
	 }
       }
     }

     r = sqrt(rsq);
     rinv = 1.0/r;
     nx = delx/r;
     ny = dely/r;

     radmin = MIN(radius[i],radius[j]);

     //L = 2*radmin*(1+(abs(radius[i] - radius[j])/r));
     L = (4.0 * radius[i] * radius[j]) / (radius[i] + radius[j]);

     // relative translational velocity
     V_FLOAT vrx = v(i,0) - v(j,0);
     V_FLOAT vry = v(i,1) - v(j,1);

     delta = radsum - r - d_firsthistory(i,size_history*jj+4);
     if (delta < 0) return ; //delta = 0;

     // Compute tangential force
     // normal component of relative translational velocity
     vnnr = vrx*nx + vry*ny;
     vnx = nx*vnnr;
     vny = ny*vnnr;

     // subtract to compute tangential component of relative translational velocity
     vtrx = vrx - vnx;
     vtry = vry - vny;

     // total relative tangential velocities at contact
     wrz = radius[i]*omega(i,2) + radius[j]*omega(j,2);
     vtx = vtrx + ny*wrz;
     vty = vtry - nx*wrz;

     vrel = vtx*vtx + vty*vty;
     vrel = sqrt(vrel);

     delta_dot = -vnnr;

     F_FLOAT contactForce;

     if (type(i) == 2 || type(j) == 2){
       F_FLOAT elasticStiffness;
       F_FLOAT elasticDamping;
       int iceIndex;
       if (type(i) == 1) iceIndex = i;
       else iceIndex = j;
       elastic_stiffness(mean_thickness(iceIndex),
			 mean_thickness(iceIndex),
			 radius(i),
			 radius(j),
			 rmass(i),
			 rmass(j),
			 L,
			 elasticStiffness,
			 elasticDamping);
       contactForce = (elasticStiffness*delta + elasticDamping*delta_dot) * L;
       kt0 = Gmod/L*mean_thickness(iceIndex);
     }
     else{
       kt0 = Gmod/L*(1/(1/mean_thickness(i) + 1/mean_thickness(j)));
       F_FLOAT particleRadius = 5000.0;

       F_FLOAT previousForce = d_firsthistory(i,size_history*jj+5);
       F_FLOAT ridgeSlip     = d_firsthistory(i,size_history*jj+6);
       F_FLOAT ridgeSlipUsed = d_firsthistory(i,size_history*jj+7);

       hopkins_ridging_model(
	   NEWTON_PAIR || j < nlocal,
	   modifyState,
           delta,
           delta_dot,
           iceConcentration(i),
           iceConcentration(j),
           mean_thickness(i),
           mean_thickness(j),
           radius(i),
           radius(j),
           rmass(i),
           rmass(j),
           ridgingIceThickness(i),
           ridgingIceThickness(j),
           ridgingIceThicknessWeight(i),
           ridgingIceThicknessWeight(j),
           netToGrossClosingRatio(i),
           netToGrossClosingRatio(j),
           changeEffectiveElementArea(i),
           changeEffectiveElementArea(j),
           particleRadius,
           plasticFrictionCoeff,
           plasticHardeningCoeff,
           exponentialIceStrengthCoeff,
           L,
           ridgeSlip,
           ridgeSlipUsed,
           previousForce,
           contactForce);
       if (HISTORYUPDATE and modifyState){
         d_firsthistory(i,size_history*jj+5) = previousForce;
         d_firsthistory(i,size_history*jj+6) = ridgeSlip;
         d_firsthistory(i,size_history*jj+7) = ridgeSlipUsed;
       }
     }

     fnx = contactForce*nx;
     fny = contactForce*ny;

     // update tangential displacement, rotate if needed
     if (HISTORYUPDATE and modifyState){
       disp_tx = d_firsthistory(i,size_history*jj+2);
       disp_ty = d_firsthistory(i,size_history*jj+3);
       ndisp = nx*disp_tx + ny*disp_ty; //Tangential displacement in normal direction
       if (fabs(ndisp) > EPSILON){
         dispmag = sqrt(disp_tx*disp_tx + disp_ty*disp_ty);
         //Remove component along normal direction
         disp_tx -= ndisp*nx;
         disp_ty -= ndisp*ny;

         //Rescale to preserve magnitude
         prjmag = sqrt(disp_tx*disp_tx + disp_ty*disp_ty);
         if (prjmag > 0) scalefac = dispmag/prjmag;
         else scalefac = 0;
         disp_tx *= scalefac;
         disp_ty *= scalefac;
         d_firsthistory(i,size_history*jj+2) = disp_tx;
         d_firsthistory(i,size_history*jj+3) = disp_ty;
       }
       d_firsthistory(i,size_history*jj+2) += vtx*update_dt;
       d_firsthistory(i,size_history*jj+3) += vty*update_dt;
     }

     disp_tx = d_firsthistory(i,size_history*jj+2);
     disp_ty = d_firsthistory(i,size_history*jj+3);
     dispmag =sqrt(disp_tx*disp_tx + disp_ty*disp_ty);

     // total tangential force
     ftx = - (kt0*disp_tx + damp_tangential*vtx);
     fty = - (kt0*disp_ty + damp_tangential*vty);

     ftmag = sqrt(ftx*ftx + fty*fty);
     ftcrit = friction_tangential*fabs(contactForce);
     if (ftmag > ftcrit){
       if (dispmag != 0){
         ftx *= ftcrit/ftmag;
         fty *= ftcrit/ftmag;
	 if (modifyState) {
	   d_firsthistory(i,size_history*jj+2) = -(ftx + damp_tangential*vtx)/kt0;
	   d_firsthistory(i,size_history*jj+3) = -(fty + damp_tangential*vty)/kt0;
	 }
       }
       else ftx = fty = 0;
     }

    //Apply forces
    fx = fnx + ftx;
    fy = fny + fty;

    //Torque
    ncrossF = nx*fty - ny*ftx;
    torque_i = -radius(i)*ncrossF;

    torque_j = 0;
    if (NEWTON_PAIR || j < nlocal){
      torque_j = -radius(j)*ncrossF;
    }

  } // rsq < radsum*radsum
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::compute_bonded_kokkos(int i,
							      int j,
							      int jj,
							      F_FLOAT &fx,
							      F_FLOAT &fy,
							      F_FLOAT &fnx,
							      F_FLOAT &fny,
							      F_FLOAT &ftx,
							      F_FLOAT &fty,
							      F_FLOAT &torque_i,
							      F_FLOAT &torque_j,
							      bool modifyState) const
{

  //See design document for definitions of these variables
  F_FLOAT s1x, s1y, s2x, s2y, mx, my, mmag, mex, mey, bex, bey, rx, ry, rxj, ryj;
  F_FLOAT An, Bn, Cn, Dn, Bt, Ct, Dt, Bnj, Cnj, Dnj, Btj;
  F_FLOAT Fnmag, Ftmag, Nn, Nt, Nnj, Ntj;

  F_FLOAT chi1 = d_firsthistory(i,size_history*jj);
  F_FLOAT chi2 = d_firsthistory(i,size_history*jj+1);
  F_FLOAT chidiff, chidiff2, chidiff3;

  F_FLOAT sig_n1, sig_s1;
  F_FLOAT chi_c, chi_t, chi_s1, chi_s2;
  F_FLOAT nprefac, sprefac;
  F_FLOAT damp_prefac, area_bond, fdampx, fdampy, torquedamp;
  F_FLOAT hmin;

  F_FLOAT dvx, dvy;

  if (HISTORYUPDATE and modifyState){

    // Update bond end points based on particle translations
    d_firsthistory(i,size_history*jj+2) += dt*v(i,0);
    d_firsthistory(i,size_history*jj+3) += dt*v(i,1);
    d_firsthistory(i,size_history*jj+4) += dt*v(i,0);
    d_firsthistory(i,size_history*jj+5) += dt*v(i,1);

    d_firsthistory(i,size_history*jj+6) += dt*v(j,0);
    d_firsthistory(i,size_history*jj+7) += dt*v(j,1);
    d_firsthistory(i,size_history*jj+8) += dt*v(j,0);
    d_firsthistory(i,size_history*jj+9) += dt*v(j,1);

    // Update bond end points based on particle rotations
    d_firsthistory(i,size_history*jj+2) += -dt*omega(i,2)*(d_firsthistory(i,size_history*jj+3)-x(i,1));
    d_firsthistory(i,size_history*jj+3) +=  dt*omega(i,2)*(d_firsthistory(i,size_history*jj+2)-x(i,0));
    d_firsthistory(i,size_history*jj+4) += -dt*omega(i,2)*(d_firsthistory(i,size_history*jj+5)-x(i,1));
    d_firsthistory(i,size_history*jj+5) +=  dt*omega(i,2)*(d_firsthistory(i,size_history*jj+4)-x(i,0));

    d_firsthistory(i,size_history*jj+6) += -dt*omega(j,2)*(d_firsthistory(i,size_history*jj+7)-x(j,1));
    d_firsthistory(i,size_history*jj+7) +=  dt*omega(j,2)*(d_firsthistory(i,size_history*jj+6)-x(j,0));
    d_firsthistory(i,size_history*jj+8) += -dt*omega(j,2)*(d_firsthistory(i,size_history*jj+9)-x(j,1));
    d_firsthistory(i,size_history*jj+9) +=  dt*omega(j,2)*(d_firsthistory(i,size_history*jj+8)-x(j,0));
  }

  //Compute s_1, s_2, m, m_e, b_e
  s1x = d_firsthistory(i,size_history*jj+6) - d_firsthistory(i,size_history*jj+2);
  s1y = d_firsthistory(i,size_history*jj+7) - d_firsthistory(i,size_history*jj+3);
  s2x = d_firsthistory(i,size_history*jj+8) - d_firsthistory(i,size_history*jj+4);
  s2y = d_firsthistory(i,size_history*jj+9) - d_firsthistory(i,size_history*jj+5);

  mx = d_firsthistory(i,size_history*jj+4) + 0.5*s2x - d_firsthistory(i,size_history*jj+2) - 0.5*s1x;
  my = d_firsthistory(i,size_history*jj+5) + 0.5*s2y - d_firsthistory(i,size_history*jj+3) - 0.5*s1y;
  mmag = sqrt(mx*mx + my*my);
  mex = mx/mmag;
  mey = my/mmag;

  bex = mey;
  bey = -mex;

  rx = d_firsthistory(i,size_history*jj+2) + 0.5*s1x - x(i,0);
  ry = d_firsthistory(i,size_history*jj+3) + 0.5*s1y - x(i,1);

  //Compute dot product 
  // If orientation of element centers is in same direction as be vector, modify 
  // signs before computing stress in update_chi
  F_FLOAT bonddirx = x(j,0) - x(i,0);
  F_FLOAT bonddiry = x(j,1) - x(i,1);
  F_FLOAT mag = sqrt(bonddirx*bonddirx + bonddiry*bonddiry);
  F_FLOAT dotprod = bonddirx/mag*bex + bonddiry/mag*bey;

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
  F_FLOAT kn0 = d_firsthistory(i,size_history*jj+11)*Emod/d_firsthistory(i,size_history*jj+10);
  F_FLOAT kt0 = d_firsthistory(i,size_history*jj+11)*Gmod/d_firsthistory(i,size_history*jj+10);
  nprefac = d_firsthistory(i,size_history*jj+10)*kn0;
  sprefac = d_firsthistory(i,size_history*jj+10)*kt0;


  Fnmag = nprefac*(Dn*chidiff + Cn*chidiff2);
  Ftmag = sprefac*(Dt*chidiff + Ct*chidiff2);
  Nn = nprefac*(An*Cn*chidiff3 + (Bn*Cn+An*Dn)*chidiff2 + Bn*Dn*chidiff);
  Nt = Bt*Ftmag;

  //Damping force
  area_bond = d_firsthistory(i,size_history*jj+10)*d_firsthistory(i,size_history*jj+11)*chidiff;
  damp_prefac = damp_bonded*area_bond;

  dvx = v(j,0) - v(i,0);
  dvy = v(j,1) - v(i,1);
  fdampx = damp_prefac*dvx;
  fdampy = damp_prefac*dvy;

  //Damping torque
  torquedamp = -damp_bonded*area_bond*area_bond*(omega(i,2)-omega(j,2));

  //Update forces and torque
  fnx = Fnmag*bex;
  fny = Fnmag*bey;
  ftx = Ftmag*mex;
  fty = Ftmag*mey;

  fx = fnx + ftx;
  fy = fny + fty;

  fx = fx+fdampx;
  fy = fy+fdampy;

  torque_i = Nn + Nt + torquedamp;

  torque_j = 0;
  if (NEWTON_PAIR || j < nlocal){

    rxj = d_firsthistory(i,size_history*jj+2) + 0.5*s1x - x(j,0);
    ryj = d_firsthistory(i,size_history*jj+3) + 0.5*s1y - x(j,1);
    Bnj = rxj*bey - ryj*bex;
    Cnj = -Cn;
    Dnj = -Dn;
    Btj = rxj*mey - ryj*mex;

    Nnj = nprefac*(An*Cnj*chidiff3 + (Bnj*Cnj+An*Dnj)*chidiff2 + Bnj*Dnj*chidiff);
    Ntj = -Btj*Ftmag;
    torque_j = Nnj + Ntj - torquedamp;
  }

  //Update chi1, chi2
  hmin = MIN(min_thickness(i), min_thickness(j));
  if (HISTORYUPDATE and modifyState){
    F_FLOAT c1, c2;
    c1 = d_firsthistory(i,size_history*jj);
    c2 = d_firsthistory(i,size_history*jj+1);
    if (dotprod > 0.0){
       Dn = -Dn;
       Cn = -Cn;
       Dt = -Dt;
       Ct = -Ct;
    }
    update_chi(kn0, kt0, Dn, Cn, Dt, Ct, hmin, d_firsthistory(i,size_history*jj+11), c1, c2);
    d_firsthistory(i,size_history*jj) = c1;
    d_firsthistory(i,size_history*jj+1) = c2;
    d_firsttouch(i,jj) = 1;
    if (d_firsthistory(i,size_history*jj) >= d_firsthistory(i,size_history*jj+1)){ //Bond just broke
        F_FLOAT dx = x(i,0) - x(j,0);
        F_FLOAT dy = x(i,1) - x(j,1);
        F_FLOAT rij = sqrt(dx*dx + dy*dy);
        F_FLOAT delta_0 = radius(i) + radius(j) - rij;
        if (delta_0 < 0) delta_0 = 0;
        for (int k = 0; k < size_history; ++k)
          d_firsthistory(i,size_history*jj+k) = 0;
        d_firsthistory(i,size_history*jj+4) = delta_0;
    }
  }
}


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::update_chi(F_FLOAT kn0,
                                                   F_FLOAT kt0,
                                                   F_FLOAT Dn,
                                                   F_FLOAT Cn,
                                                   F_FLOAT Dt,
                                                   F_FLOAT Ct,
                                                   F_FLOAT hmin,
                                                   F_FLOAT bondThickness,
                                                   F_FLOAT &chi1,
                                                   F_FLOAT &chi2) const
{ 
  F_FLOAT sig_n1 = kn0*(Dn + Cn*chi1);
  F_FLOAT sig_s1 = kt0*(Dt + Ct*chi1);
  F_FLOAT sig_n2 = kn0*(Dn + Cn*chi2);
  F_FLOAT sig_s2 = kt0*(Dt + Ct*chi2);

  // function pointers for greater efficiency?
  F_FLOAT sig_c;
  if (strcmp_sig_c0_type_constant) {
    sig_c = sig_c0 * bondThickness;
  } else if (strcmp_sig_c0_type_KovacsSodhi) {
    sig_c = sig_c0*pow(hmin,(2.0/3.0)) * 1000.0;
  } // else error case already handled previously

  F_FLOAT sig_t;
  if (strcmp_sig_t0_type_constant) {
    sig_t = sig_t0 * bondThickness;
  } else if (strcmp_sig_t0_type_multiply_sig_c0) {
    sig_t = sig_t0 * sig_c;
  } // else error case already handled previously

  F_FLOAT denom;
  sig_t = -sig_t; //Somewhat against convention, tensile load is taken to be negative

  F_FLOAT c1, c2;
  c1 = chi1;
  c2 = chi2;
  //Check for purely tensile/compressive failure
  if (sig_n1 < sig_t || sig_n2 < sig_t || sig_n1 > sig_c || sig_n2 > sig_c){
    if (sig_n1 < sig_t){
      if (Cn != 0) chi1 = (sig_t/kn0 - Dn)/Cn;
      else {
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 < sig_t) chi2 = (sig_t/kn0 - Dn)/Cn;
    //if Cn==0, then sig_n1 = sig_n2, and function would've returned above

    if (sig_n1 > sig_c){
      if (Cn != 0) chi1 = (sig_c/kn0 - Dn)/Cn;
      else{
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 > sig_c) chi2 = (sig_c/kn0 - Dn)/Cn;
    //if Cn==0, function would've returned above
  }

  //Re-compute, since stress state could still be outside of failure envelope
  sig_n1 = kn0*(Dn + Cn*chi1);
  sig_s1 = kt0*(Dt + Ct*chi1);
  sig_n2 = kn0*(Dn + Cn*chi2);
  sig_s2 = kt0*(Dt + Ct*chi2);

  //Check for 'cohesion' shear failure at chi1
  if (sig_s1 > tanphi*sig_n1 - tanphi*sig_t){
    denom = tanphi*kn0*Cn - kt0*Ct;
    if (denom != 0) chi1 = (kt0*Dt+tanphi*(sig_t-kn0*Dn))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        //Rare case of kt*Ct = -kn*tanphi*Cn, yield criterion is independent of chi,
        // therefore bond breaks.
        chi1 = chi2;
        return;
       }
       else if (Cn == 0 && Ct == 0){
        //Rare case of s1 = s2, again yield criterion is independent of chi,
        // and there is a possibility of shear failure of entire bond
        if (kt0*Dt > tanphi*(kn0*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }

  //Bottom branch of envelope
  if (sig_s1 < -tanphi*sig_n1 + tanphi*sig_t){
    denom = -kn0*tanphi*Cn - kt0*Ct;
    if (denom != 0) chi1 = (kt0*Dt+tanphi*(kn0*Dn-sig_t))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        chi1 = chi2;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt0*Dt < -tanphi*(kn0*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }
    
  //Check for 'cohesion' shear failure at chi2
  //Top branch of envelope
  if (sig_s2 > tanphi*sig_n2 - tanphi*sig_t){
    denom = tanphi*kn0*Cn - kt0*Ct;
    if (denom != 0) chi2 = (kt0*Dt+tanphi*(sig_t-kn0*Dn))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        //Rare case of kt*Ct = -kn*tanphi*Cn, yield criterion is independent of chi,
        // therefore bond breaks.
        chi2 = chi1;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt0*Dt > tanphi*(kn0*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }

  //Bottom branch of envelope
  if (sig_s2 < -tanphi*sig_n2 + tanphi*sig_t){
    denom = -kn0*tanphi*Cn - kt0*Ct;
    if (denom != 0) chi2 = (kt0*Dt+tanphi*(kn0*Dn-sig_t))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        chi1 = chi2;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt0*Dt < -tanphi*(kn0*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }
  if (chi1 < 0 || chi1 > 1 || chi2 < 0 || chi2 > 1)
    chi1 = chi2 = 0;

}


template<class DeviceType>
template<int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                                                          F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                                                          X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
  F_FLOAT v[6];

  v[0] = delx*fx;
  v[1] = dely*fy;
  v[2] = delz*fz;
  v[3] = delx*fy;
  v[4] = delx*fz;
  v[5] = dely*fz;

  if (NEWTON_PAIR) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  } else {
    if (i < nlocal) {
      ev.v[0] += 0.5*v[0];
      ev.v[1] += 0.5*v[1];
      ev.v[2] += 0.5*v[2];
      ev.v[3] += 0.5*v[3];
      ev.v[4] += 0.5*v[4];
      ev.v[5] += 0.5*v[5];
    }
    if (j < nlocal) {
      ev.v[0] += 0.5*v[0];
      ev.v[1] += 0.5*v[1];
      ev.v[2] += 0.5*v[2];
      ev.v[3] += 0.5*v[3];
      ev.v[4] += 0.5*v[4];
      ev.v[5] += 0.5*v[5];
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::ev_tally_xyz_atom(EV_FLOAT &ev, int i, int j,
                                                          F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                                                          X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
//  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  F_FLOAT v[6];

  v[0] = delx*fx;
  v[1] = dely*fy;
  v[2] = delz*fz;
  v[3] = delx*fy;
  v[4] = delx*fz;
  v[5] = dely*fz;

  if (NEWTON_PAIR || i < nlocal) {
    v_vatom(i,0) += 0.5*v[0];
    v_vatom(i,1) += 0.5*v[1];
    v_vatom(i,2) += 0.5*v[2];
    v_vatom(i,3) += 0.5*v[3];
    v_vatom(i,4) += 0.5*v[4];
    v_vatom(i,5) += 0.5*v[5];
  }
  if (NEWTON_PAIR || j < nlocal) {
    v_vatom(j,0) += 0.5*v[0];
    v_vatom(j,1) += 0.5*v[1];
    v_vatom(j,2) += 0.5*v[2];
    v_vatom(j,3) += 0.5*v[3];
    v_vatom(j,4) += 0.5*v[4];
    v_vatom(j,5) += 0.5*v[5];
  }
}



namespace LAMMPS_NS {
template class PairGranHopkinsKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairGranHopkinsKokkos<LMPHostType>;
#endif
}

