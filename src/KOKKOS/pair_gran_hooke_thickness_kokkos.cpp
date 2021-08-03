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

#include "pair_gran_hooke_thickness_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "error.h"
#include "modify.h"
#include "fix_neigh_history_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHookeThicknessKokkos<DeviceType>::PairGranHookeThicknessKokkos(LAMMPS *lmp) : PairGranHookeThickness(lmp)
{
   std::cout << "In pair gran hooke thickness kokkos... \n";
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | 
                  MASK_MASK | ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK | THICKNESS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHookeThicknessKokkos<DeviceType>::~PairGranHookeThicknessKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairGranHookeThicknessKokkos<DeviceType>::init_style()
{
  if (history && fix_history == nullptr) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",3);
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
    error->all(FLERR,"Cannot use chosen neighbor list style with gran/hooke/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGranHookeThicknessKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  if (eflag || vflag) ev_setup(eflag,vflag,0);
  else evflag = vflag_fdotr = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
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
  mean_thickness = atomKK->k_mean_thickness.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  EV_FLOAT ev;

  if (lmp->kokkos->neighflag == HALF) {
    if (force->newton_pair) {
      if (vflag_atom) {
         Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,1,2>>(0,inum),*this);
      } else if (vflag_global) {
         Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,1,1>>(0,inum),*this, ev);
      } else {
         Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,1,0>>(0,inum),*this);
      }
    } else {
      if (vflag_atom) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,0,2>>(0,inum),*this);
      } else if (vflag_global) {
	Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,0,1>>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALF,0,0>>(0,inum),*this);
      }
    }
  } else { // HALFTHREAD
    if (force->newton_pair) {
      if (vflag_atom) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,1,2>>(0,inum),*this);
      } else if (vflag_global) {
	Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,1,1>>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,1,0>>(0,inum),*this);
      }
    } else {
      if (vflag_atom) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,0,2>>(0,inum),*this);
      } else if (vflag_global) {
	Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,0,1>>(0,inum),*this, ev);
      } else {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeThicknessCompute<HALFTHREAD,0,0>>(0,inum),*this);
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
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGranHookeThicknessKokkos<DeviceType>::operator()(TagPairGranHookeThicknessCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const LMP_FLOAT irad = radius(i);
  const LMP_FLOAT imass = rmass(i);
  const int jnum = d_numneigh(i);
  //printf("gran hooke thickness %d, %f %f %f %f \n", i, imass, irad, mean_thickness(i));
  
  F_FLOAT fx_i = 0.0;
  F_FLOAT fy_i = 0.0;
  F_FLOAT fz_i = 0.0;

  F_FLOAT torquex_i = 0.0;
  F_FLOAT torquey_i = 0.0;
  F_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    
    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const LMP_FLOAT jrad = radius(j);
    const LMP_FLOAT jmass = rmass(j);
    const LMP_FLOAT radsum = irad + jrad;

    // check for touching neighbors
    if (rsq < radsum*radsum) {

      const LMP_FLOAT r = sqrt(rsq);
      const LMP_FLOAT rinv = 1.0/r;
      const LMP_FLOAT rsqinv = 1/rsq;

      // relative translational velocity
      V_FLOAT vr1 = v(i,0) - v(j,0);
      V_FLOAT vr2 = v(i,1) - v(j,1);
      V_FLOAT vr3 = v(i,2) - v(j,2);
      
      // normal component
      V_FLOAT vnnr = vr1*delx + vr2*dely + vr3*delz;
      V_FLOAT vn1 = delx*vnnr * rsqinv;
      V_FLOAT vn2 = dely*vnnr * rsqinv;
      V_FLOAT vn3 = delz*vnnr * rsqinv;
      
      // tangential component
      V_FLOAT vt1 = vr1 - vn1;
      V_FLOAT vt2 = vr2 - vn2;
      V_FLOAT vt3 = vr3 - vn3;

      // relative rotational velocity
      V_FLOAT wr1 = (irad*omega(i,0) + jrad*omega(j,0)) * rinv;
      V_FLOAT wr2 = (irad*omega(i,1) + jrad*omega(j,1)) * rinv;
      V_FLOAT wr3 = (irad*omega(i,2) + jrad*omega(j,2)) * rinv;

      // meff = effective mass of pair of particles
      // if I or J part of rigid body, use body mass
      // if I or J is frozen, meff is other particle
      LMP_FLOAT meff = imass*jmass / (imass+jmass);
      if (mask[i] & freeze_group_bit) meff = jmass;
      if (mask[j] & freeze_group_bit) meff = imass;

      F_FLOAT min_thick = MIN(mean_thickness(i),mean_thickness(j));
      //printf("gran_hooke_thick %d %d %f \n", i, j, min_thick);

      // normal forces = Hookian contact + normal velocity damping
      F_FLOAT damp = meff*gamman*vnnr*rsqinv;
      F_FLOAT ccel = min_thick*kn*(radsum-r)*rinv - damp;

      // relative velocities
      V_FLOAT vtr1 = vt1 - (delz*wr2-dely*wr3);
      V_FLOAT vtr2 = vt2 - (delx*wr3-delz*wr1);
      V_FLOAT vtr3 = vt3 - (dely*wr1-delx*wr2);
      V_FLOAT vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
      vrel = sqrt(vrel);

      // force normalization
      F_FLOAT fn = xmu*fabs(ccel*r);
      F_FLOAT fs = meff*min_thick*gammat*vrel;
      F_FLOAT ft = 0.0;
      if (vrel != 0.0) ft = MIN(fn,fs) / vrel;

      // tangential force due to tangential velocity damping
      F_FLOAT fs1 = -ft*vtr1;
      F_FLOAT fs2 = -ft*vtr2;
      F_FLOAT fs3 = -ft*vtr3;

      // forces & torques
      F_FLOAT fx = delx*ccel + fs1;
      F_FLOAT fy = dely*ccel + fs2;
      F_FLOAT fz = delz*ccel + fs3;
      fx_i += fx;
      fy_i += fy;
      fz_i += fz;

      F_FLOAT tor1 = rinv * (dely*fs3 - delz*fs2);
      F_FLOAT tor2 = rinv * (delz*fs1 - delx*fs3);
      F_FLOAT tor3 = rinv * (delx*fs2 - dely*fs1);
      torquex_i -= irad*tor1;
      torquey_i -= irad*tor2;
      torquez_i -= irad*tor3;
      
      if (NEWTON_PAIR || j < nlocal) {
        a_f(j,0) -= fx;
        a_f(j,1) -= fy;
        a_f(j,2) -= fz;
        a_torque(j,0) -= jrad*tor1;
        a_torque(j,1) -= jrad*tor2;
        a_torque(j,2) -= jrad*tor3;
      }

      if (EVFLAG == 2)
        ev_tally_xyz_atom<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, fx_i, fy_i, fz_i, delx, dely, delz);
      if (EVFLAG == 1)
        ev_tally_xyz<NEWTON_PAIR>(ev, i, j, fx_i, fy_i, fz_i, delx, dely, delz);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;

}


template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairGranHookeThicknessKokkos<DeviceType>::operator()(TagPairGranHookeThicknessCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairGranHookeThicknessCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

template<class DeviceType>
template<int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHookeThicknessKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
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
void PairGranHookeThicknessKokkos<DeviceType>::ev_tally_xyz_atom(EV_FLOAT &ev, int i, int j,
		                                                F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
					            	        X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();  

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
template class PairGranHookeThicknessKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairGranHookeThicknessKokkos<LMPHostType>;
#endif
}
