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

#include "fix_nve_sphere_demsi_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVESphereDemsiKokkos<DeviceType>::FixNVESphereDemsiKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphereDemsi(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = F_MASK | TORQUE_MASK | RMASS_MASK | RADIUS_MASK | MASK_MASK | THICKNESS_MASK;
  datamask_modify = X_MASK | V_MASK | OMEGA_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiKokkos<DeviceType>::init()
{
  FixNVESphereDemsi::init();

//  if (extra == DIPOLE) {
//    error->all(FLERR,"Fix nve/sphere/demsi/kk doesn't yet support dipole");
//  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  ice_area = atomKK->k_ice_area.view<DeviceType>();
  coriolis = atomKK->k_coriolis.view<DeviceType>();
  ocean_vel = atomKK->k_ocean_vel.view<DeviceType>();
  bvector = atomKK->k_bvector.view<DeviceType>();
  forcing = atomKK->k_forcing.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereDemsiKokkosInitialIntegrateFunctor<DeviceType> f(this); 
  Kokkos::parallel_for(nlocal,f);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereDemsiKokkos<DeviceType>::initial_integrate_item(const int i) const
{
  const double dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {
    const double vel_diff = sqrt((ocean_vel(i,0)-v(i,0))*(ocean_vel(i,0)-v(i,0)) +
         (ocean_vel(i,1)-v(i,1))*(ocean_vel(i,1)-v(i,1)));
    const double D = ice_area(i)*ocean_drag*ocean_density*vel_diff;
    const double m_prime = rmass(i)/dtf;
    const double a00 = m_prime+D;
    const double a11 = a00;
    const double a10 = rmass(i)*coriolis(i);
    const double a01 = -a10;

    const double b0 = m_prime*v(i,0) + f(i,0) + bvector(i,0) + forcing(i,0) + D*ocean_vel(i,0);
    const double b1 = m_prime*v(i,1) + f(i,1) + bvector(i,1) + forcing(i,1) + D*ocean_vel(i,1);

    const double detinv = 1.0/(a00*a11 - a01*a10);
    v(i,0) = detinv*( a11*b0 - a01*b1);
    v(i,1) = detinv*(-a10*b0 + a00*b1);

    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);

    const double dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,2) += dtirotate * torque(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  ice_area = atomKK->k_ice_area.view<DeviceType>();
  coriolis = atomKK->k_coriolis.view<DeviceType>();
  ocean_vel = atomKK->k_ocean_vel.view<DeviceType>();
  bvector = atomKK->k_bvector.view<DeviceType>();
  forcing = atomKK->k_forcing.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereDemsiKokkosFinalIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereDemsiKokkos<DeviceType>::final_integrate_item(const int i) const
{
  const double dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {

    const double vel_diff = sqrt((ocean_vel(i,0)-v(i,0))*(ocean_vel(i,0)-v(i,0)) +
         (ocean_vel(i,1)-v(i,1))*(ocean_vel(i,1)-v(i,1)));
    const double D = ice_area(i)*ocean_drag*ocean_density*vel_diff;
    const double m_prime = rmass(i)/dtf;
    const double a00 = m_prime+D;
    const double a11 = a00;
    const double a10 = rmass(i)*coriolis(i);
    const double a01 = -a10;

    const double b0 = m_prime*v(i,0) + f(i,0) + bvector(i,0) + forcing(i,0) + D*ocean_vel(i,0);
    const double b1 = m_prime*v(i,1) + f(i,1) + bvector(i,1) + forcing(i,1) + D*ocean_vel(i,1);

    const double detinv = 1.0/(a00*a11 - a01*a10);
    v(i,0) = detinv*( a11*b0 - a01*b1);
    v(i,1) = detinv*(-a10*b0 + a00*b1);

    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);

    const double dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,2) += dtirotate * torque(i,2);

  //  rke += (omega(i,0)*omega(i,0) + omega(i,1)*omega(i,1) +
  //         omega(i,2)*omega(i,2))*radius(i)*radius(i)*rmass(i);
 
  }
}

namespace LAMMPS_NS {
template class FixNVESphereDemsiKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixNVESphereDemsiKokkos<LMPHostType>;
#endif
}
