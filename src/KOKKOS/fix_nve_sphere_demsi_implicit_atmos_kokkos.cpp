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

#include "fix_nve_sphere_demsi_implicit_atmos_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::FixNVESphereDemsiImplicitAtmosKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphereDemsiImplicitAtmos(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK| F_MASK | OMEGA_MASK | TORQUE_MASK |
   RMASS_MASK | RADIUS_MASK | MASK_MASK | THICKNESS_MASK | VN_MASK; // adding vn

  datamask_modify = X_MASK | V_MASK | OMEGA_MASK |
   THICKNESS_MASK | VN_MASK; // adding vn
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
  vatom = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::init()
{
  FixNVESphereDemsiImplicitAtmos::init();

//  if (extra == DIPOLE) {
//    error->all(FLERR,"Fix nve/sphere/demsi/kk doesn't yet support dipole");
//  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::initial_integrate(int vflag)
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
  orientation = atomKK->k_orientation.view<DeviceType>();
  momentOfInertia = atomKK->k_momentOfInertia.view<DeviceType>();
  ice_area = atomKK->k_ice_area.view<DeviceType>();
  coriolis = atomKK->k_coriolis.view<DeviceType>();
  ocean_vel = atomKK->k_ocean_vel.view<DeviceType>();
  bvector = atomKK->k_bvector.view<DeviceType>();
  forcing = atomKK->k_forcing.view<DeviceType>();
  mean_thickness = atomKK->k_mean_thickness.view<DeviceType>();
  vn = atomKK->k_vn.view<DeviceType>(); // adding vn

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereDemsiImplicitAtmosKokkosInitialIntegrateFunctor<DeviceType> f(this); 
  Kokkos::parallel_for(nlocal,f);

}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::initial_integrate_item(const int i) const
{

  if (mask(i) & groupbit) {

    // acceleration with time(n) forces
    double dtm = dtf/rmass(i);
    double dvx = dtm*f(i,0);
    double dvy = dtm*f(i,1);
    double dvz = dtf*torque(i,2)/momentOfInertia(i);

    // add the forcing
    double u1 = v(i,0) + dvx + dtm*forcing(i,0);
    double v1 = v(i,1) + dvy + dtm*forcing(i,1);
    double u2 = bvector(i,0);
    double v2 = bvector(i,1);
    double u3 = ocean_vel(i,0);
    double v3 = ocean_vel(i,1);

    // 3-field finite-rate momentum exchange: field 1-ice, 2-atmosphere, 3-ocean
    double ice_Area = ice_area(i);
    double ice_Volume = fmax(1.e-99, ice_Area*mean_thickness(i));
    double ice_Density = rmass(i)/ice_Volume;
    double w12 = sqrt(pow(u2-u1,2) + pow(v2-v1,2)) + Hugoniot_Vel_Jump;
    double w13 = sqrt(pow(u3-u1,2) + pow(v3-v1,2)) + Hugoniot_Vel_Jump;
    double D12 = atmos_density*atmos_drag*w12*ice_Area;
    double D13 = ocean_density*ocean_drag*w13*ice_Area;
    double a = D12*dtf/(ice_Density*ice_Volume);
    double b = D13*dtf/(ice_Density*ice_Volume);

    // fixed u2,u3 form:
    double rden = 1./(1. + a + b);
    double up1 = (u1 + a*u2 + b*u3)*rden;
    double vp1 = (v1 + a*v2 + b*v3)*rden;
    /*
      double c = D12*dtf/(atmos_density*ice_Volume);
      double d = D13*dtf/(ocean_density*ice_Volume);
      //A = {{1 + a + b,    -a,    -b},
      //     {       -c, 1 + c,     0},
      //     {       -d,     0, 1 + d}};
      //B = Inverse[A]
      double rden = 1./(1. + a + b + c + b*c + d + a*d + c*d);
      double B11 = (1. + c + d + c*d);
      double B12 = (a + a*d);
      double B13 = (b + b*c);
      double up1 = (B11*u1 + B12*u2 + B13*u3)*rden;
      double vp1 = (B11*v1 + B12*v2 + B13*v3)*rden;
    */
    // add the coriolis acceleration
    double cdt = coriolis(i)*dtf*rden;
    double cden = 1./(1. + cdt*cdt);
    v(i,0) = (up1 + cdt*vp1)*cden;
    v(i,1) = (vp1 - cdt*up1)*cden;
    omega(i,2) += dvz;

    x(i,0) += dtv * v(i,0);  // full step with time(n + 1/2) vel
    x(i,1) += dtv * v(i,1);
    orientation(i) += dtv * omega(i,2);

    // put v, omega into vn
    vn(i,0) =     v(i,0);
    vn(i,1) =     v(i,1);
    vn(i,2) = omega(i,2);

  } // end if (mask(i) & groupbit)

} // end initial_integrate_item

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  orientation = atomKK->k_orientation.view<DeviceType>();
  momentOfInertia = atomKK->k_momentOfInertia.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  ice_area = atomKK->k_ice_area.view<DeviceType>();
  coriolis = atomKK->k_coriolis.view<DeviceType>();
  ocean_vel = atomKK->k_ocean_vel.view<DeviceType>();
  bvector = atomKK->k_bvector.view<DeviceType>();
  forcing = atomKK->k_forcing.view<DeviceType>();
  mean_thickness = atomKK->k_mean_thickness.view<DeviceType>();
  vn = atomKK->k_vn.view<DeviceType>(); // adding vn

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereDemsiImplicitAtmosKokkosFinalIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereDemsiImplicitAtmosKokkos<DeviceType>::final_integrate_item(const int i) const
{

  if (mask(i) & groupbit) {

    // get v, omega from vn
    double up1 = vn(i,0);
    double vp1 = vn(i,1);
    double wp1 = vn(i,2);

    // half step acceleration with time(n+1) forces
    double dtm = dtf/rmass(i);
    double dvx = dtm*f(i,0);
    double dvy = dtm*f(i,1);
    double dvz = dtf*torque(i,2)/momentOfInertia(i);

    // add the forcing
    double u1 = up1 + dvx + dtm*forcing(i,0);
    double v1 = vp1 + dvy + dtm*forcing(i,1);

    double u2 = bvector(i,0);
    double v2 = bvector(i,1);
    double u3 = ocean_vel(i,0);
    double v3 = ocean_vel(i,1);

    // 3-field finite-rate momentum exchange: field 1-ice, 2-atmosphere, 3-ocean
    double ice_Area = ice_area(i);
    double ice_Volume = fmax(1.e-99, ice_Area*mean_thickness(i));
    double ice_Density = rmass(i)/ice_Volume;
    double w12 = sqrt(pow(u2-u1,2) + pow(v2-v1,2)) + Hugoniot_Vel_Jump;
    double w13 = sqrt(pow(u3-u1,2) + pow(v3-v1,2)) + Hugoniot_Vel_Jump;
    double D12 = atmos_density*atmos_drag*w12*ice_Area;
    double D13 = ocean_density*ocean_drag*w13*ice_Area;
    double a = D12*dtf/(ice_Density*ice_Volume);
    double b = D13*dtf/(ice_Density*ice_Volume);

    // fixed u2,u3 form:
    double rden = 1./(1. + a + b);
    up1 = (u1 + a*u2 + b*u3)*rden;
    vp1 = (v1 + a*v2 + b*v3)*rden;
    /*
      double c = D12*dtf/(atmos_density*ice_Volume);
      double d = D13*dtf/(ocean_density*ice_Volume);
      //A = {{1 + a + b,    -a,    -b},
      //     {       -c, 1 + c,     0},
      //     {       -d,     0, 1 + d}};
      //B = Inverse[A]
      double rden = 1./(1. + a + b + c + b*c + d + a*d + c*d);
      double B11 = (1. + c + d + c*d);
      double B12 = (a + a*d);
      double B13 = (b + b*c);
      up1 = (B11*u1 + B12*u2 + B13*u3)*rden;
      vp1 = (B11*v1 + B12*v2 + B13*v3)*rden;
    */
    // add the coriolis acceleration
    double cdt = coriolis(i)*dtf*rden;
    double cden = 1./(1. + cdt*cdt);
    v(i,0) = (up1 + cdt*vp1)*cden;
    v(i,1) = (vp1 - cdt*up1)*cden;
    omega(i,2) = wp1 + dvz;

  } // end if (mask(i) & groupbit)

    // zero out diagnostics on nonmoving elements
  if(!(mask(i)&groupbit)) {
        vn(i,0) =     vn(i,1) =     vn(i,2) = 0.;
    torque(i,0) = torque(i,1) = torque(i,2) = 0.;
     omega(i,0) =  omega(i,1) =  omega(i,2) = 0.;
  }

} // end final_integrate_item

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNVESphereDemsiImplicitAtmosKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixNVESphereDemsiImplicitAtmosKokkos<LMPHostType>;
#endif
}

/* ---------------------------------------------------------------------- */
