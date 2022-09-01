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

#ifdef FIX_CLASS

FixStyle(nve/sphere/demsi_implicit_atmos/kk,FixNVESphereDemsiImplicitAtmosKokkos<LMPDeviceType>)
FixStyle(nve/sphere/demsi_implicit_atmos/kk/device,FixNVESphereDemsiImplicitAtmosKokkos<LMPDeviceType>)
FixStyle(nve/sphere/demsi_implicit_atmos/kk/host,FixNVESphereDemsiImplicitAtmosKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_NVE_SPHERE_DEMSI_IMPLICIT_ATMOS_KOKKOS_H
#define LMP_FIX_NVE_SPHERE_DEMSI_IMPLICIT_ATMOS_KOKKOS_H

#include "fix_nve_sphere_demsi_implicit_atmos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {
  
template<class DeviceType>
class FixNVESphereDemsiImplicitAtmosKokkos : public FixNVESphereDemsiImplicitAtmos {
  public:
    FixNVESphereDemsiImplicitAtmosKokkos(class LAMMPS *, int, char **);
    virtual ~FixNVESphereDemsiImplicitAtmosKokkos() {}
    void cleanup_copy();
    void init();
    void initial_integrate(int);
    void final_integrate();
  
    KOKKOS_INLINE_FUNCTION
    void initial_integrate_item(const int i) const;
    KOKKOS_INLINE_FUNCTION
    void final_integrate_item(const int i) const;

  private:
    typename ArrayTypes<DeviceType>::t_x_array x;
    typename ArrayTypes<DeviceType>::t_v_array v;
    typename ArrayTypes<DeviceType>::t_v_array omega;
    typename ArrayTypes<DeviceType>::t_f_array f;
    typename ArrayTypes<DeviceType>::t_f_array torque;
    typename ArrayTypes<DeviceType>::t_float_1d rmass;
    typename ArrayTypes<DeviceType>::t_float_1d momentOfInertia;
    typename ArrayTypes<DeviceType>::t_float_1d orientation;
    typename ArrayTypes<DeviceType>::t_float_1d radius;
    typename ArrayTypes<DeviceType>::t_int_1d mask;
    typename ArrayTypes<DeviceType>::t_float_2d forcing;
    typename ArrayTypes<DeviceType>::t_float_1d ice_area;
    typename ArrayTypes<DeviceType>::t_float_1d coriolis;
    typename ArrayTypes<DeviceType>::t_float_2d ocean_vel;
    typename ArrayTypes<DeviceType>::t_float_2d bvector;
    typename ArrayTypes<DeviceType>::t_float_1d mean_thickness;
    typename ArrayTypes<DeviceType>::t_v_array vn;
};

template <class DeviceType>
struct FixNVESphereDemsiImplicitAtmosKokkosInitialIntegrateFunctor {
  FixNVESphereDemsiImplicitAtmosKokkos<DeviceType> c;
  FixNVESphereDemsiImplicitAtmosKokkosInitialIntegrateFunctor(FixNVESphereDemsiImplicitAtmosKokkos<DeviceType> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.initial_integrate_item(i);
  }
};

template <class DeviceType>
struct FixNVESphereDemsiImplicitAtmosKokkosFinalIntegrateFunctor {
  FixNVESphereDemsiImplicitAtmosKokkos<DeviceType> c;
  FixNVESphereDemsiImplicitAtmosKokkosFinalIntegrateFunctor(FixNVESphereDemsiImplicitAtmosKokkos<DeviceType> *c_ptr): c(*c_ptr) { c.cleanup_copy(); }
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.final_integrate_item(i);
  }
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_NVE_SPHERE_DEMSI_IMPLICIT_ATMOS_KOKKOS_H
#endif // FIX_CLASS
