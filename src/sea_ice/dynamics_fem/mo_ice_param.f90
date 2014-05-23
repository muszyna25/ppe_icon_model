!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
module mo_ice_param
USE mo_kind,    ONLY: wp
  implicit none
  PUBLIC
  save
  REAL(wp), parameter  :: pi=3.141592653589793_wp, rad=pi/180.0_wp
  REAL(wp), parameter  :: omega=2.0_wp*pi/(24.0_wp*60.0_wp*60.0_wp)
  REAL(wp)             :: g=9.806_wp                            ![m/s^2]
  REAL(wp), parameter  :: r_earth =6.3675e6_wp                  ![m]
  REAL(wp), parameter  :: density_0=1027.7_wp                   ![kg/m^3]
  REAL(wp)             :: Ah0 = 1.e3_wp, Av0=1.e-4_wp
  REAL(wp)             :: Kh0 = 5.e1_wp, Kv0=1.e-5_wp
  REAL(wp)             :: dt!=3600._wp/2._wp
  REAL(wp)             :: scalevol=1.0e8_wp
  integer                  :: num_iter_solve=3
  integer(KIND=4)          :: iter
  !
  logical                  :: lfirst=.true.            ! used in solver
  logical                  :: r_restart=.false.
  logical                  :: cartesian=.false.
  !
  logical                  :: cyclic_geometry=.true.    !.true. for periodical geometry
  REAL(wp)             :: domain_length=360.0_wp*rad    !unit in rad.
  !
  !the mesh is stored there:

  character*100            :: meshpath='./mesh/'
  character*100            :: datapath='./mesh/'

end module mo_ice_param
!
