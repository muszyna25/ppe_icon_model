#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!! Type definition for the dynamical core of ICOHAM.
!!
!! @par Revision History
!! Hui Wan (MPI-M, 2010-07-06)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!     violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!     copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_icoham_dyn_types

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_hydro_atm_prog, t_hydro_atm_diag, t_hydro_atm

#ifdef HAVE_F95
  PUBLIC :: t_ptr3d
#endif

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !>
  !! Derived data type for building pointer arrays
  !!
  TYPE t_ptr3d
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr3d

  !>
  !!--------------------------------------------------------------------------
  !! Derived data type for prognostic variables. The same type is used for
  !! tendencies

  TYPE t_hydro_atm_prog

    REAL(wp), POINTER ::  &
    & pres_sfc(:,  :),  &!< surface pressure [Pa]        (nproma,     nblks_c)
    &       vn(:,:,:),  &!< normal wind [m/s]            (nproma,nlev,nblks_e)
    &     temp(:,:,:),  &!< temperature [K]              (nproma,nlev,nblks_c)
    &    theta(:,:,:),  &!< potential temperature [K]    (nproma,nlev,nblks_c)
    &   tracer(:,:,:,:)  !< tracer concentration [kg/kg] (nproma,nlev,nblks_c,ntracer)

    TYPE(t_ptr3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

  END TYPE t_hydro_atm_prog

  !>
  !!--------------------------------------------------------------------------
  !! Derived data type for diagnostic variables

  TYPE t_hydro_atm_diag

    REAL(wp), POINTER ::  &
    &          qx(:,:,:),   &!< total concentration of hydrometeors (nproma,nlev,nblks_c)
    &           u(:,:,:),   &!< zonal wind (nproma,nlev,nblks_c)
    &           v(:,:,:),   &!< meridional wind (nproma,nlev,nblks_c)
    &          vt(:,:,:),   &!< tangential wind (nproma,nlev,nblks_e)
    &    rel_vort(:,:,:),   &!< relative vorticity at dual point (nproma,nlev,nblks_v)
    &  rel_vort_e(:,:,:),   &!< needed for hexagonal model
    &  rel_vort_c(:,:,:),   &!< relative vorticity at cell centers, diagnosed for physics 
    &         div(:,:,:),   &!< wind divergence (only output) (nproma,nlev,nblks_c)
    &       e_kin(:,:,:),   &!< specific kinetic energy (nproma,nlev,nblks_c)
    &      geo_ic(:,:,:),   &!< half level geopotential (nproma,nlevp1,nblks_c)
    &      geo_mc(:,:,:),   &!< full level geopotential (nproma,nlev,nblks_c)
    &    wpres_mc(:,:,:),   &!< vert. vel. in pres. coord. at full levels (nproma,nlev,nblks_c)
    &    wpres_ic(:,:,:),   &!< vert. vel. in pres. coord. at half levels (nproma,nlevp1,nblks_c)
    &        weta(:,:,:),   &!< vert. vel. in $\eta$ coord. times dpdeta (nproma,nlevp1,nblks_c)
                             !< i.e. mass flux in $\eta$ coord. divided by delta $\eta$
    &     pres_ic(:,:,:),   &!< half level pressure (nproma,nlevp1,nblks_c)
    & pres_ic_new(:,:,:),   &!< dito, but at timestep n+1
    &     pres_mc(:,:,:),   &!< full level pressure (nproma,nlev,nblks_c)
    &       exner(:,:,:),   &!< exner function (for theta advection; nproma,nlev,nblks_c)
    &   virt_incr(:,:,:),   &!< virtual temperature increment (nproma,nlev,nblks_c)
    &     tempv(:,:,:),     &!< vertual temperature (nproma,nlev,nblks_c)
    &      delp_c(:,:,:),   &!< layer thickness at cell centers (nproma,nlev,nblks_c)
    &  delp_c_new(:,:,:),   &!< dito, but at timestep n+1
    &     rdelp_c(:,:,:),   &!< reciprocal layer thickness at cell centers (nproma,nlev,nblks_c)
    & rdelp_c_new(:,:,:),   &!< dito, but at timestep n+1
    &      delp_e(:,:,:),   &!< layer thickness at edges (nproma,nlev,nblks_e)
    &      delp_v(:,:,:),   &!< layer thickness at dual point (nproma,nlev,nblks_v)
    & hfl_tracer(:,:,:,:),  &!< horizontal tracer flux at edges (nproma,nlev,nblks_e,ntracer)
    & vfl_tracer(:,:,:,:),  &!< vertical tracer flux at cells (nproma,nlevp1,nblks_c,ntracer)
    &    rdlnpr_c(:,:,:),   &!< Rd * ln(p(k+.5)/p(k-.5)),shape:(nproma,nlev,nblks_c)
    &   rdalpha_c(:,:,:),   &!< Rd * alpha              ,shape:(nproma,nlev,nblks_c)
    &      lnp_ic(:,:,:),   &!< ln(p),shape:(nproma,nlevp1,nblks_c)
    & mass_flux_e(:,:,:)     !< mass flux at edges (nproma,nlev,nblks_e)

    TYPE(t_ptr3d),ALLOCATABLE :: hfl_tracer_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr3d),ALLOCATABLE :: vfl_tracer_ptr(:)  !< pointer array: one pointer for each tracer

  END TYPE t_hydro_atm_diag

  !>
  !!--------------------------------------------------------------------------
  !! Derived data type for the hydrostatic state vector on a single grid level.
  !! This type is in fact a wrapper which was necessary in revisions up to 1916
  !! because the actual variable used to be defined in the main program and
  !! passed through various subroutine interfaces.

  TYPE t_hydro_atm

    TYPE(t_hydro_atm_prog),ALLOCATABLE :: prog(:)  !< shape: (nTimeLevel)
    TYPE(t_hydro_atm_diag)             :: diag     !< diagnostic variables
    TYPE(t_hydro_atm_prog)             :: tend_dyn !< tendency due to dynamics
    TYPE(t_hydro_atm_prog)             :: tend_phy !< tendency due to physics

    TYPE(t_hydro_atm_prog) :: prog_out  !< for output
    TYPE(t_hydro_atm_diag) :: diag_out  !< for output

  END TYPE t_hydro_atm

END module mo_icoham_dyn_types
