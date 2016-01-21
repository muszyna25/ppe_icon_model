!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_tracer_metadata_types

  USE mo_kind,           ONLY: wp

  IMPLICIT NONE

  PRIVATE

  ! maximum string length for variable names
  INTEGER, PARAMETER :: VARNAME_LEN = 32

  ! Tracer metadata
  ! 
  ! Polymorphic type that contains tracer metadata according to the type 
  ! of tracer (aerosol, chemical, hydrometeor). 
  TYPE t_tracer_meta
    !
    LOGICAL :: lis_tracer         ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=VARNAME_LEN) :: name         ! Name of tracer
    CHARACTER(len=VARNAME_LEN) :: tracer_class ! Type of tracer
    ! Advection
    INTEGER :: ihadv_tracer       ! Method for horizontal transport
    INTEGER :: ivadv_tracer       ! Method for vertical transport
    ! Other processes covered by ICON
    LOGICAL :: lturb_tracer       ! Turbulent transport (TRUE/FALSE)
    LOGICAL :: lconv_tracer       ! Convection  (TRUE/FALSE)
    ! Processes not covered by ICON (requires ART extension)
    INTEGER :: ised_tracer        ! Sedimentation
                                  !   0 = No sedimentation
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode
                                  !   ... (hydrometeors to be added)
    LOGICAL :: ldep_tracer        ! Dry deposition  (TRUE/FALSE)
    INTEGER :: iwash_tracer       ! Washout
                                  !   0 = No washout
                                  !   1 = Monodisperse aerosol
                                  !   2 = As part of an according aerosol mode
    !
  END TYPE t_tracer_meta

  ! Aerosol-specific metadata
  TYPE, extends(t_tracer_meta) :: t_aero_meta
    !
    REAL(wp) :: solubility        ! Solubility, between 0 (insoluble) and 1 (soluble)
    REAL(wp) :: rho               ! Density [kg m-3]
    REAL(wp) :: mol_weight        ! Molar mass [g mol-1]
    !
  END TYPE

  ! Chemical tracer metadata
  TYPE, extends(t_tracer_meta) :: t_chem_meta
    !
    REAL(wp)    ::  lifetime_tracer ! Lifetime of tracer [s]
    REAL(wp)    ::  mol_weight      ! Molar mass [g mol-1]
    !
  END TYPE

  ! Hydrometeor metadata
  TYPE, extends(t_tracer_meta) :: t_hydro_meta
    ! 
  END TYPE

  PUBLIC :: t_tracer_meta
  PUBLIC :: t_aero_meta
  PUBLIC :: t_chem_meta
  PUBLIC :: t_hydro_meta

END MODULE mo_tracer_metadata_types
