!>
!! Description: 
!! This module contains types for polymorphic tracer metadata. Additionally,
!! a type bound procedure (TBP) is given as a constructor for the base type 
!! t_tracer_meta. This procedure is private and may only be adressed as TBP. 
!!
!! @author Daniel Rieger, KIT
!!
!!
!! @par Revision History
!! Initial revision by Daniel Rieger, KIT (2016-01-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
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
    CONTAINS
      procedure :: construct_base => construct_t_tracer_meta
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

CONTAINS

  SUBROUTINE construct_t_tracer_meta(meta, lis_tracer, name, ihadv_tracer, ivadv_tracer,   &
    &                                lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer, &
    &                                iwash_tracer)

    CLASS(t_tracer_meta),INTENT(OUT)     :: meta          ! Base meta container to be filled
    LOGICAL, INTENT(IN), OPTIONAL        :: lis_tracer    ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name          ! Name of tracer
    INTEGER, INTENT(IN), OPTIONAL        :: ihadv_tracer  ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL        :: ivadv_tracer  ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL        :: lturb_tracer  ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL        :: lconv_tracer  ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL        :: ised_tracer   ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL        :: ldep_tracer   ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL        :: iwash_tracer  ! Method for washout

    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      meta%lis_tracer = lis_tracer
    ELSE
      meta%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      meta%name = TRIM(name)
    ELSE
      meta%name = "unnamed"
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      meta%ihadv_tracer = ihadv_tracer
    ELSE
      meta%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      meta%ivadv_tracer = ivadv_tracer
    ELSE
      meta%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      meta%lturb_tracer = lturb_tracer
    ELSE
      meta%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      meta%lconv_tracer = lconv_tracer
    ELSE
      meta%lconv_tracer = .FALSE.
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      meta%ised_tracer = ised_tracer
    ELSE
      meta%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      meta%ldep_tracer = ldep_tracer
    ELSE
      meta%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      meta%iwash_tracer = iwash_tracer
    ELSE
      meta%iwash_tracer = 0
    ENDIF

  END SUBROUTINE construct_t_tracer_meta

END MODULE mo_tracer_metadata_types
