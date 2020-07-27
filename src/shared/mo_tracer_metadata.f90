!>
!! Description: 
!! This module contains the constructors for polymorphic tracer metadata.
!! For each extended object, a constructor for the base type is called 
!! via a type bound procedure. 
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
MODULE mo_tracer_metadata

  USE mo_kind,                  ONLY: wp
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta, t_aero_meta, &
                                  &   t_passive_meta,             &
                                  &   t_chem_meta, t_hydro_meta

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_tracer_metadata
  PUBLIC  :: create_tracer_metadata_aero
  PUBLIC  :: create_tracer_metadata_chem
  PUBLIC  :: create_tracer_metadata_passive
  PUBLIC  :: create_tracer_metadata_hydro

  !------------------------------------------------------------------------------------------------
  !>create_tracer_metadata
  ! Quasi-constructors for tracer metadata
  ! (public routine. Can be used for two things:
  ! 1.) default settings: If used without any argument, the routines return a variable
  !     of the according extended type, containing the default settings.
  ! 2.) Setting of metadata: If used with arguments, the routines return a variable
  !     of the according extended type, containing the default settings except for those components
  !     which are given in the argument list.

CONTAINS

  TYPE(t_tracer_meta) FUNCTION create_tracer_metadata(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                               lturb_tracer, lconv_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name       ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection

    ! Fill the metadata of the base type
    CALL create_tracer_metadata%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                        lturb_tracer, lconv_tracer)

  END FUNCTION create_tracer_metadata



  TYPE(t_aero_meta) FUNCTION create_tracer_metadata_aero(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                                  lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer,  &
                      &                                  iwash_tracer,                                          &
                      &                                  moment, mode, rho, mol_weight)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer      ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer    ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer    ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer    ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer    ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer     ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer     ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer    ! Method for washout
    ! Extended type (t_aero_meta) content
    INTEGER, INTENT(IN), OPTIONAL  :: moment          ! moment of distribution (e.g. 0=number, 3=proportional to mass)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: &
      &                               mode             ! Name of mode the tracer is contained in
    REAL(wp), INTENT(IN), OPTIONAL :: rho              ! Density [kg m-3]
    REAL(wp), INTENT(IN), OPTIONAL :: mol_weight       ! Molar mass [g mol-1]

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_aero%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                             lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(moment)) THEN
      create_tracer_metadata_aero%moment = moment
    ELSE
      create_tracer_metadata_aero%moment = -1
    ENDIF

    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(mode)) THEN
      create_tracer_metadata_aero%mode = TRIM(mode)
    ELSE
      create_tracer_metadata_aero%mode = 'no_mode'
    ENDIF

    IF(PRESENT(rho)) THEN
      create_tracer_metadata_aero%rho = rho
    ELSE
      create_tracer_metadata_aero%rho = -1._wp
    ENDIF

    IF(PRESENT(mol_weight)) THEN
      create_tracer_metadata_aero%mol_weight = mol_weight
    ELSE
      create_tracer_metadata_aero%mol_weight = -1._wp
    ENDIF

    ! Fill the meta of the extended type (t_chem_meta)
    IF(PRESENT(mol_weight)) THEN
      create_tracer_metadata_aero%mol_weight = mol_weight
    ELSE
      create_tracer_metadata_aero%mol_weight = -1._wp
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      create_tracer_metadata_aero%ised_tracer = ised_tracer
    ELSE
      create_tracer_metadata_aero%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      create_tracer_metadata_aero%ldep_tracer = ldep_tracer
    ELSE
      create_tracer_metadata_aero%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      create_tracer_metadata_aero%iwash_tracer = iwash_tracer
    ELSE
      create_tracer_metadata_aero%iwash_tracer = 0
    ENDIF

  END FUNCTION create_tracer_metadata_aero



  TYPE(t_chem_meta) FUNCTION create_tracer_metadata_chem(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                                  lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer,  &
                      &                                  iwash_tracer,  mol_weight)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name       ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout
    ! Extended type (t_chem_meta) content
    REAL(wp), INTENT(IN), OPTIONAL :: mol_weight       ! Molar mass [kg mol-1]

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_chem%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                             lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_chem_meta)
    IF(PRESENT(mol_weight)) THEN
      create_tracer_metadata_chem%mol_weight = mol_weight
    ELSE
      create_tracer_metadata_chem%mol_weight = -1._wp
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      create_tracer_metadata_chem%ised_tracer = ised_tracer
    ELSE
      create_tracer_metadata_chem%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      create_tracer_metadata_chem%ldep_tracer = ldep_tracer
    ELSE
      create_tracer_metadata_chem%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      create_tracer_metadata_chem%iwash_tracer = iwash_tracer
    ELSE
      create_tracer_metadata_chem%iwash_tracer = 0
    ENDIF


  END FUNCTION create_tracer_metadata_chem



  TYPE(t_passive_meta) FUNCTION create_tracer_metadata_passive(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                                  lturb_tracer, lconv_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name       ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    ! Extended type (t_passive_meta) content
    ! ...

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_passive%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                              lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_passive_meta)
    ! ...

  END FUNCTION create_tracer_metadata_passive


  TYPE(t_hydro_meta) FUNCTION create_tracer_metadata_hydro(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer, &
                      &                                  lturb_tracer, lconv_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    LOGICAL, INTENT(IN), OPTIONAL  :: lfeedback        ! feedback from child- to parent domain
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    ! Extended type (t_hydro_meta) content
    ! ...

    ! Fill the metadata of the base type
    CALL create_tracer_metadata_hydro%construct_base(lis_tracer, name, lfeedback, ihadv_tracer, ivadv_tracer,  &
      &                                              lturb_tracer, lconv_tracer)

    ! Fill the meta of the extended type (t_hydro_meta)
    ! ...

  END FUNCTION create_tracer_metadata_hydro

END MODULE mo_tracer_metadata

