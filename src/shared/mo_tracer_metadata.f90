!>
!! Description...
!!
!! @author Daniel Rieger, KIT
!!
!!
!! @par Revision History
!! Initial revision by Daniel Rieger, KIT (2016-01-21)
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
                                  &   t_chem_meta, t_hydro_meta

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_tracer_metadata
  PUBLIC  :: create_tracer_metadata_aero
  PUBLIC  :: create_tracer_metadata_chem
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

  TYPE(t_tracer_meta) FUNCTION create_tracer_metadata(lis_tracer, name, ihadv_tracer, ivadv_tracer,         &
                      &                               lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer, &
                      &                               iwash_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout
    
    ! Fill the meta of the base type (t_tracer_meta)
    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      create_tracer_metadata%lis_tracer = lis_tracer
    ELSE
      create_tracer_metadata%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      create_tracer_metadata%name = TRIM(name)
    ELSE
      create_tracer_metadata%name = "unnamed"
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      create_tracer_metadata%ihadv_tracer = ihadv_tracer
    ELSE
      create_tracer_metadata%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      create_tracer_metadata%ivadv_tracer = ivadv_tracer
    ELSE
      create_tracer_metadata%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      create_tracer_metadata%lturb_tracer = lturb_tracer
    ELSE
      create_tracer_metadata%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      create_tracer_metadata%lconv_tracer = lconv_tracer
    ELSE
      create_tracer_metadata%lconv_tracer = .FALSE.
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      create_tracer_metadata%ised_tracer = ised_tracer
    ELSE
      create_tracer_metadata%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      create_tracer_metadata%ldep_tracer = ldep_tracer
    ELSE
      create_tracer_metadata%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      create_tracer_metadata%iwash_tracer = iwash_tracer
    ELSE
      create_tracer_metadata%iwash_tracer = 0
    ENDIF
    
  END FUNCTION create_tracer_metadata



  TYPE(t_aero_meta) FUNCTION create_tracer_metadata_aero(lis_tracer, name, ihadv_tracer, ivadv_tracer,          &
                      &                                  lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer,  &
                      &                                  iwash_tracer,                                          &
                      &                                  solubility, rho, mol_weight)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout
    ! Extended type (t_aero_meta) content
    REAL(wp), INTENT(IN), OPTIONAL :: solubility        ! Solubility, between 0 (insoluble) and 1 (soluble)
    REAL(wp), INTENT(IN), OPTIONAL :: rho               ! Density [kg m-3]
    REAL(wp), INTENT(IN), OPTIONAL :: mol_weight        ! Molar mass [g mol-1]
    
    ! Fill the meta of the base type (t_tracer_meta)
    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      create_tracer_metadata_aero%lis_tracer = lis_tracer
    ELSE
      create_tracer_metadata_aero%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      create_tracer_metadata_aero%name = TRIM(name)
    ELSE
      create_tracer_metadata_aero%name = "unnamed"
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      create_tracer_metadata_aero%ihadv_tracer = ihadv_tracer
    ELSE
      create_tracer_metadata_aero%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      create_tracer_metadata_aero%ivadv_tracer = ivadv_tracer
    ELSE
      create_tracer_metadata_aero%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      create_tracer_metadata_aero%lturb_tracer = lturb_tracer
    ELSE
      create_tracer_metadata_aero%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      create_tracer_metadata_aero%lconv_tracer = lconv_tracer
    ELSE
      create_tracer_metadata_aero%lconv_tracer = .FALSE.
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
      
    ! Fill the meta of the extended type (t_aero_meta)
    IF(PRESENT(solubility)) THEN
      create_tracer_metadata_aero%solubility = solubility
    ELSE
      create_tracer_metadata_aero%solubility = -1._wp
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
    
  END FUNCTION create_tracer_metadata_aero



  TYPE(t_chem_meta) FUNCTION create_tracer_metadata_chem(lis_tracer, name, ihadv_tracer, ivadv_tracer,          &
                      &                                  lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer,  &
                      &                                  iwash_tracer,                                          &
                      &                                  lifetime_tracer, mol_weight)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout
    ! Extended type (t_chem_meta) content
    REAL(wp), INTENT(IN), OPTIONAL ::  lifetime_tracer ! Lifetime of tracer [s]
    REAL(wp), INTENT(IN), OPTIONAL ::  mol_weight      ! Molar mass [g mol-1]
    
    ! Fill the meta of the base type (t_tracer_meta)
    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      create_tracer_metadata_chem%lis_tracer = lis_tracer
    ELSE
      create_tracer_metadata_chem%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      create_tracer_metadata_chem%name = TRIM(name)
    ELSE
      create_tracer_metadata_chem%name = "unnamed"
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      create_tracer_metadata_chem%ihadv_tracer = ihadv_tracer
    ELSE
      create_tracer_metadata_chem%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      create_tracer_metadata_chem%ivadv_tracer = ivadv_tracer
    ELSE
      create_tracer_metadata_chem%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      create_tracer_metadata_chem%lturb_tracer = lturb_tracer
    ELSE
      create_tracer_metadata_chem%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      create_tracer_metadata_chem%lconv_tracer = lconv_tracer
    ELSE
      create_tracer_metadata_chem%lconv_tracer = .FALSE.
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
    
    ! Fill the meta of the extended type (t_chem_meta)
    IF(PRESENT(lifetime_tracer)) THEN
      create_tracer_metadata_chem%lifetime_tracer = lifetime_tracer
    ELSE
      create_tracer_metadata_chem%lifetime_tracer = -1._wp
    ENDIF
    
    IF(PRESENT(mol_weight)) THEN
      create_tracer_metadata_chem%mol_weight = mol_weight
    ELSE
      create_tracer_metadata_chem%mol_weight = -1._wp
    ENDIF
    
  END FUNCTION create_tracer_metadata_chem



  TYPE(t_hydro_meta) FUNCTION create_tracer_metadata_hydro(lis_tracer, name, ihadv_tracer, ivadv_tracer,       &
                      &                                  lturb_tracer, lconv_tracer, ised_tracer, ldep_tracer, &
                      &                                  iwash_tracer)
    ! Base type (t_tracer_meta) content
    LOGICAL, INTENT(IN), OPTIONAL  :: lis_tracer       ! this is a tracer field (TRUE/FALSE)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: name      ! Name of tracer
    INTEGER, INTENT(IN), OPTIONAL  :: ihadv_tracer     ! Method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL  :: ivadv_tracer     ! Method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lturb_tracer     ! Switch for turbulent transport
    LOGICAL, INTENT(IN), OPTIONAL  :: lconv_tracer     ! Switch for convection
    INTEGER, INTENT(IN), OPTIONAL  :: ised_tracer      ! Method for sedimentation
    LOGICAL, INTENT(IN), OPTIONAL  :: ldep_tracer      ! Switch for dry deposition
    INTEGER, INTENT(IN), OPTIONAL  :: iwash_tracer     ! Method for washout
    ! Extended type (t_hydro_meta) content
    ! ...
    
    ! Fill the meta of the base type (t_tracer_meta)
    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      create_tracer_metadata_hydro%lis_tracer = lis_tracer
    ELSE
      create_tracer_metadata_hydro%lis_tracer = .FALSE.
    ENDIF

    ! name
    IF ( PRESENT(name) ) THEN
      create_tracer_metadata_hydro%name = TRIM(name)
    ELSE
      create_tracer_metadata_hydro%name = "unnamed"
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      create_tracer_metadata_hydro%ihadv_tracer = ihadv_tracer
    ELSE
      create_tracer_metadata_hydro%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      create_tracer_metadata_hydro%ivadv_tracer = ivadv_tracer
    ELSE
      create_tracer_metadata_hydro%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer
    IF ( PRESENT(lturb_tracer) ) THEN
      create_tracer_metadata_hydro%lturb_tracer = lturb_tracer
    ELSE
      create_tracer_metadata_hydro%lturb_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      create_tracer_metadata_hydro%lconv_tracer = lconv_tracer
    ELSE
      create_tracer_metadata_hydro%lconv_tracer = .FALSE.
    ENDIF

    ! ised_tracer
    IF ( PRESENT(ised_tracer) ) THEN
      create_tracer_metadata_hydro%ised_tracer = ised_tracer
    ELSE
      create_tracer_metadata_hydro%ised_tracer = 0
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      create_tracer_metadata_hydro%ldep_tracer = ldep_tracer
    ELSE
      create_tracer_metadata_hydro%ldep_tracer = .FALSE.
    ENDIF

    ! iwash_tracer
    IF ( PRESENT(iwash_tracer) ) THEN
      create_tracer_metadata_hydro%iwash_tracer = iwash_tracer
    ELSE
      create_tracer_metadata_hydro%iwash_tracer = 0
    ENDIF
    
    ! Fill the meta of the extended type (t_hydro_meta)
    ! ...
    
  END FUNCTION create_tracer_metadata_hydro

END MODULE mo_tracer_metadata

