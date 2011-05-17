MODULE mo_nwp_lnd_state
!>
!!  !MODULE:  mo_nwp_phy_state\\
!!
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the NWP land physics.
!!  Constructors and destructors for these data structures as well as
!!  initialization of fields are also defined here.
!!  This module should be an analogon to 'mo_hydro_state.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2010-11-09)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
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

! !USES:

USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
!!$USE mo_global_variables,    ONLY: l_nest_rcf
USE mo_dynamics_nml,        ONLY: nsav1, nsav2
USE mo_run_nml,             ONLY: nproma
USE mo_exception,           ONLY: message, finish, message_text
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom, l_limited_area
USE mo_nonhydrostatic_nml,  ONLY: l_nest_rcf
USE mo_atm_phy_nwp_nml,     ONLY: inwp_surface

USE mo_lnd_nwp_nml,         ONLY: nlev_soil,nztlev,nlev_snow, &
                                  nsfc_subs

USE mo_linked_list,         ONLY: t_var_list
USE mo_var_list,            ONLY: default_var_list_settings, &
                                & add_var,                   &
                                & new_var_list,              &
                                & delete_var_list
USE mo_cf_convention
USE mo_grib2
USE mo_cdi_constants 

IMPLICIT NONE
PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: construct_nwp_lnd_state
PUBLIC :: destruct_nwp_lnd_state
!
!variables
PUBLIC :: t_lnd_state  !> state vector  for land scheme
PUBLIC :: t_lnd_prog   !!       for prognostic variables
PUBLIC :: t_lnd_diag   !!       for diagnostic variables
!

!PUBLIC ::  lnd_prog_nwp_list,lnd_diag_nwp_list !< variable lists


! prognostic variables state vector
  TYPE t_lnd_prog

  REAL(wp), POINTER :: & ! JH old ALLOCATABLE
         t_snow       (:,:,:)     , & ! temperature of the snow-surface               (  K  )
         t_snow_mult  (:,:,:,:,:) , & ! temperature of the snow-surface               (  K  )
         t_s          (:,:,:,:)   , & ! temperature of the ground surface             (  K  )
         t_g          (:,:)       , & ! weighted surface temperature                  (  K  )
         w_snow       (:,:,:,:)   , & ! water content of snow                         (m H2O)
         rho_snow     (:,:,:)     , & ! snow density                                  (kg/m**3)
         rho_snow_mult(:,:,:,:,:)   , & ! snow density                                  (kg/m**3)
         w_i          (:,:,:,:)     , & ! water content of interception water           (m H2O)
         t_so         (:,:,:,:,:) , & ! soil temperature (main level)                 (  K  )
         w_so         (:,:,:,:,:) , & ! total water content (ice + liquid water)       (m H20)
         w_so_ice     (:,:,:,:,:) , & ! ice content                                   (m H20)
         dzh_snow     (:,:,:,:,:)      , & ! layer thickness between half levels in snow   (  m  )
         subsfrac     (:,:,:,:)            ! 
  END TYPE t_lnd_prog

  TYPE t_lnd_diag

  REAL(wp), POINTER :: & ! JH old ALLOCATABLE
      qv_s (:,:)       , & ! specific humidity at the surface              (kg/kg)
      h_snow(:,:,:,:)    , & ! snow height                                   (  m  
      freshsnow(:,:,:) , & ! indicator for age of snow in top of snow layer(  -  )
      wliq_snow(:,:,:,:) , & ! liquid water content in the snow              (m H2O)
      wtot_snow(:,:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
      runoff_s(:,:,:)  , & ! surface water runoff; sum over forecast      (kg/m2)
      runoff_g(:,:,:)  , & ! soil water runoff; sum over forecast         (kg/m2)
      rstom(:,:,:)     , & ! stomata resistance                           ( s/m )
      lhfl_bs(:,:,:)   , & ! average latent heat flux from bare soil evap.( W/m2)
      lhfl_pl(:,:,:)   , & ! average latent heat flux from plants         ( W/m2)
      fr_seaice(:,:)        !< fraction of sea ice                [ ]   
                            !< as partition of total area of the
                            !< grid element, but set to 0 or 1
                            !< index1=1,nproma, index2=1,nblks_c
  END TYPE t_lnd_diag


! complete state vector

!!--------------------------------------------------------------------------
!!                          VARIABLE LISTS
!! has to be declared to gether in order to allocate clean for the 
!! different domains and timelevels
!!--------------------------------------------------------------------------
  TYPE t_lnd_state
    TYPE(t_lnd_prog), ALLOCATABLE :: prog_lnd(:)
    TYPE(t_var_list), ALLOCATABLE :: lnd_prog_nwp_list(:)  !<  shape: (ntimelevel)
    TYPE(t_lnd_diag), ALLOCATABLE :: diag_lnd
    TYPE(t_var_list), ALLOCATABLE :: lnd_diag_nwp_list    !< shape: (n_dom)
  END TYPE t_lnd_state
 

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! It calls constructors to
  !! single time level prognostic states, and diagnostic states.
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE construct_nwp_lnd_state(p_patch, p_lnd_state, n_timelevels)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch(n_dom) ! patch
  INTEGER, OPTIONAL, INTENT(IN)         :: n_timelevels ! number of timelevels

  TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
  INTEGER     :: ntl, &! local number of timelevels
                 ist, &! status
                 jg,  &! grid level counter
                 jt,  &! time level counter
                nblks_c ! number of blocks of cells

  CHARACTER(len=MAX_CHAR_LENGTH) :: listname
!-----------------------------------------------------------------------

     write(0,*)'begin lnd_state time level=',ntl


  DO jg = 1, n_dom

     write(0,*)'lnd_state dom=',jg

    IF(PRESENT(n_timelevels))THEN
      ntl = n_timelevels
    ELSE
      ntl = 1
    ENDIF
   
    !
     !determine size of arrays
     nblks_c = p_patch(jg)%nblks_c

!!$    ! If grid nesting is not called at every dynamics time step, an extra time
!!$    ! level is needed for full-field interpolation and boundary-tendency calculation
!!$    IF (l_nest_rcf .AND. n_dom > 1) THEN
!!$      ntl = ntl + 1
!!$      nsav1(jg) = ntl
!!$    ENDIF

    ! In the presence of grid nesting, another extra time level is needed to save
    ! the feedback increments
    ! This extra time level is also used to store the driving-model data in the
    ! limited-area mode
    IF (l_limited_area .OR. jg > 1) ntl = ntl + 1
    nsav2(jg) = ntl

    !create state arrays and lists
    ALLOCATE(p_lnd_state(jg)%prog_lnd(1:ntl), &
         &   p_lnd_state(jg)%lnd_prog_nwp_list(1:ntl),STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nwp_lnd_state:construct_lnd_state', &
           'allocation of land prognostic state array failed')
    ENDIF

    ALLOCATE(p_lnd_state(jg)%diag_lnd,          &
         &   p_lnd_state(jg)%lnd_diag_nwp_list, STAT=ist)
    IF(ist/=SUCCESS)THEN
    CALL finish ('mo_nwp_lnd_state:construct_nwp_state', &
      'allocation of diagnostic physical array and list failed')
    ENDIF

    DO jt = 1, ntl
            !construct prognostic state
      WRITE(listname,'(a,i2.2,a,i2.2)') 'lnd_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

     write(0,*)'lnd_state timelevel=',jt,'and domain= ',jg

      CALL  new_nwp_lnd_prog_list(nblks_c, TRIM(listname),&
     &   p_lnd_state(jg)%lnd_prog_nwp_list(jt), p_lnd_state(jg)%prog_lnd(jt))
    ENDDO

!
     WRITE(listname,'(a,i2.2)') 'lnd_diag_of_domain_',jg
    !construct diagnostic state
    CALL new_nwp_lnd_diag_list( nblks_c, TRIM(listname),&
         &   p_lnd_state(jg)%lnd_diag_nwp_list,p_lnd_state(jg)%diag_lnd)

  ENDDO !ndom

  CALL message ('mo_lnd_state:construct_lnd_state', &
                'land state construction completed')

  END SUBROUTINE construct_nwp_lnd_state

!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE destruct_nwp_lnd_state(p_lnd_state, n_timelevels)
   !
    INTEGER, OPTIONAL, INTENT(IN)     :: n_timelevels ! number of timelevels
    TYPE(t_lnd_state),  INTENT(INOUT) :: & 
                            &  p_lnd_state(n_dom)      ! land state at different grid levels

    INTEGER :: ntl, &! local number of timelevels
               ist, &! status
                jg, &! grid level counter
                jt   ! time level counter

!-----------------------------------------------------------------------



  DO jg = 1, n_dom

    IF(PRESENT(n_timelevels))THEN
      ntl = n_timelevels
    ELSE
      ntl = 1
    ENDIF

    IF(ntl==0)THEN
      CALL finish('mo_land_state:destruct_lnd_state', &
                  'prognostic array has no timelevels')
    ENDIF
    !destruct diagnostic state
    CALL delete_var_list( p_lnd_state(jg)%lnd_diag_nwp_list )

    DO jt = 1, ntl
      !destruct prognostic state 
       CALL delete_var_list(p_lnd_state(jg)%lnd_prog_nwp_list(jt) )
    ENDDO

    DEALLOCATE(p_lnd_state(jg)%prog_lnd,&
         &     p_lnd_state(jg)%lnd_prog_nwp_list, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_lnd_state:destruct_lnd_state', &
             'deallocation of land prognostic state array failed')
      ENDIF

    DEALLOCATE(p_lnd_state(jg)%diag_lnd,&
    &          p_lnd_state(jg)%lnd_diag_nwp_list, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_lnd_state:destruct_lnd_state', &
             'deallocation of land diagnostic state array failed')
     ENDIF
  ENDDO

  END SUBROUTINE destruct_nwp_lnd_state

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE new_nwp_lnd_prog_list( kblks,   &
                     & listname, prog_list, p_prog_lnd)

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),INTENT(INOUT) :: prog_list
    TYPE(t_lnd_prog),INTENT(INOUT) :: p_prog_lnd

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape4d_subs(4)
    INTEGER :: shape5d_snow_subs(5), shape5d_soil_subs(5)
    INTEGER :: ientr

    ientr = 16 ! "entropy" of horizontal slice

    shape2d              = (/nproma, kblks   /)
    shape3d_subs         = (/nproma, nsfc_subs, kblks  /)
    shape4d_subs         = (/nproma, nztlev, nsfc_subs,  kblks /)
    shape5d_snow_subs    = (/nproma, nlev_snow, nztlev, nsfc_subs,  kblks /)
    shape5d_soil_subs    = (/nproma, nlev_soil, nztlev, nsfc_subs,  kblks /)


    ! Register a field list and apply default settings

    CALL new_var_list( prog_list, TRIM(listname) )
    CALL default_var_list_settings( prog_list,                 &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! & p_prog_lnd%t_g(nproma,nblks_c), STAT = ist)
    cf_desc    = t_cf_var('t_g ', 'K ', 'weighted surface temperature ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'qv_s', p_prog_lnd%t_g,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape2d )    


    ! & p_prog_lnd%t_snow(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('t_snow ', 'K ', 'temperature of the snow-surface ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 't_snow', p_prog_lnd%t_snow,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )   

    ! & p_prog_lnd%t_snow_mult(nproma,nlev_snow,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('t_snow_mult ', 'K ', 'temperature of the snow-surface ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 't_snow_mult', p_prog_lnd%t_snow_mult,                    &
     & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_snow_subs )  

    ! & p_prog_lnd%t_s(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('t_s ', 'K ', 'temperature of ground surface ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 't_s', p_prog_lnd%t_s,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs )  

    ! & p_prog_lnd%w_snow(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('w_snow ', 'm H2O ', 'water content of snow ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'w_snow', p_prog_lnd%w_snow,                             &
          & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs)

    ! & p_prog_lnd%rho_snow(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('rho_snow ', 'kg/m**3 ', 'snow density ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'rho_snow', p_prog_lnd%rho_snow,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )    


    ! & p_prog_lnd%rho_snow_mult(nproma,nlev_snow,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('rho_snow_mult ', 'kg/m**3 ', 'snow density ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'rho_snow_mult', p_prog_lnd%rho_snow_mult,                  &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_snow_subs)

    ! & p_prog_lnd%w_i(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('w_i ', 'm H2O ', 'water content of interception water ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'w_i', p_prog_lnd%w_i,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs )  

    ! & p_prog_lnd%t_so(nproma,nlev_soil,nztlev,nsfc_subs,nblks_c) 
    cf_desc    = t_cf_var('t_so ', 'K ', 'soil temperature (main level) ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 't_so', p_prog_lnd%t_so,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_soil_subs )  

   ! & p_prog_lnd%w_so(nproma,nlev_soil,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('w_so ', 'm H20 ', 'total water content (ice + liquid water) ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'w_so', p_prog_lnd%w_so,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_soil_subs )

    ! & p_prog_lnd%w_so_ice(nproma,nlev_soil,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('w_so_ice ', 'm H20 ', 'ice content ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'w_so_ice', p_prog_lnd%w_so_ice,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_soil_subs )

    ! & p_prog_lnd%dzh_snow(nproma,nlev_snow,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('dzh_snow ', 'm ', 'layer thickness between half levels in snow ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'dzh_snow', p_prog_lnd%dzh_snow,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape5d_snow_subs )
 
    ! & p_prog_lnd%subsfrac(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('subsfrac ', '- ', 'subscale fraction ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, 'subsfrac', p_prog_lnd%subsfrac,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs)

    p_prog_lnd%t_g(:,:)                 = 0.0_wp
    p_prog_lnd%t_snow(:,:,:)            = 0.0_wp
    p_prog_lnd%t_snow_mult(:,:,:,:,:)   = 0.0_wp
    p_prog_lnd%t_s(:,:,:,:)             = 0.0_wp
    p_prog_lnd%w_snow(:,:,:,:)          = 0.0_wp
    p_prog_lnd%rho_snow(:,:,:)          = 0.0_wp
    p_prog_lnd%rho_snow_mult(:,:,:,:,:) = 0.0_wp
    p_prog_lnd%w_i(:,:,:,:)             = 0.0_wp
    p_prog_lnd%t_so(:,:,:,:,:)          = 0.0_wp
    p_prog_lnd%w_so(:,:,:,:,:)          = 0.0_wp
    p_prog_lnd%w_so_ice(:,:,:,:,:)      = 0.0_wp
    p_prog_lnd%dzh_snow(:,:,:,:,:)      = 0.0_wp
    p_prog_lnd%subsfrac(:,:,:,:)        = 0.0_wp

END SUBROUTINE new_nwp_lnd_prog_list

SUBROUTINE new_nwp_lnd_diag_list( kblks, listname, diag_list, p_diag_lnd)

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),INTENT(INOUT) :: diag_list
    TYPE(t_lnd_diag),INTENT(INOUT) :: p_diag_lnd

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape4d_subs(4), shapesfc(3)
    INTEGER :: ientr

    ientr = 16 ! "entropy" of horizontal slice

    shape2d      = (/nproma, kblks         /)
    shape3d_subs = (/nproma, nsfc_subs,  kblks    /)
    shape4d_subs = (/nproma, nztlev, nsfc_subs,  kblks    /)

    ! Register a field list and apply default settings

    CALL new_var_list( diag_list, TRIM(listname) )
    CALL default_var_list_settings( diag_list,                 &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! & p_diag_lnd%qv_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_s ', 'kg/kg ', 'specific humidity at the surface ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'qv_s', p_diag_lnd%qv_s,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape2d )    

    ! & p_diag_lnd%h_snow(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('h_snow ', 'm ', 'snow height ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'h_snow', p_diag_lnd%h_snow,                             &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs )  

    ! & p_diag_lnd%freshsnow(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('freshsnow ', '- ', 'indicator for age of snow in top of snow layer ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'freshsnow', p_diag_lnd%freshsnow,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%wliq_snow(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('wliq_snow ', 'm H2O ', 'liquid water content in the snow ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'wliq_snow', p_diag_lnd%wliq_snow,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs )

       ! & p_diag_lnd%wtot_snow(nproma,nztlev,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('wtot_snow ', 'm H2O ', 'total water content in the snow ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'wtot_snow', p_diag_lnd%wtot_snow,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape4d_subs )

    ! & p_diag_lnd%runoff_s(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('runoff_s ', 'kg/m2 ', 'surface water runoff; sum over forecast ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'runoff_s', p_diag_lnd%runoff_s,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%runoff_g(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('runoff_g ', 'kg/m2 ', 'soil water runoff; sum over forecast ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'runoff_g', p_diag_lnd%runoff_g,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%rstom(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('rstom ', 's/m ', 'stomata resistance ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rstom', p_diag_lnd%rstom,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%lhfl_bs(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('lhfl_bs ', 'W/m2 ', 'average latent heat flux from bare soil evap. ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_bs', p_diag_lnd%lhfl_bs,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%lhfl_pl(nproma,nsfc_subs,nblks_c)
    cf_desc    = t_cf_var('lhfl_pl ', 'W/m2 ', 'average latent heat flux from plants ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_pl', p_diag_lnd%lhfl_pl,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape3d_subs )

    ! & p_diag_lnd%fr_seaice(nproma,nblks_c)
    cf_desc    = t_cf_var(' fr_seaice ', '- ', 'fraction of sea ice ')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'fr_seaice', p_diag_lnd%fr_seaice,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc, ldims=shape2d )    

  p_diag_lnd%qv_s(:,:)          = 0._wp
  p_diag_lnd%h_snow(:,:,:,:)    = 0._wp
  p_diag_lnd%freshsnow(:,:,:)   = 0._wp
  p_diag_lnd%wliq_snow(:,:,:,:) = 0._wp
  p_diag_lnd%wtot_snow(:,:,:,:) = 0._wp
  p_diag_lnd%runoff_s(:,:,:)    = 0._wp
  p_diag_lnd%runoff_g(:,:,:)    = 0._wp
  p_diag_lnd%rstom(:,:,:)       = 0._wp
  p_diag_lnd%lhfl_bs(:,:,:)     = 0._wp
  p_diag_lnd%lhfl_pl(:,:,:)     = 0._wp
  p_diag_lnd%fr_seaice(:,:)     = 0._wp

END SUBROUTINE  new_nwp_lnd_diag_list

!
!-------------------------------------------------------------------------

END MODULE mo_nwp_lnd_state
!<
