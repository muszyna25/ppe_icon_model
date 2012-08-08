#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
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
!! Modification by Daniel Reinert, DWD (2012-04-03)
!! - encapsulated type definitions (mo_nwp_lnd_types)
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
MODULE mo_nwp_lnd_state

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_parallel_config,      ONLY: nproma
  USE mo_nwp_lnd_types,        ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag 
  USE mo_exception,            ONLY: message, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_grid_config,          ONLY: n_dom
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, ntiles_total, &
    &                                lmulti_snow, ntiles_water
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_var_list,             ONLY: default_var_list_settings, &
    &                               add_var, add_ref,           &
    &                               new_var_list,               &
    &                               delete_var_list, groups
  USE mo_var_metadata,         ONLY: t_var_metadata
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var
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
  PUBLIC :: p_lnd_state  !> state vector (variable) for land scheme


! complete state vector
!
  TYPE(t_lnd_state), TARGET, ALLOCATABLE :: p_lnd_state(:) 

  CONTAINS


!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
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
    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch(n_dom) ! patch
    INTEGER, OPTIONAL, INTENT(IN)       :: n_timelevels   ! number of timelevels

    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &        jg,  &! grid level counter
      &        jt,  &! time level counter
      &      nblks_c ! number of blocks of cells

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:construct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'land state construction started')

    DO jg = 1, n_dom


      IF(PRESENT(n_timelevels))THEN
        ntl = n_timelevels
      ELSE
        ntl = 1
      ENDIF
   
      !
      !determine size of arrays
      nblks_c = p_patch(jg)%nblks_c


      !
      !create state arrays and lists
      !
      ALLOCATE(p_lnd_state(jg)%prog_lnd(1:ntl), &
           &   p_lnd_state(jg)%lnd_prog_nwp_list(1:ntl),STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nwp_lnd_state:construct_lnd_state', &
             'allocation of land prognostic state array failed')
      ENDIF


      !
      !construct prognostic state
      !
      DO jt = 1, ntl

        WRITE(listname,'(a,i2.2,a,i2.2)') 'lnd_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        !varname_prefix = 'lnd_prog_'
        varname_prefix = ''
        CALL  new_nwp_lnd_prog_list(jg, nblks_c, TRIM(listname),             &
          &     TRIM(varname_prefix), p_lnd_state(jg)%lnd_prog_nwp_list(jt), &
          &     p_lnd_state(jg)%prog_lnd(jt), jt)

      ENDDO

      !
      !construct diagnostic state
      !
      WRITE(listname,'(a,i2.2)') 'lnd_diag_of_domain_',jg

      !varname_prefix = 'lnd_diag_'
      varname_prefix = ''
      CALL new_nwp_lnd_diag_list(jg, nblks_c, TRIM(listname),        &
        &   TRIM(varname_prefix), p_lnd_state(jg)%lnd_diag_nwp_list, &
        &   p_lnd_state(jg)%diag_lnd)

    ENDDO !ndom

    CALL message (TRIM(routine), 'land state construction completed')

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
  SUBROUTINE destruct_nwp_lnd_state(p_lnd_state)
   !
    TYPE(t_lnd_state),  INTENT(INOUT) :: & 
      &   p_lnd_state(n_dom)             ! land state at different grid levels

    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &         jg, &! grid level counter
      &         jt   ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:destruct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of nwp land state started')

    DO jg = 1, n_dom

      ntl = SIZE(p_lnd_state(jg)%prog_lnd(:))
      IF(ntl==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF


      DO jt = 1, ntl
        ! delete prognostic state list elements
        CALL delete_var_list(p_lnd_state(jg)%lnd_prog_nwp_list(jt) )
      ENDDO


      ! delete diagnostic state list elements
      CALL delete_var_list( p_lnd_state(jg)%lnd_diag_nwp_list )


      ! destruct state lists and arrays
      DEALLOCATE(p_lnd_state(jg)%prog_lnd,&
           &     p_lnd_state(jg)%lnd_prog_nwp_list, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),  &
             'deallocation of land prognostic state array failed')
      ENDIF

    ENDDO

    CALL message (TRIM(routine), 'Destruction of nwp land state completed')

  END SUBROUTINE destruct_nwp_lnd_state

!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic land state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by , MPI-M (2011-07-01)
  !!
  !!
  SUBROUTINE new_nwp_lnd_prog_list( p_jg, kblks, listname, vname_prefix, &
    &                               prog_list, p_prog_lnd, timelev )

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg  !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list),INTENT(INOUT) :: prog_list
    TYPE(t_lnd_prog),INTENT(INOUT) :: p_prog_lnd

    INTEGER, INTENT(IN) :: timelev

    ! Local variables
    !
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: shape4d_snow_subs(4), shape4d_soil_subs(4)
    INTEGER :: ibits
    INTEGER :: jsfc          !< tile counter

    CHARACTER(len=4) suffix

!-----------------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d              = (/nproma,            kblks            /)
    shape3d_subs         = (/nproma,            kblks, ntiles_total /)
    shape3d_subsw        = (/nproma,            kblks, ntiles_total+ntiles_water /)
    shape4d_snow_subs    = (/nproma, nlev_snow, kblks, ntiles_total /)
    shape4d_soil_subs    = (/nproma, nlev_soil, kblks, ntiles_total /)


    ! Suffix (mandatory for time level dependent variables)

    WRITE(suffix,'(".TL",i1)') timelev

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( prog_list, TRIM(listname), patch_id=p_jg )
    CALL default_var_list_settings( prog_list,                 &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------

    ! & p_prog_lnd%t_g(nproma,nblks_c), STAT = ist)
    cf_desc    = t_cf_var('t_g', 'K', 'weighted surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g'//suffix, p_prog_lnd%t_g,      &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
         & ldims=shape2d,                                                      &
         & tlev_source=1,                                                      &
         & in_group=groups("land_vars") ) ! for output take field from nnow_rcf slice


    ! & p_prog_lnd%w_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('w_snow_t', 'm H2O', 'water equivalent of snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 60, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_snow_t'//suffix, p_prog_lnd%w_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,            &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_snow
    ALLOCATE(p_prog_lnd%w_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_snow_t'//suffix,           &
           & vname_prefix//'w_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
           & t_cf_var('w_snow_t_'//csfc, '', '', DATATYPE_FLT32),          &
           & t_grib2_var(0, 1, 60, ibits, GRID_REFERENCE, GRID_CELL),      &
           & ldims=shape2d,                                                &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO


    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0 ) THEN

    ! & p_prog_lnd%t_g_t(nproma,nblks_c,ntiles_total+ntiles_water), STAT = ist)
    cf_desc    = t_cf_var('t_g_t', 'K', 'weighted surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g_t'//suffix, p_prog_lnd%t_g_t,  &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,        &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  


    ! fill the separate variables belonging to the container t_gt
    ALLOCATE(p_prog_lnd%t_gt_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'t_g_t'//suffix,                &
               & vname_prefix//'t_g_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
               & p_prog_lnd%t_gt_ptr(jsfc)%p_2d,                               &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('t_g_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),       &
               & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
               & ldims=shape2d,                                                &
               & tlev_source=1 ) ! for output take field from nnow_rcf slice
      ENDDO



    ! & p_prog_lnd%t_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('t_snow_t', 'K', 'temperature of the snow-surface', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_t'//suffix, p_prog_lnd%t_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,            &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. ) 

    ! fill the separate variables belonging to the container t_snow
    ALLOCATE(p_prog_lnd%t_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'t_snow_t'//suffix,             &
               & vname_prefix//'t_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
               & p_prog_lnd%t_snow_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('t_snow_t_'//csfc, '', '', DATATYPE_FLT32),          &
               & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
               & ldims=shape2d,                                                &
               & tlev_source=1 ) ! for output take field from nnow_rcf slice
      ENDDO



    ! & p_prog_lnd%t_snow_mult_t(nproma,nlev_snow+1,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('t_snow_mult_t', 'K', 'temperature of the snow-surface', &
         & DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_mult_t'//suffix, p_prog_lnd%t_snow_mult_t, &
     & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, cf_desc, grib2_desc,           &
     & ldims=(/nproma,nlev_snow+1,kblks,ntiles_total/),                         &
     & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. ) 

    IF (lmulti_snow) THEN
      ! fill the separate variables belonging to the container t_snow_mult
      !
      ALLOCATE(p_prog_lnd%t_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'t_snow_mult_t'//suffix,      &
             & vname_prefix//'t_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,  &
             & p_prog_lnd%t_snow_mult_ptr(jsfc)%p_3d,                        &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,                        &
             & t_cf_var('t_snow_mult_t_'//csfc, '', '', DATATYPE_FLT32),     &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
             & ldims=(/nproma,nlev_snow+1,kblks/), lrestart=.TRUE.,          &
             & tlev_source=1 ) ! for output take field from nnow_rcf slice 
      ENDDO
    ENDIF 



    ! & p_prog_lnd%t_s_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('t_s_t', 'K', 'temperature of ground surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_s_t'//suffix, p_prog_lnd%t_s_t,  &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,         &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container t_s
    ALLOCATE(p_prog_lnd%t_s_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_s_t'//suffix,              &
           & vname_prefix//'t_s_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
           & p_prog_lnd%t_s_ptr(jsfc)%p_2d,                                &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
           & t_cf_var('t_s_t_'//csfc, '', '', DATATYPE_FLT32),             &
           & t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL),       &
           & ldims=shape2d,                                                &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO


    ! & p_prog_lnd%rho_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rho_snow_t', 'kg/m**3', 'snow density', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'rho_snow_t'//suffix, p_prog_lnd%rho_snow_t, &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,                 &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )
    

    ! fill the separate variables belonging to the container rho_snow
    ALLOCATE(p_prog_lnd%rho_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( prog_list, vname_prefix//'rho_snow_t'//suffix,           &
           & vname_prefix//'rho_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%rho_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                          &
           & t_cf_var('rho_snow_t_'//csfc, '', '', DATATYPE_FLT32),          &
           & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
           & ldims=shape2d,                                                  &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    END DO




    ! & p_prog_lnd%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rho_snow_mult_t', 'kg/m**3', 'snow density', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,              &
         & p_prog_lnd%rho_snow_mult_t, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,      &
         & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                           &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    IF (lmulti_snow) THEN
      ! fill the separate variables belonging to the container rho_snow_mult
      !
      ALLOCATE(p_prog_lnd%rho_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,       &
             & vname_prefix//'rho_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,   &
             & p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_3d,                         &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,                           &
             & t_cf_var('rho_snow_mult_t_'//csfc, '', '', DATATYPE_FLT32),      &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),    &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,               &
             & tlev_source=1 ) ! for output take field from nnow_rcf slice 
      ENDDO
    ENDIF



    ! & p_prog_lnd%w_i_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('w_i_t', 'm H2O', 'water content of interception water', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_i_t'//suffix, p_prog_lnd%w_i_t,     &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container w_i
    ALLOCATE(p_prog_lnd%w_i_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_i_t'//suffix,                &
           & vname_prefix//'w_i_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
           & p_prog_lnd%w_i_ptr(jsfc)%p_2d,                                  &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                          &
           & t_cf_var('w_i_t_'//csfc, '', '', DATATYPE_FLT32),               &
           & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
           & ldims=shape2d,                                                  &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO




    ! & p_prog_lnd%t_so_t(nproma,nlev_soil+2,nblks_c,ntiles_total) 
    cf_desc    = t_cf_var('t_so_t', 'K', 'soil temperature (main level)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_so_t'//suffix, p_prog_lnd%t_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+2,kblks,ntiles_total/),                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container t_so
    !
    ALLOCATE(p_prog_lnd%t_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_so_t'//suffix,               &
           & vname_prefix//'t_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%t_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND,                 &
           & t_cf_var('t_so_t_'//csfc, '', '', DATATYPE_FLT32),              &
           & t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL),         &
           & ldims=(/nproma,nlev_soil+2,kblks/),                             &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO



   ! & p_prog_lnd%w_so_t(nproma,nlev_soil+1,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('w_so_t', 'm H20', 'total water content (ice + liquid water)', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 3, 20, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_t'//suffix, p_prog_lnd%w_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks,ntiles_total/),                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_so
    ALLOCATE(p_prog_lnd%w_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_t'//suffix,               &
           & vname_prefix//'w_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%w_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND,                 &
           & t_cf_var('w_so_t_'//csfc, '', '', DATATYPE_FLT32),              &
           & t_grib2_var(2, 3, 20, ibits, GRID_REFERENCE, GRID_CELL),        &
           & ldims=(/nproma,nlev_soil+1,kblks/),                             &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO



    ! & p_prog_lnd%w_so_ice_t(nproma,nlev_soil+1,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('w_so_ice_t', 'm H20', 'ice content', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 3, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_ice_t'//suffix,                   &
         & p_prog_lnd%w_so_ice_t, GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND,  &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil+1,kblks,ntiles_total/),      &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_so_ice
    ALLOCATE(p_prog_lnd%w_so_ice_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_ice_t'//suffix,           &
           & vname_prefix//'w_so_ice_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_so_ice_ptr(jsfc)%p_3d,                             &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND,                 &
           & t_cf_var('w_so_ice_t_'//csfc, '', '', DATATYPE_FLT32),          &
           & t_grib2_var(2, 3, 22, ibits, GRID_REFERENCE, GRID_CELL),        &
           & ldims=(/nproma,nlev_soil+1,kblks/),                             &
           & tlev_source=1 ) ! for output take field from nnow_rcf slice
    ENDDO


    ! & p_prog_lnd%wliq_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('wliq_snow_t', 'm H2O', 'liquid water content in snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'wliq_snow_t'//suffix,                &
         & p_prog_lnd%wliq_snow_t, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,        & 
         & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    IF (lmulti_snow) THEN
      ! fill the separate variables belonging to the container wliq_snow
      !
      ALLOCATE(p_prog_lnd%wliq_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wliq_snow_t'//suffix,          &
             & vname_prefix//'wliq_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wliq_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,                          &
             & t_cf_var('wliq_snow_t_'//csfc, '', '', DATATYPE_FLT32),         &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & tlev_source=1 ) ! for output take field from nnow_rcf slice
      ENDDO
    ENDIF




    ! & p_prog_lnd%wtot_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('wtot_snow_t', 'm H2O', 'total water content in snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'wtot_snow_t'//suffix,                &
         & p_prog_lnd%wtot_snow_t, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,        &
         & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    IF (lmulti_snow) THEN
      ! fill the separate variables belonging to the container wtot_snow
      !
      ALLOCATE(p_prog_lnd%wtot_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wtot_snow_t'//suffix,          &
             & vname_prefix//'wtot_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wtot_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,                          &
             & t_cf_var('wtot_snow_t_'//csfc, '', '', DATATYPE_FLT32),         &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & tlev_source=1 ) ! for output take field from nnow_rcf slice
      ENDDO
    ENDIF



    ! & p_prog_lnd%dzh_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('dzh_snow_t', 'm', 'layer thickness between half levels in snow', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 14, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'dzh_snow_t'//suffix,                 &
         & p_prog_lnd%dzh_snow_t, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,         &
         & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    IF (lmulti_snow) THEN
      ! fill the separate variables belonging to the container dzh_snow
      !
      ALLOCATE(p_prog_lnd%dzh_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'dzh_snow_t'//suffix,             &
               & vname_prefix//'dzh_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
               & p_prog_lnd%dzh_snow_ptr(jsfc)%p_3d,                             &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,                          &
               & t_cf_var('dzh_snow_t_'//csfc, '', '', DATATYPE_FLT32),          &
               & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
               & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
               & tlev_source=1 ) ! for output take field from nnow_rcf slice
      ENDDO
    ENDIF
 

    p_prog_lnd%t_g_t(:,:,:)   = 290.4_wp

    p_prog_lnd%t_so_t(:,1,:,:) = 290.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,2,:,:) = 290.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,3,:,:) = 290.7_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,4,:,:) = 291.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,5,:,:) = 292.7_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,6,:,:) = 293.2_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,7,:,:) = 291.1_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,8,:,:) = 283.1_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,9,:,:) = 282.1_wp !!HW: be careful about the indices!!

    p_prog_lnd%w_so_t(:,1,:,:) = 1.8E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,2,:,:) = 3.7E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,3,:,:) = 11.3E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,4,:,:) = 35.1E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,5,:,:) = 56.7E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,6,:,:) = 254.6E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,7,:,:) = 763.9E-3_wp*2._wp !! JH
    p_prog_lnd%w_so_t(:,8,:,:) = 2291.5E-3_wp*2._wp !! JH
 
    p_prog_lnd%w_so_ice_t(:,1,:,:) = 0._wp
    p_prog_lnd%w_so_ice_t(:,2,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,3,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,4,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,5,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,6,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,7,:,:) = 0._wp 
    p_prog_lnd%w_so_ice_t(:,8,:,:) = 0._wp 
    END IF !inwp_surface > 0

  END SUBROUTINE new_nwp_lnd_prog_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic land state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by , MPI-M (2011-07-01)
  !!
  !!
  SUBROUTINE new_nwp_lnd_diag_list( p_jg, kblks, listname, vname_prefix, &
    &                               diag_list, p_diag_lnd)

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list),INTENT(INOUT) :: diag_list
    TYPE(t_lnd_diag),INTENT(INOUT) :: p_diag_lnd


    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: ibits
    INTEGER :: jsfc          !< tile counter

!--------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d       = (/nproma, kblks            /)
    shape3d_subs  = (/nproma, kblks, ntiles_total /)
    shape3d_subsw = (/nproma, kblks, ntiles_total+ntiles_water /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( diag_list, TRIM(listname), patch_id=p_jg )
    CALL default_var_list_settings( diag_list,                 &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    !------------------------------

    ! & p_diag_lnd%qv_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_s', 'kg/kg', 'specific humidity at the surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'qv_s', p_diag_lnd%qv_s,        &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
           & ldims=shape2d,                                                &
           & initval_r=0.001_wp,                                           &
           & in_group=groups("land_vars") )      

    ! & p_diag_lnd%fr_seaice(nproma,nblks_c)
    ! check restart settings again, once we have a sea ice model
    cf_desc    = t_cf_var('fr_seaice', '-', 'fraction of sea ice', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'fr_seaice', p_diag_lnd%fr_seaice,  &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.FALSE.,                                  &
           & initval_r=0._wp ) 



    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0) THEN

    ! & p_diag_lnd%qv_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('qv_s_t', 'kg/kg', 'specific humidity at the surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'qv_s_t', p_diag_lnd%qv_s_t,        &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE.,         &
           & loutput=.FALSE.,                                                  &
           & initval_r=0.001_wp )

    ! fill the separate variables belonging to the container qv_s_t
    ALLOCATE(p_diag_lnd%qv_st_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'qv_s_t',                       &
               & vname_prefix//'qv_s_t_'//ADJUSTL(TRIM(csfc)),                 &
               & p_diag_lnd%qv_st_ptr(jsfc)%p_2d,                              &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('qv_s_t_'//csfc, '', '', DATATYPE_FLT32),            &
               & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
               & ldims=shape2d )
      ENDDO

!weighted variables
    ! & p_diag_lnd%w_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('w_snow', 'm H2O', 'weighted water eqivalent of snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 60, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_snow', p_diag_lnd%w_snow,          &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,          &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars") )

    ! & p_diag_lnd%t_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('t_snow', 'K', 'weighted temperature of the snow-surface', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_snow', p_diag_lnd%t_snow,          &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,          &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars", "snow_vars") )

    ! & p_diag_lnd%t_snow_mult(nproma,nlev_snow+1,nblks_c)
    cf_desc    = t_cf_var('t_snow_mult', 'K', 'weighted temperature of the snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_snow_mult', p_diag_lnd%t_snow_mult, &
     & GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC, cf_desc, grib2_desc,                &
     & ldims=(/nproma,nlev_snow+1,kblks/),                                        &
     & lrestart=.FALSE., loutput=.TRUE.,                                         &
     & in_group=groups("multisnow_vars") ) 

    ! & p_diag_lnd%t_s(nproma,nblks_c)
    cf_desc    = t_cf_var('t_s', 'K', 'weighted temperature of ground surface', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_s', p_diag_lnd%t_s,                &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

    ! & p_diag_lnd%t_seasfc(nproma,nblks_c)
    cf_desc    = t_cf_var('t_seasfc', 'K', 'sea surface temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_seasfc', p_diag_lnd%t_seasfc,     &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE. )

    ! & p_diag_lnd%rho_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('rho_snow', 'kg/m**3', 'weighted snow density', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rho_snow', p_diag_lnd%rho_snow,      &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,  cf_desc, grib2_desc,          &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars", "snow_vars") )

    ! & p_diag_lnd%rho_snow_mult(nproma,nlev_snow,nblks_c)
    cf_desc    = t_cf_var('rho_snow_mult', 'kg/m**3', 'weighted snow density', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rho_snow_mult',                      &
         & p_diag_lnd%rho_snow_mult, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,      &
         & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
         & lrestart=.FALSE., loutput=.TRUE.,                                    &
         & in_group=groups("multisnow_vars"))

    ! & p_diag_lnd%w_i(nproma,nblks_c)
    cf_desc    = t_cf_var('w_i', 'm H2O', 'weighted water content of interception water', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_i', p_diag_lnd%w_i,                &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars") )  

    ! & p_diag_lnd%t_so(nproma,nlev_soil+2,nblks_c)
    cf_desc    = t_cf_var('t_so', 'K', 'weighted soil temperature (main level)', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_so', p_diag_lnd%t_so,              &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+2,kblks/),                                   &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars") )

   ! & p_diag_lnd%w_so(nproma,nlev_soil+1,nblks_c)
    cf_desc    = t_cf_var('w_so', 'm H20', 'weighted total water content (ice + liquid water)', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 3, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_so', p_diag_lnd%w_so,              &
         & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks/),                                   &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars") )

    ! & p_diag_lnd%w_so_ice(nproma,nlev_soil+1,nblks_c)
    cf_desc    = t_cf_var('w_so_ice', 'm H20', 'weighted ice content', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)      
    CALL add_var( diag_list, vname_prefix//'w_so_ice',                           &
         & p_diag_lnd%w_so_ice, GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_LAND,  &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil+1,kblks/),              &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars") )

    ! & p_diag_lnd%wliq_snow(nproma,nlev_snow,nblks_c)
    cf_desc    = t_cf_var('wliq_snow', 'm H2O', 'weighted liquid water content in snow', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'wliq_snow',                          &
         & p_diag_lnd%wliq_snow, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,          & 
         & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
         & lrestart=.FALSE., loutput=.TRUE.,                                    &
         & in_group=groups("multisnow_vars", "snow_vars"))

    ! & p_diag_lnd%wtot_snow(nproma,nlev_snow,nblks_c)
    cf_desc    = t_cf_var('wtot_snow', 'm H2O', 'weighted total water content in snow', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'wtot_snow',                          &
         & p_diag_lnd%wtot_snow, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,          &
         & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
         & lrestart=.FALSE., loutput=.TRUE.,                                    &
         & in_group=groups("multisnow_vars", "snow_vars"))

    ! & p_diag_lnd%dzh_snow(nproma,nlev_snow,nblks_c)
    cf_desc    = t_cf_var('dzh_snow', 'm', &
         &                'weighted layer thickness between half levels in snow', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 14, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'dzh_snow',                           &
         & p_diag_lnd%dzh_snow, GRID_UNSTRUCTURED_CELL, ZAXIS_GENERIC,           &
         & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
         & lrestart=.FALSE., loutput=.TRUE.,                                    &
         & in_group=groups("multisnow_vars", "snow_vars"))

    ! & p_diag_lnd%h_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('h_snow', 'm', 'weighted snow height', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 14, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow', p_diag_lnd%h_snow,        &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )    

    ! & p_diag_lnd%freshsnow(nproma,nblks_c)
    cf_desc    = t_cf_var('freshsnow', '- ', &
           & 'weighted indicator for age of snow in top of snow layer', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow', p_diag_lnd%freshsnow,     &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

    ! & p_diag_lnd%snowfrac(nproma,nblks_c)
    cf_desc    = t_cf_var('snowfrac', '- ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac', p_diag_lnd%snowfrac,       &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("land_vars") )

    ! & p_diag_lnd%runoff_s(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_s', 'kg/m2', &
         &                'weighted surface water runoff; sum over forecast', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s', p_diag_lnd%runoff_s,         &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

    ! & p_diag_lnd%runoff_g(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_g', 'kg/m2', &
         &                'weighted soil water runoff; sum over forecast', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g', p_diag_lnd%runoff_g,         &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

!tiled variables
    ! & p_diag_lnd%h_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('h_snow_t', 'm', 'snow height', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 14, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow_t', p_diag_lnd%h_snow_t,    &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,       &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container h_snow
    ALLOCATE(p_diag_lnd%h_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'h_snow_t',                     &
               & vname_prefix//'h_snow_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%h_snow_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
               & t_cf_var('h_snow_t_'//csfc, '', '', DATATYPE_FLT32),          &
               & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
               & ldims=shape2d )
      ENDDO



    ! & p_diag_lnd%freshsnow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('freshsnow_t', '- ', &
         &                'indicator for age of snow in top of snow layer', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow_t', p_diag_lnd%freshsnow_t, &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container freshsnow
    ALLOCATE(p_diag_lnd%freshsnow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'freshsnow_t',                     &
                 & vname_prefix//'freshsnow_t_'//ADJUSTL(TRIM(csfc)),             &
                 & p_diag_lnd%freshsnow_ptr(jsfc)%p_2d,                           &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                         &
                 & t_cf_var('freshsnow_t_'//csfc, '', '', DATATYPE_FLT32),        &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &
                 & ldims=shape2d )
      END DO

    ! & p_diag_lnd%snowfrac_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_t', '- ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_t', p_diag_lnd%snowfrac_t, &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'snowfrac_t',                     &
                 & vname_prefix//'snowfrac_t_'//ADJUSTL(TRIM(csfc)),             &
                 & p_diag_lnd%snowfrac_ptr(jsfc)%p_2d,                           &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                 & t_cf_var('snowfrac_t_'//csfc, '', '', DATATYPE_FLT32),        &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                 & ldims=shape2d )
      END DO

    ! & p_diag_lnd%snowfrac_lc_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_lc_t', '- ', 'snow-cover fraction per land-cover class', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lc_t', p_diag_lnd%snowfrac_lc_t, &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_lc_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'snowfrac_lc_t',                  &
                 & vname_prefix//'snowfrac_lc_t_'//ADJUSTL(TRIM(csfc)),          &
                 & p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d,                        &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                 & t_cf_var('snowfrac_lc_t_'//csfc, '', '', DATATYPE_FLT32),     &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                 & ldims=shape2d )
      END DO

    ! & p_diag_lnd%w_snow_eff_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('w_snow_eff_t', '- ', 'effective snow-water equivalent', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_snow_eff_t', p_diag_lnd%w_snow_eff_t, &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%w_snow_eff_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'w_snow_eff_t',                  &
                 & vname_prefix//'w_snow_eff_t_'//ADJUSTL(TRIM(csfc)),          &
                 & p_diag_lnd%w_snow_eff_ptr(jsfc)%p_2d,                        &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                 & t_cf_var('w_snow_eff_t_'//csfc, '', '', DATATYPE_FLT32),     &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                 & ldims=shape2d, lrestart=.FALSE. )
      END DO

    ! & p_diag_lnd%runoff_s_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_s_t', 'kg/m2', &
         &                'surface water runoff; sum over forecast', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s_t', p_diag_lnd%runoff_s_t,     &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container runoff_s
    ALLOCATE(p_diag_lnd%runoff_s_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'runoff_s_t',                        &
                 & vname_prefix//'runoff_s_t_'//ADJUSTL(TRIM(csfc)),                &
                 & p_diag_lnd%runoff_s_ptr(jsfc)%p_2d,                              &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                           &
                 & t_cf_var('runoff_s_t_'//csfc, '', '', DATATYPE_FLT32),           &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),    &
                 & ldims=shape2d )
      END DO



    ! & p_diag_lnd%runoff_g_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_g_t', 'kg/m2', &
         &                'soil water runoff; sum over forecast', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(2, 0, 5, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g_t', p_diag_lnd%runoff_g_t,     &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container runoff_g
    ALLOCATE(p_diag_lnd%runoff_g_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'runoff_g_t',                        &
                 & vname_prefix//'runoff_g_t_'//ADJUSTL(TRIM(csfc)),                &
                 & p_diag_lnd%runoff_g_ptr(jsfc)%p_2d,                              &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                           &
                 & t_cf_var('runoff_g_t_'//csfc, '', '', DATATYPE_FLT32),           &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),    &
                 & ldims=shape2d )
      END DO


    ! & p_prog_lnd%subsfrac_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('subsfrac', '-', 'subscale fraction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'subsfrac_t', p_diag_lnd%subsfrac_t,    &
           & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,           &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,              &
           & loutput=.FALSE.,                                                      &
           & initval_r=1._wp )


    ! fill the separate variables belonging to the container subsfrac
    ALLOCATE(p_diag_lnd%subsfrac_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'subsfrac_t',                       &
                 & vname_prefix//'subsfrac_t_'//ADJUSTL(TRIM(csfc)),               &
                 & p_diag_lnd%subsfrac_ptr(jsfc)%p_2d,                             &
                 & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                          &
                 & t_cf_var('subsfrac_t_'//csfc, '', '', DATATYPE_FLT32),          &
                 & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),   &
                 & ldims=shape2d )
      END DO



    p_diag_lnd%qv_s_t(:,:,:)     = 0.001_wp
    p_diag_lnd%subsfrac_t(:,:,:)  = 1._wp

    ENDIF

    p_diag_lnd%qv_s(:,:)        = 0.001_wp
    p_diag_lnd%fr_seaice(:,:)   = 0._wp


  END SUBROUTINE  new_nwp_lnd_diag_list


END MODULE mo_nwp_lnd_state

