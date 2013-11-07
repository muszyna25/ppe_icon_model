!>
!! Contains output routines for CDI output
!! This module acts merely as a wrapper for either calling the direct
!! or the asynchronous output routines.
!!
!!
!! @par Revision History
!! Initial implementation by Rainer Johanni (2010-12-02)
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
!!
MODULE mo_output

  USE mo_exception,           ONLY: message_text, get_filename_noext !, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_units,            ONLY: filename_max
  USE mo_grid_config,         ONLY: n_dom, &
    &                               n_dom_start !, nroot, lplane
  USE mo_io_config,           ONLY: out_expname
!  USE mo_impl_constants,      ONLY: ihs_ocean !,            &
!     &                              ihs_atm_temp,         &
!     &                              ihs_atm_theta,        &
!     &                              inh_atmosphere,       &
!     &                              ishallow_water
  USE mo_dynamics_config,     ONLY: iequations, nold, nnow, nnew, nnew_rcf, nnow_rcf 
  USE mo_datetime,            ONLY: t_datetime,iso8601
  USE mo_io_restart,          ONLY: set_restart_time, set_restart_vct,         &
    &                               init_restart, open_writing_restart_files,  &
    &                               write_restart, close_writing_restart_files,&
    &                               finish_restart, set_restart_depth,         &
    &                               set_restart_depth_lnd, &  !DRset_restart_height, &
    &                               set_restart_height_snow
  USE mo_io_restart_attributes,ONLY: set_restart_attribute
  USE mo_name_list_output_init, ONLY: output_file
  USE mo_model_domain,        ONLY: t_patch,t_patch_3D, p_patch
  USE mo_intp_data_strc,      ONLY: t_lon_lat_intp
  USE mo_run_config,          ONLY: ltimer, output_mode
  USE mo_timer,               ONLY: timer_start, timer_stop,&
    &                     timer_write_restart_file, timer_write_output
#ifdef __ICON_ATMO__
  USE mo_meteogram_output,    ONLY: meteogram_flush_file
  USE mo_meteogram_config,    ONLY: meteogram_output_config
#endif

#ifdef __ICON_OCEAN__
  USE mo_oce_state,           ONLY: set_zlev
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output'

  LOGICAL :: l_omit_dom   ! flag if "'_DOM',jg" should be omitted in the filenames

  !------------------------------------------------------------------------------------------------
  !
  ! Public routines:
  PUBLIC :: create_restart_file

CONTAINS

  !-------------
  !>
  !! 
  !! Hui Wan (MPI-M, 2011-05)
  !!
  SUBROUTINE create_restart_file( patch, datetime,             &
                                & jstep,                       &
                                & opt_pvct,                    &
                                & opt_t_elapsed_phy,           &
                                & opt_lcall_phy, opt_sim_time, &
                                & opt_jstep_adv_ntsteps,       &
                                & opt_jstep_adv_marchuk_order, &
                                & opt_depth, opt_depth_lnd,    &
                                & opt_nlev_snow,               &
                                & opt_nice_class)

    TYPE(t_patch),   INTENT(IN) :: patch
    TYPE(t_datetime),INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jstep                ! simulation step

    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth_lnd   ! vertical levels soil model
    REAL(wp), INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_ntsteps
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nlev_snow
    INTEGER,  INTENT(IN), OPTIONAL :: opt_nice_class

    INTEGER :: klev, jg, kcell, kvert, kedge, icelltype 
    INTEGER :: izlev, inlev_soil, inlev_snow, i, nice_class
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)


    CHARACTER(LEN=132) :: string
    CHARACTER(len=MAX_CHAR_LENGTH) :: attname   ! attribute name
    INTEGER :: jp, jp_end   ! loop index and array size

 
    IF (ltimer) CALL timer_start(timer_write_restart_file)
    !----------------
    ! Initialization
    klev      = patch%nlev
    jg        = patch%id
    kcell     = patch%n_patch_cells_g
    kvert     = patch%n_patch_verts_g
    kedge     = patch%n_patch_edges_g
    icelltype = patch%cell_type

    CALL set_restart_attribute( 'current_caltime', datetime%caltime )
    CALL set_restart_attribute( 'current_calday' , datetime%calday )

    CALL set_restart_attribute( 'current_daysec' , datetime%daysec )

    CALL set_restart_attribute( 'nold'    , nold    (jg))
    CALL set_restart_attribute( 'nnow'    , nnow    (jg))
    CALL set_restart_attribute( 'nnew'    , nnew    (jg))
    CALL set_restart_attribute( 'nnow_rcf', nnow_rcf(jg))
    CALL set_restart_attribute( 'nnew_rcf', nnew_rcf(jg))

    ! set simulation step
    CALL set_restart_attribute( 'jstep', jstep )

    !----------------
    ! additional restart-output for nonhydrostatic model
    IF (PRESENT(opt_sim_time)) THEN
      WRITE(attname,'(a,i2.2)') 'sim_time_DOM',jg
      CALL set_restart_attribute( TRIM(attname), opt_sim_time )
    ENDIF

    !-------------------------------------------------------------
    ! DR
    ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
    ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
    ! SET_RESTART_ATTRIBUTE
    !-------------------------------------------------------------

    IF (PRESENT(opt_jstep_adv_ntsteps)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_ntsteps_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_ntsteps )
    ENDIF

    IF (PRESENT(opt_jstep_adv_marchuk_order)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_marchuk_order )
    ENDIF

    IF (PRESENT(opt_t_elapsed_phy) .AND. PRESENT(opt_lcall_phy)) THEN
      ! Inquire array size
      jp_end = SIZE(opt_t_elapsed_phy,2)
      DO jp = 1, jp_end
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_t_elapsed_phy(jg,jp) )
      ENDDO
      ! Inquire array size
      jp_end = SIZE(opt_lcall_phy,2)
      DO jp = 1, jp_end
        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_lcall_phy(jg,jp) )
      ENDDO
    ENDIF

    IF (PRESENT(opt_pvct)) CALL set_restart_vct( opt_pvct )  ! Vertical coordinate (A's and B's)
    IF (PRESENT(opt_depth_lnd)) THEN            ! geometrical depth for land module
      !This part is only called if opt_depth_lnd > 0
      IF (opt_depth_lnd > 0) THEN  
        inlev_soil = opt_depth_lnd
        ALLOCATE(zlevels_full(inlev_soil))
        ALLOCATE(zlevels_half(inlev_soil+1))
        DO i = 1, inlev_soil
          zlevels_full(i) = REAL(i,wp)
        END DO
        DO i = 1, inlev_soil+1
          zlevels_half(i) = REAL(i,wp)
        END DO
        CALL set_restart_depth_lnd(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
      ELSE
       inlev_soil = 0
      END IF
    ELSE
      inlev_soil = 0
    ENDIF
    IF (PRESENT(opt_nlev_snow)) THEN  ! number of snow levels (multi layer snow model)
      !This part is only called if opt_nlev_snow > 0
      IF (opt_nlev_snow > 0) THEN   
        inlev_snow = opt_nlev_snow
        ALLOCATE(zlevels_full(inlev_snow))
        ALLOCATE(zlevels_half(inlev_snow+1))
        DO i = 1, inlev_snow
          zlevels_full(i) = REAL(i,wp)
        END DO
        DO i = 1, inlev_snow+1
          zlevels_half(i) = REAL(i,wp)
        END DO
        CALL set_restart_height_snow(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
      ELSE
        inlev_snow = 0
      ENDIF
    ELSE
      inlev_snow = 0
    ENDIF
!DR end preliminary fix
    izlev = 0
#ifdef __ICON_OCEAN__
    IF (PRESENT(opt_depth)) THEN                              ! Ocean depth
      !This part is only called if opt_depth > 0
      IF(opt_depth>0)THEN
        izlev = opt_depth
        ALLOCATE(zlevels_full(izlev))
        ALLOCATE(zlevels_half(izlev+1))
        CALL set_zlev(zlevels_half, zlevels_full)
        CALL set_restart_depth(zlevels_half, zlevels_full)
        DEALLOCATE(zlevels_full)
        DEALLOCATE(zlevels_half)
!      ELSE
!        izlev = 0
!      END IF
!    ELSE
!      izlev = 0
    END IF
#endif

    IF (.NOT.PRESENT(opt_nice_class)) THEN
      nice_class = 1
    ELSE
      nice_class = opt_nice_class
    END IF

    CALL init_restart( TRIM(out_expname), &! exp name
                     & '1.2.2',           &! model version
                     & kcell, icelltype,  &! total # of cells, # of vertices per cell
                     & kvert, 9-icelltype,&! total # of vertices, # of vertices per dual cell
                     & kedge, 4,          &! total # of cells, shape of control volume for edge 
                     & klev,              &! total # of vertical layers
                     & izlev,             &! total # of depths below sea
                     & inlev_soil,        &! total # of depths below land (TERRA or JSBACH)
                     & inlev_snow,        &! total # of vertical snow layers (TERRA)
                     & nice_class         )! total # of ice classes (sea ice)

    CALL set_restart_time( iso8601(datetime) )  ! Time tag

    ! Open new file, write data, close and then clean-up.
    message_text = get_filename_noext(patch%grid_filename)
    WRITE(string,'(a,a)') 'restart.',TRIM(message_text)

    CALL open_writing_restart_files( TRIM(string) )

#ifdef NOMPI
    CALL write_restart
#else
    CALL write_restart( patch )
#endif

    CALL close_writing_restart_files
    CALL finish_restart

    IF (ltimer) CALL timer_stop(timer_write_restart_file)

  END SUBROUTINE create_restart_file

END MODULE mo_output
