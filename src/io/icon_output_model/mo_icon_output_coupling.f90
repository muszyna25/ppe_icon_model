!>
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

#ifdef YAC_coupling

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_icon_output_coupling

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: warning, message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_physical_constants,  ONLY: tmelt, rhoh2o
  USE mo_mpi,                 ONLY: p_pe_work
  USE mo_run_config,          ONLY: ltimer
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling, &
       &                            timer_coupling_put, timer_coupling_get,  &
       &                            timer_coupling_1stget, timer_coupling_init
  USE mo_sync,                ONLY: sync_c, sync_patch_array
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d

  USE mo_time_config,         ONLY: set_tc_current_date
  USE mtime,                  ONLY: datetime, datetimeToString, &
    &                               MAX_DATETIME_STR_LEN

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
  USE mo_yac_finterface,      ONLY: yac_finit, yac_fdef_comp, yac_fget_version,  &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize, yac_fput, yac_fget,           &
    &                               YAC_LOCATION_CELL, COUPLING, OUT_OF_BOUND
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_time_config,         ONLY: time_config 
  USE mo_hamocc_nml,          ONLY: l_cpl_co2

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_icon_output_coupling, destruct_icon_output_coupling
  PUBLIC :: couple_icon_output_tomodel

CONTAINS

  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and icon_output. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_icon_output_coupling()


  END SUBROUTINE construct_icon_output_coupling

  !--------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE destruct_icon_output_coupling()


  END SUBROUTINE destruct_icon_output_coupling

  !--------------------------------------------------------------------------

  SUBROUTINE couple_icon_output_tomodel()

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: method_name = 'couple_icon_output_tomodel'
    

  END SUBROUTINE couple_icon_output_tomodel
  !--------------------------------------------------------------------------

END MODULE mo_icon_output_coupling

#else

MODULE mo_icon_output_coupling

  USE mo_model_domain,        ONLY: t_patch_3d
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_exception,           ONLY: finish
  USE mtime,                  ONLY: datetime

  PUBLIC :: construct_icon_output_coupling, destruct_icon_output_coupling
  PUBLIC :: couple_icon_output_tomodel

CONTAINS

  SUBROUTINE construct_icon_output_coupling()


    IF ( is_coupled_run() ) THEN
       CALL finish('construct_icon_output_coupling: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE construct_icon_output_coupling

  SUBROUTINE couple_icon_output_tomodel()

    IF ( is_coupled_run() ) THEN
       CALL finish('couple_icon_output_tomodel: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE couple_icon_output_tomodel

  SUBROUTINE destruct_icon_output_coupling()

    IF ( is_coupled_run() ) THEN
       CALL finish('destruct_icon_output_coupling: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE destruct_icon_output_coupling

END MODULE mo_icon_output_coupling

#endif
