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
#ifndef __NO_ICON_ATMO__
  USE mo_meteogram_output,    ONLY: meteogram_flush_file
  USE mo_meteogram_config,    ONLY: meteogram_output_config
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output'

  LOGICAL :: l_omit_dom   ! flag if "'_DOM',jg" should be omitted in the filenames

  !------------------------------------------------------------------------------------------------
  !


END MODULE mo_output
