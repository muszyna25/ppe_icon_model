!>
!! Provide an implementation of the ocean forcing.  
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modification by Stephan Lorenz, MPI-M (2010-06):
!!   - renaming and adjustment to ocean domain and patch_oce
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!    adapted to structures discussed in 2010-01.
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_forcing
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_ocean_nml,           ONLY: itestcase_oce, iforc_oce, analyt_forc, &
  &                               wstress_coeff, iforc_stat_oce, basin_height_deg
USE mo_model_domain,        ONLY: t_patch, t_patch_3D
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_exception,           ONLY: finish, message
USE mo_math_constants,      ONLY: pi, deg2rad
USE mo_impl_constants,      ONLY: max_char_length, sea_boundary
USE mo_math_utilities,      ONLY: gvec2cvec
USE mo_sea_ice_types,       ONLY: t_sfc_flx
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

IMPLICIT NONE
PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceForcing  '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

! Public interface

PUBLIC  :: init_sfcflx

CONTAINS


  !-------------------------------------------------------------------------
  !
  !>
  !! Initialization of stationary surface fluxes for hydrostatic ocean
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2011-09)
  !
  SUBROUTINE init_sfcflx(p_patch_3D, p_sfc_flx)
  !
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: p_patch_3D
  TYPE(t_sfc_flx)                             :: p_sfc_flx

  ! Local variables
  INTEGER :: jc, jb
  INTEGER :: i_startidx_c, i_endidx_c

  REAL(wp) :: z_lat, z_lon, y_length, y_center
  REAL(wp) :: z_forc_period = 1.0_wp !=1.0: single gyre
                                     !=2.0: double gyre
                                     !=n.0: n-gyre

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_oce_forcing:init_ho_sfcflx'

  !-------------------------------------------------------------------------
  TYPE(t_subset_range), POINTER :: all_cells
  TYPE(t_patch), POINTER        :: p_patch 
  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2D(1)
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  all_cells => p_patch%cells%all

  ! analytical forcing
  IF (iforc_oce == ANALYT_FORC) THEN

    IF (itestcase_oce == 27 .OR. itestcase_oce == 29) iforc_stat_oce = 1

    SELECT CASE (iforc_stat_oce)

    CASE (0)
      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=0: no stationary wind forcing applied' )

    CASE (1)
      CALL message(TRIM(routine), 'Testcase (27,29): Apply stationary wind forcing' )
      y_length = basin_height_deg * deg2rad
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

        DO jc = i_startidx_c, i_endidx_c

          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

            z_lat = p_patch%cells%center(jc,jb)%lat
            z_lon = p_patch%cells%center(jc,jb)%lon

            p_sfc_flx%forc_wind_u(jc,jb) = wstress_coeff * &
            & cos(z_forc_period*pi*(z_lat-y_length)/y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

          !Init cartesian wind
          CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),&
                         & p_sfc_flx%forc_wind_v(jc,jb),&
                         & p_patch%cells%center(jc,jb)%lon,&
                         & p_patch%cells%center(jc,jb)%lat,&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                         & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
  !write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)

        END DO
      END DO

    CASE (2)
      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=2: stationary wind forcing - u=cos(n*(lat-lat_0)/lat_0)')

      ! Latitudes vary from -pi/2 to pi/2
      y_length = basin_height_deg * deg2rad
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lon = p_patch%cells%center(jc,jb)%lon
          IF (p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary) THEN
            p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * cos(z_forc_period*pi*(z_lat-y_length)&
              &                                / y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF
          !write(*,*)'forc',jc,jb, z_lat, z_lat*180.0_wp/pi, p_sfc_flx%forc_wind_u(jc,jb) 
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

          !Init cartesian wind
          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                           & p_sfc_flx%forc_wind_v(jc,jb),      &
                           & p_patch%cells%center(jc,jb)%lon,   &
                           & p_patch%cells%center(jc,jb)%lat,   &
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ! write(*,*)'sfc forcing', jc,jb,p_sfc_flx%forc_wind_u(jc,jb), p_sfc_flx%forc_wind_v(jc,jb)
          ELSE
            p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
          ENDIF
        END DO
      END DO

    CASE (3)
      CALL message(TRIM(routine), &
        &  'iforc_stat_oce=3: stationary wind forcing - u=cos(n*lat/lat_0)')

      ! Use here global scale:
      !y_length = basin_height_deg * deg2rad
      !y_center = basin_center_lat * deg2rad
      y_length = 180.0_wp * deg2rad
      y_center = -60.0_wp * deg2rad
      z_forc_period = 3.0_wp
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lon = p_patch%cells%center(jc,jb)%lon
          IF (p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary) THEN
            p_sfc_flx%forc_wind_u(jc,jb) =  wstress_coeff * cos(z_forc_period*pi*(z_lat-y_center)&
              &                                / y_length) 
          ELSE
            p_sfc_flx%forc_wind_u(jc,jb) = 0.0_wp
          ENDIF
          p_sfc_flx%forc_wind_v(jc,jb) = 0.0_wp

          !Init cartesian wind
          IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN
            CALL gvec2cvec(  p_sfc_flx%forc_wind_u(jc,jb),      &
                           & p_sfc_flx%forc_wind_v(jc,jb),      &
                           & p_patch%cells%center(jc,jb)%lon,   &
                           & p_patch%cells%center(jc,jb)%lat,   &
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(1),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(2),&
                           & p_sfc_flx%forc_wind_cc(jc,jb)%x(3))
          ELSE
            p_sfc_flx%forc_wind_cc(jc,jb)%x(:) = 0.0_wp
          ENDIF
        END DO
      END DO

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Stationary Analytical Forcing not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN STATIONARY FORCING OPTION NOT SUPPORTED - TERMINATE')

    END SELECT

    !---------Debug Diagnostics-------------------------------------------
    idt_src=0  ! output print level - 0: print in any case
    CALL dbg_print('analytical forcing u'      ,p_sfc_flx%forc_wind_u   ,str_module,idt_src)
    CALL dbg_print('analytical forcing v'      ,p_sfc_flx%forc_wind_u   ,str_module,idt_src)
    !---------------------------------------------------------------------

  END IF

  END SUBROUTINE init_sfcflx


END MODULE mo_oce_forcing
