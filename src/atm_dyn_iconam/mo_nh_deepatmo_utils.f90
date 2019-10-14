!>
!! This module contains subroutines related to the upper atmosphere.
!!
!! @par Revision History
!! Initial revision by the authors of the subroutines of which 
!! this module contains copies.
!! Modifications for the deep atmosphere by Sebastian Borchert, DWD (2017-06-30)
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_deepatmo_utils

  USE mo_kind,                  ONLY: wp, vp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, SUCCESS, min_rlcell, min_rlvert, &
    &                                 min_rlcell_int, min_rledge_int, min_rlvert_int
  USE mo_physical_constants,    ONLY: grav, rd, p0ref, rd_o_cpd, cpd, p0sl_bg
  USE mo_upatmo_config,         ONLY: upatmo_dyn_config, &
    &                                 imsg_thr, itmr_thr, idamtr
  USE mo_grid_config,           ONLY: grid_sphere_radius
  USE mo_run_config,            ONLY: lvert_nest, msg_level
  USE mo_nonhydrostatic_config, ONLY: lextra_diffu
  USE mo_parallel_config,       ONLY: nproma 
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_init_vgrid,            ONLY: nflatlev
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_math_divrot,           ONLY: rot_vertex_ri
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, &
    &                                 SYNC_C, SYNC_E
  USE mo_timer,                 ONLY: timer_start, timer_stop, timers_level, &
    &                                 timer_deepatmo_ztrafo, timer_solve_nh_veltend
  USE mo_icon_interpolation_scalar, ONLY: cells2verts_scalar_ri
  USE mo_util_string,           ONLY: int2string, real2string
  USE mo_util_table,            ONLY: t_table, initialize_table, add_table_column, &
    &                                 set_table_entry, print_table, finalize_table
  USE mo_mpi,                   ONLY: my_process_is_stdio

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: height_transform
  PUBLIC :: set_deepatmo_metrics
  PUBLIC :: velocity_tendencies_deepatmo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_deepatmo_utils'

  ! Unfortunately, more than one 'height_transform' is required, 
  ! because different 'intent'-attributes are required to cover 
  ! the possible input fields
  INTERFACE height_transform
    MODULE PROCEDURE height_transform_a
    MODULE PROCEDURE height_transform_b
  END INTERFACE height_transform

CONTAINS !..................................................................................

  !>
  !! Transformation of the height coordinate
  !!
  !! Variant a: in- is also out-field 
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2016-11-04)
  !!
  SUBROUTINE height_transform_a( z_inout,    &  !in/out
    &                            nblks,      &  !in
    &                            npromz,     &  !in
    &                            nlevs,      &  !in
    &                            lconstgrav, &  !in
    &                            trafo_type  )  !in

    ! In/out variables
    REAL(wp),         INTENT(INOUT) :: z_inout(:,:,:) ! Height coordinate field
    INTEGER ,         INTENT(IN)    :: nblks          ! Number of blocks
    INTEGER ,         INTENT(IN)    :: npromz         ! Length of last block
    INTEGER ,         INTENT(IN)    :: nlevs          ! Number of input levels
    LOGICAL,          INTENT(IN)    :: lconstgrav     ! .TRUE. -> const. gravitational acceleration
    CHARACTER(LEN=*), INTENT(IN)    :: trafo_type     ! Type of transformation

    ! Local variables
    REAL(wp) :: trafo_fac
    INTEGER  :: jb, jk, jc  ! (jc is habitual placeholder for jc, je, jv)
    INTEGER  :: nlen

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &    
      & routine = modname//':height_transform_a'

    !-----------------------------------------------------------------------

    ! Return, if gravitational acceleration is constant
    IF (lconstgrav) RETURN

    IF (timers_level > itmr_thr%med) CALL timer_start(timer_deepatmo_ztrafo)

    SELECT CASE(TRIM(trafo_type)) 
    CASE('zgpot2z')  
      ! Transform geopotential height z_gpot into 
      ! geometric height z by means of 
      ! z = z_gpot / ( 1 - z_gpot / a ), 
      ! where a is radius of Earth
      trafo_fac = -1._wp/grid_sphere_radius
    CASE('z2zgpot')
      ! Transform geometric height z into 
      ! geopotential height z_gpot by means of 
      ! z_gpot = z / ( 1 + z / a ), 
      ! where a is radius of Earth
      trafo_fac = 1._wp/grid_sphere_radius
    CASE DEFAULT
      CALL finish( TRIM(routine),'invalid trafo_type')
    END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z_inout(nlen+1:nproma,:,jb) = 0.0_wp
      ENDIF
      
      DO jk = 1, nlevs
        DO jc = 1, nlen
          z_inout(jc,jk,jb) = z_inout(jc,jk,jb) &
            &               / ( 1._wp + trafo_fac * z_inout(jc,jk,jb) )
        ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > itmr_thr%med) CALL timer_stop(timer_deepatmo_ztrafo)

  END SUBROUTINE height_transform_a

  !--------------------------------------------------------------------------------------------------

  !>
  !! Transformation of the height coordinate
  !!
  !! Variant b: in- and out-fields differ
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2016-11-04)
  !!
  SUBROUTINE height_transform_b( z_in,       &  !in
    &                            z_out,      &  !out
    &                            nblks,      &  !in
    &                            npromz,     &  !in
    &                            nlevs,      &  !in
    &                            lconstgrav, &  !in
    &                            trafo_type  )  !in

    ! In/out variables
    REAL(wp),         INTENT(IN)  :: z_in(:,:,:)  ! Input height coordinate field
    REAL(wp),         INTENT(OUT) :: z_out(:,:,:) ! Output height coord. field
    INTEGER ,         INTENT(IN)  :: nblks        ! Number of blocks
    INTEGER ,         INTENT(IN)  :: npromz       ! Length of last block
    INTEGER ,         INTENT(IN)  :: nlevs        ! Number of input levels
    LOGICAL,          INTENT(IN)  :: lconstgrav   ! .TRUE. -> const. gravitational acceleration
    CHARACTER(LEN=*), INTENT(IN)  :: trafo_type   ! Type of transformation

    ! Local variables
    REAL(wp) :: trafo_fac
    INTEGER  :: jb, jk, jc  ! (jc is habitual placeholder for jc, je, jv)
    INTEGER  :: nlen

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &    
      & routine = modname//':height_transform_b'

    !-----------------------------------------------------------------------

    IF (timers_level > itmr_thr%med) CALL timer_start(timer_deepatmo_ztrafo)

    IF (lconstgrav)  THEN
      trafo_fac = 0._wp
    ELSE
      SELECT CASE(TRIM(trafo_type)) 
      CASE('zgpot2z')  
        ! Transform geopotential height z_gpot into 
        ! geometric height z by means of 
        ! z = z_gpot / ( 1 - z_gpot / a ), 
        ! where a is radius of Earth
        trafo_fac = -1._wp/grid_sphere_radius
      CASE('z2zgpot')
        ! Transform geometric height z into 
        ! geopotential height z_gpot by means of 
        ! z_gpot = z / ( 1 + z / a ), 
        ! where a is radius of Earth
        trafo_fac = 1._wp/grid_sphere_radius
      CASE DEFAULT
        CALL finish( TRIM(routine),'invalid trafo_type')
      END SELECT
    ENDIF
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z_out(nlen+1:nproma,:,jb) = 0.0_wp
      ENDIF
      
      DO jk = 1, nlevs
        DO jc = 1, nlen
          z_out(jc,jk,jb) = z_in(jc,jk,jb) &
            &             / ( 1._wp + trafo_fac * z_in(jc,jk,jb) )
        ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > itmr_thr%med) CALL timer_stop(timer_deepatmo_ztrafo)

  END SUBROUTINE height_transform_b

  !--------------------------------------------------------------------------------------------------

  !>
  !! Computes metrical modification factors and geopotential heights 
  !! for the deep atmosphere
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2017-06-30)
  !!
  !!
  SUBROUTINE set_deepatmo_metrics( p_patch,       &  !inout
    &                              p_nh_metrics,  &  !inout
    &                              p_int,         &  !in
    &                              igradp_method, &  !in
    &                              is_les_phy,    &  !in
    &                              h_scal_bg,     &  !in
    &                              t0sl_bg,       &  !in
    &                              del_t_bg       )  !in

    ! In/out variables
    TYPE(t_patch),      TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_nh_metrics),         INTENT(INOUT) :: p_nh_metrics
    TYPE(t_int_state),  TARGET, INTENT(IN)    :: p_int
    INTEGER,                    INTENT(IN)    :: igradp_method  ! Type of pressure gradient discretization
    LOGICAL,                    INTENT(IN)    :: is_les_phy     ! .TRUE. -> LES mode
    REAL(wp),                   INTENT(IN)    :: h_scal_bg      ! [m] (Geopotential) scale height (parameter)
    REAL(wp),                   INTENT(IN)    :: t0sl_bg        ! [K] Sea level temperature (parameter)
    REAL(wp),                   INTENT(IN)    :: del_t_bg       ! [K] Difference between sea level
                                                                ! temperature and asymptotic
                                                                ! stratospheric temperature (parameter)

    ! Local variables
    TYPE(t_table) :: table

    REAL(wp), ALLOCATABLE :: zgpot_me(:,:,:), z_me(:,:,:), &
      & z_help(:), z_temp(:), z_aux1(:), z_aux2(:)

    REAL(wp) :: z_z_l, z_z_u, z_z_c        ! Altitudes
    REAL(wp) :: z_r_l, z_r_u, z_r_c        ! Radii (measured from center of Earth)
    REAL(wp) :: z_grad_c, z_grad_ifc       ! Modif. factors for horizontal gradients
    REAL(wp) :: z_div_c, z_div_l, z_div_u  ! Modif. factors for divergence
    REAL(wp) :: z_vol                      ! Modif. factor for grid cell volume
    REAL(wp) :: z_dPhi                     ! Geopotential difference

    INTEGER :: nlev, nlevp1
    INTEGER :: nblks_c, nblks_e, nlen, npromz_c, npromz_e
    INTEGER :: jg, jb, jk, jk1, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: istat

    LOGICAL :: lconstgrav, lcentrifugal, lnontrad

    REAL(wp), PARAMETER :: z_eps = 1.0E-25_wp

    ! Table column headings
    CHARACTER(len=14), PARAMETER :: column_jk     = "layer index jk"
    CHARACTER(len=12), PARAMETER :: column_height = "height z [m]"
    CHARACTER(len=9),  PARAMETER :: column_gradh  = "gradh [1]"
    CHARACTER(len=8),  PARAMETER :: column_divh   = "divh [1]"
    CHARACTER(len=9),  PARAMETER :: column_divzU  = "divzU [1]"
    CHARACTER(len=9),  PARAMETER :: column_divzL  = "divzL [1]"
    CHARACTER(len=7),  PARAMETER :: column_vol    = "vol [1]"
    CHARACTER(len=10), PARAMETER :: column_invr   = "invr [1/m]"
    CHARACTER(len=10), PARAMETER :: column_centri = "centri [1]"

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      & routine = modname//':set_deepatmo_metrics'

    !-------------------------------------------------------------------------------

    ! Note: apart from the deep-atmosphere modification factors, 
    ! this subroutine contains recomputations of some fields, 
    ! which have been already computed in 'src/atm_dyn_iconam/mo_vertical_grid: set_nh_metrics', 
    ! in order to incorporate deep-atmosphere modifications. 
    ! Why were these recomputations not included into 'set_nh_metrics' directly?
    ! -> An important guideline for the upper-atmosphere extension 
    !    is minimal invasivity (-> reduce bugginess potential, 
    !    prevent optical overload of existing code etc.), 
    !    which would have been violated by a direct inclusion through necessary restructuring 
    !    of existing code, the interruption of existing OpenMP-threads etc.
    ! -> The computations are performed only once during model initialization, 
    !    so that we assume the computational overhead to be bearable by today's machines
    ! -> This specific deep-atmosphere subroutine allows a somewhat freer
    !    formulation of the deep-atmosphere modifications

    ! Note: a query for 'ldeepatmo' encapsulates already the call of this subroutine
    ! (So queries for the switches in 'upatmo_dyn_config(jg)' are equivalent 
    ! to queries for the corresponding switches in 'upatmo_config(jg)%dyn'.)

    ! Domain index (domain loop is outside this subroutine)
    jg = p_patch%id

    ! Blocking
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    
    ! Number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! For convenience
    lconstgrav   = upatmo_dyn_config(jg)%lconstgrav
    lcentrifugal = upatmo_dyn_config(jg)%lcentrifugal
    lnontrad     = upatmo_dyn_config(jg)%lnontrad 

    !-----------------------------------------------------
    !     Modification factors and additional fields 
    !       for the deep-atmosphere configuration
    !-----------------------------------------------------

    ! Computation of 3d-fields related to the geopotential height, 
    ! if gravitational acceleration varies with height

    IF (.NOT. lconstgrav) THEN
      ! 1st we should make sure that the fields related 
      ! to the geopotential height are allocated 
      IF (.NOT. (ASSOCIATED(p_nh_metrics%zgpot_ifc) .AND. & 
        &        ASSOCIATED(p_nh_metrics%zgpot_mc)  .AND. & 
        &        ASSOCIATED(p_nh_metrics%dzgpot_mc)      )) THEN
        CALL finish(TRIM(routine), 'zgpot_ifc, zgpot_mc and dzgpot_mc not allocated')
      ENDIF

      ! (Note: the geopotential-height-modification depends only on the radial position, 
      ! not on the meridional position as one would expect from inclusion of 
      ! the centrifugal acceleration, this is the "spherical geopotential approximation", 
      ! see e.g. White et al. 2005 QJRMS)
      CALL height_transform( p_nh_metrics%z_ifc,     &  !in 
        &                    p_nh_metrics%zgpot_ifc, &  !out    
        &                    nblks_c,                &  !in
        &                    npromz_c,               &  !in
        &                    nlevp1,                 &  !in
        &                    lconstgrav,             &  !in
        &                    'z2zgpot'               )  !in    
      
      CALL height_transform( p_nh_metrics%z_mc,     &  !in 
        &                    p_nh_metrics%zgpot_mc, &  !out  
        &                    nblks_c,               &  !in
        &                    npromz_c,              &  !in
        &                    nlev,                  &  !in
        &                    lconstgrav,            &  !in
        &                    'z2zgpot'              )  !in 

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c
        
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        
        DO jk = 1, nlev
          DO jc = 1, nlen 
            ! Geopotential layer thickness  at full levels
            p_nh_metrics%dzgpot_mc(jc,jk,jb) =        & 
              & p_nh_metrics%zgpot_ifc(jc,jk  ,jb) -  & 
              & p_nh_metrics%zgpot_ifc(jc,jk+1,jb)
          ENDDO  !jc
        ENDDO  !jk
        
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF  !IF (.NOT. lconstgrav)

    ! Computation of metrical modification factors for deep atmosphere 
    ! (Note: for the deep-atmosphere-specific metrical modification factors 
    ! and other quantities to be 1d-fields (varying only in z-direction), 
    ! the terrain-dependence of coordinate surfaces below z = 'flat_height' 
    ! is neglected in the deep-atmosphere modifications, 
    ! otherwise the following fields would become 3d-fields which would be 
    ! too costly in terms of memory. 
    ! In addition the metrical modification factors, 
    ! e.g. for flux divergences, are relatively difficult to compute 
    ! in sperical geometry, if coordinate surfaces deviate from spherical shells, 
    ! and cell edges lose the center of Earth as curvature center, so that their 
    ! shape is no longer determined by being great circle sections)

    DO jk = 1, nlev
      
      ! Correct for differences in nested domains' vertical index
      jk1 = jk + p_patch%nshift_total
      ! Height of upper half level of full level 'jk'
      z_z_u = vct_a(jk1) 
      ! Height of lower half level of full level 'jk' 
      z_z_l = vct_a(jk1+1)
      ! Height of full level                      
      z_z_c = 0.5_wp * ( z_z_l + z_z_u ) 
      ! Radii 
      z_r_l = MAX(z_eps, grid_sphere_radius + z_z_l)
      z_r_u = MAX(z_eps, grid_sphere_radius + z_z_u)
      z_r_c = MAX(z_eps, grid_sphere_radius + z_z_c)
      
      ! Compute metrical modification factors
      
      ! Metrical modification factors ... 
      ! ... for horizontal gradients at full levels
      z_grad_c = grid_sphere_radius / z_r_c     
      
      ! ... for horizontal gradient at half levels
      z_grad_ifc = grid_sphere_radius / z_r_u
      
      ! ... for divergence: 
      ! Horizontal part (= surface of side wall / cell volume * flux denisty over side wall) 
      !                    ----------------------------------
      ! (-> modification is necessary for underlined factor)  
      ! (There is almost no difference between the magnitude of 'z_div_c'  
      ! and 'z_grad_c', but nevertheless they are not identical)
      z_div_c = z_grad_c * ( 3._wp / 4._wp ) / & 
        &       ( 1._wp - z_r_l * z_r_u / ( z_r_l + z_r_u )**2 )
      ! Vertical part
      ! 1) = surface of cell bottom / cell volume * flux density over cell bottom
      !      ------------------------------------
      z_div_l = 3._wp / ( 1._wp + z_r_u / z_r_l + ( z_r_u / z_r_l )**2 )
      ! 2) = surface of cell lid / cell volume * flux density over cell lid
      !      ---------------------------------
      z_div_u = 3._wp / ( 1._wp + z_r_l / z_r_u + ( z_r_l / z_r_u )**2 )
      
      ! Metrical modification factor for the volume of a cell
      ! (This is required e.g. for volume integrals 
      ! in 'src/atm_dyn_iconam/mo_nh_supervise/subervise_total_integrals_nh')
      z_vol = ( z_r_l**2 + z_r_l * z_r_u + z_r_u**2 ) / & 
        &     ( 3._wp * grid_sphere_radius**2 )   
      
      ! Assign computed values:
      ! Full levels
      p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh) = z_grad_c
      p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%divh)  = z_div_c
      p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)   = z_vol
      p_nh_metrics%deepatmo_t2mc(idamtr%t2mc%divzU,jk) = z_div_u
      p_nh_metrics%deepatmo_t2mc(idamtr%t2mc%divzL,jk) = z_div_l

      ! Half levels
      p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%gradh) = z_grad_ifc

      IF (lnontrad) THEN
        ! Non-traditional terms in components of momentum equation have been switched on
        ! Full levels
        p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%invr)   = 1._wp / z_r_c
        ! Half levels
        p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%invr) = 1._wp / z_r_u
        ! (The else-case should be covered by the default values in 'src/atm_dyn_iconam/mo_nonhydro_state')
      ENDIF

      IF (lcentrifugal) THEN
        ! Modification factors for centrifugal acceleration
        ! Full levels
        p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%centri)   = 1._wp / z_grad_c
        ! Half levels
        p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%centri) = 1._wp / z_grad_ifc
        ! (The else-case should be covered by the default values in 'src/atm_dyn_iconam/mo_nonhydro_state')    
      ENDIF
      
    ENDDO  !jk
    
    ! For 'nlevp1'
    p_nh_metrics%deepatmo_t1ifc(nlevp1,idamtr%t1ifc%gradh) = 1._wp ! = 'grid_sphere_radius / grid_sphere_radius'
    IF (lnontrad)     p_nh_metrics%deepatmo_t1ifc(nlevp1,idamtr%t1ifc%invr)   = 1._wp / grid_sphere_radius
    IF (lcentrifugal) p_nh_metrics%deepatmo_t1ifc(nlevp1,idamtr%t1ifc%centri) = 1._wp   

    !-----------------------------------------------------
    !                   Recomputations 
    !      for inclusion of deep-atmosphere effects      
    !-----------------------------------------------------

    ! Metrical modification of the geopotential

    IF (.NOT. lconstgrav) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, nlen, z_dPhi) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c
        
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        
        DO jk = 1, nlev   
          DO jc = 1, nlen 
            ! Geopotential on full levels: 
            ! Phi = g * z / ( 1 + z / a ), provided that 
            ! Phi(r=a)=Phi(z=0)=0 m^2/s^2 (a = radius of Earth)
            p_nh_metrics%geopot(jc,jk,jb) = grav * & 
              & p_nh_metrics%zgpot_mc(jc,jk,jb)
          ENDDO  !jc
        ENDDO  !jk
        
        ! Geopot above ground

        DO jc = 1, nlen 
          p_nh_metrics%geopot_agl_ifc(jc,nlevp1,jb) = 0._wp
        ENDDO

          DO jk = nlev, 1, -1 
            DO jc = 1, nlen              
              ! Geopotential difference between layer interfaces
              z_dPhi = grav * p_nh_metrics%dzgpot_mc(jc,jk,jb)
              
              p_nh_metrics%dgeopot_mc(jc,jk,jb) = z_dPhi 
              
              ! Geopotential (interfaces)
              p_nh_metrics%geopot_agl_ifc(jc,jk,jb) = & 
                & p_nh_metrics%geopot_agl_ifc(jc,jk+1,jb) + z_dPhi
              
              ! Unfortunately, gpot[(z1+z2)/2] /= [gpot(z1)+gpot(z2)]/2, 
              ! where gpot(z)=zgpot, so we have to compute it separately for full levels
              z_dPhi = grav * ( p_nh_metrics%zgpot_mc(jc,jk,jb) -  & 
                &               p_nh_metrics%zgpot_ifc(jc,jk+1,jb) )

              p_nh_metrics%geopot_agl(jc,jk,jb) = & 
                & p_nh_metrics%geopot_agl_ifc(jc,jk+1,jb) + z_dPhi               
            ENDDO  !jc            
          ENDDO  !jk 
        
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF  !IF (.NOT. lconstgrav)

    ! Metrical modification of interface slope

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e

      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF

      IF (is_les_phy) THEN
        DO jk = 1, nlevp1
          DO je = 1, nlen
            p_nh_metrics%ddxt_z_half_e(je,jk,jb) = p_nh_metrics%ddxt_z_half_e(je,jk,jb) * &
              & p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%gradh)
            p_nh_metrics%ddxn_z_half_e(je,jk,jb) = p_nh_metrics%ddxn_z_half_e(je,jk,jb) * &
              & p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%gradh)
          ENDDO  !je
        ENDDO  !jk
      ENDIF  !IF (is_les_phy)

      DO jk = 1, nlev
        DO je = 1, nlen
          p_nh_metrics%ddxn_z_full(je,jk,jb) = p_nh_metrics%ddxn_z_full(je,jk,jb) * &
            & p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh)          
          p_nh_metrics%ddxt_z_full(je,jk,jb) = p_nh_metrics%ddxt_z_full(je,jk,jb) * &
            & p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh)
        ENDDO  !je
      ENDDO  !jk

      ! Note: numerous quantities are determined from the slope at (or close to) the ground. 
      ! A modification of them is not necessary, since the respective modification factors 
      ! are terrain-independent and would be equal to 1 everywhere at the ground 
      ! (or at least negligibly small)
      
    ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Metrical modification of reference atmosphere

    ! The deep-atmosphere case is more or less a copy of the shallow-atmosphere case 
    ! (for computational efficiency reasons), with the geometric height z 
    ! replaced by the geopotential height zgpot (we say that the reference temperature 
    ! has the same functional dependency on zgpot as it has on z in case of the 
    ! shallow atmosphere, this way the integral of -dp/dzgpot-rho*grav=0 in the 
    ! deep-atmosphere case is formally identical to the integral of -dp/dz-rho*grav=0 
    ! in case of the shallow atmospehre). 
    ! Note: the vertical derivative in 'd2dexdz2_fac1_mc', 'd2dexdz2_fac2_mc', 
    ! and 'd_exner_dz_ref_ic' is still with respect to z, not zgpot.
    
    IF (.NOT. lconstgrav) THEN
      
      ! Auxiliary fields
      ALLOCATE( z_me(nproma,nlev,nblks_e),     &
        &       zgpot_me(nproma,nlev,nblks_e), &
        &       z_help(nproma),                &
        &       z_temp(nproma),                &
        &       z_aux1(nproma),                &
        &       z_aux2(nproma),                &
        &       STAT=istat                     )
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation of auxiliary fields failed!') 

      ! Compute geometric height at edge points
      CALL cells2edges_scalar(p_nh_metrics%z_mc, p_patch, p_int%c_lin_e, z_me)

      CALL sync_patch_array(SYNC_E, p_patch, z_me)

      ! Transform z -> zgpot 
      ! (Given an interplated height z=alpha*z1 + (1-alpha)*z2, we find unfortunately that 
      ! gpot(z) /= alpha*gpot(z1) + (1-alpha)*gpot(z2), where gpot(z) denotes the geopotential 
      ! height corresponding to z. So it seems more secure to compute gpot(z_me), instead of 
      ! interpolating it from zgpot_mc.)
      CALL height_transform( z_me,       &  !in
        &                    zgpot_me,   &  !out
        &                    nblks_e,    &  !in
        &                    npromz_e,   &  !in
        &                    nlev,       &  !in
        &                    lconstgrav, &  !in
        &                    'z2zgpot'   )  !in

!$OMP PARALLEL PRIVATE(i_startblk, i_endblk)
!$OMP DO PRIVATE(jb, jk, jc, nlen, z_help, z_temp, z_aux1, z_aux2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c
        
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        
        ! Reference surface temperature
        DO jc = 1, nlen
          p_nh_metrics%tsfc_ref(jc,jb) = ( t0sl_bg - del_t_bg ) + del_t_bg * &
            & EXP( -p_nh_metrics%zgpot_ifc(jc,nlevp1,jb) / h_scal_bg )
        ENDDO  !jc
        
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! Reference pressure, full level mass points
            z_aux1(jc) = p0sl_bg  * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
              & LOG( ( EXP( p_nh_metrics%zgpot_mc(jc,jk,jb) / h_scal_bg ) * &
              & ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )
            
            ! Reference Exner pressure, full level mass points
            p_nh_metrics%exner_ref_mc(jc,jk,jb) = ( z_aux1(jc) / p0ref )**rd_o_cpd
            
            ! Reference temperature, full level mass points
            z_temp(jc) = ( t0sl_bg - del_t_bg ) + del_t_bg * &
              & EXP( -p_nh_metrics%zgpot_mc(jc,jk,jb) / h_scal_bg )
            
            ! Reference density, full level mass points
            p_nh_metrics%rho_ref_mc(jc,jk,jb) = z_aux1(jc) / ( rd * z_temp(jc) )
            
            ! Reference potential temperature, full level mass points
            p_nh_metrics%theta_ref_mc(jc,jk,jb) = z_temp(jc) / &
              & p_nh_metrics%exner_ref_mc(jc,jk,jb)
          ENDDO  !jc
        ENDDO  !jk
        
        IF (igradp_method <= 3) THEN
          DO jk = 1, nlev
            DO jc = 1, nlen
              ! First vertical derivative of reference Exner pressure, full level mass points,
              ! divided by theta_ref
              ! Note: for computational efficiency, this field is in addition divided by
              ! the vertical layer thickness
              ! (Deep-atmosphere modification: the vertical derivative is with respect to z, 
              ! so we use dexner/dz=dexner/dzgpot*dzgpot/dz, with dzgpot/dz=(a/r)^2, 
              ! where a is the radius of Earth, and r=a+z. For the actual modification, 
              ! we can use that (a/r)^2='deepatmo_gradh_c'^2)
              p_nh_metrics%d2dexdz2_fac1_mc(jc,jk,jb) = &
                & -grav / ( cpd * p_nh_metrics%theta_ref_mc(jc,jk,jb)**2 ) * &
                & p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * & 
                & p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh)**2
              
              ! Vertical derivative of d_exner_dz/theta_ref, full level mass points
              ! (Deep-atmosphere modification: here we use that for a quantity X, 
              ! d^2X/dz^2=d^2X/dzgpot^2*dzgpot/dz*(a/r)^2 + dX/dzgpot*d(a/r)^2/dr 
              ! =d^2X/dzgpot^2*(a/r)^4 - 2*dX/dzgpot*(a/r)^2/r)
              p_nh_metrics%d2dexdz2_fac2_mc(jc,jk,jb) = &
                &  2._wp * grav / ( cpd * p_nh_metrics%theta_ref_mc(jc,jk,jb)**3 ) * ( grav / cpd &
                & - del_t_bg / h_scal_bg * EXP( -p_nh_metrics%zgpot_mc(jc,jk,jb) / h_scal_bg ) ) / &
                & p_nh_metrics%exner_ref_mc(jc,jk,jb) * &
                & p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh)**4 & 
                & - 2._wp * p_nh_metrics%d2dexdz2_fac1_mc(jc,jk,jb) * &
                & p_nh_metrics%ddqz_z_full(jc,jk,jb) * &
                & p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh) / grid_sphere_radius
            ENDDO  !jc
          ENDDO  !jk
        ENDIF
        
        DO jk = 1, nlevp1
          DO jc = 1, nlen
            ! Reference pressure, half level mass points
            z_aux1(jc) = p0sl_bg * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
              & LOG( ( EXP( p_nh_metrics%zgpot_ifc(jc,jk,jb) / h_scal_bg ) * &
              & ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )
            
            ! Reference Exner pressure, half level mass points
            z_help(jc) = ( z_aux1(jc) / p0ref )**rd_o_cpd
            
            ! Reference temperature, half level mass points
            z_temp(jc) = ( t0sl_bg - del_t_bg ) + del_t_bg * &
              & EXP( -p_nh_metrics%zgpot_ifc(jc,jk,jb) / h_scal_bg )
            
            ! Reference density, half level mass points
            z_aux2(jc) = z_aux1(jc) / ( rd * z_temp(jc) )
            
            ! Reference Potential temperature, half level mass points
            p_nh_metrics%theta_ref_ic(jc,jk,jb) = z_temp(jc) / z_help(jc)
            
            ! First vertical derivative of reference Exner pressure, half level mass points
            p_nh_metrics%d_exner_dz_ref_ic(jc,jk,jb) = &
              & -grav / cpd / p_nh_metrics%theta_ref_ic(jc,jk,jb) * & 
              & p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%gradh)**2
          ENDDO  !jc
        ENDDO  !jk
        
      ENDDO  !jb
!$OMP END DO NOWAIT

      i_startblk = p_patch%edges%start_block(2)
      i_endblk   = nblks_e

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, z_aux1, z_temp) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            
            ! Reference pressure, full level edge points
            z_aux1(je) = p0sl_bg * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
              & LOG( ( EXP( zgpot_me(je,jk,jb) / h_scal_bg ) * &
              & ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )
              
            ! Reference temperature, full level edge points
            z_temp(je) = ( t0sl_bg - del_t_bg ) + del_t_bg * &
              & EXP( -zgpot_me(je,jk,jb) / h_scal_bg )
            
            ! Reference density, full level edge points
            p_nh_metrics%rho_ref_me(je,jk,jb) = z_aux1(je) / ( rd * z_temp(je) )
            
            ! Reference potential temperature, full level edge points
            p_nh_metrics%theta_ref_me(je,jk,jb) = z_temp(je) / ( ( z_aux1(je) / p0ref )**rd_o_cpd )
          ENDDO  !je
        ENDDO  !jk
        
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      
      DEALLOCATE(z_me, zgpot_me, z_help, z_temp, z_aux1, z_aux2, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation of auxiliary fields failed!')         
      
    ENDIF  !IF (.NOT. lconstgrav)

    ! Note: quantities related to damping, diffusion and the like 
    ! are not modified for the deep atmosphere
    
    !-----------------------------------------------------
    !                  Message output       
    !-----------------------------------------------------
    
    IF (msg_level >= imsg_thr%high .AND. my_process_is_stdio() .AND. jg == 1) THEN        
      WRITE(message_text,'(a)') 'Metrical modification factors for the deep-atmosphere configuration'// &
        & ' (domain '//TRIM(int2string(jg))//'):'
      CALL message(TRIM(routine), TRIM(message_text))
      ! For full levels
      CALL message(TRIM(routine), 'At full levels:')
      ! Set up table
      CALL initialize_table(table)
      ! Set up table columns
      CALL add_table_column(table, column_jk)
      CALL add_table_column(table, column_height)
      CALL add_table_column(table, column_gradh)
      CALL add_table_column(table, column_divh)
      CALL add_table_column(table, column_divzU)
      CALL add_table_column(table, column_divzL)
      CALL add_table_column(table, column_vol)
      CALL add_table_column(table, column_invr)
      CALL add_table_column(table, column_centri)
      DO jk = 1, nlev
        ! Correct for differences in nested domains' vertical index
        ! (if message output for jg > 1 should be desired someday)
        jk1 = jk + p_patch%nshift_total
        ! Height of upper half level of full level 'jk'
        z_z_u = vct_a(jk1) 
        ! Height of lower half level of full level 'jk' 
        z_z_l = vct_a(jk1+1)
        ! Height of full level                      
        z_z_c = 0.5_wp * ( z_z_l + z_z_u )
        ! Fill the table rows
        CALL set_table_entry(table, jk, column_jk,     TRIM(int2string(jk)))
        CALL set_table_entry(table, jk, column_height, TRIM(real2string(z_z_c,'(F13.5)')))
        CALL set_table_entry(table, jk, column_gradh,  &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh))))
        CALL set_table_entry(table, jk, column_divh,   &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%divh))))
        CALL set_table_entry(table, jk, column_divzU,  &
          & TRIM(real2string(p_nh_metrics%deepatmo_t2mc(idamtr%t2mc%divzU,jk))))
        CALL set_table_entry(table, jk, column_divzL,  &
          & TRIM(real2string(p_nh_metrics%deepatmo_t2mc(idamtr%t2mc%divzL,jk))))
        CALL set_table_entry(table, jk, column_vol,    &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol))))
        CALL set_table_entry(table, jk, column_invr,   &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%invr))))
        CALL set_table_entry(table, jk, column_centri, &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1mc(jk,idamtr%t1mc%centri))))
      ENDDO  !jk
      ! Print table
      CALL print_table(table)
      ! Destruct table
      CALL finalize_table(table)
      ! Likewise for half levels
      CALL message(TRIM(routine), 'At half levels:')
      ! Set up table
      CALL initialize_table(table)
      ! Set up table columns
      CALL add_table_column(table, column_jk)
      CALL add_table_column(table, column_height)
      CALL add_table_column(table, column_gradh)
      CALL add_table_column(table, column_invr)
      CALL add_table_column(table, column_centri)
      DO jk = 1, nlevp1
        ! Correct for differences in nested domains' vertical index
        jk1 = jk + p_patch%nshift_total
        ! Height of upper half level of full level 'jk'            
        z_z_u = vct_a(jk1)     
        ! Fill the table rows
        CALL set_table_entry(table, jk, column_jk,     TRIM(int2string(jk)))
        CALL set_table_entry(table, jk, column_height, TRIM(real2string(z_z_u,'(F13.5)')))
        CALL set_table_entry(table, jk, column_gradh,  &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%gradh))))
        CALL set_table_entry(table, jk, column_invr,   &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%invr))))
        CALL set_table_entry(table, jk, column_centri, &
          & TRIM(real2string(p_nh_metrics%deepatmo_t1ifc(jk,idamtr%t1ifc%centri))))
      ENDDO  !jk
      ! Print table
      CALL print_table(table)
      ! Destruct table
      CALL finalize_table(table)   
    ENDIF  !IF (msg_level >= imsg_thr%high .AND. my_process_is_stdio() .AND. jg == 1)
    
  END SUBROUTINE set_deepatmo_metrics

  !--------------------------------------------------------------------------------------------------

  !>
  !!
  !! Discretization of nonhydrostatic momentum equation similar to hydrostatic core
  !! In particular, the Lamb transformation is applied only to the horizontal
  !! equation of motion, whereas the vertical wind equation is discretized
  !! in advective form
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl (2010-02-03)
  !! Modified by Sebastian Borchert, DWD (2017-07-24)
  !! - Upper atmosphere copy of 'src/atm_dyn_iconam/mo_velocity_advection: velocity_tendencies'
  !! - Note: any kind of damping, diffusion and the like is currently not 
  !!   modified for the deep atmosphere!
  !! - NOTE: Open-ACC parallelization of original subroutine has been removed!
  !!
  SUBROUTINE velocity_tendencies_deepatmo (p_prog, p_patch, p_int, p_metrics, p_diag, z_w_concorr_me, z_kin_hor_e, &
                                           z_vt_ie, ntnd, istep, lvn_only, dtime, opt_nrdmax)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(INOUT)       :: p_prog
    TYPE(t_nh_metrics), INTENT(INOUT)    :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)       :: p_diag

    ! Local variables from solve_nh that are passed for efficiency optimization
    REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie

    INTEGER, INTENT(IN)  :: ntnd     ! time level of ddt_adv fields used to store tendencies
    INTEGER, INTENT(IN)  :: istep    ! 1: predictor step, 2: corrector step
    LOGICAL, INTENT(IN)  :: lvn_only ! true: compute only vn tendency
    REAL(wp),INTENT(IN)  :: dtime    ! time step
    ! 'nrdmax' cannot be included from 'mo_vertical_grid' within this module, 
    ! because of circular dependencies
    INTEGER, INTENT(IN), OPTIONAL :: opt_nrdmax(:) 

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_startblk_2, i_endblk_2, i_startidx_2, i_endidx_2
    INTEGER :: rl_start, rl_end, rl_start_2, rl_end_2
    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
    REAL(vp):: z_w_concorr_mc(nproma,p_patch%nlev)
    REAL(vp):: z_w_con_c(nproma,p_patch%nlevp1)
    REAL(vp):: z_w_con_c_full(nproma,p_patch%nlev,p_patch%nblks_c)
    ! These fields in addition have reversed index order (vertical first) for optimization
#ifdef __LOOP_EXCHANGE
    REAL(vp):: z_v_grad_w(p_patch%nlev,nproma,p_patch%nblks_e)
    REAL(vp):: z_w_v(p_patch%nlevp1,nproma,p_patch%nblks_v)
    REAL(vp):: zeta(p_patch%nlev,nproma,p_patch%nblks_v)
    REAL(vp):: z_ekinh(p_patch%nlev,nproma,p_patch%nblks_c)
#else
    REAL(vp):: z_v_grad_w(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(vp):: z_w_v(nproma,p_patch%nlevp1,p_patch%nblks_v)
    REAL(vp):: zeta(nproma,p_patch%nlev,p_patch%nblks_v)
    REAL(vp):: z_ekinh(nproma,p_patch%nlev,p_patch%nblks_c)
#endif

    ! Pointers
    INTEGER, DIMENSION(:,:,:), POINTER   &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS                       &
#endif
      ::                                 &
      icidx, icblk, ieidx, ieblk, iqidx, iqblk, ividx, ivblk, incidx, incblk

    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    ! Local control variable for vertical nesting
    LOGICAL :: l_vert_nested

    INTEGER :: jg

    ! Variables for conditional additional diffusion for vertical advection
    REAL(vp) :: cfl_w_limit, vcfl, maxvcfl, vcflmax(p_patch%nblks_c)
    REAL(wp) :: w_con_e, scalfac_exdiff, difcoef, max_vcfl_dyn
                
    INTEGER  :: ie, nrdmax_jg, nflatlev_jg
    LOGICAL  :: levmask(p_patch%nblks_c,p_patch%nlev),levelmask(p_patch%nlev)
    LOGICAL  :: cfl_clipping(nproma,p_patch%nlevp1)   ! CFL > 0.85

    ! (upper-atmosphere/deep-atmosphere-related variables, 
    ! unfortunately, for full and half levels separately, 
    ! since it might be bad to change the pointer within an open-mp thread)
    REAL(wp), DIMENSION(:), POINTER :: deepatmo_gradh_mc, deepatmo_gradh_ifc, &
      &                                deepatmo_invr_mc, deepatmo_invr_ifc

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &    
      routine = modname//':velocity_tendencies_deepatmo'

    !--------------------------------------------------------------------------

    IF (timers_level > 5) CALL timer_start(timer_solve_nh_veltend)

    IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF

    ! in this temporary, pure deep-atmosphere version of the subroutine, 
    ! 'opt_nrdmax' should be input
    IF (.NOT. PRESENT(opt_nrdmax)) CALL finish( TRIM(routine),'opt_nrdmax should be present')

    !Get patch id
    jg = p_patch%id
    nrdmax_jg     = opt_nrdmax(jg)
    nflatlev_jg   = nflatlev(jg)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells/edges/vertices
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    incidx => p_patch%cells%neighbor_idx
    incblk => p_patch%cells%neighbor_blk

    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    ! (deep atmosphere: associate pointers for metrical modification factors)
    deepatmo_gradh_mc  => p_metrics%deepatmo_t1mc(:,idamtr%t1mc%gradh)
    deepatmo_invr_mc   => p_metrics%deepatmo_t1mc(:,idamtr%t1mc%invr)
    deepatmo_gradh_ifc => p_metrics%deepatmo_t1ifc(:,idamtr%t1ifc%gradh)
    deepatmo_invr_ifc  => p_metrics%deepatmo_t1ifc(:,idamtr%t1ifc%invr)

    ! Limit on vertical CFL number for applying extra diffusion
    IF (lextra_diffu) THEN
      cfl_w_limit = 0.65_wp/dtime   ! this means 65% of the nominal CFL stability limit

      ! Scaling factor for extra diffusion
      scalfac_exdiff = 0.05_wp / ( dtime*(0.85_wp - cfl_w_limit*dtime) )
    ELSE
      cfl_w_limit = 0.85_wp/dtime   ! this means 65% of the nominal CFL stability limit
      scalfac_exdiff = 0._wp
    ENDIF

    ! Compute w at vertices
    IF (.NOT. lvn_only) CALL cells2verts_scalar_ri(p_prog%w, p_patch, &
      p_int%cells_aw_verts, z_w_v, opt_rlend=min_rlvert_int-1)

    ! Compute vertical vorticity component at vertices
    ! (deep atmosphere: 'rot_vertex_ri' is not modified for spherical geometry, 
    ! but rather the output field is multiplied by a modification factor later on)    
    CALL rot_vertex_ri (p_prog%vn, p_patch, p_int, zeta, opt_rlend=min_rlvert_int-1)

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk, rl_start_2, rl_end_2, i_startblk_2, i_endblk_2)

    IF (istep == 1) THEN ! Computations of velocity-derived quantities that come from solve_nh in istep=2

      rl_start = 5
      rl_end = min_rledge_int - 2

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! RBF reconstruction of tangential wind component
            p_diag%vt(je,jk,jb) = &
              p_int%rbf_vec_coeff_e(1,je,jb) * p_prog%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) + &
              p_int%rbf_vec_coeff_e(2,je,jb) * p_prog%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) + &
              p_int%rbf_vec_coeff_e(3,je,jb) * p_prog%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) + &
              p_int%rbf_vec_coeff_e(4,je,jb) * p_prog%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
          ENDDO
        ENDDO

        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        DO jk = 2, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,jk,jb) =                                    &
              p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
            z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)**2 + p_diag%vt(je,jk,jb)**2)
          ENDDO
        ENDDO

        IF (.NOT. lvn_only) THEN ! Interpolate also vt to interface levels
          DO jk = 2, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              z_vt_ie(je,jk,jb) =                                         &
                p_metrics%wgtfac_e(je,jk,jb)*p_diag%vt(je,jk,jb) +        &
               (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_diag%vt(je,jk-1,jb)
            ENDDO
          ENDDO
        ENDIF

        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)
        ! (deep atmosphere: 'ddx(n)(t)_z_full' are modified for spherical geometry 
        ! in 'src/atm_dyn_iconam/mo_nh_deepatmo_init_utils: set_deepatmo_metrics')
        DO jk = nflatlev_jg, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                              &
              p_prog%vn(je,jk,jb)*p_metrics%ddxn_z_full(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO

        IF (.NOT. l_vert_nested) THEN
          ! Top and bottom levels
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            ! Quadratic extrapolation at the top turned out to cause numerical instability in pathological cases,
            ! thus we use a no-gradient condition in the upper half layer
            p_diag%vn_ie(je,1,jb) = p_prog%vn(je,1,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
        ELSE
          ! vn_ie(jk=1) is extrapolated using parent domain information in this case
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,1,jb) = p_diag%vn_ie(je,2,jb) + p_diag%dvn_ie_ubc(je,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF ! istep = 1


    rl_start = 7
    rl_end = min_rledge_int - 1

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    IF (.NOT. lvn_only) THEN
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Compute v*grad w on edges (level nlevp1 is not needed because w(nlevp1) is diagnostic)
        ! Note: this implicitly includes a minus sign for the gradients, which is needed later on
        ! (deep atmosphere: horizontal gradient (dw/dn,dw/dt) is multiplied by a modification factor 
        ! to account for the full spherical geometry, in addition metrical terms and 
        ! Coriolis acceleration have been added: 
        ! vn * dw/dn + vt * dw/dt -> vn * ( dw/dn - vn / r + ft ) + vt * ( dw/dt - vt / r - fn ), 
        ! with horizontal Coriolis parameters fn and ft, and radius r)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
            z_v_grad_w(jk,je,jb) = p_diag%vn_ie(je,jk,jb) *                                               &
             (                                                                                            &
              p_patch%edges%inv_dual_edge_length(je,jb) *                                                 & 
              (p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) - p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2))) * & 
              deepatmo_gradh_ifc(jk) + p_diag%vn_ie(je,jk,jb) * deepatmo_invr_ifc(jk)                     &
              - p_patch%edges%ft_e(je,jb)                                                                 &
             )                                                                                            &                                          
             + z_vt_ie(je,jk,jb) *                                                                        &
             (                                                                                            &
              p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%tangent_orientation(je,jb) *    &
              (z_w_v(jk,ividx(je,jb,1),ivblk(je,jb,1)) - z_w_v(jk,ividx(je,jb,2),ivblk(je,jb,2))) *       &
              deepatmo_gradh_ifc(jk) + z_vt_ie(je,jk,jb) * deepatmo_invr_ifc(jk)                          &
              + p_patch%edges%fn_e(je,jb)                                                                 &
             )
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_v_grad_w(je,jk,jb) = p_diag%vn_ie(je,jk,jb) *                                               &
             (                                                                                            &
              p_patch%edges%inv_dual_edge_length(je,jb) *                                                 &
              (p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) - p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2))) * &
              deepatmo_gradh_ifc(jk) + p_diag%vn_ie(je,jk,jb) * deepatmo_invr_ifc(jk)                     &
              - p_patch%edges%ft_e(je,jb)                                                                 & 
             )                                                                                            &
             + z_vt_ie(je,jk,jb) *                                                                        &
             (                                                                                            &
              p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%tangent_orientation(je,jb) *    &
              (z_w_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) - z_w_v(ividx(je,jb,2),jk,ivblk(je,jb,2))) *       &
              deepatmo_gradh_ifc(jk) + z_vt_ie(je,jk,jb) * deepatmo_invr_ifc(jk)                          &
              + p_patch%edges%fn_e(je,jb)                                                                 &
             )                    
#endif

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF

    rl_start = 4
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    rl_start_2 = grf_bdywidth_c+1
    rl_end_2   = min_rlcell_int

    i_startblk_2 = p_patch%cells%start_block(rl_start_2)
    i_endblk_2   = p_patch%cells%end_block(rl_end_2)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, i_startidx_2, i_endidx_2, z_w_con_c, &
!$OMP            z_w_concorr_mc, difcoef, vcfl, maxvcfl, cfl_clipping) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate horizontal kinetic energy to cell centers
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev
        z_ekinh(jk,jc,jb) =  &
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
        z_ekinh(jc,jk,jb) =  &
#endif
          p_int%e_bln_c_s(jc,1,jb)*z_kin_hor_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
          p_int%e_bln_c_s(jc,2,jb)*z_kin_hor_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
          p_int%e_bln_c_s(jc,3,jb)*z_kin_hor_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

        ENDDO
      ENDDO

      IF (istep == 1) THEN

        ! Interpolate contravariant correction to cell centers ...
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = nflatlev_jg, nlev
#else
        DO jk = nflatlev_jg, nlev
          DO jc = i_startidx, i_endidx
#endif

            z_w_concorr_mc(jc,jk) =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

          ENDDO
        ENDDO

        ! ... and to interface levels
        ! Remark: computation of w_concorr_c at nlevp1 is needed in solve_nh only
        ! because this serves solely for setting the lower boundary condition for w
        DO jk = nflatlev_jg+1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_diag%w_concorr_c(jc,jk,jb) =                                &
              p_metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._vp - p_metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1) 
          ENDDO
        ENDDO

      ENDIF

      z_w_con_c(:,1:nlev) = p_prog%w(:,1:nlev,jb)
      z_w_con_c(:,nlevp1) = 0._wp

      ! Contravariant vertical velocity on w points and interpolation to full levels
      DO jk = nlev, nflatlev_jg+1, -1
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_w_con_c(jc,jk) = z_w_con_c(jc,jk) - p_diag%w_concorr_c(jc,jk,jb)
        ENDDO
      ENDDO

      ! Search for grid points for which w_con is close to or above the CFL stability limit
      ! At these points, additional diffusion is applied in order to prevent numerical 
      ! instability if lextra_diffu = .TRUE.
      DO jk = MAX(3,nrdmax_jg-2), nlev-3
        levmask(jb,jk) = .FALSE.
      ENDDO

      maxvcfl = 0
      DO jk = MAX(3,nrdmax_jg-2), nlev-3
        DO jc = i_startidx, i_endidx
          cfl_clipping(jc,jk) = (ABS(z_w_con_c(jc,jk)) > cfl_w_limit*p_metrics%ddqz_z_half(jc,jk,jb))
          IF ( cfl_clipping(jc,jk) ) THEN
            levmask(jb,jk) = .TRUE.
            vcfl = z_w_con_c(jc,jk)*dtime/p_metrics%ddqz_z_half(jc,jk,jb)
            maxvcfl = MAX( maxvcfl, ABS( vcfl ) )
            !
            ! limit w_con to 85% of the nominal CFL stability threshold
            IF (vcfl < -0.85_vp) THEN
              z_w_con_c(jc,jk)           = -0.85_vp*p_metrics%ddqz_z_half(jc,jk,jb)/dtime
            ELSE IF (vcfl > 0.85_vp) THEN
              z_w_con_c(jc,jk)           = 0.85_vp*p_metrics%ddqz_z_half(jc,jk,jb)/dtime
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      vcflmax(jb) = maxvcfl

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          z_w_con_c_full(jc,jk,jb) = 0.5_vp*(z_w_con_c(jc,jk)+z_w_con_c(jc,jk+1))
        ENDDO
      ENDDO

      ! The remaining computations are not needed in vn_only mode and only on prognostic grid points
      IF (lvn_only) CYCLE
      IF (jb < i_startblk_2 .OR. jb > i_endblk_2) CYCLE

      CALL get_indices_c(p_patch, jb, i_startblk_2, i_endblk_2, &
                         i_startidx_2, i_endidx_2, rl_start_2, rl_end_2)


      ! Compute vertical derivative terms of vertical wind advection
      DO jk = 2, nlev
!DIR$ IVDEP
        DO jc = i_startidx_2, i_endidx_2
          p_diag%ddt_w_adv(jc,jk,jb,ntnd) =  - z_w_con_c(jc,jk)   *                                 &
            (p_prog%w(jc,jk-1,jb)*p_metrics%coeff1_dwdz(jc,jk,jb) -                                 &
             p_prog%w(jc,jk+1,jb)*p_metrics%coeff2_dwdz(jc,jk,jb) +                                 &
             p_prog%w(jc,jk,jb)*(p_metrics%coeff2_dwdz(jc,jk,jb) - p_metrics%coeff1_dwdz(jc,jk,jb)) )
        ENDDO
      ENDDO

      ! Interpolate horizontal advection of w from edges to cells and add to advective tendency
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx_2, i_endidx_2
!DIR$ IVDEP
        DO jk = 2, nlev
          p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd)       + &
            p_int%e_bln_c_s(jc,1,jb)*z_v_grad_w(jk,ieidx(jc,jb,1),ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_v_grad_w(jk,ieidx(jc,jb,2),ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_v_grad_w(jk,ieidx(jc,jb,3),ieblk(jc,jb,3))
#else
      DO jk = 2, nlev
        DO jc = i_startidx_2, i_endidx_2
          p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd)       + &
            p_int%e_bln_c_s(jc,1,jb)*z_v_grad_w(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_v_grad_w(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_v_grad_w(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))
#endif
        ENDDO
      ENDDO

      IF (lextra_diffu) THEN
        ! Apply extra diffusion at grid points where w_con is close to or above the CFL stability limit
        DO jk = MAX(3,nrdmax_jg-2), nlev-3
          IF (levmask(jb,jk)) THEN
            DO jc = i_startidx_2, i_endidx_2
              IF (cfl_clipping(jc,jk) .AND. p_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
                difcoef = scalfac_exdiff * MIN(0.85_wp - cfl_w_limit*dtime,                       &
                  ABS(z_w_con_c(jc,jk))*dtime/p_metrics%ddqz_z_half(jc,jk,jb) - cfl_w_limit*dtime )

                ! nabla2 diffusion on w
                p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd)        + &
                  difcoef * p_patch%cells%area(jc,jb) * (                                  &
                  p_prog%w(jc,jk,jb)                          *p_int%geofac_n2s(jc,1,jb) + &
                  p_prog%w(incidx(jc,jb,1),jk,incblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
                  p_prog%w(incidx(jc,jb,2),jk,incblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
                  p_prog%w(incidx(jc,jb,3),jk,incblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb)   )

              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jk)
    DO jk = MAX(3,nrdmax_jg-2), nlev-3
      levelmask(jk) = ANY(levmask(i_startblk:i_endblk,jk))
    ENDDO
!$OMP END DO

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, ie, w_con_e, difcoef) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection: grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
      ! (deep atmosphere: grad(Ekin_h) is multiplied by metrical modification factor to account 
      ! for spherical geometry, in addition metrical terms and contribution of vertical wind to 
      ! Coriolis acceleration have been added: wcon_e * dvn/dz -> wcon_e * ( dvn/dz + vn / r - ft ), 
      ! where r is radius and ft is tangential component of horizontal Coriolis parameter. 
      ! The vorticity 'zeta' has to be multiplied by 'deepatmo_gradh_mc(jk)', 
      ! because the subroutine 'rot_vertex_ri', which computes 'zeta', is itself not 
      ! modified for spherical geometry.)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
        DO jk = 1, nlev
          p_diag%ddt_vn_adv(je,jk,jb,ntnd) = - (                                                & 
             (                                                                                  & 
              z_kin_hor_e(je,jk,jb) *                                                           &
              (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +         &
              p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(jk,icidx(je,jb,2),icblk(je,jb,2)) -     &
              p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(jk,icidx(je,jb,1),icblk(je,jb,1))       &
             ) * deepatmo_gradh_mc(jk)                                                          &
             + p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                       &
             (zeta(jk,ividx(je,jb,1),ivblk(je,jb,1)) + zeta(jk,ividx(je,jb,2),ivblk(je,jb,2)))  &
             * deepatmo_gradh_mc(jk) )                                                          &
             + (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +       &
             p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))) *         & 
             (                                                                                  &
              (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb) &
              + p_prog%vn(je,jk,jb) * deepatmo_invr_mc(jk)                                      & 
              - p_patch%edges%ft_e(je,jb)                                                       &
             )                                                                                  &
            )
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          p_diag%ddt_vn_adv(je,jk,jb,ntnd) = - (                                                &
             (                                                                                  &
              z_kin_hor_e(je,jk,jb) *                                                           &
              (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +         &
              p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -     &
              p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(icidx(je,jb,1),jk,icblk(je,jb,1))       &         
             ) * deepatmo_gradh_mc(jk)                                                          &
             + p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                       &
             (zeta(ividx(je,jb,1),jk,ivblk(je,jb,1)) + zeta(ividx(je,jb,2),jk,ivblk(je,jb,2)))  &
             * deepatmo_gradh_mc(jk) )                                                          &
             + (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +       &
             p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))) *         &
             (                                                                                  &
             (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb) & 
             + p_prog%vn(je,jk,jb) * deepatmo_invr_mc(jk)                                       &
             - p_patch%edges%ft_e(je,jb)                                                        &
             )                                                                                  &
            )
#endif

        ENDDO
      ENDDO


      IF (lextra_diffu) THEN
        ! Search for grid points for which w_con is close to or above the CFL stability limit
        ! At these points, additional diffusion is applied in order to prevent numerical instability
        ie = 0

        DO jk = MAX(3,nrdmax_jg-2), nlev-4
          IF (levelmask(jk) .OR. levelmask(jk+1)) THEN
            DO je = i_startidx, i_endidx
              w_con_e = p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) + &
                        p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))
              IF (ABS(w_con_e) > cfl_w_limit*p_metrics%ddqz_z_full_e(je,jk,jb)) THEN
                difcoef = scalfac_exdiff * MIN(0.85_wp - cfl_w_limit*dtime,                &
                  ABS(w_con_e)*dtime/p_metrics%ddqz_z_full_e(je,jk,jb) - cfl_w_limit*dtime )

                p_diag%ddt_vn_adv(je,jk,jb,ntnd) = p_diag%ddt_vn_adv(je,jk,jb,ntnd)   +                 &
                  difcoef * p_patch%edges%area_edge(je,jb) * (                                          &
                  p_int%geofac_grdiv(je,1,jb)*p_prog%vn(je,jk,jb)                         +             &
                  p_int%geofac_grdiv(je,2,jb)*p_prog%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) +             &
                  p_int%geofac_grdiv(je,3,jb)*p_prog%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) +             &
                  p_int%geofac_grdiv(je,4,jb)*p_prog%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) +             &
                  p_int%geofac_grdiv(je,5,jb)*p_prog%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4)) +             &
                  p_patch%edges%tangent_orientation(je,jb)*p_patch%edges%inv_primal_edge_length(je,jb) * &
#ifdef __LOOP_EXCHANGE
                  (zeta(jk,ividx(je,jb,2),ivblk(je,jb,2)) - zeta(jk,ividx(je,jb,1),ivblk(je,jb,1))) )
#else
                  (zeta(ividx(je,jb,2),jk,ivblk(je,jb,2)) - zeta(ividx(je,jb,1),jk,ivblk(je,jb,1))) )
#endif
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO 
!$OMP END PARALLEL

    ! Save maximum vertical CFL number for substep number adaptation
    i_startblk = p_patch%cells%start_block(grf_bdywidth_c)
    i_endblk   = p_patch%cells%end_block(min_rlcell_int)

    max_vcfl_dyn = MAX(p_diag%max_vcfl_dyn,MAXVAL(vcflmax(i_startblk:i_endblk)))
    p_diag%max_vcfl_dyn = max_vcfl_dyn

    NULLIFY( deepatmo_gradh_mc, deepatmo_invr_mc,  &
      &      deepatmo_gradh_ifc, deepatmo_invr_ifc )

    IF (timers_level > 5) CALL timer_stop(timer_solve_nh_veltend)

  END SUBROUTINE velocity_tendencies_deepatmo

END MODULE mo_nh_deepatmo_utils
