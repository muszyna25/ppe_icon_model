!>
!! Subroutine to initialize testcase 12 (Hadley-like meridional circulation) proposed 
!! for the DCMIP summer school
!!
!!
!!
!! @par Revision History
!! - initial revision by Daniel Reinert, DWD, (2012-05-30)
!!
!! @par Literature
!! - Dynamical Core Model Intercomparison Project (DCMIP) 
!!   Test Case Document (P. Ullrich et al, 2012)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nh_dcmip_hadley

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, grav, p0ref, cpd, cvd_o_rd
   USE mo_math_constants,       ONLY: pi
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge
   USE mo_parallel_config,      ONLY: nproma
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
   USE mo_model_domain,         ONLY: t_patch
   USE mo_grid_config,          ONLY: grid_sphere_radius
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_sync,                 ONLY: sync_patch_array, SYNC_E

   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_dcmip_hadley
   PUBLIC :: set_nh_velocity_hadley


   ! test case parameters
   !
   REAL(wp), PARAMETER :: t0    = 300.0_wp    !< temperature           [K]
   REAL(wp), PARAMETER :: scale_hgt = rd*t0/grav !< scale height       [m] 
   REAL(wp), PARAMETER :: u0    = 40.0_wp     !< maximum amplitude 
                                              !< of the zonal wind     [m s^-1]
   REAL(wp), PARAMETER :: w0    = 0.15_wp     !< maximum amplitude 
                                              !< of the vertical wind  [m s^-1]
   REAL(wp), PARAMETER :: z1    = 2000.0_wp   !< lower tracer bound    [m]
   REAL(wp), PARAMETER :: z2    = 5000.0_wp   !< upper tracer bound    [m]
   REAL(wp), PARAMETER :: z_mid = 0.5_wp*(z1+z2) !< midpoint              [m]
   REAL(wp), PARAMETER :: ztop  = 12000._wp   !< model top             [m]
   REAL(wp), PARAMETER :: tau   = 1.0_wp * 86400.0_wp ! period of motion 1 day [s]
   INTEGER , PARAMETER :: nhadley = 5         !< number of overturning hadley cells

!--------------------------------------------------------------------

CONTAINS
!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the DCMIP Hadley-like 
  !! meridional circulation test 
  !!
  !! @par Revision History
  !! - initial revision by Daniel Reinert, DWD (2012-05-30)
  !!
  SUBROUTINE init_nh_dcmip_hadley( p_patch, p_nh_prog, p_nh_diag, &
    &                         p_int, p_metrics )

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state vector
      &  p_int

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics


    INTEGER  :: jc, jk, jb                    !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1                  !< number of full/half levels

 !--------------------------------------------------------------------


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    i_rlstart = 1
    i_rlend   = min_rlcell

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)



    !
    ! Init prognostic variables
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! init tracer field
          !
          IF (p_metrics%z_mc(jc,jk,jb) < z2 .AND. p_metrics%z_mc(jc,jk,jb) > z1) THEN
            p_nh_prog%tracer(jc,jk,jb,1) = 0.5_wp * (1.0_wp + cos( 2.0_wp*pi      &
              &                          *(p_metrics%z_mc(jc,jk,jb)-z_mid)/(z2-z1) ) )
          ELSE
            p_nh_prog%tracer(jc,jk,jb,1) = 0.0_wp
          ENDIF

! checks whether the specified velocity field is non-divergent
!Test         p_nh_prog%tracer(jc,jk,jb,2) = 1._wp


          ! temperature is constant
          !
          p_nh_diag%temp(jc,jk,jb) = t0


          ! init pressure field
          !
          p_nh_diag%pres(jc,jk,jb) = p0ref * exp(-p_metrics%z_mc(jc,jk,jb)/scale_hgt)


          ! init exner pressure
          !
          p_nh_prog%exner(jc,jk,jb) = (p_nh_diag%pres(jc,jk,jb)/p0ref)**(rd/cpd)


          ! init virtual potential temperature
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb)/p_nh_prog%exner(jc,jk,jb)


          ! init density of moist air
          !
          p_nh_prog%rho(jc,jk,jb) = p_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref  &
            &                     /rd/p_nh_prog%theta_v(jc,jk,jb)

        ENDDO  ! jc
      ENDDO  ! jk



      ! half level initialization
      !
      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          ! init pressure field (half levels)
          !
          p_nh_diag%pres_ifc(jc,jk,jb) = p0ref * exp(-p_metrics%z_ifc(jc,jk,jb)/scale_hgt)

          ! init density field (half levels)
          !
          p_nh_diag%rho_ic(jc,jk,jb)   = p_nh_diag%pres_ifc(jc,jk,jb) / (rd * t0)

        ENDDO  ! jc
      ENDDO  ! jk

    ENDDO ! jb
!$OMP ENDDO
!$OMP END PARALLEL



   ! set initial velocity field for t=0s
   CALL set_nh_velocity_hadley( p_patch, p_nh_prog, p_nh_diag, p_int, &
     &                          p_metrics, 0.0_wp )


  END SUBROUTINE init_nh_dcmip_hadley

!--------------------------------------------------------------------


!--------------------------------------------------------------------

  !>
  !! Initialization of horizontal and vertical velocity field for 
  !! the DCMIP Hadley-like meridional circulation test 
  !!
  !! @par Revision History
  !! - initial revision by Daniel Reinert, DWD (2012-05-30)
  !!
  SUBROUTINE set_nh_velocity_hadley( p_patch, p_nh_prog, p_nh_diag, p_int,  &
    &                                p_metrics, time )

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state vector
      &  p_int

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    REAL(wp), INTENT(IN) :: time              !< simulation time   [s]

    REAL(wp) ::           &                   !< zonal and meridional velocity
      & zu(nproma), zv(nproma)

    REAL(wp) ::        &                      !< geometric height at edge points
      &  z_me(nproma,p_patch%nlev,p_patch%nblks_e), &
      &  z_ife(nproma,p_patch%nlevp1,p_patch%nblks_e)

    REAL(wp) ::        &                      !< density at edge points
      &  z_rho_e(nproma,p_patch%nlev,p_patch%nblks_e), &
      &  z_rho_ie(nproma,p_patch%nlevp1,p_patch%nblks_e)

    REAL(wp) :: z_lat                         !< geographical latitude
    REAL(wp) :: ztop                          !< model top
    INTEGER  :: jc, je, jk, jb                !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1                  !< number of full and half levels

    !-----------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! Compute geometric height of full levels at edge midpoints
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me)

    ! Compute geometric height of half levels at edge midpoints
    CALL cells2edges_scalar(p_metrics%z_ifc, p_patch, p_int%c_lin_e, z_ife)


    ! Compute rho at full level edge midpoints
    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, z_rho_e)

    ! Compute rho at half level at edge midpoints
    CALL cells2edges_scalar(p_nh_diag%rho_ic, p_patch, p_int%c_lin_e, z_rho_ie)

    ! syncs
    CALL sync_patch_array(SYNC_E, p_patch, z_me)
    CALL sync_patch_array(SYNC_E, p_patch, z_ife)
    CALL sync_patch_array(SYNC_E, p_patch, z_rho_e)
    CALL sync_patch_array(SYNC_E, p_patch, z_rho_ie)

    !
    ! set normal velocity field
    !

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_lat,zu,zv,ztop)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO je = i_startidx, i_endidx


          z_lat  = p_patch%edges%center(je,jb)%lat

          zu(je) = u0*cos(z_lat)

          ztop   = z_ife(je,1,jb)

          zv(je) = -(z_rho_ie(je,nlevp1,jb)/z_rho_e(je,jk,jb))               &
            &    * (grid_sphere_radius*w0*pi) /(nhadley*ztop)                &
            &    * cos(z_lat)*sin(nhadley*z_lat)*cos(pi*z_me(je,jk,jb)/ztop) &
            &    * cos(pi*time/tau)

          p_nh_prog%vn(je,jk,jb) = zu(je) * p_patch%edges%primal_normal(je,jb)%v1  &
            &                    + zv(je) * p_patch%edges%primal_normal(je,jb)%v2

        ENDDO  !je
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO


    !
    ! set vertical velocity field
    !

    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          z_lat = p_patch%cells%center(jc,jb)%lat

          p_nh_prog%w(jc,jk,jb) = (w0/nhadley)                                    &
            &       * (p_nh_diag%rho_ic(jc,nlevp1,jb)/p_nh_diag%rho_ic(jc,jk,jb)) & 
            &       * (-2.0_wp*SIN(nhadley*z_lat)*SIN(z_lat)                      &
            &       + nhadley*COS(z_lat)*COS(nhadley*z_lat))                      &
            &       * SIN(pi*p_metrics%z_ifc(jc,jk,jb)/p_metrics%z_ifc(jc,1,jb))  &
            &       * COS(pi*time/tau)

        ENDDO  !jc
      ENDDO  !jk

      p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb) = 0._wp

    ENDDO  !jb
!$OMP ENDDO
!$OMP END PARALLEL


  END SUBROUTINE set_nh_velocity_hadley


END MODULE mo_nh_dcmip_hadley
