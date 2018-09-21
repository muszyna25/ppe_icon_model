#if defined(__COMPAD_ADJLOOP__)
  !-------------------------------------------------------------------------
  SUBROUTINE restore_checkpoint(cp,jg,ocean_state, sea_ice_state, p_oce_sfc, atmos_fluxes, &
          &  current_time, time_config, output_file)
    INTEGER, INTENT(IN) :: cp   !   number of checkpoint
    INTEGER, INTENT(IN) :: jg
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_sea_ice) :: sea_ice_state
    TYPE(t_ocean_surface)  :: p_oce_sfc
    TYPE(t_atmos_fluxes)   :: atmos_fluxes
    TYPE(datetime), intent(inout) :: current_time
    TYPE(t_time_config) :: time_config
    TYPE(t_output_file) :: output_file(:)

    current_time                                    = ocean_state(jg)%p_check(cp)%datetime
    CALL set_OceanCurrentTime( current_time )
#ifdef __COMPAD_CPTEST__
    PRINT *, '  RESTORE_CP  ::::  ', cp, current_time
#endif /* __COMPAD_CPTEST__ */

    ocean_state(jg)%p_prog(nold(1))%h(:,:)          = ocean_state(jg)%p_check(cp)%h (:,:)  
    ocean_state(jg)%p_prog(nold(1))%vn(:,:,:)       = ocean_state(jg)%p_check(cp)%vn(:,:,:)
    ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1) = ocean_state(jg)%p_check(cp)%t (:,:,:)
    ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2) = ocean_state(jg)%p_check(cp)%s (:,:,:)

    ocean_state(jg)%p_prog(nnew(1))%h(:,:)          = ocean_state(jg)%p_check(cp)%h0 (:,:)  
    ocean_state(jg)%p_prog(nnew(1))%vn(:,:,:)       = ocean_state(jg)%p_check(cp)%vn0(:,:,:)
    ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,1) = ocean_state(jg)%p_check(cp)%t0 (:,:,:)
    ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2) = ocean_state(jg)%p_check(cp)%s0 (:,:,:)

    ocean_state(jg)%p_diag%w                        = ocean_state(jg)%p_check(cp)%w
    ocean_state(jg)%p_diag%zgrad_rho                = ocean_state(jg)%p_check(cp)%zgrad_rho

    ocean_state(jg)%p_aux%g_nm1                     = ocean_state(jg)%p_check(cp)%g_nm1
 
    current_time                                    = ocean_state(jg)%p_check(cp)%datetime
    time_config                                     = ocean_state(jg)%p_check(cp)%time_config

    ! sea ice variables
    sea_ice_state%hi(:,:,:)      = ocean_state(jg)%p_check(cp)%hi(:,:,:)
    sea_ice_state%hs(:,:,:)      = ocean_state(jg)%p_check(cp)%hs(:,:,:)
    sea_ice_state%hiold(:,:,:)   = ocean_state(jg)%p_check(cp)%hiold(:,:,:)
    sea_ice_state%hsold(:,:,:)   = ocean_state(jg)%p_check(cp)%hsold(:,:,:)
    sea_ice_state%conc(:,:,:)    = ocean_state(jg)%p_check(cp)%conc(:,:,:)
    sea_ice_state%zUnderIce(:,:) = ocean_state(jg)%p_check(cp)%zUnderIce(:,:)
    sea_ice_state%Tsurf(:,:,:)   = ocean_state(jg)%p_check(cp)%Tsurf(:,:,:)
    sea_ice_state%T1(:,:,:)      = ocean_state(jg)%p_check(cp)%T1(:,:,:)
    sea_ice_state%T2(:,:,:)      = ocean_state(jg)%p_check(cp)%T2(:,:,:)
    sea_ice_state%concSum        = ocean_state(jg)%p_check(cp)%concSum
    sea_ice_state%vol            = ocean_state(jg)%p_check(cp)%vol
    sea_ice_state%vols           = ocean_state(jg)%p_check(cp)%vols
    sea_ice_state%Tfw            = ocean_state(jg)%p_check(cp)%Tfw
    sea_ice_state%draft          = ocean_state(jg)%p_check(cp)%draft
    sea_ice_state%draftave       = ocean_state(jg)%p_check(cp)%draftave
    sea_ice_state%E1             = ocean_state(jg)%p_check(cp)%E1
    sea_ice_state%E2             = ocean_state(jg)%p_check(cp)%E2
    sea_ice_state%alb            = ocean_state(jg)%p_check(cp)%alb
    sea_ice_state%Qtop           = ocean_state(jg)%p_check(cp)%Qtop
    sea_ice_state%Qbot           = ocean_state(jg)%p_check(cp)%Qbot
    sea_ice_state%zHeatOceI      = ocean_state(jg)%p_check(cp)%zHeatOceI
    sea_ice_state%HeatOceI       = ocean_state(jg)%p_check(cp)%HeatOceI
    sea_ice_state%HeatOceW       = ocean_state(jg)%p_check(cp)%HeatOceW
    sea_ice_state%newice         = ocean_state(jg)%p_check(cp)%newice

    !   sea_ice_state%u_prog     = ocean_state(jg)%p_check(cp)%u_prog
    !   sea_ice_state%v_prog     = ocean_state(jg)%p_check(cp)%v_prog
    !   sea_ice_state%u          = ocean_state(jg)%p_check(cp)%u
    !   sea_ice_state%v          = ocean_state(jg)%p_check(cp)%v
    !   sea_ice_state%vn_e       = ocean_state(jg)%p_check(cp)%vn_e
    !   atmos_fluxes%stress_x    = ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_x
    !   atmos_fluxes%stress_y    = ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_y
    !   atmos_fluxes%stress_xw   = ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_xw
    !   atmos_fluxes%stress_yw   = ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_yw


    p_oce_sfc%SaltFlux_Relax     = ocean_state(jg)%p_check(cp)%SaltFlux_Relax
    p_oce_sfc%FrshFlux_Relax     = ocean_state(jg)%p_check(cp)%FrshFlux_Relax
    p_oce_sfc%HeatFlux_Relax     = ocean_state(jg)%p_check(cp)%HeatFlux_Relax
    p_oce_sfc%TempFlux_Relax     = ocean_state(jg)%p_check(cp)%TempFlux_Relax
    atmos_fluxes%SWnet           = ocean_state(jg)%p_check(cp)%  SWnet
    atmos_fluxes%lat             = ocean_state(jg)%p_check(cp)% lat
    atmos_fluxes%sens            = ocean_state(jg)%p_check(cp)% sens
    atmos_fluxes%LWnet           = ocean_state(jg)%p_check(cp)% LWnet
    atmos_fluxes%dlatdT          = ocean_state(jg)%p_check(cp)% dlatdT
    atmos_fluxes%dsensdT         = ocean_state(jg)%p_check(cp)% dsensdT
    atmos_fluxes%dLWdT           = ocean_state(jg)%p_check(cp)%dLWdT
    atmos_fluxes%albvisdir       = ocean_state(jg)%p_check(cp)%albvisdir
    atmos_fluxes%albvisdif       = ocean_state(jg)%p_check(cp)%albvisdif
    atmos_fluxes%albnirdir       = ocean_state(jg)%p_check(cp)%albnirdir
    atmos_fluxes%albnirdif       = ocean_state(jg)%p_check(cp)%albnirdif

  END SUBROUTINE restore_checkpoint


  
  SUBROUTINE store_checkpoints(cp,jg, ocean_state, sea_ice_state, p_oce_sfc, atmos_fluxes, &
          &  current_time,time_config, output_file)
    INTEGER, INTENT(IN) :: cp   !   number of checkpoint
    INTEGER, INTENT(IN) :: jg
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in) :: ocean_state(n_dom)
    TYPE (t_sea_ice),       INTENT(IN)  :: sea_ice_state
    TYPE(t_ocean_surface),  INTENT(IN)  :: p_oce_sfc
    TYPE(t_atmos_fluxes)    :: atmos_fluxes
    TYPE(datetime)                      :: current_time
    TYPE(t_time_config)                 :: time_config
    TYPE(t_output_file)                 :: output_file(:)

    integer :: status

#ifdef __COMPAD_CPTEST__
    PRINT *, '  STORE_CP    ::::  ', cp, current_time
#endif /* __COMPAD_CPTEST__ */
    ocean_state(jg)%p_check(cp)%h (:,:)                    = ocean_state(jg)%p_prog(nold(1))%h(:,:)
    ocean_state(jg)%p_check(cp)%vn(:,:,:)                  = ocean_state(jg)%p_prog(nold(1))%vn(:,:,:)
    ocean_state(jg)%p_check(cp)%t (:,:,:)                  = ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1)
    ocean_state(jg)%p_check(cp)%s (:,:,:)                  = ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2)

    ocean_state(jg)%p_check(cp)%h0 (:,:)                   = ocean_state(jg)%p_prog(nnew(1))%h(:,:)
    ocean_state(jg)%p_check(cp)%vn0(:,:,:)                 = ocean_state(jg)%p_prog(nnew(1))%vn(:,:,:)
    ocean_state(jg)%p_check(cp)%t0 (:,:,:)                 = ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,1)
    ocean_state(jg)%p_check(cp)%s0 (:,:,:)                 = ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2)

    ocean_state(jg)%p_check(cp)%w                          = ocean_state(jg)%p_diag%w
    ocean_state(jg)%p_check(cp)%zgrad_rho                  = ocean_state(jg)%p_diag%zgrad_rho
    ocean_state(jg)%p_check(cp)%g_nm1                      = ocean_state(jg)%p_aux%g_nm1

    !    ocean_state(jg)%p_check(cp)%datetime              = current_time
    ocean_state(jg)%p_check(cp)%datetime                   = get_OceanCurrentTime()
    ocean_state(jg)%p_check(cp)%time_config                = time_config

    ! sea ice variables
    ocean_state(jg)%p_check(cp)%hi(:,:,:)                  = sea_ice_state%hi(:,:,:)
    ocean_state(jg)%p_check(cp)%hs(:,:,:)                  = sea_ice_state%hs(:,:,:)
    ocean_state(jg)%p_check(cp)%hiold(:,:,:)               = sea_ice_state%hiold(:,:,:)
    ocean_state(jg)%p_check(cp)%hsold(:,:,:)               = sea_ice_state%hsold(:,:,:)
    ocean_state(jg)%p_check(cp)%conc(:,:,:)                = sea_ice_state%conc(:,:,:)
    ocean_state(jg)%p_check(cp)%zUnderIce(:,:)             = sea_ice_state%zUnderIce(:,:)
    ocean_state(jg)%p_check(cp)%Tsurf(:,:,:)               = sea_ice_state%Tsurf(:,:,:)
    ocean_state(jg)%p_check(cp)%T1(:,:,:)                  = sea_ice_state%T1(:,:,:)
    ocean_state(jg)%p_check(cp)%T2(:,:,:)                  = sea_ice_state%T2(:,:,:)

    ocean_state(jg)%p_check(cp)%concSum                    = sea_ice_state%concSum
    ocean_state(jg)%p_check(cp)%vol                        = sea_ice_state%vol
    ocean_state(jg)%p_check(cp)%vols                       = sea_ice_state%vols
    ocean_state(jg)%p_check(cp)%Tfw                        = sea_ice_state%Tfw
    ocean_state(jg)%p_check(cp)%draft                      = sea_ice_state%draft
    ocean_state(jg)%p_check(cp)%draftave                   = sea_ice_state%draftave
    ocean_state(jg)%p_check(cp)%E1                         = sea_ice_state%E1
    ocean_state(jg)%p_check(cp)%E2                         = sea_ice_state%E2
    ocean_state(jg)%p_check(cp)%alb                        = sea_ice_state%alb
    ocean_state(jg)%p_check(cp)%Qtop                       = sea_ice_state%Qtop
    ocean_state(jg)%p_check(cp)%Qbot                       = sea_ice_state%Qbot

    !   ocean_state(jg)%p_check(cp)%u_prog                 = sea_ice_state%u_prog
    !   ocean_state(jg)%p_check(cp)%v_prog                 = sea_ice_state%v_prog
    !   ocean_state(jg)%p_check(cp)%u                      = sea_ice_state%u
    !   ocean_state(jg)%p_check(cp)%v                      = sea_ice_state%v
    !   ocean_state(jg)%p_check(cp)%vn_e                   = sea_ice_state%vn_e
    !   ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_x  = atmos_fluxes%stress_x
    !   ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_y  = atmos_fluxes%stress_y
    !   ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_xw = atmos_fluxes%stress_xw
    !   ocean_state(jg)%p_check(cp)%atmos_fluxes_stress_yw = atmos_fluxes%stress_yw


    ! ocean surface stuff
    ocean_state(jg)%p_check(cp)%SaltFlux_Relax             = p_oce_sfc%SaltFlux_Relax
    ocean_state(jg)%p_check(cp)%FrshFlux_Relax             = p_oce_sfc%FrshFlux_Relax
    ocean_state(jg)%p_check(cp)%HeatFlux_Relax             = p_oce_sfc%HeatFlux_Relax
    ocean_state(jg)%p_check(cp)%TempFlux_Relax             = p_oce_sfc%TempFlux_Relax

    ocean_state(jg)%p_check(cp)%  SWnet                    = atmos_fluxes%SWnet
    ocean_state(jg)%p_check(cp)% lat                       = atmos_fluxes%lat
    ocean_state(jg)%p_check(cp)% sens                      = atmos_fluxes%sens
    ocean_state(jg)%p_check(cp)% LWnet                     = atmos_fluxes%LWnet
    ocean_state(jg)%p_check(cp)% dlatdT                    = atmos_fluxes%dlatdT
    ocean_state(jg)%p_check(cp)% dsensdT                   = atmos_fluxes%dsensdT
    ocean_state(jg)%p_check(cp)%dLWdT                      = atmos_fluxes%dLWdT
    ocean_state(jg)%p_check(cp)%albvisdir                  = atmos_fluxes%albvisdir
    ocean_state(jg)%p_check(cp)%albvisdif                  = atmos_fluxes%albvisdif
    ocean_state(jg)%p_check(cp)%albnirdir                  = atmos_fluxes%albnirdir
    ocean_state(jg)%p_check(cp)%albnirdif                  = atmos_fluxes%albnirdif

  END SUBROUTINE store_checkpoints

  SUBROUTINE compare_checkpoints(jg,cp1,cp2,strt,ocean_state)
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: cp1,cp2,strt   !   number of checkpoint
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    CHARACTER(*), PARAMETER        :: sfmt = '(A,I5,2x,I5,2x,I5,A,E30.17)'
!! COMPAD-DCO  BEGIN  JR2018
!!!! deactivate I/O
#ifdef __COMPAD__
    WRITE(0,sfmt) ' CPTEST ', strt, cp1,cp2, '  H     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%h%value  - ocean_state(jg)%p_check(cp2)%h%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt, cp1,cp2, '  VN    ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vn%value - ocean_state(jg)%p_check(cp2)%vn%value)**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  T     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%t%value  - ocean_state(jg)%p_check(cp2)%t%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  S     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%s%value  - ocean_state(jg)%p_check(cp2)%s%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  W     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%w%value  - ocean_state(jg)%p_check(cp2)%w%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  G_NM1 ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%g_nm1%value-ocean_state(jg)%p_check(cp2)%g_nm1%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  HI    ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%hi%value-ocean_state(jg)%p_check(cp2)%hi%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  HS    ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%hs%value-ocean_state(jg)%p_check(cp2)%hs%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  CONC  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%conc%value-ocean_state(jg)%p_check(cp2)%conc%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, 'ZUNDERIC', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%zUnderIce%value-ocean_state(jg)%p_check(cp2)%zUnderIce%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  ALB   ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%alb%value-ocean_state(jg)%p_check(cp2)%alb%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, 'SFL_REL ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%SaltFlux_Relax%value-ocean_state(jg)%p_check(cp2)%SaltFlux_Relax%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  TSURF ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Tsurf%value-ocean_state(jg)%p_check(cp2)%Tsurf%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  CSUM  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%concSum%value-ocean_state(jg)%p_check(cp2)%concSum%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  QTOP  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Qtop%value-ocean_state(jg)%p_check(cp2)%Qtop%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  QBOT  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Qbot%value-ocean_state(jg)%p_check(cp2)%Qbot%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '   VOL  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vol%value-ocean_state(jg)%p_check(cp2)%vol%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  VOLS  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vols%value-ocean_state(jg)%p_check(cp2)%vols%value )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  DRAFT ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%draft%value-ocean_state(jg)%p_check(cp2)%draft%value )**2))
#else /* __COMPAD__ */
    WRITE(0,sfmt) ' CPTEST ', strt, cp1,cp2, '  H     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%h  - ocean_state(jg)%p_check(cp2)%h )**2))
    WRITE(0,sfmt) ' CPTEST ', strt, cp1,cp2, '  VN    ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vn - ocean_state(jg)%p_check(cp2)%vn)**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  T     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%t  - ocean_state(jg)%p_check(cp2)%t )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  S     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%s  - ocean_state(jg)%p_check(cp2)%s )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  W     ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%w  - ocean_state(jg)%p_check(cp2)%w )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  G_NM1 ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%g_nm1-ocean_state(jg)%p_check(cp2)%g_nm1 )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '    HI  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%hi-ocean_state(jg)%p_check(cp2)%hi )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '    HS  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%hs-ocean_state(jg)%p_check(cp2)%hs )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  CONC  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%conc-ocean_state(jg)%p_check(cp2)%conc )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, 'ZUNDERIC', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%zUnderIce-ocean_state(jg)%p_check(cp2)%zUnderIce )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  ALB   ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%alb-ocean_state(jg)%p_check(cp2)%alb )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, 'SFL_REL ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%SaltFlux_Relax-ocean_state(jg)%p_check(cp2)%SaltFlux_Relax )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  TSURF ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Tsurf-ocean_state(jg)%p_check(cp2)%Tsurf )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  CSUM  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%concSum-ocean_state(jg)%p_check(cp2)%concSum )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  QTOP  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Qtop-ocean_state(jg)%p_check(cp2)%Qtop )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  QBOT  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%Qbot-ocean_state(jg)%p_check(cp2)%Qbot )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  VOL   ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vol-ocean_state(jg)%p_check(cp2)%vol )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  VOLS  ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%vols-ocean_state(jg)%p_check(cp2)%vols )**2))
    WRITE(0,sfmt) ' CPTEST ', strt,  cp1,cp2, '  DRAFT ', &
         &  SQRT(SUM((ocean_state(jg)%p_check(cp1)%draft-ocean_state(jg)%p_check(cp2)%draft )**2))
#endif /* __COMPAD__ */
!! COMPAD-DCO  END    JR2018 
    WRITE(0,*) ' -----------------------    CPTEST ', strt,  cp1,cp2

  END SUBROUTINE compare_checkpoints
#endif /*  __COMPAD_ADJLOOP__  */
