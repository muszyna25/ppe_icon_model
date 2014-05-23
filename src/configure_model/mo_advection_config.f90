!>
!! @brief configuration setup for tracer transport
!!
!! configuration setup for tracer transport
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-07-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_advection_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom,  &
    &                              MIURA, MIURA3, FFSL, FFSL_HYB, MCYCL,   &
    &                              MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL,  &
    &                              FFSL_HYB_MCYCL, ippm_vcfl, ippm_v,      &
    &                              ino_flx, izero_grad, iparent_flx, inwp

  IMPLICIT NONE
  PUBLIC


  CHARACTER(len=*), PARAMETER, PRIVATE :: &
    &  version = '$Id$'



  ! Derived type to allow for the onetime computation and cleanup 
  ! of tracer independent parts
  !
  TYPE t_compute                                                               
    LOGICAL :: ppm_v     (MAX_NTRACER)                                           
    LOGICAL :: miura_h   (MAX_NTRACER)                                           
    LOGICAL :: miura3_h  (MAX_NTRACER)
    LOGICAL :: ffsl_h    (MAX_NTRACER)
    LOGICAL :: ffsl_hyb_h(MAX_NTRACER)
    LOGICAL :: mcycl_h   (MAX_NTRACER)                            
  END TYPE t_compute                                                           


  TYPE t_scheme
    INTEGER :: iadv_min_slev     !< scheme dependent minimum vertical start level
                                 !< required for tracer-independent computations
    !INTEGER :: iadv_max_elev
    !LOGICAL :: lcompute (MAX_NTRACER)
    !LOGICAL :: lcleanup (MAX_NTRACER)
  END TYPE t_scheme


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for tracer advection
  !!--------------------------------------------------------------------------
  TYPE :: t_advection_config

    ! namelist variables

    CHARACTER(len=MAX_CHAR_LENGTH) :: &  !< list of tracers to initialize        
      &  ctracer_list                                                            

    INTEGER :: &                    !< selects horizontal transport scheme       
      &  ihadv_tracer(MAX_NTRACER)  !< 0:  no horizontal advection                
                                    !< 1:  1st order upwind                       
                                    !< 2:  2nd order miura                        
                                    !< 3:  3rd order miura with quadr./cubic reconstr.
                                    !< 4:  Flux form semi lagrange (FFSL)
                                    !< 5:  hybrid FFSL Miura3
                                    !< 6:  3rd or 4th order upstream (on hexagons only)
                                    !< 20: subcycling version of miura
                                    !< 22: 2nd order miura and miura_cycl
                                    !< 32: 3rd order miura with miura_cycl
                                    !< 42: FFSL with miura_cycl
                                    !< 52: FFSL_HYB with miura_cycl

    INTEGER :: &                    !< selects vertical transport scheme         
      &  ivadv_tracer(MAX_NTRACER)  !< 0 : no vertical advection                 
                                    !< 1 : 1st order upwind                      
                                    !< 3 : 3rd order PPM for CFL>                         
                                    !< 30: 3rd order PPM               

    INTEGER :: &                    !< advection of TKE
      &  iadv_tke                   !< 0 : none
                                    !< 1 : vertical advection only
                                    !< 2 : vertical and horizontal advection

    LOGICAL :: lvadv_tracer         !< if .TRUE., calculate vertical tracer advection
    LOGICAL :: lclip_tracer         !< if .TRUE., clip negative tracer values    
    LOGICAL :: lstrang              !< if .TRUE., use complete Strang splitting  
                                    !< (\Delta t/2 vert)+(\Delta t hor)+(\Delta t/2 vert)  
                                                   
    LOGICAL :: llsq_svd             !< least squares reconstruction with         
                                    !< singular value decomposition (TRUE) or    
                                    !< QR decomposition (FALSE) of design matrix A
    INTEGER :: &                    !< parameter used to select the limiter      
      &  itype_vlimit(MAX_NTRACER)  !< for vertical transport                    

    INTEGER :: &                    !< parameter used to select the limiter
      &  itype_hlimit(MAX_NTRACER)  !< for horizontal transport                  
                                    !< 0: no limiter                             
                                    !< 1: semi-monotonous slope limiter          
                                    !< 2: monotonous slope limiter               
                                    !< 3: monotonous flux limiter                

    INTEGER :: niter_fct            !< number of iterations for monotone
                                    !< flux correction procedure

    REAL(wp):: beta_fct             !< factor of allowed over-/undershooting in monotonous limiter

    INTEGER :: iord_backtraj        !< parameter to select the spacial order     
                                    !< of accuracy for the backward trajectory   

    INTEGER :: igrad_c_miura        !< parameter used to select the gradient     
                                    !< reconstruction method at cell center      
                                    !< for second order miura scheme             

    INTEGER :: ivcfl_max            !< determines stability range of vertical    
                                    !< ppm-scheme (approximate allowable maximum 
                                    !< CFL-number)                               

    REAL(wp) :: upstr_beta_adv      !< later, it should be combined with         
                                    !< upstr_beta in non-hydrostatic namelist


    ! derived variables

    REAL(wp) :: cSTR             !< if complete Strang-splitting is used,        
                                 !< this constant adapts the time step           
                                                                                 
    INTEGER  :: iubc_adv         !< selects upper boundary condition             
                                 !< for tracer transport                         
                                 !< 0: no flux                                   
                                 !< 1: zero gradient                             
                                 !< 2: interpolated flux from parent grid        
                                                                                 
    INTEGER ::  &                !< selects vertical start level for each patch  
      &  iadv_slev(MAX_NTRACER)  !< and each tracer.

    INTEGER ::  &                !< vertical end level down to which qv is 
      &  iadv_qvsubstep_elev     !< advected with internal substepping (to 
                                 !< circumvent CFL instability in the 
                                 !< stratopause region).
                                                                                 
    REAL(wp) :: coeff_grid       !< parameter which is used to make the vertical 
                                 !< advection scheme applicable to a height      
                                 !< based coordinate system (coeff_grid=-1)      

    ! scheme specific derived variables
    !
    TYPE(t_scheme) :: ppm_v      !< vertical PPM scheme
    TYPE(t_scheme) :: miura_h    !< horizontal miura scheme (linear reconstr.)
    TYPE(t_scheme) :: miura3_h   !< horizontal miura scheme (higher order reconstr.)
    TYPE(t_scheme) :: ffsl_h     !< horizontal FFSL scheme
    TYPE(t_scheme) :: ffsl_hyb_h !< horizontal hybrid FFSL/miura3 scheme
    TYPE(t_scheme) :: mcycl_h    !< horizontal miura scheme (linear) with subcycling
    TYPE(t_scheme) ::  &         !< combination of horizontal miura scheme and miura 
      &  miura_mcycl_h           !< scheme with subcycling
    TYPE(t_scheme) ::  &         !< combination of horizontal miura3 scheme and miura
      &  miura3_mcycl_h          !< scheme with subcycling
    TYPE(t_scheme) ::  &         !< combination of horizontal ffsl scheme and miura
      &  ffsl_mcycl_h            !< scheme with subcycling
    TYPE(t_scheme) ::  &         !< combination of horizontal hybrid ffsl scheme and miura
      &  ffsl_hyb_mcycl_h        !< scheme with subcycling 
  END TYPE t_advection_config

  !>
  !!
  TYPE(t_advection_config), TARGET :: advection_config(0:max_dom)

!DR For the time being lcompute and lcleanup are not added to the 
!DR advection_config state
  TYPE(t_compute)  :: lcompute
  TYPE(t_compute)  :: lcleanup


  ! for first order Gauss-Legendre quadrature
  !
  REAL(wp) :: shape_func_l(4)  !< shape functions for mapping the FFSL departure
                               !< region onto the standard rectangle
                                                                                 
  REAL(wp) :: zeta_l, eta_l    !< Gauss quadrature point in \zeta-\eta space                        
                                                                                 
  REAL(wp) :: wgt_zeta_l       !< Gauss quadrature weights for zeta and eta    
  REAL(wp) :: wgt_eta_l        !< points


  ! for second order Gauss-Legendre quadrature
  !
  REAL(wp) :: shape_func(4,4)  !< shape functions for mapping the FFSL departure
                               !< region onto the standard rectangle (miura3 only)
                                                                                 
  REAL(wp) :: zeta(4), eta(4)  !< Gauss quadrature points in \zeta-\eta space  
                               !< (miura3 only)                                
                                                                                 
  REAL(wp) :: wgt_zeta(4)      !< Gauss quadrature weights for zeta and eta    
  REAL(wp) :: wgt_eta(4)       !< points (miura3 only) 



CONTAINS

  !>
  !! setup components of the transport scheme depending on this namelist
  !!
  !! Setup of additional transport control variables depending on the 
  !! transport-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-04-20)
  !!
  SUBROUTINE configure_advection( jg, num_lev, num_lev_1, iequations, iforcing,        &
    &                            iqc, iqt,                                             &
    &                            kstart_moist, kend_qvsubstep, lvert_nest, l_open_ubc, &
    &                            ntracer )
  !
    INTEGER, INTENT(IN) :: jg           !< patch 
    INTEGER, INTENT(IN) :: num_lev      !< number of vertical levels
    INTEGER, INTENT(IN) :: num_lev_1    !< vertical levels of global patch
    INTEGER, INTENT(IN) :: iequations
    INTEGER, INTENT(IN) :: iforcing
    INTEGER, INTENT(IN) :: iqc, iqt     !< hydrometeor indices
    INTEGER, INTENT(IN) :: kstart_moist
    INTEGER, INTENT(IN) :: kend_qvsubstep
    INTEGER, INTENT(IN) :: ntracer
    LOGICAL, INTENT(IN) :: lvert_nest
    LOGICAL, INTENT(IN) :: l_open_ubc

    INTEGER :: jt          !< tracer loop index
    INTEGER :: jm          !< loop index for shape functions
    INTEGER :: ivadv_tracer(MAX_NTRACER)
    INTEGER :: ihadv_tracer(MAX_NTRACER)
    !-----------------------------------------------------------------------

    !
    ! set dependent transport variables/model components, depending on 
    ! the transport namelist and potentially other namelsists.
    !

    ! check, whether Strang-splitting has been chosen and adapt cSTR accordingly
    IF ( advection_config(jg)%lstrang ) THEN
      advection_config(jg)%cSTR = 0.5_wp
    ELSE
      advection_config(jg)%cSTR = 1._wp
    ENDIF


    ! Set grid-coefficient according to the applied vertical grid.
    !
    ! coeff_grid=1   : pressure based vertical coordinate system
    ! coeff_grid=-1  : height based vertical coordinate system
    !
    IF (iequations == 3) THEN  ! non-hydrostatic equation-set
      advection_config(jg)%coeff_grid = -1._wp
    ELSE
      advection_config(jg)%coeff_grid = 1._wp
    ENDIF


    !
    ! set vertical start level for each patch and each tracer
    !
    advection_config(jg)%iadv_slev(:) = 1
    advection_config(jg)%iadv_qvsubstep_elev = 1
    IF (iforcing == inwp) THEN
      ! Set iadv_slev to kstart_moist for all moisture fields but QV
      ! note: iqt denotes the first tracer index not related to moisture
      advection_config(jg)%iadv_slev(iqc:iqt-1) = kstart_moist
      advection_config(jg)%iadv_qvsubstep_elev = kend_qvsubstep
    ENDIF



    ! set boundary condition for vertical transport
    !
    IF (iequations == 3) THEN  ! non-hydrostatic equation-set

      IF (.NOT. lvert_nest ) THEN ! no vertical nesting

        IF (l_open_ubc) THEN
          advection_config(jg)%iubc_adv = izero_grad ! zero gradient ubc
        ELSE
          advection_config(jg)%iubc_adv = ino_flx    ! no flux ubc
        ENDIF

      ELSE ! vertical nesting

        IF (num_lev < num_lev_1) THEN
          advection_config(jg)%iubc_adv = iparent_flx
        ELSE IF ( (num_lev >= num_lev_1) .AND. l_open_ubc) THEN
          advection_config(jg)%iubc_adv = izero_grad
        ELSE IF ( (num_lev >= num_lev_1) .AND. .NOT. l_open_ubc) THEN
          advection_config(jg)%iubc_adv = ino_flx
        ENDIF
      ENDIF

    ELSE ! hydrostatic or shallow water equation set
      advection_config(jg)%iubc_adv = ino_flx    ! no flux ubc
    ENDIF




    ! to save some paperwork
    ivadv_tracer(:) = advection_config(1)%ivadv_tracer(:)
    ihadv_tracer(:) = advection_config(1)%ihadv_tracer(:)


    ! PPM_V[CFL] specific settings (vertical transport)
    !
    lcompute%ppm_v(:)   = .FALSE.
    lcleanup%ppm_v(:)   = .FALSE.

    advection_config(jg)%ppm_v%iadv_min_slev = HUGE(1)

    IF ( ANY(ivadv_tracer == ippm_v) .OR. ANY(ivadv_tracer == ippm_vcfl)  ) THEN

      ! compute minimum required slev for this group of tracers
      DO jt=1,ntracer
        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
          advection_config(jg)%ppm_v%iadv_min_slev =                           &
            &                  MIN( advection_config(jg)%ppm_v%iadv_min_slev,  &
            &                        advection_config(jg)%iadv_slev(jt) )
        ENDIF
      ENDDO

      ! Search for the first tracer jt for which vertical advection of
      ! type PPM has been selected.
      DO jt=1,ntracer
        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
          lcompute%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which vertical advection of
      ! type PPM has been selected.
      DO jt=ntracer,1,-1
        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
          lcleanup%ppm_v(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    !
    ! MIURA specific settings (horizontal transport)
    !
    lcompute%miura_h(:) = .FALSE.
    lcleanup%miura_h(:) = .FALSE.

    advection_config(jg)%miura_h%iadv_min_slev = HUGE(1)

    IF ( ANY(ihadv_tracer(1:ntracer) == MIURA) ) THEN

      ! compute minimum required slev for this group of tracers
      DO jt=1,ntracer
        IF ( ANY( (/MIURA, MIURA_MCYCL/) == ihadv_tracer(jt) ) ) THEN
          advection_config(jg)%miura_h%iadv_min_slev =  &
            &                  MIN( advection_config(jg)%miura_h%iadv_min_slev, &
            &                       advection_config(jg)%iadv_slev(jt) )
        ENDIF
      ENDDO

      ! Search for the first tracer jt for which horizontal advection of
      ! type MIURA has been selected.
      DO jt=1,ntracer
        IF ( ANY( (/MIURA, MIURA_MCYCL/) == ihadv_tracer(jt) ) ) THEN
          lcompute%miura_h(jt) = .TRUE.
          exit
        ENDIF
      ENDDO

      ! Search for the last tracer jt for which horizontal advection of
      ! type MIURA has been selected.
      DO jt=ntracer,1,-1
        IF ( ANY( (/MIURA, MIURA_MCYCL/) == ihadv_tracer(jt) ) ) THEN
          lcleanup%miura_h(jt) = .TRUE.
          exit
        ENDIF
      ENDDO
    END IF


    !
    ! MIURA3 specific settings (horizontal transport)
    !
    lcompute%miura3_h(:) = .FALSE.
    lcleanup%miura3_h(:) = .FALSE.

    advection_config(jg)%miura3_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%miura3_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%miura3_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type MIURA3 has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/MIURA3, MIURA3_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%miura3_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! FFSL specific settings (horizontal transport)
    !
    lcompute%ffsl_h(:) = .FALSE.
    lcleanup%ffsl_h(:) = .FALSE.

    advection_config(jg)%ffsl_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%ffsl_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%ffsl_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type FFSL has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%ffsl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type FFSL has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/FFSL, FFSL_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%ffsl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO



    !
    ! FFSL_HYB specific settings (horizontal transport)
    !
    lcompute%ffsl_hyb_h(:) = .FALSE.
    lcleanup%ffsl_hyb_h(:) = .FALSE.

    advection_config(jg)%ffsl_hyb_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%ffsl_hyb_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%ffsl_hyb_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type FFSL_HYB has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%ffsl_hyb_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type FFSL_HYB has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/FFSL_HYB, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%ffsl_hyb_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! MCYCL specific settings (horizontal transport)
    !
    lcompute%mcycl_h(:) = .FALSE.
    lcleanup%mcycl_h(:) = .FALSE.

    advection_config(jg)%mcycl_h%iadv_min_slev = HUGE(1)

    ! compute minimum required slev for this group of tracers
    DO jt=1,ntracer
      IF ( ANY( (/MCYCL, MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        advection_config(jg)%mcycl_h%iadv_min_slev =   &
          &                  MIN( advection_config(jg)%mcycl_h%iadv_min_slev, &
          &                       advection_config(jg)%iadv_slev(jt) )
      ENDIF
    ENDDO

    ! Search for the first tracer jt for which horizontal advection of
    ! type MCYCL has been selected.
    DO jt=1,ntracer
      IF ( ANY( (/MCYCL, MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcompute%mcycl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO

    ! Search for the last tracer jt for which horizontal advection of
    ! type MCYCL has been selected.
    DO jt=ntracer,1,-1
      IF ( ANY( (/MCYCL, MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL, FFSL_HYB_MCYCL/) == ihadv_tracer(jt) ) ) THEN
        lcleanup%mcycl_h(jt) = .TRUE.
        exit
      ENDIF
    ENDDO


    !
    ! Compute shape functions for mapping the departure region onto the
    ! standard rectangle. Integration points, shape functions and quadrature 
    ! weights are provided for a first and second order Gauss-Legendre
    ! quadrature.
    !

    IF (jg == 1) THEN

      !
      ! First order
      !

      ! Coordinates of integration points (in \zeta,\eta-System)
      !
      zeta_l = 0._wp
      eta_l  = 0._wp

      ! shape function for mapping
      shape_func_l(1) = 0.25_wp * (1._wp-zeta_l)*(1._wp-eta_l)
      shape_func_l(2) = 0.25_wp * (1._wp+zeta_l)*(1._wp-eta_l)
      shape_func_l(3) = 0.25_wp * (1._wp+zeta_l)*(1._wp+eta_l)
      shape_func_l(4) = 0.25_wp * (1._wp-zeta_l)*(1._wp+eta_l)


      ! Gauss quadrature weights
      !
      wgt_zeta_l = 2._wp
      wgt_eta_l  = 2._wp


      !
      ! Second order
      !

      ! Coordinates of integration points (in \zeta,\eta-System)
      !
      zeta(1) = -1._wp/SQRT(3._wp)
      zeta(2) =  1._wp/SQRT(3._wp)
      zeta(3) =  1._wp/SQRT(3._wp)
      zeta(4) = -1._wp/SQRT(3._wp)

      eta(1)  = -1._wp/SQRT(3._wp)
      eta(2)  = -1._wp/SQRT(3._wp)
      eta(3)  =  1._wp/SQRT(3._wp)
      eta(4)  =  1._wp/SQRT(3._wp)

      ! shape functions for mapping
      !
      DO jm = 1,4
        shape_func(1,jm) = 0.25_wp * (1._wp-zeta(jm))*(1._wp-eta(jm))
        shape_func(2,jm) = 0.25_wp * (1._wp+zeta(jm))*(1._wp-eta(jm))
        shape_func(3,jm) = 0.25_wp * (1._wp+zeta(jm))*(1._wp+eta(jm))
        shape_func(4,jm) = 0.25_wp * (1._wp-zeta(jm))*(1._wp+eta(jm))
      END DO

      ! Gauss quadrature weights
      !
      wgt_zeta(1) = 1._wp
      wgt_zeta(2) = 1._wp
      wgt_zeta(3) = 1._wp
      wgt_zeta(4) = 1._wp

      wgt_eta(1)  = 1._wp
      wgt_eta(2)  = 1._wp
      wgt_eta(3)  = 1._wp
      wgt_eta(4)  = 1._wp
    END IF

  END SUBROUTINE configure_advection


END MODULE mo_advection_config
