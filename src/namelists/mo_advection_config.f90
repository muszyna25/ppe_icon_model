!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_advection_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom,  &
    &                              imiura, imiura3, ippm_vcfl, ippm_v,     &
    &                              ino_flx, izero_grad, iparent_flx, inwp

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'



!!$  ! Derived types to allow for the onetime computation of tracer independent parts
!!$  !
!!$  TYPE t_compute                                                               
!!$    LOGICAL :: muscl_v (MAX_NTRACER)                                           
!!$    LOGICAL :: ppm_v   (MAX_NTRACER)                                           
!!$    LOGICAL :: miura_h (MAX_NTRACER)                                           
!!$    LOGICAL :: miura3_h(MAX_NTRACER)                                           
!!$  END TYPE t_compute                                                           
!!$                                                                               
!!$  TYPE t_cleanup                                                              
!!$    LOGICAL :: muscl_v (MAX_NTRACER)                                           
!!$    LOGICAL :: ppm_v   (MAX_NTRACER)                                           
!!$    LOGICAL :: miura_h (MAX_NTRACER)                                           
!!$    LOGICAL :: miura3_h(MAX_NTRACER)                                           
!!$  END TYPE t_cleanup



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for tracer advection
  !!--------------------------------------------------------------------------
  TYPE :: t_advection_config

    ! namelist variables

    CHARACTER(len=MAX_CHAR_LENGTH) :: &  !< list of tracers to initialize        
      &  ctracer_list                                                            

    INTEGER :: &                    !< selects horizontal transport scheme       
      &  ihadv_tracer(MAX_NTRACER)  !< 0: no horizontal advection                
                                    !< 1: 1st order upwind                       
                                    !< 2: 2nd order muscl                        
                                    !< 3: 2nd order miura                        
                                    !< 4: 3rd order miura with quadr. reconstr.  

    INTEGER :: &                    !< selects vertical transport scheme         
      &  ivadv_tracer(MAX_NTRACER)  !< 0 : no vertical advection                 
                                    !< 1 : 1st order upwind                      
                                    !< 2 : 2nd order muscl                       
                                    !< 20: 2nd order muscl for CFL>1             
                                    !< 3 : 3rd order PPM                         
                                    !< 30: 3rd order PPM for CFL>1               

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
                                                                                 
    REAL(wp) :: coeff_grid       !< parameter which is used to make the vertical 
                                 !< advection scheme applicable to a height      
                                 !< based coordinate system (coeff_grid=-1)      
                   
!!$    TYPE(t_compute) :: lcompute
!!$    TYPE(t_cleanup) :: lcleanup
                                           
  END TYPE t_advection_config

  !>
  !!
  TYPE(t_advection_config), TARGET :: advection_config(max_dom)



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
  SUBROUTINE configure_advection( jg, num_lev, num_lev_1, iequations,    &
    &                            iforcing, iqv, kstart_moist, kstart_qv, &
    &                            lvert_nest, l_open_ubc, ntracer )
  !
    INTEGER, INTENT(IN) :: jg           !< patch 
    INTEGER, INTENT(IN) :: num_lev      !< number of vertical levels
    INTEGER, INTENT(IN) :: num_lev_1    !< vertical levels of global patch
    INTEGER, INTENT(IN) :: iequations
    INTEGER, INTENT(IN) :: iforcing
    INTEGER, INTENT(IN) :: iqv
    INTEGER, INTENT(IN) :: kstart_moist
    INTEGER, INTENT(IN) :: kstart_qv
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
    IF (iforcing == inwp) THEN
      ! Set iadv_slev to kstart_moist for all tracers but QV
      advection_config(jg)%iadv_slev(:)   = kstart_moist
      advection_config(jg)%iadv_slev(iqv) = kstart_qv
    ELSE
      advection_config(jg)%iadv_slev(:) = 1
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




!!$    ! to save some paperwork
!!$    ivadv_tracer(:) = advection_config(jg)%ivadv_tracer(:)
!!$    ihadv_tracer(:) = advection_config(jg)%ihadv_tracer(:)
!!$
!!$    ! MUSCL_V[CFL] specific settings (vertical transport)
!!$    !
!!$    advection_config(jg)%lcompute%muscl_v(:) = .FALSE.
!!$    advection_config(jg)%lcleanup%muscl_v(:) = .FALSE.
!!$
!!$    IF ( ANY(ivadv_tracer == imuscl_v) .OR. ANY(ivadv_tracer == imuscl_vcfl)  ) THEN
!!$      ! Search for the first tracer jt for which vertical advection of
!!$      ! type MUSCL has been selected.
!!$      DO jt=1,ntracer
!!$        IF ( ivadv_tracer(jt) == imuscl_v .OR. ivadv_tracer(jt) == imuscl_vcfl ) THEN
!!$          advection_config(jg)%lcompute%muscl_v(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$
!!$      ! Search for the last tracer jt for which vertical advection of
!!$      ! type MUSCL has been selected.
!!$      DO jt=ntracer,1,-1
!!$        IF ( ivadv_tracer(jt) == imuscl_v .OR. ivadv_tracer(jt) == imuscl_vcfl ) THEN
!!$          advection_config(jg)%lcleanup%muscl_v(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$    END IF


!!$    ! PPM_V[CFL] specific settings (vertical transport)
!!$    !
!!$    advection_config(jg)%lcompute%ppm_v(:)   = .FALSE.
!!$    advection_config(jg)%lcleanup%ppm_v(:)   = .FALSE.
!!$
!!$    IF ( ANY(ivadv_tracer == ippm_v) .OR. ANY(ivadv_tracer == ippm_vcfl)  ) THEN
!!$      ! Search for the first tracer jt for which vertical advection of
!!$      ! type PPM has been selected.
!!$      DO jt=1,ntracer
!!$        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
!!$          advection_config(jg)%lcompute%ppm_v(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$
!!$      ! Search for the last tracer jt for which vertical advection of
!!$      ! type PPM has been selected.
!!$      DO jt=ntracer,1,-1
!!$        IF ( ivadv_tracer(jt) == ippm_v .OR. ivadv_tracer(jt) == ippm_vcfl ) THEN
!!$          advection_config(jg)%lcleanup%ppm_v(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$    END IF


!!$    !
!!$    ! MIURA specific settings
!!$    !
!!$    advection_config(jg)%lcompute%miura_h(:) = .FALSE.
!!$    advection_config(jg)%lcleanup%miura_h(:) = .FALSE.
!!$
!!$    IF ( ANY(ihadv_tracer(1:ntracer) == imiura) ) THEN
!!$      ! Search for the first tracer jt for which horizontal advection of
!!$      ! type MIURA has been selected.
!!$      DO jt=1,ntracer
!!$        IF ( ihadv_tracer(jt) == imiura ) THEN
!!$          advection_config(jg)%lcompute%miura_h(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$
!!$      ! Search for the last tracer jt for which horizontal advection of
!!$      ! type MIURA has been selected.
!!$      DO jt=ntracer,1,-1
!!$        IF ( ihadv_tracer(jt) == imiura ) THEN
!!$          advection_config(jg)%lcleanup%miura_h(jt) = .TRUE.
!!$          exit
!!$        ENDIF
!!$      ENDDO
!!$    END IF


!!$    !
!!$    ! MIURA3 specific settings
!!$    !
!!$    advection_config(jg)%lcompute%miura3_h(:) = .FALSE.
!!$    advection_config(jg)%lcleanup%miura3_h(:) = .FALSE.
!!$
!!$    ! Search for the first tracer jt for which horizontal advection of
!!$    ! type MIURA3 has been selected.
!!$    DO jt=1,ntracer
!!$      IF ( ihadv_tracer(jt) == imiura3 ) THEN
!!$        advection_config(jg)%lcompute%miura3_h(jt) = .TRUE.
!!$        exit
!!$      ENDIF
!!$    ENDDO
!!$
!!$    ! Search for the last tracer jt for which horizontal advection of
!!$    ! type MIURA3 has been selected.
!!$    DO jt=ntracer,1,-1
!!$      IF ( ihadv_tracer(jt) == imiura3 ) THEN
!!$        advection_config(jg)%lcleanup%miura3_h(jt) = .TRUE.
!!$        exit
!!$      ENDIF
!!$    ENDDO


    !
    ! Compute shape functions for mapping the departure region onto the
    ! standard rectangle. It is assumed that a second order Gauss-Legendre
    ! quadrature is applied
    !

    IF (jg == 1) THEN

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
