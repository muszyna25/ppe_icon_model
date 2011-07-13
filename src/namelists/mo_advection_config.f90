!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
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

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  ! Derived types to allow for the onetime computation of tracer independent parts

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



  PUBLIC

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
                                                                                 
    REAL(wp) :: shape_func(4,4)  !< shape functions for mapping the FFSL departure
                                 !< region onto the standard rectangle (miura3 only)
                                                                                 
    REAL(wp) :: zeta(4), eta(4)  !< Gauss quadrature points in \zeta-\eta space  
                                 !< (miura3 only)                                
                                                                                 
    REAL(wp) :: wgt_zeta(4)      !< Gauss quadrature weights for zeta and eta    
    REAL(wp) :: wgt_eta(4)       !< points (miura3 only) 

!!$    TYPE(t_compute) :: lcompute                                                  
!!$    TYPE(t_cleanup) :: lcleanup

    REAL(wp), POINTER :: ptr_delp_mc_now(:,:,:)  !< pointer to old layer thickness      
                                                 !< at cell center                      
    REAL(wp), POINTER :: ptr_delp_mc_new(:,:,:)  !< pointer to new layer thickness      
                                                 !< at cell center
  END TYPE t_advection_config
  !>
  !!
 TYPE(t_advection_config), TARGET :: advection_config(max_dom)


!SUBROUTINE config_advection !(iequations)
!END SUBROUTINE config_advection


END MODULE mo_advection_config
