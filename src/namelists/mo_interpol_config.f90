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
MODULE mo_interpol_config

  USE mo_kind,                ONLY: wp
  USE mo_intp_data_strc,      ONLY: t_lsq_set
  USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE

  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !!
  TYPE :: t_interpol_config

    ! namelist variables
    LOGICAL  :: llsq_high_consv     ! flag to determine whether the high order least 
                                    ! squares reconstruction should be conservative
                                                                                 
    INTEGER  :: lsq_high_ord        ! specific order for higher order lsq        
                                                                                 
    INTEGER  :: rbf_vec_kern_c,   & ! parameter determining the type             
       &        rbf_vec_kern_v,   & ! of vector rbf kernel                       
       &        rbf_vec_kern_e                                                   
                                                                                 
    ! Parameter fields determining the scale factor used by the vector rbf       
    ! interpolator.                                                              
    ! Note: these fields are defined on each grid level; to allow the namelist input
    ! going from 1 to depth (rather than from start_lev to end_lev), the namelist input         
    ! fields defined here differ from those used in the model                    
  
    REAL(wp) :: rbf_vec_scale_c
    REAL(wp) :: rbf_vec_scale_e
    REAL(wp) :: rbf_vec_scale_v
                                                                                 
    INTEGER  :: i_cori_method       ! Identifier for the method with wich the tangential        
                                    ! wind reconstruction in Coriolis force is computed,        
                                    ! if the Thuburn method is used. (To be      
                                    ! implemented for triangles, currently only for
                                    ! hexagons)                                  
                                    ! i_cori_method = 1 : Almut's method for reconstruction    
                                    !                     but TRSK method for PV 
                                    ! i_cori_method = 2 : Thuburn/Ringler/Skamarock/Klemp      
                                    ! i_cori_method = 3 : Almut's method for reconstruction    
                                    !                     Almut's method also for PV
                                    ! i_cori_method = 4 : Almut's method for reconstruction,   
                                    !                     but PV on averaged on vertices       

    ! Namelist variables setting up the lateral boundary nudging (applicable to limited-area   
    ! runs and one-way nesting). The nudging coefficients start with nudge_max_coeff in        
    ! the cell row bordering to the boundary interpolation zone, and decay exponentially       
    ! with nudge_efold_width (in units of cell rows)
  
    REAL(wp) :: nudge_max_coeff, nudge_efold_width
    INTEGER  :: nudge_zone_width    ! total width of nudging zone in units of cell rows
                                                                                 
    LOGICAL :: l_corner_vort        ! yields for i_cori_method>=3
                                    ! Decision wheter the hexagon vector reconstruction is
                                    ! combined with either of the two vorticities :
                                    ! .TRUE. : Three rhombi are combined to the corner
                                    !          and afterwards averaged to the hexagon center
                                    ! .FALSE.: 6 rhombi are directly averaged to the
                                    !          hexagon center (original method).
                                    ! After the writing of the paper to be published in JCP
                                    ! it seems that l_corner_vort=.TRUE. should be the right way.
  
    INTEGER  ::  rbf_vec_dim_c,    & ! parameter determining the size            
       &         rbf_vec_dim_v,    & ! of vector rbf stencil                     
       &         rbf_vec_dim_e,    & !                                           
       &         rbf_c2grad_dim      ! ... and for cell-to-gradient reconstruction
                                                                                 
                                                                                 
    TYPE(t_lsq_set) :: lsq_lin_set, &! Parameter settings for linear and higher order  
      &                lsq_high_set  ! least squares
  
  END TYPE t_interpol_config
  !>
  !!
  TYPE(t_interpol_config) :: interpol_config(max_dom)

END MODULE mo_interpol_config
