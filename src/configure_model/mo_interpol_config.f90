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
  USE mo_intp_data_strc,      ONLY: t_lsq_set, sick_a, sick_o
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: llsq_high_consv, lsq_high_ord                                !< variables
  PUBLIC :: rbf_vec_kern_c, rbf_vec_kern_v, rbf_vec_kern_e               !< variables
  PUBLIC :: rbf_vec_scale_c, rbf_vec_scale_e, rbf_vec_scale_v            !< variables
  PUBLIC :: i_cori_method, l_corner_vort                                 !< variables
  PUBLIC :: nudge_max_coeff, nudge_efold_width, nudge_zone_width         !< variables
  PUBLIC :: rbf_vec_dim_c, rbf_vec_dim_v, rbf_vec_dim_e, rbf_c2grad_dim  !< variables
  PUBLIC :: lsq_lin_set, lsq_high_set                                    !< variables
  PUBLIC :: configure_interpolation                                      !< subroutine

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'
  !>
  !!
  !TYPE :: t_interpol_config

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
  
    REAL(wp) :: rbf_vec_scale_c(max_dom)
    REAL(wp) :: rbf_vec_scale_e(max_dom)
    REAL(wp) :: rbf_vec_scale_v(max_dom)
                                                                                 
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
  
  !END TYPE t_interpol_config
  !>
  !!
  !TYPE(t_interpol_config) :: interpol_config(max_dom)

CONTAINS
  !>
  !!
  SUBROUTINE configure_interpolation( global_cell_type, n_dom, grid_level )

    INTEGER,INTENT(IN) :: global_cell_type
    INTEGER,INTENT(IN) :: n_dom
    INTEGER,INTENT(IN) :: grid_level(n_dom)
    
    INTEGER :: jg, jlev
    CHARACTER(len=*),PARAMETER :: routine = 'mo_interpol_config:configure_interpol'

    !-----------------------------------------------------------------------
    ! Set stencil size for RBF vector reconstruction
    !-----------------------------------------------------------------------
    rbf_vec_dim_c   = 9    ! use 2nd order reconstruction for cell centers
    rbf_vec_dim_v   = 6    ! use 6-point reconstruction for vertices
    rbf_vec_dim_e   = 4    ! use 4-point reconstruction for edge midpoints
    rbf_c2grad_dim  = 10   ! for reconstructing gradient at cell midpoints

    !-----------------------------------------------------------------------
    ! If RBF scaling factors are not supplied by the namelist, they are now
    ! initialized with meaningful values depending on the grid level and the
    ! stencil size
    ! Please note: RBF scaling factors for p_patch(0) (if it exists)
    ! are not set here - they are taken from p_patch(1) in the setup routines
    !-----------------------------------------------------------------------
    ! rbf_vec_scale_c 
    !-----------------
    ! - values are specified for Gaussian kernel
    ! (need to be smaller for inv. multiquadric)

    DO jg = 1,n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_c(jg) > 0.0_wp) CYCLE

      jlev = grid_level(jg)
      IF      (jlev <= 9 ) THEN ; rbf_vec_scale_c(jg) = 0.5_wp
      ELSE IF (jlev == 10) THEN ; rbf_vec_scale_c(jg) = 0.45_wp
      ELSE IF (jlev == 11) THEN ; rbf_vec_scale_c(jg) = 0.3_wp
      ELSE IF (jlev == 12) THEN ; rbf_vec_scale_c(jg) = 0.1_wp
      ELSE IF (jlev == 13) THEN ; rbf_vec_scale_c(jg) = 0.03_wp
      ELSE                      ; rbf_vec_scale_c(jg) = 0.01_wp
      ENDIF
    ENDDO

    !-----------------
    ! rbf_vec_scale_v
    !-----------------
    DO jg = 1, n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_v(jg) > 0.0_wp) CYCLE

      jlev = grid_level(jg)
      IF      (jlev <= 10) THEN ; rbf_vec_scale_v(jg) = 0.5_wp
      ELSE IF (jlev == 11) THEN ; rbf_vec_scale_v(jg) = 0.4_wp
      ELSE IF (jlev == 12) THEN ; rbf_vec_scale_v(jg) = 0.25_wp
      ELSE IF (jlev == 13) THEN ; rbf_vec_scale_v(jg) = 0.07_wp
      ELSE                      ; rbf_vec_scale_v(jg) = 0.02_wp
      ENDIF
    ENDDO

    !-----------------
    ! rbf_vec_scale_e 
    !-----------------
    ! - values are specified for inverse multiquadric kernel

    DO jg = 1, n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_e(jg) > 0.0_wp) CYCLE

      jlev = grid_level(jg)
      IF      (jlev <= 10) THEN ; rbf_vec_scale_e(jg) = 0.5_wp
      ELSE IF (jlev == 11) THEN ; rbf_vec_scale_e(jg) = 0.45_wp
      ELSE IF (jlev == 12) THEN ; rbf_vec_scale_e(jg) = 0.37_wp
      ELSE IF (jlev == 13) THEN ; rbf_vec_scale_e(jg) = 0.25_wp
      ELSE                      ; rbf_vec_scale_e(jg) = 0.1_wp
      ENDIF
    ENDDO

    !-----------------------------------------------------------------------
    ! Now check the RBF scaling factors
    !-----------------------------------------------------------------------
    DO jg = 1, n_dom
      IF (rbf_vec_scale_c(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_c')
      ELSE IF ((rbf_vec_scale_c(jg) < 0.01_wp).OR.(rbf_vec_scale_c(jg) > 0.6_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_c')
        CALL message('','! is 0.01 <= rbf_vec_scale_c <= 0.6')
      ENDIF
    ENDDO

    DO jg = 1, n_dom
      IF (rbf_vec_scale_v(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_v')
      ELSE IF ((rbf_vec_scale_v(jg) < 0.02_wp).OR.(rbf_vec_scale_v(jg) > 1.0_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_v')
        CALL message('','! is 0.02 <= rbf_vec_scale_v <= 1.0')
      ENDIF
    ENDDO
   
    DO jg = 1, n_dom
      IF (rbf_vec_scale_e(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_e')
      ELSE IF ((rbf_vec_scale_e(jg) < 0.05_wp).OR.(rbf_vec_scale_e(jg) > 0.6_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_e')
        CALL message('','! is 0.05 <= rbf_vec_scale_e <= 0.6')
      ENDIF
    ENDDO

    !-----------------------------------------------------------------------
    ! set the number of unknowns in the least squares reconstruction, and the
    ! stencil size, depending on the chosen polynomial order (lsq_high_ord).
    ! Default value for lsq_high_ord is 3. The user might have made a different
    ! choise. The possibilities are:
    !  lsq_high_ord=2 : quadratic polynomial        : 5 unknowns with a 9-point stencil
    !  lsq_high_ord=30: poor man's cubic polynomial : 7 unknowns with a 9-point stencil
    !  lsq_high_ord=3 : full cubic polynomial       : 9 unknowns with a 9-point stencil
    !-----------------------------------------------------------------------
    IF (global_cell_type == 3) THEN

      ! Settings for linear lsq reconstruction

      lsq_lin_set%l_consv = .FALSE.
      lsq_lin_set%dim_c   = 3
      lsq_lin_set%dim_unk = 2
      lsq_lin_set%wgt_exp = 2

      ! Settings for high order lsq reconstruction

      lsq_high_set%l_consv = llsq_high_consv
      
      IF (lsq_high_ord == 2) THEN
        lsq_high_set%dim_c   = 9
        lsq_high_set%dim_unk = 5
        lsq_high_set%wgt_exp = 3
      ELSE IF (lsq_high_ord == 30) THEN
        lsq_high_set%dim_c   = 9
        lsq_high_set%dim_unk = 7
        lsq_high_set%wgt_exp = 2
      ELSE IF (lsq_high_ord == 3) THEN
        lsq_high_set%dim_c   = 9
        lsq_high_set%dim_unk = 9
        lsq_high_set%wgt_exp = 0
      ELSE
        CALL finish( TRIM(routine),'wrong value of lsq_high_ord, must be 2,30 or 3')
      ENDIF

      ! triangular grid: just avoid that thickness is averaged at edges as it is
      ! needed in the hexagonal grid
      i_cori_method=1
      
    ENDIF

    ! In case of a hexagonal model, we perform a quadratic reconstruction, and check
    ! for i_cori_method
    IF (global_cell_type == 6) THEN

      ! ... quadratic reconstruction
      lsq_high_set%dim_c   = 6
      lsq_high_set%dim_unk = 5
      lsq_high_set%wgt_exp = 2

      ! ... linear reconstruction is not used in the hexagonal model.
      ! However, if we do not initialize these variables, they will get 
      ! value zero, and cause problem in the dump/restore functionality.
      lsq_lin_set%dim_c   = 3
      lsq_lin_set%dim_unk = 2
      
      sick_a=0.0_wp
      SELECT CASE (i_cori_method)
      CASE (2) 
        sick_a=0.375_wp
      CASE (3)
        sick_a=0.75_wp
      END SELECT
      sick_o = 1.0_wp-sick_a
      
    ENDIF

  END SUBROUTINE configure_interpolation

END MODULE mo_interpol_config
