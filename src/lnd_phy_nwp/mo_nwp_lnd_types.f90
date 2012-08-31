#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!!  !MODULE:  mo_nwp_lnd_types\\
!!
!! Description:  Data type definition for land surface scheme (TERRA)
!
!!
!! @author Daniel Reinert, DWD
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Initial by Kristina Froehlich (2010-11-09)
!! Modification by Daniel Reinert, DWD (2012-04-03)
!! - encapsulated type definitions (mo_nwp_lnd_types)
!!
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
MODULE mo_nwp_lnd_types

  USE mo_kind,                 ONLY: wp
  USE mo_linked_list,          ONLY: t_var_list


  IMPLICIT NONE
  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !
  !variables
  PUBLIC :: t_lnd_state  !> state vector for land scheme
  PUBLIC :: t_lnd_prog   !!  for prognostic variables
  PUBLIC :: t_wtr_prog   !!  for prognostic variables related to lake and sea ice models
  PUBLIC :: t_lnd_diag   !!  for diagnostic variables
  PUBLIC :: t_ptr_lnd

  TYPE t_ptr_lnd
    REAL(wp), POINTER :: p_3d(:,:,:) ! pointer to 3D (spatial) array
    REAL(wp), POINTER :: p_2d(:,:)   ! pointer to 2D (spatial) array
  END TYPE t_ptr_lnd

  !
  ! prognostic variables state vector (land)
  !
  TYPE t_lnd_prog

    REAL(wp), POINTER :: &
    &  t_snow_t       (:,:,:)   , & ! temperature of the snow-surface               (  K  )
    &  t_snow_mult_t  (:,:,:,:) , & ! temperature of snow                           (  K  )
    &  t_s_t          (:,:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_g            (:,:)     , & ! weighted surface temperature                  (  K  )
    &  t_g_t          (:,:,:)   , & ! weighted surface temperature on tiles         (  K  )
    &  w_snow_t       (:,:,:)   , & ! water content of snow                         (m H2O)
    &  rho_snow_t     (:,:,:)   , & ! snow density                                  (kg/m**3)
    &  rho_snow_mult_t(:,:,:,:) , & ! snow density                                  (kg/m**3)
    &  w_i_t          (:,:,:)   , & ! water content of interception water           (m H2O)
    &  t_so_t         (:,:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so_t         (:,:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice_t     (:,:,:,:) , & ! ice content                                   (m H20)
    &  wliq_snow_t    (:,:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  wtot_snow_t    (:,:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  dzh_snow_t     (:,:,:,:)     ! layer thickness between half levels in snow   (  m  )

    TYPE(t_ptr_lnd), ALLOCATABLE :: t_snow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: t_snow_mult_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: t_s_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: t_gt_ptr(:) 
    TYPE(t_ptr_lnd), ALLOCATABLE :: w_snow_ptr(:) 
    TYPE(t_ptr_lnd), ALLOCATABLE :: rho_snow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: rho_snow_mult_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: w_i_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: t_so_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: w_so_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: w_so_ice_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: wliq_snow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: wtot_snow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: dzh_snow_ptr(:)

  END TYPE t_lnd_prog


  !
  ! prognostic variables state vector (sea)
  !
  TYPE t_wtr_prog

    REAL(wp), POINTER :: &
    &  t_ice        (:,:)   , & ! temperature of the sea ice          (  K  )
    &  h_ice        (:,:)   , & ! depth of the sea ice                (  m  )
    &  t_snow_si    (:,:)   , & ! temperature of the snow on sea ice  (  K  )
    &  h_snow_si    (:,:)       ! depth of the snow on sea ice        (  m  )

  END TYPE t_wtr_prog

  !
  ! diagnostic variables state vector
  !
  TYPE t_lnd_diag

    REAL(wp), POINTER ::   &
    &  qv_s         (:,:)   , & ! specific humidity at the surface              (kg/kg)
    &  t_snow       (:,:)   , & ! temperature of the snow-surface               (  K  )
    &  t_snow_mult  (:,:,:) , & ! temperature of snow                           (  K  )
    &  t_s          (:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_seasfc     (:,:)   , & ! temperature of the sea surface                (  K  )
    &  w_snow       (:,:)   , & ! snow water equivalent                         (m H2O)
    &  w_snow_eff_t (:,:,:) , & ! snow water equivalent / snow-cover fraction   (m H2O)
    &  rho_snow     (:,:)   , & ! snow density                                  (kg/m**3)
    &  rho_snow_mult(:,:,:) , & ! snow density                                  (kg/m**3)
    &  w_i          (:,:)   , & ! water content of interception water           (m H2O)
    &  t_so         (:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so         (:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice     (:,:,:) , & ! ice content                                   (m H20)
    &  wliq_snow    (:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  wtot_snow    (:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  dzh_snow     (:,:,:) , & ! layer thickness between half levels in snow   (  m  )
    &  h_snow       (:,:)   , & ! snow height                                   (  m  ) 
    &  freshsnow    (:,:)   , & ! indicator for age of snow in top of snow layer(  -  )
    &  snowfrac     (:,:)   , & ! snow-cover fraction                           (  -  )
    &  runoff_s     (:,:)   , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g     (:,:)   , & ! soil water runoff; sum over forecast          (kg/m2)
    &  fr_seaice    (:,:)   , & !< fraction of sea ice                          ( )   
                                !< as partition of total area of the
                                !< grid element, but set to 0 or 1
                                !< index1=1,nproma, index2=1,nblks_c
    &  qv_s_t       (:,:,:) , & ! specific humidity at the surface              (kg/kg)
    &  h_snow_t     (:,:,:) , & ! snow height                                   (  m  ) 
    &  freshsnow_t  (:,:,:) , & ! indicator for age of snow in top of snow layer(  -  )
    &  snowfrac_lc_t(:,:,:) , & ! snow-cover fraction per land-cover class      (  -  )
    &  snowfrac_t   (:,:,:) , & ! snow-cover fraction                           (  -  )
    &  runoff_s_t   (:,:,:) , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g_t   (:,:,:) , & ! soil water runoff; sum over forecast          (kg/m2)
    &  subsfrac_t   (:,:,:)

    TYPE(t_ptr_lnd), ALLOCATABLE :: qv_st_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: h_snow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: freshsnow_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: snowfrac_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: snowfrac_lc_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: w_snow_eff_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: runoff_s_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: runoff_g_ptr(:)
    TYPE(t_ptr_lnd), ALLOCATABLE :: subsfrac_ptr(:)

  END TYPE t_lnd_diag


! complete state vector type
!
  TYPE t_lnd_state
    TYPE(t_lnd_prog), ALLOCATABLE  :: prog_lnd(:)          ! number of time levels
    TYPE(t_wtr_prog), ALLOCATABLE  :: prog_wtr(:)          ! number of time levels
    TYPE(t_var_list), ALLOCATABLE  :: lnd_prog_nwp_list(:) ! number of time levels
    TYPE(t_var_list), ALLOCATABLE  :: wtr_prog_nwp_list(:) ! number of time levels
    TYPE(t_lnd_diag)               :: diag_lnd
    TYPE(t_var_list)               :: lnd_diag_nwp_list
  END TYPE t_lnd_state
 

END MODULE mo_nwp_lnd_types
