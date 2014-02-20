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
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
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


  !
  ! prognostic variables state vector (land)
  !
  TYPE t_lnd_prog

    REAL(wp), POINTER             &
#ifdef _CRAYFTN
    , CONTIGUOUS                  &
#endif
    &  ::                         &
    &  t_s_t          (:,:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_g            (:,:)     , & ! weighted surface temperature                  (  K  )
    &  t_g_t          (:,:,:)   , & ! weighted surface temperature on tiles         (  K  )
    &  w_i_t          (:,:,:)   , & ! water content of interception water           (m H2O)
    &  w_p_t          (:,:,:)   , & ! water content of pond interception water      (m H2O)
    &  w_s_t          (:,:,:)   , & ! water content of interception snow            (m H2O)
    &  t_so_t         (:,:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so_t         (:,:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice_t     (:,:,:,:) , & ! ice content                                   (m H20)
    &  t_snow_t       (:,:,:)   , & ! temperature of the snow-surface               (  K  )
    &  w_snow_t       (:,:,:)   , & ! water content of snow                         (m H2O)
    &  rho_snow_t     (:,:,:)   , & ! snow density                                  (kg/m**3)
    &  t_snow_mult_t  (:,:,:,:) , & ! temperature of snow                           (  K  )
    &  wtot_snow_t    (:,:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  wliq_snow_t    (:,:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  rho_snow_mult_t(:,:,:,:) , & ! snow density                                  (kg/m**3)
    &  dzh_snow_t     (:,:,:,:)     ! layer thickness between half levels in snow   (  m  )

    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_snow_mult_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_gt_ptr(:) 
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_snow_ptr(:) 
    TYPE(t_ptr_2d3d), ALLOCATABLE :: rho_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: rho_snow_mult_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_i_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_p_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_so_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_so_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_so_ice_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wliq_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wtot_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: dzh_snow_ptr(:)

  END TYPE t_lnd_prog


  !
  ! prognostic variables state vector (sea)
  !
  TYPE t_wtr_prog

    REAL(wp), POINTER         &
#ifdef _CRAYFTN
    , CONTIGUOUS              &
#endif
    &  ::                     &
    &  t_ice        (:,:)   , & ! temperature of the sea ice             (  K  )
    &  h_ice        (:,:)   , & ! depth of the sea ice                   (  m  )
    &  t_snow_si    (:,:)   , & ! temperature of the snow on sea ice     (  K  )
    &  h_snow_si    (:,:)   , & ! depth of the snow on sea ice           (  m  )
    &  t_snow_lk    (:,:)   , & ! temperature of snow on lake ice        (  K  )
    &  h_snow_lk    (:,:)   , & ! depth of snow on lake ice              (  K  )
    &  t_mnw_lk     (:,:)   , & ! mean temperature of the water column   (  K  )
    &  t_wml_lk     (:,:)   , & ! temperature of the mixed-layer         (  K  )
    &  h_ml_lk      (:,:)   , & ! thickness of the mixed-layer           (  m  )
    &  t_bot_lk     (:,:)   , & ! temperature at the water-bottom sediment interface  (  K  ) 
    &  c_t_lk       (:,:)   , & ! shape factor with respect to the       (  -  )
                                ! temperature profile in lake thermocline
    &  t_b1_lk      (:,:)   , & ! temperature at the bottom of the       (  K  ) 
                                ! upper layer of the sediments 
    &  h_b1_lk      (:,:)       ! thickness of the upper layer of the    (  K  )
                                ! sediments
  END TYPE t_wtr_prog

  !
  ! diagnostic variables state vector
  !
  TYPE t_lnd_diag

    REAL(wp), POINTER         &
#ifdef _CRAYFTN
    , CONTIGUOUS              &
#endif
    &  ::                     &
    &  qv_s         (:,:)   , & ! specific humidity at the surface              (kg/kg)
    &  t_s          (:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_seasfc     (:,:)   , & ! temperature of the sea surface                (  K  )
    &  w_i          (:,:)   , & ! water content of interception water           (m H2O)
    &  w_p          (:,:)   , & ! water content of pond interception water      (m H2O)
    &  w_s          (:,:)   , & ! water content of interception snow            (m H2O)
    &  t_so         (:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so         (:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice     (:,:,:) , & ! ice content                                   (m H20)
    &  runoff_s     (:,:)   , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g     (:,:)   , & ! soil water runoff; sum over forecast          (kg/m2)
    &  fr_seaice    (:,:)   , & !< fraction of sea ice                          ( )   
                                !< as partition of total area of the
                                !< grid element, but set to 0 or 1
                                !< index1=1,nproma, index2=1,nblks_c
    &  qv_s_t       (:,:,:) , & ! specific humidity at the surface              (kg/kg)
    &  runoff_s_t   (:,:,:) , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g_t   (:,:,:) , & ! soil water runoff; sum over forecast          (kg/m2)
    &  rstom        (:,:)   , & ! stomatal resistance                           ( s/m )
    &  rstom_t      (:,:,:) , & ! tile based stomatal resistance                ( s/m )
    &  t_snow       (:,:)   , & ! temperature of the snow-surface               (  K  )
    &  rho_snow     (:,:)   , & ! snow density                                  (kg/m**3)
    &  w_snow       (:,:)   , & ! snow water equivalent                         (m H2O)
    &  h_snow       (:,:)   , & ! snow height                                   (  m  ) 
    &  h_snow_t     (:,:,:) , & ! snow height                                   (  m  ) 
    &  freshsnow    (:,:)   , & ! indicator for age of snow in top of snow layer(  -  )
    &  freshsnow_t  (:,:,:) , & ! indicator for age of snow in top of snow layer(  -  )
    &  snowfrac     (:,:)   , & ! snow-cover fraction                           (  -  )
    &  snowfrac_t   (:,:,:) , & ! snow-cover fraction                           (  -  )
    &  snowfrac_lc_t(:,:,:) , & ! snow-cover fraction per land-cover class      (  -  )
    &  t_snow_mult  (:,:,:) , & ! temperature of snow                           (  K  )
    &  rho_snow_mult(:,:,:) , & ! snow density                                  (kg/m**3)
    &  wliq_snow    (:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  wtot_snow    (:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  dzh_snow     (:,:,:) , & ! layer thickness between half levels in snow   (  m  )
    &  w_snow_eff_t (:,:,:)     ! snow water equivalent / snow-cover fraction   !! OBSOLETE ? !! (m H2O)

    TYPE(t_ptr_2d3d), ALLOCATABLE :: qv_st_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: h_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: freshsnow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snowfrac_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snowfrac_lc_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_snow_eff_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: runoff_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: runoff_g_ptr(:)

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
