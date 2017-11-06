!>
!!        Contains the variables to set up the ocean model.
!=============================================================================================
!!
!! @par Revision History
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
MODULE mo_hamocc_types

  USE mo_kind,                ONLY: wp, sp
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    & success, max_char_length, min_dolic,               &
    & full_coriolis, beta_plane_coriolis,                &
    & f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates,      &
    & t_geographical_coordinates
  USE mo_linked_list,        ONLY: t_var_list


  PUBLIC :: t_hamocc_diag
  PUBLIC :: t_hamocc_monitor
  PUBLIC :: t_hamocc_state
  PUBLIC :: t_hamocc_tend
  PUBLIC :: t_hamocc_sed
  PUBLIC :: t_hamocc_acc
  PUBLIC :: t_hamocc_bcond


  !
  
  TYPE t_hamocc_monitor
    REAL(wp), POINTER :: phosy(:)
    REAL(wp), POINTER :: phosy_cya(:)
    REAL(wp), POINTER :: grazing(:)
    REAL(wp), POINTER :: omex90(:)
    REAL(wp), POINTER :: calex90(:)
    REAL(wp), POINTER :: opex90(:)
    REAL(wp), POINTER :: omex1000(:)
    REAL(wp), POINTER :: calex1000(:)
    REAL(wp), POINTER :: opex1000(:)
    REAL(wp), POINTER :: omex2000(:)
    REAL(wp), POINTER :: calex2000(:)
    REAL(wp), POINTER :: opex2000(:)
    REAL(wp), POINTER :: net_co2_flux(:)
    REAL(wp), POINTER :: sfalk(:)
    REAL(wp), POINTER :: sfdic(:)
    REAL(wp), POINTER :: sfphos(:)
    REAL(wp), POINTER :: sfsil(:)
    REAL(wp), POINTER :: bacfra(:)
    REAL(wp), POINTER :: phymor(:)
    REAL(wp), POINTER :: zoomor(:)
    REAL(wp), POINTER :: exudz(:)
    REAL(wp), POINTER :: graton(:)
    REAL(wp), POINTER :: exud(:)
    REAL(wp), POINTER :: n2fix(:)
    REAL(wp), POINTER :: wcdenit(:)
    REAL(wp), POINTER :: remina(:)
    REAL(wp), POINTER :: remins(:)
    REAL(wp), POINTER :: seddenit(:)
    REAL(wp), POINTER :: cyaldet(:)
    REAL(wp), POINTER :: cyaldoc(:)
    REAL(wp), POINTER :: delsil(:)
    REAL(wp), POINTER :: delcar(:)
    REAL(wp), POINTER :: zalkn2(:)
  END TYPE t_hamocc_monitor

  TYPE t_hamocc_diag
    !--------------------------------------------
    ! dimension: (nproma,n_zlev, nblks_e)
    REAL(wp), POINTER ::  phy(:,:,:)           
    REAL(wp), POINTER ::  zoo(:,:,:)       
    REAL(wp), POINTER ::  doc(:,:,:)         
    REAL(wp), POINTER ::  hi(:,:,:)          
    REAL(wp), POINTER ::  co3(:,:,:)     
    REAL(wp), POINTER ::  cya(:,:,:)         
    REAL(wp), POINTER ::  det(:,:,:)         
    REAL(wp), POINTER ::  dic(:,:,:)         
    REAL(wp), POINTER ::  alk(:,:,:)         
    REAL(wp), POINTER ::  no3(:,:,:)         
    REAL(wp), POINTER ::  po4(:,:,:)         
    REAL(wp), POINTER ::  n2(:,:,:)         
    REAL(wp), POINTER ::  o2(:,:,:)         
    REAL(wp), POINTER ::  si(:,:,:)         
    REAL(wp), POINTER ::  iron(:,:,:)           
    REAL(wp), POINTER ::  n2o(:,:,:)         
    REAL(wp), POINTER ::  dms(:,:,:)         
    REAL(wp), POINTER ::  h2s(:,:,:)         
    REAL(wp), POINTER ::  calc(:,:,:)         
    REAL(wp), POINTER ::  opal(:,:,:)         
    REAL(wp), POINTER ::  dust(:,:,:)         
    !--------------------------------------------
  END TYPE t_hamocc_diag
  !

  TYPE t_hamocc_sed
    !--------------------------------------------
    ! dimension: (nproma, ks, nblks_e)
    REAL(wp), POINTER ::  so12(:,:,:)       
    REAL(wp), POINTER ::  sc12(:,:,:)       
    REAL(wp), POINTER ::  ssil(:,:,:)       
    REAL(wp), POINTER ::  ster(:,:,:)       
    REAL(wp), POINTER ::  pwic(:,:,:)       
    REAL(wp), POINTER ::  pwal(:,:,:)       
    REAL(wp), POINTER ::  pwph(:,:,:)       
    REAL(wp), POINTER ::  pwox(:,:,:)       
    REAL(wp), POINTER ::  pwsi(:,:,:)       
    REAL(wp), POINTER ::  pwfe(:,:,:)       
    REAL(wp), POINTER ::  pwh2s(:,:,:)       
    REAL(wp), POINTER ::  pwn2(:,:,:)       
    REAL(wp), POINTER ::  pwno3(:,:,:)       
    REAL(wp), POINTER ::  sedhi(:,:,:)       
    REAL(wp), POINTER ::  pwh2ob(:,:,:)       
    REAL(wp), POINTER ::  pwn2b(:,:,:)       
    REAL(wp), POINTER ::  bo12(:,:)       
    REAL(wp), POINTER ::  bc12(:,:)       
    REAL(wp), POINTER ::  bsil(:,:)       
    REAL(wp), POINTER ::  bter(:,:)       
    INTEGER, POINTER ::  kbo(:,:)       
    REAL(wp), POINTER ::  bolay(:,:)       
  END TYPE t_hamocc_sed

  TYPE t_hamocc_acc
    !--------------------------------------------
    ! dimension: (nproma,n_zlev, nblks_e)
    REAL(wp), POINTER ::  npp(:,:,:)           
    REAL(wp), POINTER ::  phoc(:,:,:)           
    REAL(wp), POINTER ::  cyloss(:,:,:)           
    REAL(wp), POINTER ::  graz(:,:,:)       
    REAL(wp), POINTER ::  graton(:,:,:)       
    REAL(wp), POINTER ::  exud(:,:,:)       
    REAL(wp), POINTER ::  exudz(:,:,:)       
    REAL(wp), POINTER ::  zoomor(:,:,:)       
    REAL(wp), POINTER ::  phymor(:,:,:)       
    REAL(wp), POINTER ::  nfix(:,:,:)       
    REAL(wp), POINTER ::  remina(:,:,:)       
    REAL(wp), POINTER ::  remins(:,:,:)       
    REAL(wp), POINTER ::  reminn(:,:,:)       
    REAL(wp), POINTER ::  bacfra(:,:,:)       
    REAL(wp), POINTER ::  co2mr(:,:)       
    REAL(wp), POINTER ::  cflux(:,:)       
    REAL(wp), POINTER ::  nflux(:,:)       
    REAL(wp), POINTER ::  n2oflux(:,:)       
    REAL(wp), POINTER ::  dmsflux(:,:)       
    REAL(wp), POINTER ::  oflux(:,:)       
    REAL(wp), POINTER ::  nfixd(:,:)       
    REAL(wp), POINTER ::  delsil(:,:,:)       
    REAL(wp), POINTER ::  delcar(:,:,:)       
    REAL(wp), POINTER ::  prcaca(:,:)       
    REAL(wp), POINTER ::  produs(:,:)       
    REAL(wp), POINTER ::  prorca(:,:)       
    REAL(wp), POINTER ::  silpro(:,:)       
    REAL(wp), POINTER ::  coex90(:,:)       
    REAL(wp), POINTER ::  opex90(:,:)       
    REAL(wp), POINTER ::  calex90(:,:)       
    REAL(wp), POINTER ::  coex1000(:,:)       
    REAL(wp), POINTER ::  opex1000(:,:)       
    REAL(wp), POINTER ::  calex1000(:,:)       
    REAL(wp), POINTER ::  coex2000(:,:)       
    REAL(wp), POINTER ::  opex2000(:,:)       
    REAL(wp), POINTER ::  calex2000(:,:)       
    REAL(wp), POINTER ::  orginp(:,:)   ! for now constant value read via nml    
    REAL(wp), POINTER ::  silinp(:,:)   ! later riverine input possible     
    REAL(wp), POINTER ::  calinp(:,:)   ! via mo_bgc_bcond    
    REAL(wp), POINTER ::  h2obudget(:,:,:)       
    REAL(wp), POINTER ::  n2budget(:,:,:)       
    REAL(wp), POINTER ::  sedflic(:,:)       
    REAL(wp), POINTER ::  sedflal(:,:)       
    REAL(wp), POINTER ::  sedflph(:,:)       
    REAL(wp), POINTER ::  sedflox(:,:)       
    REAL(wp), POINTER ::  sedflsi(:,:)       
    REAL(wp), POINTER ::  sedflfe(:,:)       
    REAL(wp), POINTER ::  sedfln2(:,:)       
    REAL(wp), POINTER ::  sedflno3(:,:)       
    REAL(wp), POINTER ::  sedflh2s(:,:)       
    REAL(wp), POINTER ::  sedro2(:,:,:)       
    REAL(wp), POINTER ::  sedrn(:,:,:)       
    REAL(wp), POINTER ::  sedrs(:,:,:)       
    REAL(wp), POINTER ::  dmsprod(:,:,:)       
    REAL(wp), POINTER ::  dmsbac(:,:,:)       
    REAL(wp), POINTER ::  dmsuv(:,:,:)       
    REAL(wp), POINTER ::  euexp(:,:,:)       
    REAL(wp), POINTER ::  plim(:,:,:)       
    REAL(wp), POINTER ::  flim(:,:,:)       
    REAL(wp), POINTER ::  nlim(:,:,:)       
    REAL(wp), POINTER ::  aksp(:,:,:)       
    REAL(wp), POINTER ::  akb(:,:,:)       
    REAL(wp), POINTER ::  akw(:,:,:)       
    REAL(wp), POINTER ::  ak1(:,:,:)       
    REAL(wp), POINTER ::  ak2(:,:,:)       
    REAL(wp), POINTER ::  satoxy(:,:,:)       
    REAL(wp), POINTER ::  satn2(:,:)       
    REAL(wp), POINTER ::  satn2o(:,:)       
    REAL(wp), POINTER ::  solco2(:,:)       
    REAL(wp), POINTER ::  aou(:,:,:)       
    REAL(wp), POINTER ::  cTlim(:,:,:)       
    REAL(wp), POINTER ::  cLlim(:,:,:)       
    REAL(wp), POINTER ::  cPlim(:,:,:)       
    REAL(wp), POINTER ::  cFlim(:,:,:)       
    REAL(wp), POINTER ::  o2min(:,:)       
    REAL(wp), POINTER ::  zo2min(:,:)       
    REAL(wp), POINTER ::  h2sprod(:,:,:)       
    REAL(wp), POINTER ::  h2sloss(:,:,:)       
  END TYPE t_hamocc_acc

  TYPE, EXTENDS(t_hamocc_acc):: t_hamocc_tend
    ! all tendencies (not all are accumulated)
    !--------------------------------------------

    TYPE(t_hamocc_monitor) :: monitor
  END TYPE t_hamocc_tend

! array of states
   TYPE t_hamocc_state

    TYPE(t_hamocc_diag) :: p_diag
    TYPE(t_hamocc_tend) :: p_tend
    TYPE(t_hamocc_sed)  :: p_sed
    TYPE(t_hamocc_acc)  :: p_acc

  END TYPE t_hamocc_state

  TYPE t_hamocc_bcond
   REAL(wp), POINTER:: dusty(:,:)       !  index1=1,nproma, nblks_e
   REAL(wp), POINTER:: nitro(:,:)       !  index1=1,nproma, nblks_e
  END TYPE t_hamocc_bcond

END MODULE 


