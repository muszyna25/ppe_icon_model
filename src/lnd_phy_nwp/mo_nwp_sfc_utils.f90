!>
!! Utility routines related to the TERRA surface model
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_sfc_utils

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, zml_soil
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_parallel_config,     ONLY: nproma
  USe mo_extpar_config,       ONLY: itopo
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, nsfc_subs, t_tiles,  &
    &                               lseaice, llake, lmulti_snow, idiag_snowfrac
  USE mo_soil_ml,             ONLY: terra_multlay_init
  USE mo_phyparam_soil,       ONLY: cf_snow     ! soil and vegetation parameters for TILES
  USE mo_physical_constants,  ONLY: rdocp => rd_o_cpd  ! r_d / cp_d
!  USE mo_aggregate_surface,   ONLY: subsmean,subs_disaggregate_radflux,subsmean_albedo
  
  IMPLICIT NONE 

  PRIVATE

#ifdef __SX__
! parameters for loop unrolling
INTEGER, PARAMETER :: nlsoil= 7
INTEGER, PARAMETER :: nlsnow= 2
#endif

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: nwp_surface_init
  PUBLIC :: diag_snowfrac_tg
  PUBLIC :: aggregate_landvars

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init surface model TERRA
  !!
  !! Init surface model TERRA.
  !!
  !! @par Revision History
  !! Initial revision by Ekaterina Machulskaya, DWD (2011-07-??)
  !! Modification by Daniel Reienrt, DWD (2011-07-29)
  !! - initialize climatological layer t_so(nlev_soil+2)
  !!
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now, &
    &                          p_prog_lnd_new, p_lnd_diag )
 
                             

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(IN)    :: ext_data
!    TYPE(t_tiles)        , INTENT(INOUT) :: p_tiles(:)
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag
!!$    REAL(wp)             , INTENT(IN)   :: subsfrac(nproma,1,nsfc_subs)
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jc,jb,isubs,jk


    REAL(wp) :: t_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_snow_mult_now_t(nproma, 1:nlev_snow+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_s_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_g_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_s_new_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rho_snow_now_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rho_snow_mult_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_so_now_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_so_new_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_ice_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: w_so_ice_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wliq_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: wtot_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: dzh_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, nsfc_subs)

    REAL(wp) ::          freshsnow_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          snowfrac_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) ::          sso_sigma_t(nproma,  p_patch%nblks_c)
    INTEGER  ::          lc_class_t (nproma,  p_patch%nblks_c, nsfc_subs)

!    INTEGER  :: i_tile(nproma,nsfc_subs),lu_subs
    INTEGER  :: i_count, ic

  !-------------------------------------------------------------------------


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,ic,jk), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF (itopo == 1) THEN
        DO isubs = 1, nsfc_subs
          DO jc = i_startidx, i_endidx

            ! initialize climatological layer (deepest layer of t_so)
            p_prog_lnd_now%t_so_t(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)
            p_prog_lnd_new%t_so_t(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)

            p_prog_lnd_now%t_g_t(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)

          END DO 
        END DO
      ENDIF



      IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point

!---------- Copy input fields for each tile

!----------------------------------
      DO isubs = 1,nsfc_subs
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs) 

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          t_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%t_snow_t(jc,jb,isubs) 
          t_s_now_t(ic,jb,isubs)             =  p_prog_lnd_now%t_s_t(jc,jb,isubs)   
          t_s_new_t(ic,jb,isubs)             =  p_prog_lnd_new%t_s_t(jc,jb,isubs)   
          w_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%w_snow_t(jc,jb,isubs)  
          rho_snow_now_t(ic,jb,isubs)        =  p_prog_lnd_now%rho_snow_t(jc,jb,isubs)

          sso_sigma_t(ic,jb)                 = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic,jb,isubs)            = ext_data%atm%lc_class_t(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           = p_lnd_diag%freshsnow_t(jc,jb,isubs)
        ENDDO


        IMSNOWI: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
              t_snow_mult_now_t(ic,jk,jb,isubs)   =  p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) 
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
              rho_snow_mult_now_t(ic,jk,jb,isubs) =  p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs)
              wliq_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wliq_snow_t    (jc,jk,jb,isubs) 
              wtot_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wtot_snow_t    (jc,jk,jb,isubs) 
              dzh_snow_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%dzh_snow_t     (jc,jk,jb,isubs) 
            ENDDO
          ENDDO

        END IF  IMSNOWI

!CDIR UNROLL=nlsoil+2
        DO jk=1,nlev_soil+2
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            t_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) 
            t_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            w_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) 
            w_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) 
            w_so_ice_new_t(ic,jk,jb,isubs)      =  p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

         
        CALL terra_multlay_init(                                  &
        &  ie=nproma,                                             & ! array dimensions
        &  istartpar=1, iendpar= i_count,                         & ! optional start/end indicies
!       &  ke=nlev, &! nsubs0=1, nsubs1=nsfc_subs               , & ! nsfc_subs
        &  ke_soil=nlev_soil, ke_snow=nlev_snow                 , &
        &  czmls=zml_soil                                       , & ! processing soil level structure 
        &  soiltyp_subs  =  ext_data%atm%soiltyp_t(:,jb,isubs)  , & ! type of the soil (keys 0-9)  --
        &  rootdp        =  ext_data%atm%rootdp_t(:,jb,isubs)   , & ! depth of the roots                ( m  )
        &  t_snow_now    =  t_snow_now_t(:,jb,isubs)            , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_now   = t_snow_mult_now_t(:,:,jb,isubs)  , & ! temperature of the snow-surface   (  K  )
        &  t_s_now           = t_s_now_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  w_snow_now        = w_snow_now_t(:,jb,isubs)         , & ! water content of snow             (m H2O)
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)       , & ! snow density                      (kg/m**3)
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs), & ! snow density                      (kg/m**3)
        &  t_so_now          = t_so_now_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  t_so_new          = t_so_new_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  w_so_now          = w_so_now_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_new          = w_so_new_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_ice_now      = w_so_ice_now_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  w_so_ice_new      = w_so_ice_new_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  wliq_snow_now     = wliq_snow_now_t(:,:,jb,isubs)    , & ! liquid water content in the snow  (m H2O)
        &  wtot_snow_now     = wtot_snow_now_t(:,:,jb,isubs)    , & ! total (liquid + solid) water content of snow (m H2O)
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs)       & ! layer thickness between half levels in snow  (  m  )
                                                      )


       IF (lmulti_snow) THEN
!CDIR NOIEXPAND
          CALL diag_snowfrac_tg(                           &
            &  istart = 1, iend = i_count                , & ! start/end indices
            &  z0_lcc    = ext_data%atm%z0_lcc(:)        , & ! roughness length
            &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
            &  t_snow    = t_snow_mult_now_t (:,2,jb,isubs), & ! snow temp
            &  t_soiltop = t_s_now_t         (:,jb,isubs), & ! soil top temp
            &  w_snow    = w_snow_now_t      (:,jb,isubs), & ! snow WE
            &  rho_snow  = rho_snow_now_t    (:,jb,isubs), & ! snow depth
            &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma = sso_sigma_t       (:,jb),       & ! sso stdev
            &  tai       = ext_data%atm%tai_t(:,jb,isubs), & ! effective leaf area index
            &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp
        ELSE
!CDIR NOIEXPAND
          CALL diag_snowfrac_tg(                           &
            &  istart = 1, iend = i_count                , & ! start/end indices
            &  z0_lcc    = ext_data%atm%z0_lcc(:)        , & ! roughness length
            &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
            &  t_snow    = t_snow_now_t      (:,jb,isubs), & ! snow temp
            &  t_soiltop = t_s_now_t         (:,jb,isubs), & ! soil top temp
            &  w_snow    = w_snow_now_t      (:,jb,isubs), & ! snow WE
            &  rho_snow  = rho_snow_now_t    (:,jb,isubs), & ! snow depth
            &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma = sso_sigma_t       (:,jb),       & ! sso stdev
            &  tai       = ext_data%atm%tai_t(:,jb,isubs), & ! effective leaf area index
            &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp
        ENDIF

!  Recover fields from index list
!
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          p_prog_lnd_now%t_snow_t(jc,jb,isubs)   = t_snow_now_t(ic,jb,isubs)
          p_prog_lnd_now%t_s_t(jc,jb,isubs)      = t_s_now_t(ic,jb,isubs)  
          p_prog_lnd_new%t_s_t(jc,jb,isubs)      = t_s_new_t(ic,jb,isubs) 
          p_prog_lnd_now%w_snow_t(jc,jb,isubs)   = w_snow_now_t(ic,jb,isubs) 
          p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = rho_snow_now_t(ic,jb,isubs)
          p_lnd_diag %snowfrac_t(jc,jb,isubs)    = snowfrac_t(ic,jb,isubs)
          p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          p_prog_lnd_new%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
        ENDDO

        IMSNOWO: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
              p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) =  t_snow_mult_now_t(ic,jk,jb,isubs)   
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
              p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs) = rho_snow_mult_now_t(ic,jk,jb,isubs) 
              p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)   
              p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)
              p_prog_lnd_now%dzh_snow_t(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)    
            ENDDO
          ENDDO

        END IF  IMSNOWO

!CDIR UNROLL=nlsoil+2
        DO jk=1,nlev_soil+2
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) = t_so_now_t(ic,jk,jb,isubs)          
            p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)          
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) = w_so_now_t(ic,jk,jb,isubs)        
            p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)        
            p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO

      END DO ! isubs

    ENDDO  ! jb loop
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_surface_init


!-------------------------------------------------------------------------

  SUBROUTINE aggregate_landvars( p_patch, ext_data, lnd_prog, lnd_diag)

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(in)   :: ext_data      !< external data
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog      !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc


    ! Local scalars:
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    INTEGER :: jc,jb,isubs,jk
    INTEGER :: i_count, ic

    REAL(wp) :: tilefrac ! fractional area covered by tile


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,ic,jk,tilefrac) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      i_count = ext_data%atm%lp_count(jb) ! for aggregation process land points only

      IF (nsfc_subs == 1) THEN  ! just copy prognostic variables from tile 1 to diagnostic aggregated variable
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          lnd_diag%t_snow(jc,jb)       = lnd_prog%t_snow_t(jc,jb,1)
          lnd_diag%t_s(jc,jb)          = lnd_prog%t_s_t(jc,jb,1)
          lnd_diag%w_snow(jc,jb)       = lnd_prog%w_snow_t(jc,jb,1)
          lnd_diag%rho_snow(jc,jb)     = lnd_prog%rho_snow_t(jc,jb,1)
          lnd_diag%w_i(jc,jb)          = lnd_prog%w_i_t(jc,jb,1)
          lnd_diag%h_snow(jc,jb)       = lnd_diag%h_snow_t(jc,jb,1)
          lnd_diag%freshsnow(jc,jb)    = lnd_diag%freshsnow_t(jc,jb,1)
          lnd_diag%snowfrac(jc,jb)     = lnd_diag%snowfrac_t(jc,jb,1)
          lnd_diag%runoff_s(jc,jb)     = lnd_diag%runoff_s_t(jc,jb,1)
          lnd_diag%runoff_g(jc,jb)     = lnd_diag%runoff_g_t(jc,jb,1)
          lnd_diag%t_so(jc,nlev_soil+2,jb) = lnd_prog%t_so_t(jc,nlev_soil+2,jb,1)

          IF(lmulti_snow) THEN
            lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,1)
          ENDIF
        ENDDO

        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            lnd_diag%t_so(jc,jk,jb)      = lnd_prog%t_so_t(jc,jk,jb,1)
            lnd_diag%w_so(jc,jk,jb)      = lnd_prog%w_so_t(jc,jk,jb,1)
            lnd_diag%w_so_ice(jc,jk,jb)  = lnd_prog%w_so_ice_t(jc,jk,jb,1)
          ENDDO
        ENDDO

        IF (lmulti_snow) THEN
          DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp(ic,jb)
              lnd_diag%t_snow_mult(jc,jk,jb)   = lnd_prog%t_snow_mult_t(jc,jk,jb,1)
              lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_prog%rho_snow_mult_t(jc,jk,jb,1)
              lnd_diag%wliq_snow(jc,jk,jb)     = lnd_prog%wliq_snow_t(jc,jk,jb,1)
              lnd_diag%wtot_snow(jc,jk,jb)     = lnd_prog%wtot_snow_t(jc,jk,jb,1)
              lnd_diag%dzh_snow(jc,jk,jb)      = lnd_prog%dzh_snow_t(jc,jk,jb,1)
            ENDDO
          ENDDO
        ENDIF

      ELSE ! aggregate fields over tiles

        ! First initialize fields to zero in order to prepare subsequent summation over the tiles

        lnd_diag%t_snow(:,jb)      = 0._wp
        lnd_diag%t_s(:,jb)         = 0._wp
        lnd_diag%w_snow(:,jb)      = 0._wp
        lnd_diag%rho_snow(:,jb)    = 0._wp
        lnd_diag%w_i(:,jb)         = 0._wp
        lnd_diag%h_snow(:,jb)      = 0._wp
        lnd_diag%freshsnow(:,jb)   = 0._wp
        lnd_diag%snowfrac(:,jb)    = 0._wp
        lnd_diag%runoff_s(:,jb)    = 0._wp
        lnd_diag%runoff_g(:,jb)    = 0._wp

        lnd_diag%t_so(:,:,jb)      = 0._wp
        lnd_diag%w_so(:,:,jb)      = 0._wp
        lnd_diag%w_so_ice(:,:,jb)  = 0._wp
 
        IF (lmulti_snow) THEN
          lnd_diag%t_snow_mult(:,:,jb)   = 0._wp
          lnd_diag%rho_snow_mult(:,:,jb) = 0._wp
          lnd_diag%wliq_snow(:,:,jb)     = 0._wp
          lnd_diag%wtot_snow(:,:,jb)     = 0._wp
          lnd_diag%dzh_snow(:,:,jb)      = 0._wp
        ENDIF

        DO isubs = 1,nsfc_subs
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            tilefrac = ext_data%atm%lc_frac_t(jc,jb,isubs)
            lnd_diag%t_snow(jc,jb)       = lnd_diag%t_snow(jc,jb) + tilefrac * &
                                           lnd_prog%t_snow_t(jc,jb,isubs)
            lnd_diag%t_s(jc,jb)          = lnd_diag%t_s(jc,jb) + tilefrac * &
                                           lnd_prog%t_s_t(jc,jb,isubs)
            lnd_diag%w_snow(jc,jb)       = lnd_diag%w_snow(jc,jb) + tilefrac * &
                                           lnd_prog%w_snow_t(jc,jb,isubs)
            lnd_diag%rho_snow(jc,jb)     = lnd_diag%rho_snow(jc,jb) + tilefrac * &
                                           lnd_prog%rho_snow_t(jc,jb,isubs)
            lnd_diag%w_i(jc,jb)          = lnd_diag%w_i(jc,jb) + tilefrac * &
                                           lnd_prog%w_i_t(jc,jb,isubs)
            lnd_diag%h_snow(jc,jb)       = lnd_diag%h_snow(jc,jb) + tilefrac * &
                                           lnd_diag%h_snow_t(jc,jb,isubs)
            lnd_diag%freshsnow(jc,jb)    = lnd_diag%freshsnow(jc,jb) + tilefrac * &
                                           lnd_diag%freshsnow_t(jc,jb,isubs)
            lnd_diag%snowfrac(jc,jb)     = lnd_diag%snowfrac(jc,jb) + tilefrac * &
                                           lnd_diag%snowfrac_t(jc,jb,isubs)
            lnd_diag%runoff_s(jc,jb)     = lnd_diag%runoff_s(jc,jb) + tilefrac * &
                                           lnd_diag%runoff_s_t(jc,jb,isubs)
            lnd_diag%runoff_g(jc,jb)     = lnd_diag%runoff_g(jc,jb) + tilefrac * &
                                           lnd_diag%runoff_g_t(jc,jb,isubs)
            lnd_diag%t_so(jc,nlev_soil+2,jb) = lnd_diag%t_so(jc,nlev_soil+2,jb) + tilefrac * &
                                               lnd_prog%t_so_t(jc,nlev_soil+2,jb,isubs)

            IF(lmulti_snow) THEN
              lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_diag%t_snow_mult(jc,nlev_snow+1,jb)+ &
                                        tilefrac * lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,isubs)
            ENDIF
          ENDDO

          DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp(ic,jb)
              tilefrac = ext_data%atm%lc_frac_t(jc,jb,isubs)
              lnd_diag%t_so(jc,jk,jb)      = lnd_diag%t_so(jc,jk,jb) + tilefrac* &
                                             lnd_prog%t_so_t(jc,jk,jb,isubs)
              lnd_diag%w_so(jc,jk,jb)      = lnd_diag%w_so(jc,jk,jb) + tilefrac * &
                                             lnd_prog%w_so_t(jc,jk,jb,isubs)
              lnd_diag%w_so_ice(jc,jk,jb)  = lnd_diag%w_so_ice(jc,jk,jb) + tilefrac * &
                                             lnd_prog%w_so_ice_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

          IF (lmulti_snow) THEN
            DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
              DO ic = 1, i_count
                jc = ext_data%atm%idx_lst_lp(ic,jb)
                tilefrac = ext_data%atm%lc_frac_t(jc,jb,isubs)
                lnd_diag%t_snow_mult(jc,jk,jb)   = lnd_diag%t_snow_mult(jc,jk,jb) + tilefrac * &
                                                   lnd_prog%t_snow_mult_t(jc,jk,jb,isubs)
                lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_diag%rho_snow_mult(jc,jk,jb) + tilefrac * &
                                                   lnd_prog%rho_snow_mult_t(jc,jk,jb,isubs)
                lnd_diag%wliq_snow(jc,jk,jb)     = lnd_diag%wliq_snow(jc,jk,jb) + tilefrac * &
                                                   lnd_prog%wliq_snow_t(jc,jk,jb,isubs)
                lnd_diag%wtot_snow(jc,jk,jb)     = lnd_diag%wtot_snow(jc,jk,jb) + tilefrac * &
                                                   lnd_prog%wtot_snow_t(jc,jk,jb,isubs)
                lnd_diag%dzh_snow(jc,jk,jb)      = lnd_diag%dzh_snow(jc,jk,jb) + tilefrac * &
                                                   lnd_prog%dzh_snow_t(jc,jk,jb,isubs)
              ENDDO
            ENDDO
          ENDIF

        ENDDO
      ENDIF

     i_count = ext_data%atm%sp_count(jb) ! second step: just copy variables for sea points 

!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_sp(ic,jb)
        lnd_diag%t_snow(jc,jb)       = lnd_prog%t_snow_t(jc,jb,1)
        lnd_diag%t_s(jc,jb)          = lnd_prog%t_s_t(jc,jb,1)
        lnd_diag%w_snow(jc,jb)       = lnd_prog%w_snow_t(jc,jb,1)
        lnd_diag%rho_snow(jc,jb)     = lnd_prog%rho_snow_t(jc,jb,1)
        lnd_diag%w_i(jc,jb)          = lnd_prog%w_i_t(jc,jb,1)
        lnd_diag%h_snow(jc,jb)       = lnd_diag%h_snow_t(jc,jb,1)
        lnd_diag%freshsnow(jc,jb)    = lnd_diag%freshsnow_t(jc,jb,1)
        lnd_diag%snowfrac(jc,jb)     = lnd_diag%snowfrac_t(jc,jb,1)
        lnd_diag%runoff_s(jc,jb)     = lnd_diag%runoff_s_t(jc,jb,1)
        lnd_diag%runoff_g(jc,jb)     = lnd_diag%runoff_g_t(jc,jb,1)
        lnd_diag%t_so(jc,nlev_soil+2,jb) = lnd_prog%t_so_t(jc,nlev_soil+2,jb,1)

        IF(lmulti_snow) THEN
          lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,1)
        ENDIF
      ENDDO

      DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_sp(ic,jb)
          lnd_diag%t_so(jc,jk,jb)      = lnd_prog%t_so_t(jc,jk,jb,1)
          lnd_diag%w_so(jc,jk,jb)      = lnd_prog%w_so_t(jc,jk,jb,1)
          lnd_diag%w_so_ice(jc,jk,jb)  = lnd_prog%w_so_ice_t(jc,jk,jb,1)
        ENDDO
      ENDDO

      IF (lmulti_snow) THEN
        DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_sp(ic,jb)
            lnd_diag%t_snow_mult(jc,jk,jb)   = lnd_prog%t_snow_mult_t(jc,jk,jb,1)
            lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_prog%rho_snow_mult_t(jc,jk,jb,1)
            lnd_diag%wliq_snow(jc,jk,jb)     = lnd_prog%wliq_snow_t(jc,jk,jb,1)
            lnd_diag%wtot_snow(jc,jk,jb)     = lnd_prog%wtot_snow_t(jc,jk,jb,1)
            lnd_diag%dzh_snow(jc,jk,jb)      = lnd_prog%dzh_snow_t(jc,jk,jb,1)
          ENDDO
        ENDDO
      ENDIF

    ENDDO    
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE aggregate_landvars


  SUBROUTINE subsmean(i_count, t_g_s, qv_s_s, t_s_s, t_snow_s, t_snow_mult_s, &
        &  w_snow_s, wtot_snow_s, rho_snow_s, rho_snow_mult_s, h_snow_s, freshsnow_s,   &
        &  w_i_s, t_so_s, w_so_s, w_so_ice_s, runoff_s_s, runoff_g_s, tch_s, tfv_s,     &
        &  t_2m_s, qv_2m_s, td_2m_s, rh_2m_s, u_10m_s, v_10m_s, shfl_s_s, lhfl_s_s,     &
        &  ie, nsubs0, nsubs1, ke_soil, ke_snow, idx_lst_lndp,      &
        &  subsfrac, t, qv, p0, ps, t_g_now, t_g_new, qv_s, t_s, t_snow, t_snow_mult,   &
        &  w_snow, wtot_snow, rho_snow, rho_snow_mult, h_snow, freshsnow, w_i, t_so,    &
        &  w_so, w_so_ice, runoff_s, runoff_g, tch, tfv, t_2m, qv_2m, td_2m, rh_2m,     &
        &  u_10m, v_10m, shfl_s, lhfl_s)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ie, i_count, nsubs0, nsubs1, ke_soil, ke_snow 
  INTEGER, DIMENSION(i_count), INTENT(IN) ::                     &
                  idx_lst_lndp         ! index list for land points                    (  -  )
  REAL(wp), DIMENSION(ie), INTENT(IN) ::                    &
                  t,                 & ! temperature                                   (  K  )
                  qv,                & ! humidity                                      (kg/kg)
                  p0,                & ! pressure                                      ( Pa  )
                  ps                   ! surface pressure                              ( Pa  )
  REAL(wp), DIMENSION(ie,nsubs1), INTENT(IN) ::             &
                  subsfrac,          & ! fractions of tiles                            (  -  )
                  t_g_now,           & ! surface temperature                           (  K  )
                  t_g_new,           & ! surface temperature                           (  K  )
                  qv_s,              & ! specific humidity at the surface              (kg/kg)
                  t_s,               & ! temperature of the ground surface             (  K  )
                  t_snow,            & ! temperature of the snow                       (  K  )
                  w_snow,            & ! water content of snow                         (m H2O)
                  rho_snow,          & ! snow density                                  (kg/m**3)
                  h_snow,            & ! snow height
                  w_i,               & ! water content of interception water           (m H2O)
                  freshsnow,         & ! indicator for age of snow in top of snow layer(  -  )
                  runoff_s,          & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g,          & ! soil water runoff; sum over forecast          (kg/m2)
                  tch,               & ! turbulent transfer coefficient for heat       ( -- )
                  tfv,               & ! laminar reduction factor for evaporation      ( -- )
                  t_2m,              & ! temperature in 2m                             (  K  )
                  qv_2m,             & ! humidity in 2m                                (kg/kg)
                  td_2m,             & ! dew point in 2m                               (  K  )
                  rh_2m,             & ! relative humidity in 2m                       (  K  )
                  u_10m,             & ! zonal wind in 10m                             ( m/s )
                  v_10m,             & ! meridional wind in 10m                        ( m/s )
                  shfl_s,            & ! sensible heat flux                            (W/m2 )
                  lhfl_s               ! latent heat flux                              (W/m2 )
  REAL(wp), DIMENSION(ie,0:ke_snow,nsubs1), INTENT(IN) ::   &
                  t_snow_mult          ! temperature of the snow-surface               (  K  )
  REAL(wp), DIMENSION(ie,1:ke_snow,nsubs1), INTENT(IN) ::   &
                  rho_snow_mult,     & ! snow density                                  (kg/m**3)
                  wtot_snow            ! total water content of snow layers            (m H2O)  
  REAL(wp), DIMENSION(ie,0:ke_soil+1,nsubs1), INTENT(IN) :: &
                  t_so                 ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(ie,ke_soil+1,nsubs1), INTENT(IN) ::   &
                  w_so,              & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice             ! ice content       
!mean over tiles:
  REAL(wp), DIMENSION(ie), INTENT(OUT) ::                   &
                  t_g_s,             & ! surface temperature                           (  K  )
                  qv_s_s,            & ! specific humidity at the surface              (kg/kg)
                  t_s_s,             & ! temperature of the ground surface             (  K  )
                  t_snow_s,          & ! temperature of the snow                       (  K  )
                  w_snow_s,          & ! water content of snow                         (m H2O)
                  rho_snow_s,        & ! snow density                                  (kg/m**3)
                  h_snow_s,          & ! snow height
                  w_i_s,             & ! water content of interception water           (m H2O)
                  freshsnow_s,       & ! indicator for age of snow in top of snow layer(  -  )
                  runoff_s_s,        & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g_s,        & ! soil water runoff; sum over forecast          (kg/m2)
                  tch_s,             & ! turbulent transfer coefficient for heat       ( -- )
                  tfv_s,             & ! laminar reduction factor for evaporation      ( -- )
                  t_2m_s,            & ! temperature in 2m                             (  K  )
                  qv_2m_s,           & ! humidity in 2m                                (kg/kg)
                  td_2m_s,           & ! dew point in 2m                               (  K  )
                  rh_2m_s,           & ! relative humidity in 2m                       (  K  )
                  u_10m_s,           & ! zonal wind in 10m                             ( m/s )
                  v_10m_s,           & ! meridional wind in 10m                        ( m/s )
                  shfl_s_s,          & ! sensible heat flux                            (W/m2 )
                  lhfl_s_s             ! latent heat flux                              (W/m2 )
  REAL(wp), DIMENSION(ie,0:ke_snow), INTENT(OUT) ::         &
                  t_snow_mult_s        ! temperature of the snow-surface               (  K  )
  REAL(wp), DIMENSION(ie,1:ke_snow), INTENT(OUT) ::         &
                  rho_snow_mult_s,   & ! snow density                                  (kg/m**3)
                  wtot_snow_s          ! total water content of snow layers            (m H2O)  
  REAL(wp), DIMENSION(ie,0:ke_soil+1), INTENT(OUT) ::       &
                  t_so_s               ! soil temperature (main level)                 (  K  )
  REAL(wp), DIMENSION(ie,ke_soil+1), INTENT(OUT) ::         &
                  w_so_s,            & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_s           ! ice content       

  ! Local scalars:

  INTEGER  :: ic,jc,jk,ns,kso,ksn      !loop indices
  INTEGER  :: counter
  REAL(wp), PARAMETER :: small = 1.0E-07_wp

  ! Local arrays:

  REAL(wp) :: tmp(i_count)
  LOGICAL  :: llandmask(i_count,nsubs1)
  REAL(wp) :: t_g_slp(i_count), qv_s_slp(i_count), t_s_slp(i_count), t_snow_slp(i_count),        &
    &         t_snow_mult_slp(i_count, nlev_snow), w_snow_slp(i_count),                          &
    &         wtot_snow_slp(i_count, nlev_snow),                                                 &
    &         rho_snow_slp(i_count), rho_snow_mult_slp(i_count, nlev_snow), h_snow_slp(i_count), &
    &         freshsnow_slp(i_count), w_i_slp(i_count), t_so_slp(i_count, nlev_soil+2),          &
    &         w_so_slp(i_count, nlev_soil+1), w_so_ice_slp(i_count, nlev_soil+1),                &
    &         runoff_s_slp(i_count),                                                             &
    &         runoff_g_slp(i_count), tch_slp(i_count), tfv_slp(i_count), t_2m_slp(i_count),      &
    &         qv_2m_slp(i_count), td_2m_slp(i_count), rh_2m_slp(i_count), u_10m_slp(i_count),    &
    &         v_10m_slp(i_count), shfl_s_slp(i_count), lhfl_s_slp(i_count)
  REAL(wp) :: t_lp(i_count), qv_lp(i_count), p0_lp(i_count), ps_lp(i_count),                     &
    &         t_g_now_tlp(i_count,nsubs1),                                                       &
    &         t_g_new_tlp(i_count,nsubs1), qv_s_tlp(i_count,nsubs1), t_s_tlp(i_count,nsubs1),    &
    &         t_snow_tlp(i_count,nsubs1), t_snow_mult_tlp(i_count,nlev_snow,nsubs1),             &
    &         w_snow_tlp(i_count,nsubs1), wtot_snow_tlp(i_count,nlev_snow,nsubs1),               &
    &         rho_snow_tlp(i_count,nsubs1), rho_snow_mult_tlp(i_count,nlev_snow,nsubs1),         &
    &         h_snow_tlp(i_count,nsubs1), freshsnow_tlp(i_count,nsubs1), w_i_tlp(i_count,nsubs1),&
    &         t_so_tlp(i_count,nlev_soil+2,nsubs1), w_so_tlp(i_count,nlev_soil+1,nsubs1),        &
    &         w_so_ice_tlp(i_count,nlev_soil+1,nsubs1), runoff_s_tlp(i_count,nsubs1),            &
    &         runoff_g_tlp(i_count,nsubs1), tch_tlp(i_count,nsubs1), tfv_tlp(i_count,nsubs1),    &
    &         t_2m_tlp(i_count,nsubs1), qv_2m_tlp(i_count,nsubs1), td_2m_tlp(i_count,nsubs1),    &
    &         rh_2m_tlp(i_count,nsubs1), u_10m_tlp(i_count,nsubs1), v_10m_tlp(i_count,nsubs1),   &
    &         shfl_s_tlp(i_count,nsubs1), lhfl_s_tlp(i_count,nsubs1)
!-------------------------------------------------------------------------

  tmp(:)    = 0._wp
  tch_s(:)  = 0._wp
  qv_s_s(:) = 0._wp

!Set llandmask           !EM: it is a constant field, should be set once, not at each time step!
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
      llandmask(ic,ns) = .TRUE.
    END DO
  END DO
  DO ic = 1, i_count
    llandmask(ic,20) = .FALSE.    !water
  END DO

  DO ic = 1, i_count       
    jc = idx_lst_lndp(ic)
    t_lp(ic) = t(jc)
    qv_lp(ic) = qv(jc)
    p0_lp(ic) = p0(jc)
    ps_lp(ic) = ps(jc)
  END DO

  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count       
      jc = idx_lst_lndp(ic)
      t_g_now_tlp  (ic,ns) = t_g_now  (jc,ns)
      t_g_new_tlp  (ic,ns) = t_g_new  (jc,ns)
      qv_s_tlp     (ic,ns) = qv_s     (jc,ns)
      t_s_tlp      (ic,ns) = t_s      (jc,ns)
      t_snow_tlp   (ic,ns) = t_snow   (jc,ns)
      w_snow_tlp   (ic,ns) = w_snow   (jc,ns)
      rho_snow_tlp (ic,ns) = rho_snow (jc,ns)
      h_snow_tlp   (ic,ns) = h_snow   (jc,ns)
      freshsnow_tlp(ic,ns) = freshsnow(jc,ns)
      w_i_tlp      (ic,ns) = w_i      (jc,ns)
      runoff_s_tlp (ic,ns) = runoff_s (jc,ns)
      runoff_g_tlp (ic,ns) = runoff_g (jc,ns)
      tch_tlp      (ic,ns) = tch      (jc,ns)
      tfv_tlp      (ic,ns) = tfv      (jc,ns)
      t_2m_tlp     (ic,ns) = t_2m     (jc,ns)
      qv_2m_tlp    (ic,ns) = qv_2m    (jc,ns)
      td_2m_tlp    (ic,ns) = td_2m    (jc,ns)
      rh_2m_tlp    (ic,ns) = rh_2m    (jc,ns)
      u_10m_tlp    (ic,ns) = u_10m    (jc,ns)
      v_10m_tlp    (ic,ns) = v_10m    (jc,ns)
      shfl_s_tlp   (ic,ns) = shfl_s   (jc,ns)
      lhfl_s_tlp   (ic,ns) = lhfl_s   (jc,ns)
      DO jk = 0, nlev_soil+1                         !EM order of loops? 
        t_so_tlp    (ic,jk,ns) = t_so    (jc,jk,ns)
      END DO 
      DO jk = 1, nlev_soil+1
        w_so_tlp    (ic,jk,ns) = w_so    (jc,jk,ns)
        w_so_ice_tlp(ic,jk,ns) = w_so_ice(jc,jk,ns)
      END DO 
      DO jk = 0, nlev_snow+1
        t_snow_mult_tlp  (ic,jk,ns) = t_snow_mult  (jc,jk,ns)
      END DO 
      DO jk = 1, nlev_snow+1
        wtot_snow_tlp    (ic,jk,ns) = wtot_snow    (jc,jk,ns)
        rho_snow_mult_tlp(ic,jk,ns) = rho_snow_mult(jc,jk,ns)
      END DO 
    END DO 
  END DO 

!mean (effective) transfer coefficient for scalars
  CALL subsmean_power4(i_count, nsubs0, nsubs1, t_g_now_tlp, t_g_slp, subsfrac)
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
!intermediate storage for mean sensible heat flux
      tmp(ic) = tmp(ic) + tch_tlp(ic,ns)*(t_g_now_tlp(ic,ns) - &
                t_lp(ic) * (ps_lp(ic)/p0_lp(ic))**rdocp)*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    tch_slp(ic) = tmp(ic)/SIGN(MAX(ABS(t_g_slp(ic) - t_lp(ic)                 &
      &           *(ps_lp(ic)/p0_lp(ic))**rdocp),small), t_g_slp(ic)-t_lp(ic) &
      &           *(ps_lp(ic)/p0_lp(ic))**rdocp)
    tch_slp(ic) = SIGN(MAX(ABS(tch_slp(ic)),small),tch_slp(ic))
  END DO

!mean (effective) surface humidity
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
!intermediate storage for mean latent heat flux
      tmp(ic) = tmp(ic) + tch_tlp(ic,ns)*(qv_s_tlp(ic,ns)-qv_lp(ic))*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    qv_s_slp(ic) = MAX(tmp(ic)/tch_slp(ic) + qv_lp(ic),small)
  END DO

  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_s_tlp, t_s_slp, subsfrac)
  CALL subsmean_power4    (i_count, nsubs0, nsubs1, t_g_new_tlp, t_g_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_snow_tlp, t_snow_slp, subsfrac, &
    &                      w_snow_tlp>0.0_wp)
  DO ksn = 0, ke_snow
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_snow_mult_tlp(:,:,ksn),   &
      &                      t_snow_mult_slp(:,ksn), subsfrac, llandmask)
  END DO
  DO ic = 1, i_count
    counter = COUNT(w_snow_tlp(ic,:) .GT. small)
    IF(counter .EQ. 0) t_snow_slp(ic) = t_s_slp(ic)
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_snow_tlp, w_snow_slp, subsfrac, llandmask)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, rho_snow_tlp,  rho_snow_slp, subsfrac,  &
    &                      w_snow_tlp>0.0_wp)
  DO ksn = 1, ke_snow
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, rho_snow_mult_tlp(:,:,ksn), &
                             rho_snow_mult_slp(:,ksn), subsfrac, w_snow_tlp>0.0_wp)
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, wtot_snow_tlp(:,:,ksn),     &
      &                      wtot_snow_slp(:,ksn), subsfrac, llandmask)
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, h_snow_tlp, h_snow_slp, subsfrac, llandmask)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, freshsnow_tlp, freshsnow_slp, subsfrac,  &
    &                      w_snow_tlp>0.0_wp)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_i_tlp, w_i_slp, subsfrac, llandmask)
  DO kso = 0, ke_soil+1
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, t_so_tlp(:,:,kso), t_so_slp(:,kso),    &
      &                      subsfrac, llandmask)
  END DO
  DO kso = 1, ke_soil+1
!    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_tlp(:,:,kso), w_so_slp(:,kso), subsfrac, xlsmmask)
!    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_ice_tlp(:,:,kso), w_so_ice_slp(:,kso), subsfrac, xlsmmask)
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_tlp(:,:,kso), &
      w_so_slp(:,kso), subsfrac, llandmask)           !temporary; should not be averaged over rock
    CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, w_so_ice_tlp(:,:,kso), &
      w_so_ice_slp(:,kso), subsfrac, llandmask)   !temporary; should not be averaged over rock   
  END DO
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, runoff_s_tlp, runoff_s_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, runoff_g_tlp, runoff_g_slp, subsfrac)

  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, tfv_tlp, tfv_slp, subsfrac)

! These variables are averaged inside the radiation module:
! alb_rad, sobs, thbs, pabs

!diagnostics:
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,   t_2m_tlp,   t_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  qv_2m_tlp,  qv_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  td_2m_tlp,  td_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  rh_2m_tlp,  rh_2m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  u_10m_tlp,  u_10m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1,  v_10m_tlp,  v_10m_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, shfl_s_tlp, shfl_s_slp, subsfrac)
  CALL subsmean_arithmetic(i_count, nsubs0, nsubs1, lhfl_s_tlp, lhfl_s_slp, subsfrac)

  DO ic = 1, i_count                                  
    jc = idx_lst_lndp(ic)
    t_g_s      (jc) = t_g_slp      (ic)               
    qv_s_s     (jc) = qv_s_slp     (ic)
    t_s_s      (jc) = t_s_slp      (ic)
    t_snow_s   (jc) = t_snow_slp   (ic)
    w_snow_s   (jc) = w_snow_slp   (ic)
    rho_snow_s (jc) = rho_snow_slp (ic)
    h_snow_s   (jc) = h_snow_slp   (ic)
    freshsnow_s(jc) = freshsnow_slp(ic)
    w_i_s      (jc) = w_i_slp      (ic)
    runoff_s_s (jc) = runoff_s_slp (ic)
    runoff_g_s (jc) = runoff_g_slp (ic)
    tch_s      (jc) = tch_slp      (ic)
    tfv_s      (jc) = tfv_slp      (ic)
    t_2m_s     (jc) = t_2m_slp     (ic)
    qv_2m_s    (jc) = qv_2m_slp    (ic)
    td_2m_s    (jc) = td_2m_slp    (ic)
    rh_2m_s    (jc) = rh_2m_slp    (ic)
    u_10m_s    (jc) = u_10m_slp    (ic)
    v_10m_s    (jc) = v_10m_slp    (ic)
    shfl_s_s   (jc) = shfl_s_slp   (ic)
    lhfl_s_s   (jc) = lhfl_s_slp   (ic)
    DO jk = 0, nlev_soil+1                         !EM order of loops?
      t_so_s    (jc,jk) = t_so_slp    (ic,jk)
    END DO
    DO jk = 1, nlev_soil+1
      w_so_s    (jc,jk) = w_so_slp    (ic,jk)
      w_so_ice_s(jc,jk) = w_so_ice_slp(ic,jk)
    END DO
    DO jk = 0, nlev_snow+1
      t_snow_mult_s  (jc,jk) = t_snow_mult_slp  (ic,jk)
    END DO
    DO jk = 1, nlev_snow+1
      wtot_snow_s    (jc,jk) = wtot_snow_slp    (ic,jk)
      rho_snow_mult_s(jc,jk) = rho_snow_mult_slp(ic,jk)
    END DO
  ENDDO

  END SUBROUTINE subsmean

!-------------------------------------------------------------------------

  SUBROUTINE subsmean_arithmetic(i_count, nsubs0, nsubs1, field, field_s, subsfrac, mask)

  INTEGER , INTENT(IN) :: i_count, nsubs0, nsubs1
  REAL(wp), DIMENSION(i_count, nsubs1), INTENT(IN) :: field, subsfrac
  LOGICAL , DIMENSION(i_count, nsubs1), INTENT(IN), OPTIONAL :: mask
  REAL(wp), DIMENSION(i_count), INTENT(OUT) :: field_s 

  ! Local scalars:
  INTEGER :: ic, ns

  ! Local arrays:
  REAL(wp):: tmp(i_count)       ! total fraction of tiles with mask=true
  LOGICAL :: mask_tmp(i_count)

  IF(PRESENT(mask)) THEN
    DO ic = 1, i_count
      mask_tmp(ic) = .FALSE.
    END DO 
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        IF(mask(ic,ns)) mask_tmp(ic) = .TRUE.
      END DO
    END DO 
    DO ic = 1, i_count
      IF(mask_tmp(ic)) THEN
        field_s(ic) = 0._wp
        tmp(ic)     = 0._wp
      END IF
    END DO
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        IF(mask(ic,ns)) THEN
          field_s(ic) = field_s(ic) + field(ic,ns)*subsfrac(ic,ns)
          tmp(ic)     = tmp(ic) + subsfrac(ic,ns)
        END IF
      END DO
    END DO
    DO ic = 1, i_count
      IF(mask_tmp(ic)) THEN
        field_s(ic) = field_s(ic)/MAX(tmp(ic),1.E-07_wp)
      END IF
    END DO
  ELSE  ! mask is not present
    DO ic = 1, i_count
      field_s(ic) = 0._wp
    END DO 
    DO ns = nsubs0, nsubs1
      DO ic = 1, i_count
        field_s(ic) = field_s(ic) + field(ic,ns)*subsfrac(ic,ns)
      END DO
    END DO
  END IF

  END SUBROUTINE subsmean_arithmetic

!-------------------------------------------------------------------------

  SUBROUTINE subsmean_power4(i_count, nsubs0, nsubs1, field, field_s, subsfrac)

  INTEGER , INTENT(IN) :: i_count, nsubs0, nsubs1
  REAL(wp), DIMENSION(i_count, nsubs1), INTENT(IN) :: field, subsfrac
  REAL(wp), DIMENSION(i_count), INTENT(OUT) :: field_s 

  ! Local scalars:
  INTEGER :: ic, ns

!Without any mask because this kind of averaging (power4) is being done only for t_g

  DO ic = 1, i_count
    field_s(ic) = 0._wp
  END DO
  DO ns = nsubs0, nsubs1
    DO ic = 1, i_count
      field_s(ic) = field_s(ic) + field(ic,ns)**4*subsfrac(ic,ns)
    END DO
  END DO
  DO ic = 1, i_count
    field_s(ic) = SQRT(SQRT(field_s(ic)))
  END DO

  END SUBROUTINE subsmean_power4


  SUBROUTINE diag_snowfrac_tg(istart, iend, z0_lcc, lc_class, t_snow, t_soiltop, w_snow, &
    & rho_snow, freshsnow, sso_sigma, tai, snowfrac, t_g)

    INTEGER, INTENT (IN) :: istart, iend ! start and end-indices of the computation

    INTEGER, INTENT (IN) :: lc_class(:)  ! list of land-cover classes
    REAL(wp), DIMENSION(:), INTENT(IN) :: z0_lcc(:)    ! roughness length
    REAL(wp), DIMENSION(:), INTENT(IN) :: t_snow, t_soiltop, w_snow, rho_snow, &
      freshsnow, sso_sigma, tai

    REAL(wp), DIMENSION(:), INTENT(INOUT) :: snowfrac, t_g

    INTEGER  :: ic
    REAL(wp) :: h_snow, snowdepth_fac, sso_fac, z0_fac, z0_limit, lc_limit

    IF (idiag_snowfrac == 1) THEN
      DO ic = istart, iend
        snowfrac(ic) = MIN(1.0_wp, w_snow(ic)/cf_snow)
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    ELSE
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.04_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(20._wp*freshsnow(ic)+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          z0_fac   = MAX(1._wp,SQRT(20._wp*z0_lcc(MAX(1,lc_class(ic)))))
          z0_limit = MIN(1._wp,SQRT(SQRT(1.5_wp/z0_fac)))
          lc_limit = MIN(1._wp,1._wp/SQRT(MAX(0.1_wp,tai(ic))))
          snowfrac(ic) = MIN(lc_limit,z0_limit,snowdepth_fac/z0_fac)
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO

    ENDIF

  END SUBROUTINE diag_snowfrac_tg


END MODULE mo_nwp_sfc_utils

