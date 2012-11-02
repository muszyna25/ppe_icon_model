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
  USE mo_exception,           ONLY: message, message_text
  USE mo_exception,           ONLY: finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, zml_soil
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,     ONLY: nproma
  USe mo_extpar_config,       ONLY: itopo
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
    &                               lseaice, llake, lmulti_snow, idiag_snowfrac, ntiles_lnd, &
    &                               lsnowtile, isub_water, isub_seaice
  USE mo_soil_ml,             ONLY: terra_multlay_init
  USE mo_seaice_nwp,          ONLY: seaice_init_nwp, hice_min
  USE mo_phyparam_soil,       ONLY: cf_snow     ! soil and vegetation parameters for TILES
  USE mo_physical_constants,  ONLY: tmelt, rdocp => rd_o_cpd  ! r_d / cp_d
  USE mo_seaice_nwp,          ONLY: frsi_min
  USE mo_sync,                ONLY: global_sum_array
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
  PUBLIC :: update_idx_lists_lnd
  PUBLIC :: update_idx_lists_sea
  PUBLIC :: init_snowtile_lists
  PUBLIC :: init_seaice_lists

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init surface model TERRA and seaice model
  !!
  !! Init surface model TERRA and seaice model
  !!
  !! @par Revision History
  !! Initial revision by Ekaterina Machulskaya, DWD (2011-07-??)
  !! Modification by Daniel Reinert, DWD (2011-07-29)
  !! - initialize climatological layer t_so(nlev_soil+2)
  !! Modification by Daniel Reinert, DWD (2012-08-31)
  !! - initialize seaice model
  !!
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now, &
    &                          p_prog_lnd_new, p_prog_wtr_now,    &
    &                          p_prog_wtr_new, p_lnd_diag )
 
                             

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_wtr_prog)     , INTENT(INOUT) :: p_prog_wtr_now, p_prog_wtr_new
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jc,jb,isubs,jk


    REAL(wp) :: t_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_snow_mult_now_t(nproma, 1:nlev_snow+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_g_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_new_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_mult_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_now_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_new_t(nproma, 1:nlev_soil+2, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wliq_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wtot_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: dzh_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)

    INTEGER  ::          soiltyp_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) ::          rootdp_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) ::          tai_t     (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) ::          freshsnow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) ::          snowfrac_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) ::          sso_sigma_t(nproma,  p_patch%nblks_c, ntiles_total)
    INTEGER  ::          lc_class_t (nproma,  p_patch%nblks_c, ntiles_total)

    ! local fields for sea ice model
    !
    REAL(wp) :: frsi     (nproma)   ! sea ice fraction
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time 
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time 
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level
    INTEGER  :: icount_ice          ! total number of sea-ice points per block

    INTEGER  :: i_count, ic, i_count_snow, isubs_snow
    REAL(wp), PARAMETER :: small = 1.E-06_wp
    REAL(wp) :: t_g_s(nproma)

  !-------------------------------------------------------------------------


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,i_count,i_count_snow,icount_ice,    &
!$OMP            ic,jk,isubs_snow,t_g_s,frsi,tice_now,hice_now,tsnow_now,hsnow_now), &
!$OMP            SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF (itopo == 1) THEN
        DO isubs = 1, ntiles_total
          DO jc = i_startidx, i_endidx

            ! initialize climatological layer (deepest layer of t_so)
            p_prog_lnd_now%t_so_t(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)
            p_prog_lnd_new%t_so_t(jc,nlev_soil+2,jb,isubs) = ext_data%atm%t_cl(jc,jb)

            p_prog_lnd_now%t_g_t(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,isubs) = p_prog_lnd_now%t_g(jc,jb)

          END DO 
        END DO
      ENDIF

!      IF (ext_data%atm%lp_count(jb) == 0) CYCLE ! skip loop if there is no land point

!---------- Copy input fields for each tile

!----------------------------------
      DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          t_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%t_snow_t(jc,jb,isubs) 
          t_s_now_t(ic,jb,isubs)             =  p_prog_lnd_now%t_s_t(jc,jb,isubs)   
          t_s_new_t(ic,jb,isubs)             =  p_prog_lnd_new%t_s_t(jc,jb,isubs)   
          w_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%w_snow_t(jc,jb,isubs)  
          rho_snow_now_t(ic,jb,isubs)        =  p_prog_lnd_now%rho_snow_t(jc,jb,isubs)

          sso_sigma_t(ic,jb,isubs)           = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic,jb,isubs)            = ext_data%atm%lc_class_t(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           = p_lnd_diag%freshsnow_t(jc,jb,isubs)

          soiltyp_t(ic,jb,isubs)             =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          rootdp_t(ic,jb,isubs)              =  ext_data%atm%rootdp_t(jc,jb,isubs)
          tai_t(ic,jb,isubs)                 =  ext_data%atm%tai_t(jc,jb,isubs)
        ENDDO


        IMSNOWI: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              t_snow_mult_now_t(ic,jk,jb,isubs)   =  p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) 
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
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
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            t_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) 
            t_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            w_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) 
            w_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) 
            w_so_ice_new_t(ic,jk,jb,isubs)      =  p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) 
          ENDDO
        ENDDO

        ! Preliminary diagnosis of snow-cover fraction for initialization of splitted tile index list

! Remark(GZ): this directive is needed because of a NEC compiler bug - OpenMP parallelization causes
! a segmentation fault otherwise
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
          &  sso_sigma = sso_sigma_t       (:,jb,isubs), & ! sso stdev
          &  tai       = tai_t             (:,jb,isubs), & ! effective leaf area index
          &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
          &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp   

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
          p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
        ENDDO

      ENDDO

      ! create index lists for snow tiles (first call)
      IF(lsnowtile) THEN      ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd
          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                   lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                   idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                   gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                   idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                   gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                   lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                   partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                   partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                   snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                   snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

          ! Set w_snow to zero for grid points having a snow-tile counterpart
          DO ic = 1, ext_data%atm%lp_count_t(jb,isubs)
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) > 0 ) p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
          END DO

        END DO
      END IF


      DO isubs = 1,ntiles_total

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        CALL terra_multlay_init(                                  &
        &  ie=nproma,                                             & ! array dimensions
        &  istartpar=1, iendpar= i_count,                         & ! optional start/end indicies
        &  ke_soil=nlev_soil, ke_snow=nlev_snow                 , &
        &  czmls=zml_soil                                       , & ! processing soil level structure 
        &  soiltyp_subs      = soiltyp_t(:,jb,isubs)            , & ! type of the soil (keys 0-9)  --
        &  rootdp            = rootdp_t(:,jb,isubs)             , & ! depth of the roots                ( m  )
        &  t_snow_now        = t_snow_now_t(:,jb,isubs)         , & ! temperature of the snow-surface   (  K  )
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
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs)       ) ! layer thickness between half levels in snow  (  m  )

! Remark(GZ): this directive is needed because of a NEC compiler bug - OpenMP parallelization causes
! a segmentation fault otherwise
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
          &  sso_sigma = sso_sigma_t       (:,jb,isubs), & ! sso stdev
          &  tai       = tai_t             (:,jb,isubs), & ! effective leaf area index
          &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
          &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp

!  Recover fields from index list
!
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          p_prog_lnd_now%t_snow_t(jc,jb,isubs)   = t_snow_now_t(ic,jb,isubs)
          p_prog_lnd_now%t_s_t(jc,jb,isubs)      = t_s_now_t(ic,jb,isubs)  
          p_prog_lnd_new%t_s_t(jc,jb,isubs)      = t_s_new_t(ic,jb,isubs) 
          p_prog_lnd_now%w_snow_t(jc,jb,isubs)   = w_snow_now_t(ic,jb,isubs) 
          p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = rho_snow_now_t(ic,jb,isubs)
          p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
          p_lnd_diag%snowfrac_t(jc,jb,isubs)     = snowfrac_t(ic,jb,isubs)
          p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          p_prog_lnd_new%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
        ENDDO

        IF (lsnowtile .AND. isubs > ntiles_lnd) THEN ! copy snowfrac_t to snow-free tile
!CDIR NODEP,VOVERTAKE,VOB                            ! (needed for index list computation)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs-ntiles_lnd) = p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)
          ENDDO
        ENDIF

        IMSNOWO: IF(lmulti_snow) THEN

!CDIR UNROLL=nlsnow+1
          DO jk=1,nlev_snow+1
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) =  t_snow_mult_now_t(ic,jk,jb,isubs)   
            ENDDO
          ENDDO

!CDIR UNROLL=nlsnow
          DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
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
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) = t_so_now_t(ic,jk,jb,isubs)          
            p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)          
          ENDDO
        ENDDO

!CDIR UNROLL=nlsoil+1
        DO jk=1,nlev_soil+1
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) = w_so_now_t(ic,jk,jb,isubs)        
            p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)        
            p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO
      END DO ! isubs

      i_count = ext_data%atm%lp_count(jb)

      ! Rediagnose t_g - needed for initial calculation of transfer coefficients

      IF (ntiles_total == 1) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          p_prog_lnd_now%t_g(jc,jb) = p_prog_lnd_now%t_g_t(jc,jb,1)
          p_prog_lnd_new%t_g(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
        ENDDO
      ELSE ! aggregate fields over tiles
        t_g_s(:)  =  0._wp
        DO isubs = 1,ntiles_total+ntiles_water
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            t_g_s(jc) = t_g_s(jc) + ext_data%atm%frac_t(jc,jb,isubs)* &
              p_prog_lnd_now%t_g_t(jc,jb,isubs)**4
          ENDDO
        ENDDO
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          p_prog_lnd_now%t_g(jc,jb)  = SQRT(SQRT(t_g_s(jc)))
          p_prog_lnd_new%t_g(jc,jb)  = p_prog_lnd_now%t_g(jc,jb)
        ENDDO
      END IF

      IF(lsnowtile) THEN      ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd 

          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                   lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                   idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                   gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                   idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                   gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                   lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                   partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                   partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                   snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                   snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

          i_count = ext_data%atm%gp_count_t(jb,isubs)
          i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)

          DO ic = 1, i_count_snow
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

            ! snow depth per surface unit -> snow depth per fraction
            p_lnd_diag%w_snow_eff_t(jc,jb,isubs_snow) = &
              p_prog_lnd_now%w_snow_t(jc,jb,isubs_snow)/MAX(p_lnd_diag%snowfrac_lc_t(jc,jb,isubs_snow),small)

            ! reset field for actual snow-cover for grid points / land-cover classes for which there
            ! are seperate snow-free and snow-covered tiles 
            p_lnd_diag%snowfrac_t(jc,jb,isubs)   = 0._wp
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
            p_prog_lnd_now%t_snow_t(jc,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)    = p_prog_lnd_now%t_s_t(jc,jb,isubs)

            ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
            ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
            ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
            ! time steps when associating traces of snow with a snow-cover fraction of 1)
            p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) = MIN(1._wp,p_lnd_diag%h_snow_t(jc,jb,isubs_snow)/0.01_wp)

            ! Rediagnose t_g according to the modified snow-cover fraction
            p_prog_lnd_now%t_g_t(jc,jb,isubs_snow) =  &
              p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) * p_prog_lnd_now%t_snow_t(jc,jb,isubs_snow) + &
              (1._wp-p_lnd_diag%snowfrac_t(jc,jb,isubs_snow))*p_prog_lnd_now%t_s_t(jc,jb,isubs_snow)

             IF (lmulti_snow) THEN
               p_prog_lnd_now%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
             ENDIF

          END DO

          IF (lmulti_snow) THEN
!CDIR UNROLL=nlsnow
            DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
              DO ic = 1, i_count_snow
                jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%dzh_snow_t (jc,jk,jb,isubs) = 0._wp
              ENDDO
            ENDDO
          ENDIF

        END DO
      END IF


      !
      ! Init sea ice parameterization
      ! 
      IF ( lseaice ) THEN


        icount_ice = ext_data%atm%spi_count(jb) ! number of sea-ice points in block jb

!IF (icount_ice > 0) THEN
!write(0,*) "nwp_surface_init:ext_data%atm%spi_count(jb), jb: ", ext_data%atm%spi_count(jb), jb
!ENDIF

        DO ic = 1, icount_ice

          jc = ext_data%atm%idx_lst_spi(ic,jb)

          frsi     (ic) = p_lnd_diag%fr_seaice(jc,jb)
          tice_now (ic) = p_prog_wtr_now%t_ice(jc,jb)
          hice_now (ic) = p_prog_wtr_now%h_ice(jc,jb)
          tsnow_now(ic) = p_prog_wtr_now%t_snow_si(jc,jb)
          hsnow_now(ic) = p_prog_wtr_now%h_snow_si(jc,jb)
        ENDDO  ! jc


        CALL seaice_init_nwp ( icount_ice, frsi,                         & ! in
          &                    tice_now, hice_now, tsnow_now, hsnow_now, & ! inout
          &                    tice_new, hice_new, tsnow_new, hsnow_new  ) ! inout


        !  Recover fields from index list
        !
        DO ic = 1, icount_ice
          jc = ext_data%atm%idx_lst_spi(ic,jb)

          ! fields at time level now may have changed, potentially!
          p_prog_wtr_now%t_ice(jc,jb)    = tice_now(ic)
          p_prog_wtr_now%h_ice(jc,jb)    = hice_now(ic)
IF (p_prog_wtr_now%t_ice(jc,jb) < 1._wp) THEN
write(0,*) "after init tice_now(ic), hice_now(ic), lat/lon: ", &
  &  p_prog_wtr_now%t_ice(jc,jb), p_prog_wtr_now%h_ice(jc,jb), &
  &  p_patch%cells%center(jc,jb)%lat, p_patch%cells%center(jc,jb)%lon
ENDIF


          p_prog_wtr_now%t_snow_si(jc,jb)= tsnow_now(ic)
          p_prog_wtr_now%h_snow_si(jc,jb)= hsnow_now(ic)

          p_prog_wtr_new%t_ice(jc,jb)    = tice_new(ic)
          p_prog_wtr_new%h_ice(jc,jb)    = hice_new(ic)
          p_prog_wtr_new%t_snow_si(jc,jb)= tsnow_new(ic)
          p_prog_wtr_new%h_snow_si(jc,jb)= hsnow_new(ic)

        ENDDO  ! ic

      ENDIF  ! lseaice

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

      IF (ntiles_total == 1) THEN  ! just copy prognostic variables from tile 1 to diagnostic aggregated variable
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

        DO isubs = 1,ntiles_total
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp(ic,jb)
            tilefrac = ext_data%atm%frac_t(jc,jb,isubs)
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
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)
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
                tilefrac = ext_data%atm%frac_t(jc,jb,isubs)
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

     i_count = ext_data%atm%fp_count(jb) ! third step: copy variables also for lake points 

!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_fp(ic,jb)
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
          jc = ext_data%atm%idx_lst_fp(ic,jb)
          lnd_diag%t_so(jc,jk,jb)      = lnd_prog%t_so_t(jc,jk,jb,1)
          lnd_diag%w_so(jc,jk,jb)      = lnd_prog%w_so_t(jc,jk,jb,1)
          lnd_diag%w_so_ice(jc,jk,jb)  = lnd_prog%w_so_ice_t(jc,jk,jb,1)
        ENDDO
      ENDDO

      IF (lmulti_snow) THEN
        DO jk=1,nlev_snow
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_fp(ic,jb)
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



  !! Driver routine to (re-)initialize the snowtile index lists in the case of a restart
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD (2012-08-03)
  !!
  SUBROUTINE init_snowtile_lists(p_patch, ext_data, p_lnd_diag)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jb,isubs,i_count,isubs_snow


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,isubs,i_count,isubs_snow), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO isubs = 1, ntiles_lnd
        isubs_snow = isubs + ntiles_lnd

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                 lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                 idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                 gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                 idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                 gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                 lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                 partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                 partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                 snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                 snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                 snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

      ENDDO

    ENDDO    
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_snowtile_lists



  !-------------------------------------------------------------------------
  !>
  !! Initialize the sea-ice and open water index lists for restart and 
  !! non-restart runs. Seaice and open water points sre distinguished 
  !! based on fr_seaice which is provided by analysis.
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD (2012-08-03)
  !!
  SUBROUTINE init_seaice_lists(p_patch, ext_data, p_lnd_diag, lseaice)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch        !< grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag)     , INTENT(IN)    :: p_lnd_diag     !< diag vars for sfc
    LOGICAL              , INTENT(IN)    :: lseaice        !< seaice model on/off
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc
    INTEGER :: jg
    INTEGER :: i_count_sea, i_count_ice, i_count_water
    INTEGER :: npoints_ice, npoints_wtr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_seaice'
!-------------------------------------------------------------------------


    ! patch ID
    jg = p_patch%id

    i_nchdom = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    IF (lseaice) THEN

    ! generate sea-ice and open-water index list
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,i_count_sea,i_count_ice,i_count_water), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk



       ! Init sub-index lists for sea points. We distinguish between sea-water 
       ! (i.e. ice free) points and sea-ice points. diag_lnd%fr_seaice is used 
       ! to distinguish between sea-water and sea-ice points.
       ! Needed by sea-ice model
       !
       i_count_sea   = ext_data%atm%sp_count(jb)
       i_count_ice   = 0
       i_count_water = 0

       IF (i_count_sea == 0) CYCLE ! skip loop if the index list for the given block is empty

!CDIR NODEP,VOVERTAKE,VOB
       DO ic = 1, i_count_sea

         jc = ext_data%atm%idx_lst_sp(ic,jb)


         ! set sea-ice area fraction (static)
         ! simply copy from water tile. time-dependent fraction will be set lateron
         ext_data%atm%lc_frac_t(jc,jb,isub_seaice)= ext_data%atm%lc_frac_t(jc,jb,isub_water)

         !
         ! seaice point
         !
         IF ( p_lnd_diag%fr_seaice(jc,jb) >= frsi_min ) THEN
           i_count_ice = i_count_ice + 1
           ext_data%atm%idx_lst_spi(i_count_ice,jb) = jc
           ext_data%atm%spi_count(jb)               = i_count_ice
           ! set land-cover class
           ext_data%atm%lc_class_t(jc,jb,isub_seaice)= ext_data%atm%i_lc_snow_ice
         ENDIF

         !
         ! water point: all sea points with fr_ice < 1
         !
         IF ( p_lnd_diag%fr_seaice(jc,jb) < 1._wp ) THEN
           i_count_water = i_count_water + 1
           ext_data%atm%idx_lst_spw(i_count_water,jb) = jc
           ext_data%atm%spw_count(jb)                 = i_count_water
         ENDIF

         ! Initialize frac_t for seaice and update frac_t for water tile
         ext_data%atm%frac_t(jc,jb,isub_water)  = ext_data%atm%lc_frac_t(jc,jb,isub_water)  &
           &                                    * (1._wp - p_lnd_diag%fr_seaice(jc,jb))
         ext_data%atm%frac_t(jc,jb,isub_seaice) = ext_data%atm%lc_frac_t(jc,jb,isub_seaice) &
           &                                    * p_lnd_diag%fr_seaice(jc,jb)

       ENDDO  ! ic 

    ENDDO  ! jb    
!$OMP END DO
!$OMP END PARALLEL


      ! Some diagnostics
      npoints_ice = SUM(ext_data%atm%spi_count(i_startblk:i_endblk))
      npoints_ice = global_sum_array(npoints_ice)
      npoints_wtr = SUM(ext_data%atm%spw_count(i_startblk:i_endblk))
      npoints_wtr = global_sum_array(npoints_wtr)
      WRITE(message_text,'(a,i3,a,i10)') 'Number of seaice points in domain',jg, &
        &  ':',npoints_ice
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,a,i10)') 'Number of water points in domain',jg, &
        &  ':',npoints_wtr
      CALL message('', TRIM(message_text))


    ELSE   ! seaice model switched off

      ! copy sea points index list to open-water index list
      !
      ext_data%atm%spw_count(i_startblk:i_endblk)     =     &
        &                     ext_data%atm%sp_count(i_startblk:i_endblk)
      ext_data%atm%idx_lst_spw(:,i_startblk:i_endblk) =     &
        &                     ext_data%atm%idx_lst_sp(:,i_startblk:i_endblk)

      ext_data%atm%spi_count(i_startblk:i_endblk) = 0
 
    ENDIF  ! lseaice


  END SUBROUTINE init_seaice_lists

  !-------------------------------------------------------------------------



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
    ELSE IF (idiag_snowfrac == 2) THEN
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
    ELSE    ! idiag_snowfrac = 3 - similar to option 2, but different tuning
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.025_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(17.5_wp*freshsnow(ic)+5._wp+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          z0_fac   = MAX(1._wp,SQRT(12.5_wp*z0_lcc(MAX(1,lc_class(ic)))))
          z0_limit = MIN(1._wp,SQRT(SQRT(2.5_wp/z0_fac)))
          lc_limit = MIN(1._wp,1.75_wp/SQRT(MAX(0.1_wp,tai(ic))))
          snowfrac(ic) = MIN(lc_limit,z0_limit,snowdepth_fac/z0_fac)
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    ENDIF

  END SUBROUTINE diag_snowfrac_tg




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Updating of dynamic index lists
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic 
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD (2012-07-01)
  !!
  SUBROUTINE update_idx_lists_lnd (idx_lst_lp, lp_count, idx_lst, gp_count, idx_lst_snow, &
    &             gp_count_snow, lc_frac, partial_frac, partial_frac_snow, snowtile_flag, &
    &             snowtile_flag_snow, snowfrac)


    INTEGER ,    INTENT(   IN) ::  &   !< static list of all land points of a tile index
      &  idx_lst_lp(:), lp_count       !< and corresponding grid point counts

 
    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-free or mixed land points
      &  idx_lst(:), gp_count          !< of a tile index
                                       !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-covered points for land-cover 
      &  idx_lst_snow(:), gp_count_snow !< classes eligible for separate treatment
                                       !< and corresponding grid point counts

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-free or mixed points
      &  snowtile_flag(:)              !< -1: no separation between snow tile and snow-free tile
                                       !< inactive
                                       !< active
                                       !< newly activated; initialization from corresponding tile required

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-covered points
      &  snowtile_flag_snow(:)         !< -1: no separation between snow tile and snow-free tile
                                       !< inactive
                                       !< active
                                       !< newly activated; initialization from corresponding tile required

    REAL(wp),    INTENT(   IN) ::  &   !< area fraction of land-cover class
      &  lc_frac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-cover fraction
      &  snowfrac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-free and snow-covered sub-tile area fraction
      &  partial_frac(:), partial_frac_snow(:)


    ! Local variables
    INTEGER  :: ic, jc, icount, icount_snow
    REAL(wp) :: eps = 1.e-6_wp

    !-------------------------------------------------------------------------

    icount = 0
    icount_snow = 0

!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, lp_count
      jc = idx_lst_lp(ic)

      IF (snowtile_flag(jc) == -1) THEN
        icount = icount + 1
        idx_lst(icount) = jc
        partial_frac(jc) = lc_frac(jc)
      ELSE
        ! Reset snowfrac to 0/1 in case of very small deviations (just to be safe)
        IF (snowfrac(jc) < eps) snowfrac(jc) = 0._wp
        IF (1._wp - snowfrac(jc) < eps) snowfrac(jc) = 1._wp
#ifdef __SX__
      ENDIF
    ENDDO

!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, lp_count
      jc = idx_lst_lp(ic)

      IF (snowtile_flag(jc) /= -1) THEN
#endif
        IF (snowfrac(jc) > 0._wp) THEN ! snow tile is active
          icount_snow = icount_snow + 1
          idx_lst_snow(icount_snow) = jc
          partial_frac_snow(jc) = lc_frac(jc)*snowfrac(jc)
          IF (snowtile_flag_snow(jc) == 0) THEN
            snowtile_flag_snow(jc) = 2 ! newly activated, initialization needed
          ELSE
            snowtile_flag_snow(jc) = 1
          ENDIF
        ELSE
          snowtile_flag_snow(jc) = 0
          partial_frac_snow(jc) = 0._wp
        ENDIF
        IF (snowfrac(jc) < 1._wp) THEN ! snow-free tile is active
          icount = icount + 1
          idx_lst(icount) = jc
          partial_frac(jc) = lc_frac(jc)*(1._wp-snowfrac(jc))
          IF (snowtile_flag(jc) == 0) THEN
            snowtile_flag(jc) = 2 ! newly activated, initialization needed
          ELSE
            snowtile_flag(jc) = 1
          ENDIF
        ELSE
          snowtile_flag(jc) = 0
          partial_frac(jc) = 0._wp
        ENDIF
      ENDIF
    ENDDO

    gp_count = icount
    gp_count_snow = icount_snow

  END SUBROUTINE update_idx_lists_lnd




  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Updating of dynamic index lists for sea points
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic sea-water and sea-ice point index list, which are sub-index lists 
  !!   of the static sea point index list idx_lst_sp 
  !!
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2012-08-31)
  !!
  SUBROUTINE update_idx_lists_sea (hice_n, idx_lst_spw, spw_count, idx_lst_spi,    &
    &                              spi_count, partial_frac_ice, partial_frac_water,&
    &                              fr_seaice, t_g_t_n )


    REAL(wp),    INTENT(IN)    ::  &   !< sea ice depth at new time level  [m]
      &  hice_n(:)                     !< dim: (nproma)

    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic sea water point index list 
      &  idx_lst_spw(:), spw_count     !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic sea ice point index list 
      &  idx_lst_spi(:), spi_count     !< and corresponding grid point counts

    REAL(wp),    INTENT(INOUT) ::  &   !< ice-covered and ice-free ocean area fraction
      &  partial_frac_ice(:), partial_frac_water(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< seaice fraction
      &  fr_seaice(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< temperature of water tile (new)
      &  t_g_t_n(:)

    ! Local variables
    INTEGER, DIMENSION(SIZE(idx_lst_spi,1)) :: &
      &   idx_lst_spi_old       !< seaice index list local copy               
    INTEGER  :: ic, jc          !< loop indices
    INTEGER  :: spi_count_old   !< current seaice grid point count

    ! should be moved to mo_physical_constants
    REAL (wp), PARAMETER ::                             &
      &   tf_salt      = 271.45_wp      !< salt-water freezing point [K] 

    !-------------------------------------------------------------------------


    ! save old seaice index list and grid point count
    idx_lst_spi_old(:) = idx_lst_spi(:)
    spi_count_old      = spi_count

    ! re-initialize seaice index list and grid point count
    idx_lst_spi(:) = 0
    spi_count      = 0



    !
    ! update index list for sea-ice and open water
    !
    ! Since the current seaice model does not allow for new seaice points to be 
    ! created during model integration, the number of seaice points can only 
    ! decrease with time. I.e. seaice points may be converted into water points, 
    ! but not vice versa. 
    !
    ! Loop over old seaice-points, only
!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, spi_count_old
      jc = idx_lst_spi_old(ic)

      IF ( hice_n(jc) >= hice_min )  THEN ! still seaice point
        spi_count = spi_count + 1
        idx_lst_spi(spi_count) = jc
        ! seaice fraction remains unchanged
      ELSE                          ! seaice point has turned into water point
        ! Check whether we need to initialize a new water tile, or whether a water tile 
        ! already exists for the given point:
        IF ( fr_seaice(jc) == 1._wp ) THEN   ! water tile does not exist for given point
          ! add new water tile to water-points index list and initialize
          spw_count = spw_count + 1
          idx_lst_spw(spw_count) = jc
          ! Initialize temperature of water tile with salt water freezing point
          t_g_t_n(jc) = tf_salt   ! if the SST analysis contains a meaningful water 
                                  ! temperature for this point, one may also take 
                                  ! the latter
          ! Initialize saturation specific humidity for new water tile
          ! TO BE CODED
          !qv_s_t = spec_humi(sat_pres_water(t_g_t_n(jc)),p_diag%pres_sfc(jc,jb))
        ENDIF

        ! re-set partial fractions of water and seaice
        !
        ! new water area fraction is the sum of the current water and seaice area fractions
        partial_frac_water(jc)= partial_frac_water(jc) + partial_frac_ice(jc)
        ! since seaice melted away, the seaice fraction is re-set to 0
        fr_seaice(jc)         = 0._wp
        partial_frac_ice(jc)  = 0._wp

!DR Debug output: will be removed lateron
        write(0,*) "ice->water: partial_frac_water(jc), partial_frac_ice(jc): ", &
          & partial_frac_water(jc), partial_frac_ice(jc), jc
!DR END DEBUG

      ENDIF

    ENDDO  ! ic


  END SUBROUTINE update_idx_lists_sea

END MODULE mo_nwp_sfc_utils

