#ifndef __NO_ICON_OCEAN__

      MODULE mo_hamocc_output

! icon specific routines for output and restart

      USE mo_control_bgc,          ONLY: dtb


      USE mo_exception, ONLY      : message, finish

      USE mo_model_domain,   ONLY: t_patch_3D, t_patch

      USE mo_ocean_nml,  ONLY:  n_zlev

      USE mo_kind,     ONLY: wp

      USE mo_impl_constants,      ONLY: max_char_length

      USE mo_var_list,            ONLY: add_var,                  &
    &                               new_var_list,             &
    &                               delete_var_list,          &
    &                               default_var_list_settings,&
    &                               add_ref

      USE mo_grid_config,         ONLY: n_dom


      USE mo_hamocc_types,       ONLY: t_hamocc_diag, t_hamocc_state, &
    &                                  t_hamocc_sed, t_hamocc_tend,   &
    &                                  t_hamocc_acc, t_hamocc_monitor                    
   
      USE mo_cdi_constants

      USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, &
   &                                    DATATYPE_PACK16, DATATYPE_INT8, &
   &                                    GRID_LONLAT

      USE mo_cf_convention

      USE mo_linked_list,            ONLY: t_var_list


      USE mo_grib2,               ONLY:  grib2_var
       
      USE mo_parallel_config,     ONLY: nproma

      USE mo_hamocc_nml,         ONLY: io_stdo_bgc

      USE mo_var_metadata,       ONLY: groups, post_op

      USE mo_var_metadata_types, ONLY: POST_OP_SCALE

      USE mo_sedmnt,          ONLY: ks
  
      USE mo_bgc_constants,    ONLY: kilo, s2year, n2tgn, c2gtc

     

      IMPLICIT NONE

      PUBLIC

      TYPE(t_var_list)                              :: hamocc_default_list ! for output
      TYPE(t_var_list)                              :: hamocc_restart_list ! for hi, co3
      TYPE(t_var_list)                              :: hamocc_tendency_list ! for NPP etc
      TYPE(t_var_list)                              :: hamocc_acc_list ! for NPP etc
      TYPE(t_var_list)                              :: hamocc_sediment_list ! for sediment outout

      CONTAINS

!================================================================================== 
    SUBROUTINE construct_hamocc_state( patch_2d, hamocc_state )
    
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d(n_dom)
    TYPE(t_hamocc_state), TARGET :: hamocc_state!(n_dom)
    
    ! local variables
    INTEGER :: jg
    
    INTEGER :: i_status, jp, prlength ! local prognostic array length
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_state'
    
    CALL message(TRIM(routine), 'start to construct hamocc state' )
    
    !create state array for each domain
    DO jg = 1, n_dom
      
      !CALL construct_hamocc_diag(patch_2d(jg), hamocc_state(jg)%p_diag)
      !CALL construct_hamocc_tend(patch_2d(jg), hamocc_state(jg)%p_tend)
      !CALL construct_hamocc_acc(patch_2d(jg), hamocc_state(jg)%p_acc)
      !CALL construct_hamocc_sed(patch_2d(jg), hamocc_state(jg)%p_sed)
      CALL construct_hamocc_diag(patch_2d(jg), hamocc_state%p_diag)
    CALL message(TRIM(routine), 'start to construct hamocc state: tend' )
      CALL construct_hamocc_tend(patch_2d(jg), hamocc_state%p_tend)
      CALL construct_hamocc_moni(patch_2d(jg), hamocc_state%p_tend%monitor)
   
    CALL message(TRIM(routine), 'start to construct hamocc state: acc' )
      CALL construct_hamocc_acc(patch_2d(jg), hamocc_state%p_acc)
    CALL message(TRIM(routine), 'start to construct hamocc state: sed' )
      CALL construct_hamocc_sed(patch_2d(jg), hamocc_state%p_sed)
      
      CALL message(TRIM(routine),'construction of hamocc state finished')
      
    END DO
    
  END SUBROUTINE 

!================================================================================== 
  SUBROUTINE construct_hamocc_diag(patch_2d, hamocc_state_diag)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hamocc_diag), INTENT(inout)         :: hamocc_state_diag
    
    ! local variables
    
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_diag'

     ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v



    !-----DIAG W/O restart-----------------------------------------------------------------
    ! for tracers restart is handled by ICON
    CALL message(TRIM(routine), 'start to construct diagnostic hamocc state')
      CALL add_var(hamocc_default_list, 'phy',hamocc_state_diag%phy,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('phy','kmolP m-3','phytoplankton concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

      CALL add_var(hamocc_default_list, 'zoo',hamocc_state_diag%zoo,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('zoo','kmolP m-3','zooplankton concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)
    
     CALL add_var(hamocc_default_list, 'doc',hamocc_state_diag%doc,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('doc','kmolP m-3','DOC concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)
    
     CALL add_var(hamocc_default_list, 'det',hamocc_state_diag%det,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('det','kmolP m-3','POC concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_default_list, 'fdust',hamocc_state_diag%dust,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('fdust','kmol m-3','Dust concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)


      CALL add_var(hamocc_default_list, 'dic',hamocc_state_diag%dic,    &
       & grid_unstructured_cell, za_depth_below_sea,&
       & t_cf_var('dic','kmolP m-3','DIC concentration', DATATYPE_FLT32), &
       & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
       & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
       & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'cya',hamocc_state_diag%cya,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cya','kmolP m-3','cyanobacteria concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)


     CALL add_var(hamocc_default_list, 'alk',hamocc_state_diag%alk,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('alk','kmol m-3','alkalinity', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)
     
     CALL add_var(hamocc_default_list, 'no3',hamocc_state_diag%no3,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('no3','kmol N m-3','Nitrate concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'po4',hamocc_state_diag%po4,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('po4','kmolP m-3','phosphate concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'n2',hamocc_state_diag%n2,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('n2','kmol N2 m-3','gaseous nitrogen concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'o2',hamocc_state_diag%o2,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('o2','kmolP O2 m-3','oxygen concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'si',hamocc_state_diag%si,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('si','kmolP m-3','silicate concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

      CALL add_var(hamocc_default_list, 'iron',hamocc_state_diag%iron,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('iron','kmolFe m-3','dissolved iron concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

      CALL add_var(hamocc_default_list, 'dms',hamocc_state_diag%dms,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms','kmol S m-3','dimethylsulfide concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 218, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

      CALL add_var(hamocc_default_list, 'h2s',hamocc_state_diag%h2s,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2s','kmol m-3','H2S alkalinity concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'n2o',hamocc_state_diag%n2o,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('n2o','kmol  m-3','N2O concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'calc',hamocc_state_diag%calc,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('calc','kmol C m-3','calcium carbonate shells', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

     CALL add_var(hamocc_default_list, 'opal',hamocc_state_diag%opal,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('opal','kmol Si m-3','opal shells', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL message(TRIM(routine), 'construct diagnostic hamocc end base')

! !-----------------DIAG WITH RESTART---------------------------------------------------------------
     CALL add_var(hamocc_restart_list, 'hion',hamocc_state_diag%hi,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('hion','kmol m-3','hydrogen ion concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

     CALL add_var(hamocc_restart_list, 'co3',hamocc_state_diag%co3,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('co3','kmol m-3','carbonate ion concentration', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_BASE"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL message(TRIM(routine), 'construct diagnostic hamocc end')
  END SUBROUTINE 
!================================================================================== 
  SUBROUTINE construct_hamocc_acc(patch_2d, hamocc_state_acc)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hamocc_acc), INTENT(inout)         :: hamocc_state_acc
    
    ! local variables
    
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_acc'

     ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v



    !-----DIAG W/O restart-----------------------------------------------------------------
    ! for tracers restart is handled by ICON
    CALL message(TRIM(routine), 'start to construct hamocc acc state')
   
    CALL add_var(hamocc_acc_list, 'NPP',hamocc_state_acc%npp,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('NPP','kmolP m-3 s-1','net primary production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 200, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)
  
   CALL add_var(hamocc_acc_list, 'BACFRA',hamocc_state_acc%bacfra,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('BACFRA','kmolP m-3 s-1','bacterial decomposition of DOC', DATATYPE_FLT32), &
      & grib2_var(255, 255, 293, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)


    CALL add_var(hamocc_acc_list, 'SRED',hamocc_state_acc%remins,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('SRED','kmolP m-3 s-1','sulfate reduction', DATATYPE_FLT32), &
      & grib2_var(255, 255, 296, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)
      
    CALL add_var(hamocc_acc_list, 'REMIN',hamocc_state_acc%remina,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('REMIN','kmolP m-3 s-1','aerob detritus remineralization', DATATYPE_FLT32), &
      & grib2_var(255, 255, 294, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'DENIT',hamocc_state_acc%reminn,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('DENIT','kmolP m-3 s-1','denitrification', DATATYPE_FLT32), &
      & grib2_var(255, 255, 159, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'nfix',hamocc_state_acc%nfix,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('nfix','kmol N2 m-3 s-1','n fixation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 166, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'dms_prod',hamocc_state_acc%dmsprod,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_prod','kmol S m-3 s-1','DMS production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 69, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'h2s_prod',hamocc_state_acc%h2sprod,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2s_prod','kmol m-3 s-1','H2S production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'h2s_loss',hamocc_state_acc%h2sloss,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2s_loss','kmol m-3 s-1','H2S oxidation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'dms_bac',hamocc_state_acc%dmsbac,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_bac','kmol S m-3 s-1','DMS microbial consumption', DATATYPE_FLT32), &
      & grib2_var(255, 255, 70, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'dms_uv',hamocc_state_acc%dmsuv,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_uv','kmol S m-3 s-1','DMS photolysis', DATATYPE_FLT32), &
      & grib2_var(255, 255, 71, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'eu_export',hamocc_state_acc%euexp,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('eu_export','kmol P m-3 s-1','_export_ variable (ecan*zoomor etc)', DATATYPE_FLT32), &
      & grib2_var(255, 255, 71, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'pho_cya',hamocc_state_acc%phoc,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('pho_cya','kmol P m-3 s-1','total cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 141, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cya_loss',hamocc_state_acc%cyloss,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cya_loss','kmol P m-3 s-1','cyanobacteria dying', DATATYPE_FLT32), &
      & grib2_var(255, 255, 142, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'grazing',hamocc_state_acc%graz,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('graz','kmol P m-3 s-1','zooplankton grazing', DATATYPE_FLT32), &
      & grib2_var(255, 255, 101, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'remin_via_grazer',hamocc_state_acc%graton,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('graz','kmol P m-3 s-1','remineralization_via_grazer', DATATYPE_FLT32), &
      & grib2_var(255, 255, 292, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'exudation_phy',hamocc_state_acc%exud,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('exudation_phy','kmol P m-3 s-1','phytoplankton exudation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'exudation_zoo',hamocc_state_acc%exudz,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('exudation_zoo','kmol P m-3 s-1','zooplankton exudation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'zoomor',hamocc_state_acc%zoomor,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('zoomor','kmol P m-3 s-1','zooplankton mortality', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'phymor',hamocc_state_acc%phymor,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('phymor','kmol P m-3 s-1','phytoplankton mortality', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'co2flux',hamocc_state_acc%cflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('co2flux','kmol C m-2 s-1','co2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 92, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'orginp',hamocc_state_acc%orginp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('orginp','check','accumulated organic input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 204, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_acc_list, 'dmsflux',hamocc_state_acc%dmsflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('dmsflux','kmol S m-2 s-1','co2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 68, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'silinp',hamocc_state_acc%silinp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('silinp','check','accumulated silicate input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 205, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_acc_list, 'calinp',hamocc_state_acc%calinp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calinp','check','accumulated calc input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 203, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_acc_list, 'delsil',hamocc_state_acc%delsil,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('delsil','kmol Si m-3 s-1','opal production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 86, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'delcar',hamocc_state_acc%delcar,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('delcar','kmol C m-3 s-1','calcium carbonate production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 85, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'o2flux',hamocc_state_acc%oflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('o2flux','kmol O m-2 s-1','o2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 72, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=kilo,new_cf=t_cf_var('o2flux','mol O m-2 s-1','o2 flux', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_acc_list, 'n2flux',hamocc_state_acc%nflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('n2flux','kmol N m-2 s-1','n2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 74, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'n2oflux',hamocc_state_acc%n2oflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('n2oflux','kmol N m-2 s-1','n2o flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 93, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)


    CALL add_var(hamocc_acc_list, 'nfix_diag',hamocc_state_acc%nfixd,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('nfix_diag','kmol N2 m-2 s-1','diagnostic n fixation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 157, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'silpro',hamocc_state_acc%silpro,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('silpro','kmol Si m-2 s-1','Si flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 96, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'produs',hamocc_state_acc%produs,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('produs','kmol m-2 s-1','dust flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 97, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'prcaca',hamocc_state_acc%prcaca,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('prcaca','kmol C m-2 s-1','C inorg flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 95, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'prorca',hamocc_state_acc%prorca,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('prorca','kmol P m-2 s-1','C org flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 94, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'coex90',hamocc_state_acc%coex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex90','kmol P m-2 s-1','OM flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 81, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'calex90',hamocc_state_acc%calex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex90','kmol C m-2 s-1','Calc flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 78, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'opex90',hamocc_state_acc%opex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex90','kmol Si m-2 s-1','Opal flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 75, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'coex1000',hamocc_state_acc%coex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex1000','kmol P m-2 s-1','OM flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 81, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'calex1000',hamocc_state_acc%calex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex1000','kmol C m-2 s-1','Calc flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 78, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'opex1000',hamocc_state_acc%opex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex1000','kmol Si m-2 s-1','Opal flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 75, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'coex2000',hamocc_state_acc%coex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex2000','kmol P m-2 s-1','OM flux at 2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 81, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'calex2000',hamocc_state_acc%calex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex2000','kmol C m-2 s-1','Calc flux at 2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 78, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'opex2000',hamocc_state_acc%opex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex2000','kmol Si m-2 s-1','Opal flux at 2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 75, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_dic',hamocc_state_acc%sedflic,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_dic','kmol C m-2 s-1 ','sediment-ocean flux DIC', DATATYPE_FLT32), &
      & grib2_var(255, 255, 280, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_alk',hamocc_state_acc%sedflal,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_alk','kmol m-2 s-1','sediment-ocean flux alkalinity', DATATYPE_FLT32), &
      & grib2_var(255, 255, 281, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_po4',hamocc_state_acc%sedflph,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_po4','kmol P m-2 s-1','sediment-ocean flux phosphate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 282, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_o2',hamocc_state_acc%sedflox,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_po2','kmol O2 m-2 s-1','sediment-ocean flux oxygen', DATATYPE_FLT32), &
      & grib2_var(255, 255, 283, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_n2',hamocc_state_acc%sedfln2,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_n2','kmol N m-2 s-1','sediment-ocean flux nitrogen', DATATYPE_FLT32), &
      & grib2_var(255, 255, 284, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_no3',hamocc_state_acc%sedflno3,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_no3','kmol N m-2 s-1','sediment-ocean flux nitrate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 285, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_si',hamocc_state_acc%sedflsi,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_si','kmol Si m-2 s-1','sediment-ocean flux nitrate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 286, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_fe',hamocc_state_acc%sedflfe,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_fe','kmol Fe m-2 s-1','sediment-ocean flux iron', DATATYPE_FLT32), &
      & grib2_var(255, 255, 287, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_remino2',hamocc_state_acc%sedro2,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_remino2','kmol P m-3 s-1','sediment aerob remineralization', DATATYPE_FLT32), &
      & grib2_var(255, 255, 300, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_reminn',hamocc_state_acc%sedrn,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_reminn','kmol P m-3 s-1','sediment denitrification', DATATYPE_FLT32), &
      & grib2_var(255, 255, 301, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_remins',hamocc_state_acc%sedrs,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_remins','kmol P m-3 s-1','sediment sulfate reduction', DATATYPE_FLT32), &
      & grib2_var(255, 255, 302, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'flim',hamocc_state_acc%flim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('flim',' ','Fe limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 138, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'plim',hamocc_state_acc%plim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('plim',' ','P_limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 139, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'nlim',hamocc_state_acc%nlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('nlim',' ','N_limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 140, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'aou',hamocc_state_acc%aou,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('aou','kmol O2 m-3','apparent O2 utilisation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 143, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'cTlim',hamocc_state_acc%cTlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cTlim',' ','Temp_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'cLlim',hamocc_state_acc%cLlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cLlim',' ','Light_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'cPlim',hamocc_state_acc%cPlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cPlim',' ','P_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_acc_list, 'cFlim',hamocc_state_acc%cFlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cFlim',' ','Fe_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'o2min',hamocc_state_acc%o2min,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('o2min','mol m-3','mole_concentration_of_dissolved_molecular_&
      &oxygen_in_sea_water_at_shallowest_local_minimum_in_vertical_profile', DATATYPE_FLT32), &
      & grib2_var(255, 255, 162, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'zo2min',hamocc_state_acc%zo2min,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('zo2min','m','depth_at_shallowest_local_minimum_&
      &in_vertical_profile_of_mole_concentration_of_dissolved_molecular&
      &_oxygen_in_sea_water', DATATYPE_FLT32), &
      & grib2_var(255, 255, 163, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL message(TRIM(routine), 'construct hamocc acc end')
  END SUBROUTINE 

!================================================================================== 
  SUBROUTINE construct_hamocc_moni(patch_2d, hamocc_state_moni)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hamocc_monitor), INTENT(inout)      :: hamocc_state_moni
    
    ! local variables
    
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_tend'

     ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v



    !-----DIAG W/O restart-----------------------------------------------------------------
    ! for tracers restart is handled by ICON
    CALL message(TRIM(routine), 'start to construct hamocc monitoring')
 ! 

    ! add monitoring 
   CALL add_var(hamocc_tendency_list, 'global_primary_production', hamocc_state_moni%phosy , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_primary_production', 'GtC s-1', 'global net primary production', DATATYPE_FLT32),&
      & grib2_var(255, 255, 500, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_primary_production','GtC yr-1','global net_npp', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_zooplankton_grazing', hamocc_state_moni%grazing , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_zooplankton_grazing', 'GtC s-1', 'global zooplankton grazing', DATATYPE_FLT32),&
      & grib2_var(255, 255, 501, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_zooplankton_grazing','GtC yr-1','global zooplankton grazing', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_OM_export_at_90m', hamocc_state_moni%omex90, &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_OM_export_at_90m', 'GtC s-1', 'global_om_export_at_90m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 502, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_OM_export_at_90m','GtCyr-1','global_om_export_at_90m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_calc_export_at_90m', hamocc_state_moni%calex90 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_calc_export_at_90m', 'GtC s-1', 'global_calc_export_at_90m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 503, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_calc_export_at_90m','GtCyr-1','global_calc_export_at_90m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_opal_export_at_90m', hamocc_state_moni%opex90 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_opal_export_at_90m', 'kmol Si s-1', 'global_opal_export_at_90m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 504, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_opal_production', hamocc_state_moni%delsil , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_opal_production', 'kmol Si s-1', 'global opal production', DATATYPE_FLT32),&
      & grib2_var(255, 255, 505, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_caco3_production', hamocc_state_moni%delcar , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_caco3_production', 'GtC s-1', 'global caco3 production', DATATYPE_FLT32),&
      & grib2_var(255, 255, 506, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_caco3_production','GtC yr-1','global caco3 production', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_net_co2_flux', hamocc_state_moni%net_co2_flux , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_net_co2_flux', 'GtC s-1', 'global net co2 flux', DATATYPE_FLT32),&
      & grib2_var(255, 255, 507, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_net_co2_flux','GtC yr-1','global net_co2_flux', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_surface_alk', hamocc_state_moni%sfalk , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_surface_alk', 'kmol  m-3', 'global_surface_alkalinity', DATATYPE_FLT32),&
      & grib2_var(255, 255, 508, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_surface_dic', hamocc_state_moni%sfdic , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_surface_dic', 'kmol C m-3', 'global_surface_dissolved_inorganic_carbon', DATATYPE_FLT32),&
      & grib2_var(255, 255, 509, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'global_remin_via_grazer', hamocc_state_moni%graton , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_remin_via_grazer', 'GtC s-1', 'global_remineralization_via_grazer', DATATYPE_FLT32),&
      & grib2_var(255, 255, 511, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_remin_via_grazer','GtC yr-1','global_remineraliation_via_grazer', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_exudation_phytoplankton', hamocc_state_moni%exud , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_exudation_phytoplankton', 'GtC s-1', 'global_phytoplankton_exudation', DATATYPE_FLT32),&
      & grib2_var(255, 255, 512, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_exudation_phytoplankton','GtC yr-1','global_phytoplankton_exudation', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_phytoplankton_dying', hamocc_state_moni%phymor , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_phytoplankton_dying', 'GtC s-1', 'global_phytoplankton_dying', DATATYPE_FLT32),&
      & grib2_var(255, 255, 513, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_phytoplankton_dying','GtC yr-1','global_phytoplankton_dying', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'EU_bacterial_activity', hamocc_state_moni%bacfra , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('EU_bacterial_activity', 'GtC s-1', 'EU_bacterial_activity', DATATYPE_FLT32),&
      & grib2_var(255, 255, 514, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('EU_bacterial_activity','GtC yr-1','EU_bacterial_activity', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_exudation_zooplankton', hamocc_state_moni%exudz , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_exudation_zooplankton', 'GtC s-1', 'global_zooplankton_exudation', DATATYPE_FLT32),&
      & grib2_var(255, 255, 515, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_exudation_zooplankton','GtC yr-1','global_zooplankton_exudation', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_zooplankton_dying', hamocc_state_moni%zoomor , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_zooplankton_dying', 'GtC s-1', 'global_zooplankton_dying', DATATYPE_FLT32),&
      & grib2_var(255, 255, 516, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_zooplankton_dying','GtC yr-1','global_zooplankton_dying', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'N2_fixation', hamocc_state_moni%n2fix , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('N2_fixation', 'TgN s-1', 'global N2 fixation', DATATYPE_FLT32),&
      & grib2_var(255, 255, 517, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('N2_fixation','TgN yr-1','global N2 fixation', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'WC_denit', hamocc_state_moni%wcdenit , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('WC_denit', 'TgN s-1', 'global water column denitrification', DATATYPE_FLT32),&
      & grib2_var(255, 255, 519, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('WC_denit','TgN yr-1','global watern column denitrification', &
      &  DATATYPE_FLT32)))

    CALL add_var(hamocc_tendency_list, 'global_npp_cya', hamocc_state_moni%phosy_cya , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_npp_cya', 'GtC s-1', 'annual net primary production of cyanobacteria', DATATYPE_FLT32),&
      & grib2_var(255, 255, 520, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_primary_production_cya','GtC yr-1',&
      & 'global annual primary production of cyanobacteria', &
      &  DATATYPE_FLT32)))


   CALL add_var(hamocc_tendency_list, 'Aerob_remin_of_detritus', hamocc_state_moni%remina , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('Aerob_remin_of_detritus', 'GtC s-1', 'Aerob_remineralization_of_detritus', DATATYPE_FLT32),&
      & grib2_var(255, 255, 522, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('Aerob_remin_of_detritus','GtC yr-1','Aerob_remineralization_of_detritus', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'SED_denit', hamocc_state_moni%seddenit , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('SED_denit', 'TgN s-1', 'Sediment_denitrification', DATATYPE_FLT32),&
      & grib2_var(255, 255, 523, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('SED_denit','TgN yr-1','Sediment_denitrification', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_cya_loss_det', hamocc_state_moni%cyaldet , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_cya_loss_det', 'GtC s-1', 'Global_cyanobacteria_loss_to_detritus', DATATYPE_FLT32),&
      & grib2_var(255, 255, 524, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_cya_loss_det','GtC yr-1','Global_cyanobacteria_loss_to_detritus', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_cya_loss_doc', hamocc_state_moni%cyaldoc , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_cya_loss_doc', 'GtC s-1', 'Global_cyanobacteria_loss_to_DOC', DATATYPE_FLT32),&
      & grib2_var(255, 255, 525, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_cya_loss_doc','GtC yr-1','Global_cyanobacteria_loss_to_DOC', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_surface_phosphate', hamocc_state_moni%sfphos , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_surface_phoshate', 'kmol  m-3', 'global_surface_phosphate', DATATYPE_FLT32),&
      & grib2_var(255, 255, 536, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_surface_silicate', hamocc_state_moni%sfsil , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_surface_silicate', 'kmol  m-3', 'global_surface_silicate', DATATYPE_FLT32),&
      & grib2_var(255, 255, 537, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)


   CALL add_var(hamocc_tendency_list, 'zalkn2', hamocc_state_moni%zalkn2 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('zalkn2', 'kmol ', 'global_H+_compensation_for_N2_production/fixation', DATATYPE_FLT32),&
      & grib2_var(255, 255, 538, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_OM_export_at_1000m', hamocc_state_moni%omex1000, &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_OM_export_at_1000m', 'GtC s-1', 'global_om_export_at_1000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 539, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_OM_export_at_1000m','GtCyr-1','global_om_export_at_1000m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_calc_export_at_1000m', hamocc_state_moni%calex1000 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_calc_export_at_1000m', 'GtC s-1', 'global_calc_export_at_1000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 540, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_calc_export_at_1000m','GtCyr-1','global_calc_export_at_1000m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_opal_export_at_1000m', hamocc_state_moni%opex1000 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_opal_export_at_1000m', 'kmol Si s-1', 'global_opal_export_at_1000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 541, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'global_OM_export_at_2000m', hamocc_state_moni%omex2000, &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_OM_export_at_2000m', 'GtC s-1', 'global_om_export_at_2000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 542, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_OM_export_at_2000m','GtCyr-1','global_om_export_at_2000m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_calc_export_at_2000m', hamocc_state_moni%calex2000 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_calc_export_at_2000m', 'GtC s-1', 'global_calc_export_at_2000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 543, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE., post_op=post_op(POST_OP_SCALE,&
      &  arg1=s2year,new_cf=t_cf_var('global_calc_export_at_2000m','GtCyr-1','global_calc_export_at_2000m', &
      &  DATATYPE_FLT32)))

   CALL add_var(hamocc_tendency_list, 'global_opal_export_at_2000m', hamocc_state_moni%opex2000 , &
      & GRID_LONLAT, za_surface,    &
      & t_cf_var('global_opal_export_at_2000m', 'kmol Si s-1', 'global_opal_export_at_2000m', DATATYPE_FLT32),&
      & grib2_var(255, 255, 544, DATATYPE_PACK16, grid_reference, grid_lonlat),&
      & in_group=groups("HAMOCC_MONI"),ldims=(/1/), &
      & loutput=.TRUE., lrestart=.FALSE.)

  END SUBROUTINE 




!================================================================================== 
  SUBROUTINE construct_hamocc_tend(patch_2d, hamocc_state_tend)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hamocc_tend), INTENT(inout)         :: hamocc_state_tend
    
    ! local variables
    
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_moni'

     ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v



    !-----DIAG W/O restart-----------------------------------------------------------------
    ! for tracers restart is handled by ICON
    CALL message(TRIM(routine), 'start to construct hamocc tendency state')
   
    CALL add_var(hamocc_tendency_list, 'NPP',hamocc_state_tend%npp,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('NPP','kmolP m-3 s-1','net primary production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'BACFRA',hamocc_state_tend%bacfra,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('BACFRA','kmolP m-3 s-1','bacterial decomposition of DOC', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)
   

    CALL add_var(hamocc_tendency_list, 'SRED',hamocc_state_tend%remins,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('SRED','kmolP m-3 s-1','sulfate reduction', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)
      
    CALL add_var(hamocc_tendency_list, 'REMIN',hamocc_state_tend%remina,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('REMIN','kmolP m-3 s-1','aerob detritus remineralization', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'DENIT',hamocc_state_tend%reminn,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('DENIT','kmolP m-3 s-1','denitrification', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'nfix',hamocc_state_tend%nfix,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('nfix','kmol N2 m-3 s-1','n fixation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'pho_cya',hamocc_state_tend%phoc,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('pho_cya','kmol P m-3 s-1','total cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cya_loss',hamocc_state_tend%cyloss,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cya_loss','kmol P m-3 s-1','cyanobacteria dying', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'h2obudget',hamocc_state_tend%h2obudget,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2obudget','kmol','h2o budget', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'n2budget',hamocc_state_tend%n2budget,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('n2budget','kmol','N2 budget', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_tendency_list, 'grazing',hamocc_state_tend%graz,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('graz','kmol P m-3 s-1','zooplankton grazing', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'remin_via_grazer',hamocc_state_tend%graton,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('remin_via_grazer','kmol P m-3 s-1','remineralization_via_grazer', DATATYPE_FLT32), &
      & grib2_var(255, 255, 292, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'exudation_phy',hamocc_state_tend%exud,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('exudation_phy','kmol P m-3 s-1','phytoplankton exudation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'exudation_zoo',hamocc_state_tend%exudz,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('exudation_zoo','kmol P m-3 s-1','zooplankton exudation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'zoomor',hamocc_state_tend%zoomor,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('zoomor','kmol P m-3 s-1','zooplankton mortality', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'phymor',hamocc_state_tend%phymor,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('phymor','kmol P m-3 s-1','phytoplankton mortality', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'orginp',hamocc_state_tend%orginp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('orginp','check','organic input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'silinp',hamocc_state_tend%silinp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('silinp','check','silicate input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'calinp',hamocc_state_tend%calinp,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calinp','check','calc input', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'dms_prod',hamocc_state_tend%dmsprod,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_prod','kmol S m-3 s-1','DMS production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'h2s_prod',hamocc_state_tend%h2sprod,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2s_prod','kmol m-3 s-1','H2S production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'h2s_loss',hamocc_state_tend%h2sloss,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('h2s_loss','kmol m-3 s-1','H2S oxidation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'dms_bac',hamocc_state_tend%dmsbac,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_bac','kmol S m-3 s-1','DMS microbial consumption', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'dms_uv',hamocc_state_tend%dmsuv,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('dms_uv','kmol S m-3 s-1','DMS photolysis', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'eu_export',hamocc_state_tend%euexp,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('eu_export','kmol P m-3 s-1','_export_ variable (ecan*zoomor etc)', DATATYPE_FLT32), &
      & grib2_var(255, 255, 71, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'flim',hamocc_state_tend%flim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('flim','','Fe limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 138, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'plim',hamocc_state_tend%plim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('plim','','P_limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 139, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'nlim',hamocc_state_tend%nlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('nlim','','N_limitation of PP', DATATYPE_FLT32), &
      & grib2_var(255, 255, 140, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cTlim',hamocc_state_tend%cTlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cTlim','','Temp_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cLlim',hamocc_state_tend%cLlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cLlim','','Light_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cPlim',hamocc_state_tend%cPlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cPlim','','P_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'cFlim',hamocc_state_tend%cFlim,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('cFlim','','Fe_limitation of cyanobacteria growth', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'o2min',hamocc_state_tend%o2min,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('o2min','mol m-3','mole_concentration_of_dissolved_molecular_&
      &oxygen_in_sea_water_at_shallowest_local_minimum_in_vertical_profile', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_acc_list, 'zo2min',hamocc_state_tend%zo2min,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('zo2min','m','depth_at_shallowest_local_minimum_&
      &in_vertical_profile_of_mole_concentration_of_dissolved_molecular&
      &_oxygen_in_sea_water', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'dmsflux',hamocc_state_tend%dmsflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('dmsflux','kmol S m-2 s-1','co2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'co2flux',hamocc_state_tend%cflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('co2flux','kmol C m-2 s-1','co2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'o2flux',hamocc_state_tend%oflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('o2flux','kmol O m-2 s-1','o2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'n2flux',hamocc_state_tend%nflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('n2flux','kmol N m-2 s-1','n2 flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'n2oflux',hamocc_state_tend%n2oflux,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('n2oflux','kmol N m-2 s-1','n2o flux', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'nfix_diag',hamocc_state_tend%nfixd,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('nfix_diag','kmol N m-2 s-1','diagnostic n fixation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'silpro',hamocc_state_tend%silpro,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('silpro','kmol Si m-2 s-1','Si flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

   CALL add_var(hamocc_tendency_list, 'produs',hamocc_state_tend%produs,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('produs','kmol m-2 s-1','Dust flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)


    CALL add_var(hamocc_tendency_list, 'prcaca',hamocc_state_tend%prcaca,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('prcaca','kmol C m-2 s-1','C inorg flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'prorca',hamocc_state_tend%prorca,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('prorca','kmol P m-2 s-1','C org flux to sediment', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'coex90',hamocc_state_tend%coex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex90','kmol P m-2 s-1','OM flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 81, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'calex90',hamocc_state_tend%calex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex90','kmol C m-2 s-1','Calc flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 78, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'opex90',hamocc_state_tend%opex90,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex90','kmol Si m-2 s-1','Opal flux at 90 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 75, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'coex1000',hamocc_state_tend%coex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex1000','kmol P m-2 s-1','OM flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 82, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'calex1000',hamocc_state_tend%calex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex1000','kmol C m-2 s-1','Calc flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 79, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'opex1000',hamocc_state_tend%opex1000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex1000','kmol Si m-2 s-1','Opal flux at 1000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 76, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'coex2000',hamocc_state_tend%coex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('coex2000','kmol P m-2 s-1','OM flux at  2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 83, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'calex2000',hamocc_state_tend%calex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('calex2000','kmol C m-2 s-1','Calc flux at  2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 80, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'opex2000',hamocc_state_tend%opex2000,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('opex2000','kmol Si m-2 s-1','Opal flux at  2000 m', DATATYPE_FLT32), &
      & grib2_var(255, 255, 77, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'delsil',hamocc_state_tend%delsil,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('delsil','kmol Si m-3 s-1','opal production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 86, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'delcar',hamocc_state_tend%delcar,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('delcar','kmol C m-3 s-1','calcium carbonate production', DATATYPE_FLT32), &
      & grib2_var(255, 255, 85, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_dic',hamocc_state_tend%sedflic,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_dic','kmol C m-2 s-1 ','sediment-ocean flux DIC', DATATYPE_FLT32), &
      & grib2_var(255, 255, 280, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_alk',hamocc_state_tend%sedflal,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_alk','kmol m-2 s-1','sediment-ocean flux alkalinity', DATATYPE_FLT32), &
      & grib2_var(255, 255, 281, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_po4',hamocc_state_tend%sedflph,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_po4','kmol P m-2 s-1','sediment-ocean flux phosphate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 282, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_o2',hamocc_state_tend%sedflox,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_po2','kmol O2 m-2 s-1','sediment-ocean flux oxygen', DATATYPE_FLT32), &
      & grib2_var(255, 255, 283, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_n2',hamocc_state_tend%sedfln2,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_n2','kmol N m-2 s-1','sediment-ocean flux nitrogen', DATATYPE_FLT32), &
      & grib2_var(255, 255, 284, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_no3',hamocc_state_tend%sedflno3,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_no3','kmol N m-2 s-1','sediment-ocean flux nitrate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 285, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_si',hamocc_state_tend%sedflsi,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_si','kmol Si m-2 s-1','sediment-ocean flux nitrate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 286, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_fe',hamocc_state_tend%sedflfe,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_fe','kmol Fe m-2 s-1','sediment-ocean flux iron', DATATYPE_FLT32), &
      & grib2_var(255, 255, 287, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sedflux_h2s',hamocc_state_tend%sedflh2s,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sedflux_h2s','kmol m-2 s-1','sediment-ocean flux H2S', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_remino2',hamocc_state_tend%sedro2,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_remino2','kmol P m-3 s-1','sediment aerob remineralization', DATATYPE_FLT32), &
      & grib2_var(255, 255, 300, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_reminn',hamocc_state_tend%sedrn,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_reminn','kmol P m-3 s-1','sediment denitrification', DATATYPE_FLT32), &
      & grib2_var(255, 255, 301, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_sediment_list, 'sed_remins',hamocc_state_tend%sedrs,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('sed_remins','kmol P m-3 s-1','sediment sulfate reduction', DATATYPE_FLT32), &
      & grib2_var(255, 255, 302, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'akb',hamocc_state_tend%akb,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('akb','','dissociation constant boric acid', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'akw',hamocc_state_tend%akw,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('akw','','dissociation constant water', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'ak1',hamocc_state_tend%ak1,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('ak1','','dissociation constant k1', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'ak2',hamocc_state_tend%ak2,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('ak2','','dissociation constant k2', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'aksp',hamocc_state_tend%aksp,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('aksp','','apparent solubility product', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'satoxy',hamocc_state_tend%satoxy,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('satoxy','','O2 at saturation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'satn2',hamocc_state_tend%satn2,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('satn2','','N2 at saturation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'satn2o',hamocc_state_tend%satn2o,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('satn2o','','N2O at saturation', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'solco2',hamocc_state_tend%solco2,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('solco2','','CO2 solubility', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL add_var(hamocc_tendency_list, 'aou',hamocc_state_tend%aou,    &
      & grid_unstructured_cell, za_depth_below_sea,&
      & t_cf_var('aou','','apparent O2 utilisation ', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,n_zlev,alloc_cell_blocks/),in_group=groups("HAMOCC_TEND"),&
      & loutput=.FALSE., lrestart=.FALSE.)

    CALL message(TRIM(routine), 'construct hamocc tend end')

  END SUBROUTINE 

!================================================================================== 
  SUBROUTINE construct_hamocc_sed(patch_2d, hamocc_state_sed)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: patch_2d
    TYPE(t_hamocc_sed), INTENT(inout)          :: hamocc_state_sed
    
    ! local variables
    
    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:construct_hamocc_sed'

     ! determine size of arrays
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    nblks_v = patch_2d%nblks_v


    CALL message(TRIM(routine), 'start to construct hamocc sed state')
  
    CALL add_var(hamocc_sediment_list, 'SED_C12org',hamocc_state_sed%so12,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('SED_C12org','kmol P m-3 ','solid sediment C org', DATATYPE_FLT64), &
      & grib2_var(255, 255, 38, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'SED_C12',hamocc_state_sed%sc12,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('SED_C12','kmol C m-3 ','solid sediment calcium carbonate', DATATYPE_FLT64), &
      & grib2_var(255, 255, 41, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'SED_Si',hamocc_state_sed%ssil,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('SED_Si','kmol Si m-3','solid sediment opal', DATATYPE_FLT64), &
      & grib2_var(255, 255, 44, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'SED_clay',hamocc_state_sed%ster,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('SED_clay','kmol m-3','solid sediment clay', DATATYPE_FLT64), &
      & grib2_var(255, 255, 45, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_DIC',hamocc_state_sed%pwic,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_DIC','kmol C m -3','sediment pore water DIC', DATATYPE_FLT64), &
      & grib2_var(255, 255, 51, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_alk',hamocc_state_sed%pwal,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_alk','kmol m-3','sediment pore water alkalinity', DATATYPE_FLT64), &
      & grib2_var(255, 255, 54, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_phos',hamocc_state_sed%pwph,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_phos','kmol P m-3','sediment pore water phosphate', DATATYPE_FLT64), &
      & grib2_var(255, 255, 55, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_o2',hamocc_state_sed%pwox,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_o2','kmol O2 m-3','sediment pore water oxygen', DATATYPE_FLT64), &
      & grib2_var(255, 255, 56, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_si',hamocc_state_sed%pwsi,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_si','kmol Si m-3','sediment pore water silicate', DATATYPE_FLT64), &
      & grib2_var(255, 255, 59, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_fe',hamocc_state_sed%pwfe,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_fe','kmol Fe m-3','sediment pore water iron', DATATYPE_FLT64), &
      & grib2_var(255, 255, 60, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_n2',hamocc_state_sed%pwn2,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_n2','kmoll N2 m-3','sediment pore water N2', DATATYPE_FLT64), &
      & grib2_var(255, 255, 57, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_no3',hamocc_state_sed%pwno3,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_no3','kmol N m-3','sediment pore water nitrate', DATATYPE_FLT64), &
      & grib2_var(255, 255, 58, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_h2s',hamocc_state_sed%pwh2s,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_h2s','kmol m-3','sediment pore water H2S alkalinity', DATATYPE_FLT64), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_n2b',hamocc_state_sed%pwn2b,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_n2b','kmol N m-3','sediment pore water n2 budget', DATATYPE_FLT64), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'POW_h2ob',hamocc_state_sed%pwh2ob,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('POW_h2ob','kmol N m-3','sediment pore water h2o budget', DATATYPE_FLT64), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.FALSE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'sedhpl',hamocc_state_sed%sedhi,    &
      & grid_unstructured_cell, za_ocean_sediment,&
      & t_cf_var('sedhpl','kmol m-3','sediment hydrogen ion concentration', DATATYPE_FLT64), &
      & grib2_var(255, 255, 50, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,ks,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'BUR_C12org',hamocc_state_sed%bo12,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('BUR_C12org','kmol P m-2 ','sediment burial C org', DATATYPE_FLT64), &
      & grib2_var(255, 255, 46, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'BUR_C12',hamocc_state_sed%bc12,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('BUR_C12','kmol C m-2 ','sediment burial calcium carbonate', DATATYPE_FLT32), &
      & grib2_var(255, 255, 47, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'BUR_clay',hamocc_state_sed%bter,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('BUR_clay','kmol  m-2 ','sediment burial clay', DATATYPE_FLT32), &
      & grib2_var(255, 255, 49, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)

    CALL add_var(hamocc_sediment_list, 'BUR_si',hamocc_state_sed%bsil,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('BUR_si','kmol Si m-2 ','sediment burial opal', DATATYPE_FLT32), &
      & grib2_var(255, 255, 48, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.TRUE.,lrestart_cont=.TRUE.)
  
   CALL add_var(hamocc_sediment_list, 'bolay',hamocc_state_sed%bolay,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('bolay','m ','bottom layer thickness', DATATYPE_FLT32), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

   CALL add_var(hamocc_sediment_list, 'kbo',hamocc_state_sed%kbo,    &
      & grid_unstructured_cell, za_surface,&
      & t_cf_var('kbo',' ','index of bottom layer', DATATYPE_INT8), &
      & grib2_var(255, 255, 255, DATATYPE_PACK16, grid_reference, grid_cell),&
      & ldims=(/nproma,alloc_cell_blocks/),in_group=groups("HAMOCC_SED"),&
      & loutput=.TRUE., lrestart=.FALSE.)

    CALL message(TRIM(routine), 'construct hamocc sed end')

  END SUBROUTINE 

!================================================================================== 
  SUBROUTINE destruct_hamocc_state(hamocc_state)
    TYPE(t_hamocc_state), TARGET,INTENT(inout)   :: hamocc_state!(n_dom)
    
    ! local variables
    
    INTEGER :: jg
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_bgc_icon_comm:destruct_hydro_ocean_state'
    
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start to destruct hamocc state ')
    
    
    CALL delete_var_list(hamocc_restart_list)
    CALL delete_var_list(hamocc_default_list)
    CALL delete_var_list(hamocc_tendency_list)
    CALL delete_var_list(hamocc_acc_list)
    CALL delete_var_list(hamocc_sediment_list)
    
!     DO jg = 1, n_dom
!       CALL destruct_hamocc_diag(hamocc_state(jg)%p_diag)
!       
!       
!     END DO
!     
    CALL message(TRIM(routine),'destruction of hamocc state finished')
    CALL close_bgcout 

   
  END SUBROUTINE 

!================================================================================== 
 
!  SUBROUTINE destruct_hamocc_diag(hamocc_diag)
!     
!     TYPE(t_hamocc_diag), INTENT(inout) :: hamocc_diag
!     
!     ! local variables
!     
!     INTEGER :: ist
!     
!     CHARACTER(LEN=max_char_length), PARAMETER :: &
!       & routine = 'mo_bgc_icon_comm:destruct_hamocc_diag'
!     
!     DEALLOCATE(hamocc_diag%p_vn, stat=ist)
!     IF (ist/=success) THEN
!       CALL finish(TRIM(routine), 'deallocation for hamocc_diag failed')
!     END IF
!     
!   END SUBROUTINE
!================================================================================== 

    SUBROUTINE construct_hamocc_var_lists(patch_2d)

    TYPE(t_patch), TARGET, INTENT(in) :: patch_2d
    
    CHARACTER(LEN=max_char_length) :: listname
    
    WRITE(listname,'(a)')  'hamocc_restart_list'
    CALL new_var_list(hamocc_restart_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( hamocc_restart_list,             &
      & lrestart=.TRUE.,loutput=.TRUE.,&
      & model_type='oce' )

    WRITE(listname,'(a)')  'hamocc_default_list'
    CALL new_var_list(hamocc_default_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( hamocc_default_list,            &
      & lrestart=.FALSE.,model_type='oce',loutput=.TRUE. )

    WRITE(listname,'(a)')  'hamocc_acc_list'
    CALL new_var_list(hamocc_acc_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( hamocc_acc_list,            &
      & lrestart=.TRUE.,model_type='oce',loutput=.TRUE. )

    WRITE(listname,'(a)')  'hamocc_tendency_list'
    CALL new_var_list(hamocc_tendency_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( hamocc_tendency_list,            &
      & lrestart=.TRUE.,model_type='oce',loutput=.TRUE. )

    WRITE(listname,'(a)')  'hamocc_sediment_list'
    CALL new_var_list(hamocc_sediment_list, listname, patch_id=patch_2d%id)
    CALL default_var_list_settings( hamocc_sediment_list,            &
      & lrestart=.TRUE.,model_type='oce',loutput=.TRUE. )

    END SUBROUTINE construct_hamocc_var_lists
!================================================================================== 
  SUBROUTINE close_bgcout
!

    INTEGER :: istat

!-------------------------------------------------------------------------

    CLOSE (io_stdo_bgc, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish ('close_bgcout','Could not close bgcout')
    END IF

  END SUBROUTINE close_bgcout

 END MODULE

#endif
