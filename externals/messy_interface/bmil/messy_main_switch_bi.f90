! ***************************************************************************
MODULE messy_main_switch_bi
! ***************************************************************************

  ! MODULE FOR SWITCHING SUBMODEL ON/OFF
  ! Authors:
  !   Patrick Joeckel, MPICH, September 2002

  ! NOTE: TO ADD A NEW SWITCH LOOK FOR
  !       '### ADD HERE'
  !       AND ADD IT ALSO TO messy_main_switch.f90

  ! MESSy
  USE messy_main_switch
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: main_switch_setup
  
CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE main_switch_setup

    ! ECHAM5/MESSy
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_channel_bi, ONLY: channel_halt
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_switch_setup'
    INTEGER     :: iou    ! I/O unit
    LOGICAL     :: lex    ! file exists ?
    INTEGER     :: fstat  ! file status
    INTEGER     :: status

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL messy_main_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    !
    CALL p_bcast(L_TIME_INFO, p_io)
    !
    CALL p_bcast(USE_AEROPT,  p_io)
    CALL p_bcast(USE_AIRSEA,  p_io)
    CALL p_bcast(USE_BUFLY,   p_io)
    CALL p_bcast(USE_CLOUD,   p_io)
    CALL p_bcast(USE_CONVECT, p_io)
    CALL p_bcast(USE_CVTRANS, p_io)
    CALL p_bcast(USE_D14CO,   p_io)
    CALL p_bcast(USE_DDEP,    p_io)
    CALL p_bcast(USE_DRADON,  p_io)
    CALL p_bcast(USE_E4CHEM,  p_io)
    CALL p_bcast(USE_GMXE,    p_io)
    CALL p_bcast(USE_GWAVE,   p_io)  
    CALL p_bcast(USE_H2O,     p_io)
    CALL p_bcast(USE_HETCHEM, p_io)
    CALL p_bcast(USE_JVAL,    p_io)
    CALL p_bcast(USE_LNOX,    p_io)
    CALL p_bcast(USE_M7,      p_io)
    CALL p_bcast(USE_MADE,    p_io)
    CALL p_bcast(USE_MECCA1,  p_io)
    CALL p_bcast(USE_MECCA,   p_io)
    CALL p_bcast(USE_MEGAN,   p_io)
    CALL p_bcast(USE_MLOCEAN, p_io)
    CALL p_bcast(USE_MMFORCE, p_io)
    CALL p_bcast(USE_MSBM,    p_io)
    CALL p_bcast(USE_O3ORIG,  p_io)
    CALL p_bcast(USE_OFFEMIS, p_io)
    CALL p_bcast(USE_ONEMIS,  p_io)
    CALL p_bcast(USE_ORBIT,   p_io)
    CALL p_bcast(USE_PHOTO,   p_io)
    CALL p_bcast(USE_PLUMEGAS,p_io)
    CALL p_bcast(USE_PSC,     p_io)
    CALL p_bcast(USE_PTRAC,   p_io)
    CALL p_bcast(USE_QBO,     p_io)
    CALL p_bcast(USE_RAD4ALL, p_io)
    CALL p_bcast(USE_SATSIMS, p_io)
    CALL p_bcast(USE_SCALC,   p_io)
    CALL p_bcast(USE_SCAV,    p_io)
    CALL p_bcast(USE_SCOUT,   p_io)
    CALL p_bcast(USE_SEDI,    p_io)
    CALL p_bcast(USE_S4D,     p_io)
    CALL p_bcast(USE_SORBIT,  p_io)
    CALL p_bcast(USE_SPE,     p_io)
    CALL p_bcast(USE_SPACENOX,p_io)
    CALL p_bcast(USE_TIMEPOS, p_io)
    CALL p_bcast(USE_TNUDGE,  p_io)
    CALL p_bcast(USE_TREXP,   p_io)
    CALL p_bcast(USE_TROPOP,  p_io)
    CALL p_bcast(USE_VAHR,    p_io)
    CALL p_bcast(USE_VISO,    p_io)
    ! ### example submodels
    CALL p_bcast(USE_SUBMOD1, p_io)
    CALL p_bcast(USE_SUBMOD2, p_io)
    CALL p_bcast(USE_SUBMOD3, p_io)
    ! ### ADD HERE

    CALL switch_init(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_switch_setup
  ! ------------------------------------------------------------------------

  ! op_pj_20110621+
  ! temprarily moved from messy_main_switch.f90 to reduce dependencies
  ! (via USE_MLOCEAN) of BML to SMCL
  ! -------------------------------------------------------------
  SUBROUTINE switch_init(status)

    ! MESSy
    USE messy_main_channel,       ONLY: new_attribute
    USE messy_main_constants_mem, ONLY: modstr_MESSy => modstr &
                                      , modver_MESSy => modver
    USE messy_main_tracer,        ONLY: modstr_tracer => modstr &
                                      , modver_tracer => modver
    USE messy_main_channel,       ONLY: modstr_channel => modstr &
                                      , modver_channel => modver
    USE messy_main_timer,         ONLY: modstr_timer => modstr &
                                      , modver_timer => modver
    ! op_bk_20130820+
#ifndef __ICON__
    USE messy_main_qtimer,        ONLY: modstr_qtimer => modstr &
                                      , modver_qtimer => modver
    USE messy_main_import,        ONLY: modstr_import => modstr &
                                       , modver_import => modver
#endif
    ! op_bk_20130820-

#if defined(ECHAM5) || defined(COSMO)
    ! MESSY SUBMODELS
    USE messy_aeropt,     ONLY: modstr_aeropt=>modstr,  modver_aeropt=>modver
    USE messy_airsea,     ONLY: modstr_airsea=>modstr,  modver_airsea=>modver
    USE messy_bufly,      ONLY: modstr_bufly=>modstr,   modver_bufly=>modver
    USE messy_cloud,      ONLY: modstr_cloud=>modstr,   modver_cloud=>modver
    USE messy_convect,    ONLY: modstr_convect=>modstr, modver_convect=>modver
    USE messy_cvtrans,    ONLY: modstr_cvtrans=>modstr, modver_cvtrans=>modver
    USE messy_d14co,      ONLY: modstr_d14co=>modstr,   modver_d14co=>modver
    USE messy_ddep,       ONLY: modstr_ddep=>modstr,    modver_ddep=>modver
    USE messy_dradon,     ONLY: modstr_dradon=>modstr,  modver_dradon=>modver
    USE messy_e4chem,     ONLY: modstr_e4chem=>modstr,  modver_e4chem=>modver
    USE messy_gmxe,       ONLY: modstr_gmxe=>modstr,    modver_gmxe=>modver
    USE messy_gwave,      ONLY: modstr_gwave=>modstr,   modver_gwave=>modver
    USE messy_h2o,        ONLY: modstr_h2o=>modstr,     modver_h2o=>modver
    USE messy_hetchem,    ONLY: modstr_hetchem=>modstr, modver_hetchem=>modver
    USE messy_jval,       ONLY: modstr_jval=>modstr,    modver_jval=>modver
    USE messy_lnox,       ONLY: modstr_lnox=>modstr,    modver_lnox=>modver
    USE messy_m7,         ONLY: modstr_m7=>modstr,      modver_m7=>modver
    USE messy_made,       ONLY: modstr_made=>modstr,    modver_made=>modver
    USE messy_mecca1,     ONLY: modstr_mecca1=>modstr,  modver_mecca1=>modver
    USE messy_mecca,      ONLY: modstr_mecca=>modstr,   modver_mecca=>modver
    USE messy_megan,      ONLY: modstr_megan=>modstr,   modver_megan=>modver
    USE messy_mlocean,    ONLY: modstr_mlocean=>modstr, modver_mlocean=>modver
    USE messy_mmforce,    ONLY: modstr_mmforce=>modstr, modver_mmforce=>modver
    USE messy_msbm,       ONLY: modstr_msbm=>modstr,    modver_msbm=>modver
    USE messy_o3orig,     ONLY: modstr_o3orig=>modstr,  modver_o3orig=>modver
    USE messy_offemis,    ONLY: modstr_offemis=>modstr, modver_offemis=>modver
    USE messy_onemis,     ONLY: modstr_onemis=>modstr,  modver_onemis=>modver
    USE messy_orbit,      ONLY: modstr_orbit=>modstr,   modver_orbit=>modver
    USE messy_psc,        ONLY: modstr_psc=>modstr,     modver_psc=>modver
    USE messy_photo,      ONLY: modstr_photo=>modstr,   modver_photo=>modver
    USE messy_plumegas,   ONLY: modstr_plumegas=>modstr,modver_plumegas=>modver
    USE messy_ptrac,      ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
    USE messy_qbo,        ONLY: modstr_qbo=>modstr,     modver_qbo=>modver
    USE messy_rad4all,    ONLY: modstr_rad4all=>modstr, modver_rad4all=>modver
    USE messy_satsims,    ONLY: modstr_satsims=>modstr, modver_satsims=>modver
    USE messy_scalc,      ONLY: modstr_scalc=>modstr,   modver_scalc=>modver
    USE messy_scav,       ONLY: modstr_scav=>modstr,    modver_scav=>modver
    USE messy_scout,      ONLY: modstr_scout=>modstr,   modver_scout=>modver
    USE messy_sedi,       ONLY: modstr_sedi=>modstr,    modver_sedi=>modver
    USE messy_s4d,        ONLY: modstr_s4d=>modstr,     modver_s4d=>modver
    USE messy_sorbit,     ONLY: modstr_sorbit=>modstr,  modver_sorbit=>modver
    USE messy_spacenox,   ONLY: modstr_spacenox=>modstr,modver_spacenox=>modver
    USE messy_spe,        ONLY: modstr_spe=>modstr,     modver_spe=>modver
    USE messy_timepos,    ONLY: modstr_timepos=>modstr, modver_timepos=>modver
    USE messy_tnudge,     ONLY: modstr_tnudge=>modstr,  modver_tnudge=>modver
    USE messy_trexp,      ONLY: modstr_trexp=>modstr,   modver_trexp=>modver
    USE messy_tropop,     ONLY: modstr_tropop=>modstr,  modver_tropop=>modver
    USE messy_vahr,       ONLY: modstr_vahr=>modstr,    modver_vahr=>modver
    USE messy_viso,       ONLY: modstr_viso=>modstr,    modver_viso=>modver
#endif
#ifdef BLANK
    USE messy_ptrac,      ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
    USE messy_submod1,    ONLY: modstr_submod1=>modstr, modver_submod1=>modver
    USE messy_submod2,    ONLY: modstr_submod2=>modstr, modver_submod2=>modver
    USE messy_submod3,    ONLY: modstr_submod3=>modstr, modver_submod3=>modver
#endif
#ifdef MBM_CMAT
#endif
    ! ### ADD HERE

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! put MESSy information to output file
    CALL new_attribute(status, 'MESSy' &
         , c = modstr_MESSy//' version '//modver_MESSy//&
         &', http://www.messy-interface.org' )
    IF (status /= 0) RETURN
    !
    CALL put_submodel_att (status, .TRUE., modstr, modver) ! switch
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_channel, modver_channel)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_tracer,  modver_tracer)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_timer,   modver_timer)
    IF (status /= 0) RETURN
    ! op_bk_20130820+
#ifndef __ICON__
    CALL put_submodel_att (status, .TRUE., modstr_qtimer,  modver_qtimer)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_import,  modver_import)
    IF (status /= 0) RETURN
#endif
    ! op_bk_20130820-

#if defined(ECHAM5) || defined(COSMO)
    ! put submodel version into output file
    ! MESSy SUBMODELS
    CALL put_submodel_att (status, USE_AEROPT,   modstr_aeropt,   modver_aeropt)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AIRSEA,   modstr_airsea,   modver_airsea)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_BUFLY,    modstr_bufly,    modver_bufly)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLOUD,    modstr_cloud,    modver_cloud)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CONVECT,  modstr_convect,  modver_convect)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CVTRANS,  modstr_cvtrans,  modver_cvtrans)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_D14CO,    modstr_d14co,    modver_d14co)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DDEP,     modstr_ddep,     modver_ddep)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DRADON,   modstr_dradon,   modver_dradon)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_GMXE,     modstr_gmxe,     modver_gmxe)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_GWAVE,    modstr_gwave,    modver_gwave)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_H2O,      modstr_h2o,      modver_h2o)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_HETCHEM,  modstr_hetchem,  modver_hetchem)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_JVAL,     modstr_jval,     modver_jval)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_LNOX,     modstr_lnox,     modver_lnox)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_M7,       modstr_m7,       modver_m7)       ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MADE,     modstr_made,     modver_made)       ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MECCA1,   modstr_mecca1,   modver_mecca1)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MECCA,    modstr_mecca,    modver_mecca)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MEGAN,    modstr_megan,    modver_megan)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MLOCEAN,  modstr_mlocean,  modver_mlocean)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MMFORCE,  modstr_mmforce,  modver_mmforce)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MSBM,     modstr_msbm,     modver_msbm)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_O3ORIG,   modstr_o3orig,   modver_o3orig)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OFFEMIS,  modstr_offemis,  modver_offemis)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ONEMIS,   modstr_onemis,   modver_onemis)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ORBIT,    modstr_orbit,    modver_orbit)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PHOTO,    modstr_photo,    modver_photo)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PLUMEGAS, modstr_plumegas, modver_plumegas)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PSC,      modstr_psc,      modver_psc)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_QBO,      modstr_qbo,      modver_qbo)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_RAD4ALL,  modstr_rad4all,  modver_rad4all)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SATSIMS,  modstr_satsims,  modver_satsims)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCALC,    modstr_scalc,    modver_scalc)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCAV,     modstr_scav,     modver_scav)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCOUT,    modstr_scout,    modver_scout)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SEDI,     modstr_sedi,     modver_sedi)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_S4D,      modstr_s4d,      modver_s4d)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SORBIT,   modstr_sorbit,   modver_sorbit)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SPACENOX, modstr_spacenox, modver_spacenox) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_E4CHEM,   modstr_e4chem,   modver_e4chem)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SPE,      modstr_spe,      modver_spe)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TIMEPOS,  modstr_timepos,  modver_timepos)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TNUDGE,   modstr_tnudge,   modver_tnudge)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TREXP,    modstr_trexp,    modver_trexp)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TROPOP,   modstr_tropop,   modver_tropop)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VAHR,     modstr_vahr,     modver_vahr)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VISO,     modstr_viso,     modver_viso)     ; IF (status/=0) RETURN
#endif
#ifdef BLANK
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD1,  modstr_submod1,  modver_submod1)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD2,  modstr_submod2,  modver_submod2)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD3,  modstr_submod3,  modver_submod3)  ; IF (status/=0) RETURN
#endif
#ifdef MBM_CMAT
#endif
    ! ### ADD HERE

    status = 0

  END SUBROUTINE switch_init
  ! -------------------------------------------------------------

  ! -------------------------------------------------------------
  SUBROUTINE put_submodel_att(status, LUSE, modstr, modver)
    
    ! MESSy
    USE messy_main_channel,       ONLY: new_attribute
    
    ! I/O
    INTEGER,      INTENT(OUT) :: status
    LOGICAL,      INTENT(in)  :: LUSE
    CHARACTER(*), INTENT(in)  :: modstr
    CHARACTER(*), INTENT(in)  :: modver
    
    status = 0

    IF (LUSE) THEN
       CALL new_attribute(status,  &
            'MESSy_'//TRIM(modstr), c= 'version '//TRIM(modver))
    END IF

  END SUBROUTINE put_submodel_att
  ! -------------------------------------------------------------
  ! op_pj_20110621-

! ***************************************************************************
END MODULE messy_main_switch_bi
! ***************************************************************************
