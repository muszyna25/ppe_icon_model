!>
!! Namelist variables for assimilation schemes
!!        
!! @author Klaus Stephan, DWD
!!
!! @par Revision History
!! Initial revision by Klaus Stephan , DWD (2014-12-18)
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
MODULE mo_assimilation_nml

  USE mo_assimilation_config,     ONLY: assimilation_config

  USE mo_kind,                ONLY: wp,i4
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 

  USE mo_impl_constants,      ONLY: max_dom
  
  USE mo_master_config,       ONLY: isRestart
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,   &
                            &       open_and_restore_namelist, close_tmpfile
  USE mo_assimilation_config, ONLY: assimilation_config !,n_noobs
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_run_config,          ONLY: ldass_lhn

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_assimilation_namelist

  !---------------------------------------------------------------
  ! Namelist variables 
  !---------------------------------------------------------------
  
  LOGICAL                          ::           &
    llhn(max_dom)             ,& ! on/off switch for latent heat nudging (lhn)
    llhnverif(max_dom)        ,& ! on/off switch for verification against radar
    lhn_artif        ,& ! apply artificial profile
    lhn_filt         ,& ! vertical filtering of lhn t-increments
    lhn_relax        ,& ! horizontal filtering of lhn t-increments
    lhn_limit        ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim)
    lhn_limitp       ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim), restrict increments to Gaussian profile
    lhn_no_ttend     ,& ! do not apply the t increments (except for artif points), but only the q increments
    lhn_hum_adj      ,& ! apply a humidity adjustment along with t-increments
    lhn_incloud      ,& ! apply the LHN-scaling in cloudy layers only
    lhn_spqual       ,& ! switch for use of a spatial quality function
    lhn_black        ,& ! use blacklist for radar data
    lhn_qrs          ,& ! calculate the integrated precipitation flux
    lhn_logscale     ,& ! apply logarithmic scaling factors
    lhn_wweight      ,& ! apply a weighting with respect to the mean horizontal wind
    lhn_height       ,& ! use height infos for radar data
    lhn_bright       ,& ! apply bright band detection
    lhn_diag         ,& ! produce more detailed diagnostic output during lhn
    lhn_artif_only      ! apply only artificial temperature profile instead of applying modelled tt_lheat profile

  INTEGER ::  &
    nlhn_start       ,& ! start of latent heat nudging period in timesteps
    nlhn_end         ,& ! end of latent heat nudging period in timesteps
    nlhnverif_start  ,& ! start of latent heat nudging period in timesteps
    nlhnverif_end    ,& ! end of latent heat nudging period in timesteps
    nlhn_relax       ,& ! number of interations of horizontal filtering
    nradar(max_dom)     ! max. number of radar stations within input data

  REAL (KIND=wp)                   ::           &
    lhn_coef          ,& ! factor for reduction of lhn t-increments
    lhn_dt_obs        ,& ! time step of input data in minutes
    abs_lhn_lim       ,& ! absolute limit for lhn t-increments (used if lhn_limit or lhn_limitp)
    fac_lhn_artif     ,& ! factor when artificial profile will applied
    fac_lhn_up_in     ,& ! limiting factor for upscaling of model heating profile incoming
    fac_lhn_up        ,& ! limiting factor for upscaling of model heating profile
    fac_lhn_down_in   ,& ! limiting factor for downscaling model heating profile incoming
    fac_lhn_down      ,& ! limiting factor for downscaling model heating profile
    thres_lhn         ,& ! threshold of rain rates to be consinderd within lhn approach
    rqrsgmax          ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
    tt_artif_max      ,& ! maximum of latent heat used as artificial profile
    zlev_artif_max    ,& ! altidude of maximum of artificial profile
    std_artif_max     ,& ! parameter to define vertical width of artifical temperature profile
    start_fadeout         ! time relative to lhn_end, when the lhn coefficient in decreased toward 0

 CHARACTER (LEN=100)              ::           &
    radar_in             ,& ! directory for reading radar-files
    radardata_file(max_dom)       ,& ! filename of radar data
    blacklist_file(max_dom)       ,& ! filename of blacklist for radar data
    height_file(max_dom)             ! filename of radar beam heights

! CHARACTER (LEN=12)               ::           &
!    noobs_date (n_noobs)    ! array of missing observations

!  NAMELIST/assimilation_nml/  llhn         ,llhnverif                  ,           &
  NAMELIST/assimilation_nml/  nlhn_start   ,nlhn_end                   ,           &
                              nlhnverif_start ,nlhnverif_end           ,           &
                              lhn_coef, fac_lhn_up  ,fac_lhn_down      ,           &
                              thres_lhn    ,                                       &  ! noobs_date
                              rqrsgmax                   ,           &
                              radar_in     ,                                       &
                              lhn_black    ,blacklist_file             ,           &
                              lhn_artif    ,fac_lhn_artif              ,           &
                              lhn_filt     ,lhn_hum_adj, lhn_no_ttend  ,           &
                              lhn_limit    ,lhn_limitp  ,abs_lhn_lim   ,           &
                              lhn_relax    ,nlhn_relax                 ,           &
                              lhn_incloud  ,lhn_diag, lhn_qrs          ,           &
                              lhn_logscale ,lhn_wweight                ,           &
                              lhn_bright   ,lhn_height                 ,           &
                              lhn_artif_only   ,                                   &
                              height_file  ,lhn_spqual                 ,           &
                              lhn_dt_obs   ,nradar, radardata_file     ,           &
                              tt_artif_max  ,zlev_artif_max, std_artif_max,        &
                              start_fadeout 
CONTAINS
  !>
  !!
  SUBROUTINE read_assimilation_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    INTEGER :: jg
    CHARACTER(LEN=*),PARAMETER :: routine='mo_assimilation_nml:read_assimilation_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    llhn(:)               = ldass_lhn
    llhnverif(:)          = .TRUE.
    lhn_artif             = .TRUE.
    lhn_filt           = .TRUE.
    lhn_relax          = .FALSE.
    lhn_limit          = .TRUE.
    lhn_limitp         = .FALSE.
    lhn_hum_adj        = .TRUE.
    lhn_no_ttend       = .FALSE.
    lhn_spqual         = .FALSE.
    lhn_black          = .FALSE.
    lhn_incloud        = .TRUE.
    lhn_qrs            = .TRUE.
    lhn_logscale       = .TRUE.
    lhn_wweight        = .FALSE.
    lhn_diag           = .FALSE.
    lhn_artif_only     = .FALSE.
    lhn_height         = .FALSE.
    lhn_bright         = .FALSE.
    nlhn_start         = -9999
    nlhn_end           = -9999
    nlhnverif_start    = -9999
    nlhnverif_end      = -9999
    nradar(:)          = 200
    nlhn_relax         = 2_i4
    lhn_dt_obs         = 300.0_wp
    lhn_coef           = 1.0_wp
    abs_lhn_lim        = 50._wp / 3600._wp    ! max. change in heating: 4K/h
    fac_lhn_artif      = 5._wp
    fac_lhn_up         = 2.0_wp
    fac_lhn_down       = 1.0_wp / 2.0_wp
    thres_lhn          = 0.1_wp / 3600._wp
    radar_in           = './'
    radardata_file(:)     = 'radardata.g2'
    blacklist_file(:)     = 'blacklist_dx.g2'
    height_file(:)        = 'height_dx.g2'
!    noobs_date  (:)    = '            '
    rqrsgmax           = 1.0_wp
    tt_artif_max       = 0.0015
    zlev_artif_max     = 1000.
    std_artif_max      = 4.
    start_fadeout      = 1.0

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('assimilation_nml')
      READ(funit,NML=assimilation_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('assimilation_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, assimilation_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, assimilation_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, assimilation_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=assimilation_nml)
      CALL store_and_close_namelist(funit, 'assimilation_nml')
    ENDIF
    
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=assimilation_nml)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    DO jg= 1,max_dom
        assimilation_config(jg)%llhn            = llhn(jg)
        assimilation_config(jg)%llhnverif       = llhnverif(jg)
        assimilation_config(jg)%lhn_artif       = lhn_artif
        assimilation_config(jg)%lhn_filt        = lhn_filt
        assimilation_config(jg)%lhn_relax       = lhn_relax
        assimilation_config(jg)%lhn_limit       = lhn_limit
        assimilation_config(jg)%lhn_limitp      = lhn_limitp
        assimilation_config(jg)%lhn_hum_adj     = lhn_hum_adj
        assimilation_config(jg)%lhn_no_ttend    = lhn_no_ttend
        assimilation_config(jg)%lhn_diag        = lhn_diag
        assimilation_config(jg)%lhn_qrs         = lhn_qrs
        assimilation_config(jg)%lhn_logscale    = lhn_logscale
        assimilation_config(jg)%lhn_wweight     = lhn_wweight
        assimilation_config(jg)%lhn_spqual      = lhn_spqual
        assimilation_config(jg)%lhn_black       = lhn_black
        assimilation_config(jg)%lhn_incloud     = lhn_incloud
        assimilation_config(jg)%lhn_artif_only  = lhn_artif_only
        assimilation_config(jg)%lhn_height      = lhn_height
        assimilation_config(jg)%lhn_bright      = lhn_bright
        assimilation_config(jg)%nlhn_start      = nlhn_start
        assimilation_config(jg)%nlhn_end        = nlhn_end
        assimilation_config(jg)%nlhnverif_start = nlhnverif_start
        assimilation_config(jg)%nlhnverif_end   = nlhnverif_end
        assimilation_config(jg)%nradar          = nradar(jg)
        assimilation_config(jg)%nlhn_relax      = nlhn_relax
        assimilation_config(jg)%lhn_dt_obs      = lhn_dt_obs/60.
        assimilation_config(jg)%lhn_coef        = lhn_coef
        assimilation_config(jg)%abs_lhn_lim     = abs_lhn_lim
        assimilation_config(jg)%fac_lhn_artif   = fac_lhn_artif
        assimilation_config(jg)%fac_lhn_up      = fac_lhn_up
        assimilation_config(jg)%fac_lhn_down    = fac_lhn_down
        assimilation_config(jg)%thres_lhn       = thres_lhn
        assimilation_config(jg)%rqrsgmax        = rqrsgmax
        assimilation_config(jg)%tt_artif_max    = tt_artif_max
        assimilation_config(jg)%zlev_artif_max  = zlev_artif_max
        assimilation_config(jg)%std_artif_max   = std_artif_max
        assimilation_config(jg)%start_fadeout   = start_fadeout
        assimilation_config(jg)%radar_in        = radar_in
        assimilation_config(jg)%radardata_file  = radardata_file(jg)
        assimilation_config(jg)%blacklist_file  = blacklist_file(jg)
        assimilation_config(jg)%height_file     = height_file(jg)
    ENDDO 

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=assimilation_nml)                    
      CALL store_and_close_namelist(funit, 'assimilation_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=assimilation_nml)


  END SUBROUTINE read_assimilation_namelist
  !-------------

END MODULE mo_assimilation_nml
