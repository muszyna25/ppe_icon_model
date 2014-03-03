!>
!! Contains the setup of variables related to large eddy simulation setup
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_les_nml

  USE mo_les_config,          ONLY: les_config
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_les_namelist
  PUBLIC :: turb_profile_list, turb_tseries_list

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  REAL(wp) :: sst        ! prescribed SST
  REAL(wp) :: shflx      ! prescribed sensible heat flux (Km/s)
  REAL(wp) :: lhflx      ! prescribed latent heat flux   (Km/s)
  INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux, 3=fixed buyancy flux

  REAL(wp) :: ufric      ! friction velocity
 
  LOGICAL  :: is_dry_cbl  !special case for CBL testcase
 
  !For isrf_type==3
  REAL(wp) :: bflux      !Buoyancy flux
  REAL(wp) :: tran_coeff !Surface transfer coefficient in units of velocity (m/s)

  !Some parameters
  REAL(wp) :: smag_constant
  REAL(wp) :: turb_prandtl 

  !Scheme for vertical discretization
  INTEGER :: vert_scheme_type !1=explicit, 2=implicit

  !Parameters for additional diagnostic output
  LOGICAL  :: ldiag_les_out                    !.TRUE. to turn it on
  REAL(wp) :: avg_interval_sec, sampl_freq_sec !averaging and sampling time 
  CHARACTER(LEN=7) :: turb_tseries_list(9), turb_profile_list(28) !list of variables  
  CHARACTER(MAX_CHAR_LENGTH) :: expname        !name of experiment for naming the file

  NAMELIST/les_nml/ sst, shflx, lhflx, isrfc_type, ufric, is_dry_cbl, &
                    smag_constant, turb_prandtl, bflux, tran_coeff,   &
                    vert_scheme_type, avg_interval_sec, sampl_freq_sec,  &
                    expname, ldiag_les_out

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for LES
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  !! @par Revision History
  !!  by Anurag Dipankar, MPIM (2013-04)
  !!
  SUBROUTINE read_les_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, jg
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_les_nml: read_les_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    sst          = 300._wp
    shflx        = -999._wp 
    lhflx        = -999._wp 
    isrfc_type   = 1 
    ufric        = -999._wp 

    is_dry_cbl   = .FALSE.

    !parameters
    smag_constant    = 0.23_wp
    turb_prandtl     = 0.33333333333_wp

    bflux       = -999._wp
    tran_coeff  = -999._wp

    vert_scheme_type = 2 !implicit

    !output parameters
    ldiag_les_out = .TRUE. 
    expname  = 'ICOLES'
    avg_interval_sec = 900._wp
    sampl_freq_sec   = 60._wp

    turb_profile_list = (/                                                     &
      'u      ','v      ','w      ','th     ','exner  ','rho    ','qv     ',   & !1-7
      'qc     ','wu     ','wv     ','wth    ','wqv    ','wqc    ','ww     ',   & !8-14
      'thth   ','qvqv   ','qcqc   ','uu     ','vv     ','kh     ','km     ',   & !15-21
      'thv    ','wthv   ','wqvd   ','wthd   ','wqcd   ','bynprd ','mechprd' /)   !22-28

    turb_tseries_list = (/                                          &
      'ccover ','shflx  ','lhflx  ','ustress','vstress','tsfc   ',  & !1-6
      'qsfc   ','hbl    ','psfc   '  /)                               !7-9


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('les_nml')
      READ(funit,NML=les_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('les_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, les_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, les_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, les_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1 , max_dom
      les_config(jg)% sst          =  sst
      les_config(jg)% shflx        =  shflx
      les_config(jg)% lhflx        =  lhflx
      les_config(jg)% isrfc_type   =  isrfc_type
      les_config(jg)% ufric        =  ufric
      les_config(jg)% is_dry_cbl   =  is_dry_cbl
      les_config(jg)% smag_constant     =  smag_constant
      les_config(jg)% turb_prandtl      =  turb_prandtl
      les_config(jg)% rturb_prandtl     =  1._wp/turb_prandtl
      les_config(jg)% bflux             =  bflux
      les_config(jg)% tran_coeff        =  tran_coeff
      les_config(jg)% vert_scheme_type  =  vert_scheme_type
      les_config(jg)% ldiag_les_out     =  ldiag_les_out
      les_config(jg)% expname           =  expname
      les_config(jg)% avg_interval_sec  =  avg_interval_sec
      les_config(jg)% sampl_freq_sec    =  sampl_freq_sec
    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=les_nml)                    
      CALL store_and_close_namelist(funit,'les_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=les_nml)

  END SUBROUTINE read_les_namelist

END MODULE mo_les_nml
