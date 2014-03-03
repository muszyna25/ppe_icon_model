!>
!! Contains the setup of variables related to large eddy simulation setup
!!
!! @Anurag Dipankar, MPIM (2013-04)
!!
!!
!! @par Revision History
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
MODULE mo_les_config

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: max_dom, MAX_CHAR_LENGTH
  USE mo_grid_config,         ONLY: is_plane_torus
  USE mo_math_constants,      ONLY: dbl_eps

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_les_config, les_config, configure_les

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for LES with or without TORUS grid
  !--------------------------------------------------------------------------
  TYPE t_les_config

    ! variables from namelist
    REAL(wp) :: sst        ! prescribed SST
    REAL(wp) :: shflx      ! prescribed sensible heat flux (W/m2)
    REAL(wp) :: lhflx      ! prescribed latent heat flux   (W/m2)
    INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux

    REAL(wp) :: ufric      ! friction velocity
 
    LOGICAL  :: is_dry_cbl  !special case for CBL testcase

    !For isrf_type==3
    REAL(wp) :: bflux      !Buoyancy flux
    REAL(wp) :: tran_coeff !Surface transfer coefficient in units of velocity (m/s)

    !Some parameters
    REAL(wp) :: smag_constant
    REAL(wp) :: turb_prandtl 
    REAL(wp) :: rturb_prandtl     !inverse turbulent prandtl number

    !Scheme for vertical discretization
    INTEGER :: vert_scheme_type !1=explicit, 2=implicit

    !Parameters for additional diagnostic output
    LOGICAL  :: ldiag_les_out                    !.TRUE. to turn it on
    REAL(wp) :: avg_interval_sec, sampl_freq_sec !averaging and sampling time 
    CHARACTER(MAX_CHAR_LENGTH) :: expname        !name of experiment for naming the file

  END TYPE t_les_config
  !>
  !!
  TYPE(t_les_config), TARGET :: les_config(max_dom)

  CONTAINS

  SUBROUTINE configure_les(jg, dtime_adv)
  !--------------------------------------------------------------------------------------
  !  Set up LES parameters 
  !--------------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: jg !patch id
    REAL(wp),INTENT(IN) :: dtime_adv

    CHARACTER(*), PARAMETER :: routine = "mo_les_config:configure_les:"

    !----------------------------------------------------
    ! Sanity check and Prints
    !----------------------------------------------------

    IF(les_config(jg)%isrfc_type==1)THEN

       les_config(jg)%shflx = 0._wp   
       les_config(jg)%lhflx = 0._wp   

       WRITE(message_text,'(a,e14.6)')'LES with surface scheme TERRA'

       CALL message(TRIM(routine),message_text)

    ELSEIF(les_config(jg)%isrfc_type==2)THEN

       WRITE(message_text,'(a,e14.6,e14.6)')'LES with fixed fluxes=', &
           les_config(jg)%shflx,les_config(jg)%lhflx

       CALL message(TRIM(routine),message_text)

       IF(les_config(jg)%shflx==-999._wp .OR. les_config(jg)%lhflx==-999._wp) &
          CALL finish(TRIM(routine),'Wrong input for irsfc_type=2')

    ELSEIF(les_config(jg)%isrfc_type==3)THEN

       WRITE(message_text,'(a,e14.6,e14.6)') 'LES with fixed Buoyancy flux and tran coeff=', &
               les_config(jg)%bflux,les_config(jg)%tran_coeff

       CALL message(TRIM(routine),message_text)

       IF(les_config(jg)%bflux==-999._wp .OR. les_config(jg)%tran_coeff==-999._wp) &
          CALL finish(TRIM(routine),'Wrong input for irsfc_type=3')

    END IF
  
    IF(les_config(jg)%is_dry_cbl)THEN
       les_config(jg)%lhflx = 0._wp
    END IF

    !Adjust sampling and frequency time
    IF(les_config(jg)%ldiag_les_out)THEN

     IF( MOD(les_config(jg)%avg_interval_sec, dtime_adv) > 10._wp*dbl_eps )  THEN
      WRITE(message_text,'(a,2F13.4)') &
        &'WARNING: averaging interval and advective timesteps not synchronized: ', &
        & les_config(jg)%avg_interval_sec, dtime_adv
      CALL message(TRIM(routine), TRIM(message_text))

      les_config(jg)%avg_interval_sec =   &
           REAL((FLOOR(les_config(jg)%avg_interval_sec/dtime_adv) + 1),wp) * dtime_adv

      WRITE(message_text,'(a,2F13.4)') &
          'implicit synchronization in configure_les: avg_interval_sec =', &
          les_config(jg)%avg_interval_sec
      CALL message(TRIM(routine), TRIM(message_text))
     END IF

     IF( MOD(les_config(jg)%sampl_freq_sec, dtime_adv) > 10._wp*dbl_eps )  THEN
      WRITE(message_text,'(a,2F13.4)') &
        &'WARNING: sampling frequency and advective timesteps not synchronized: ', &
        & les_config(jg)%sampl_freq_sec, dtime_adv
      CALL message(TRIM(routine), TRIM(message_text))

      les_config(jg)%sampl_freq_sec =   &
           REAL((FLOOR(les_config(jg)%sampl_freq_sec/dtime_adv) + 1),wp) * dtime_adv

      WRITE(message_text,'(a,2F13.4)') &
         'implicit synchronization in configure_les: sampl_freq_sec =', &
         les_config(jg)%sampl_freq_sec
      CALL message(TRIM(routine), TRIM(message_text))
     END IF

    ENDIF
     
 
  END SUBROUTINE configure_les

END MODULE mo_les_config
