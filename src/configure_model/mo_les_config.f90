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
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_grid_config,         ONLY: is_plane_torus

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: les_config, configure_les  

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

    REAL(wp) :: ugeo(2)    ! ugeo(1)=constant, ugeo(2)=gradient
    REAL(wp) :: vgeo(2)    ! vgeo(1)=constant, vgeo(2)=gradient
    REAL(wp) :: ufric      ! friction velocity
 
    LOGICAL  :: is_dry_cbl  !special case for CBL testcase
    LOGICAL  :: set_geowind !TRUE is geostrophic wind is set

    !For isrf_type==3
    REAL(wp) :: bflux      !Buoyancy flux
    REAL(wp) :: tran_coeff !Surface transfer coefficient in units of velocity (m/s)

    !Some parameters
    REAL(wp) :: karman_constant
    REAL(wp) :: rkarman_constant  !inverse karman constant
    REAL(wp) :: smag_constant
    REAL(wp) :: turb_prandtl 
    REAL(wp) :: rturb_prandtl     !inverse turbulent prandtl number

  END TYPE t_les_config
  !>
  !!
  TYPE(t_les_config), TARGET :: les_config(max_dom)

  CONTAINS

  SUBROUTINE configure_les(jg)
  !--------------------------------------------------------------------------------------
  !  Set up LES parameters 
  !--------------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: jg !patch id

    CHARACTER(*), PARAMETER :: routine = "mo_les_config:configure_les:"

    !----------------------------------------------------
    ! Sanity check and Prints
    !----------------------------------------------------

    IF(les_config(jg)%isrfc_type==1)THEN

       les_config(jg)%shflx = 0._wp   
       les_config(jg)%lhflx = 0._wp   

       WRITE(message_text,'(a,e14.6)')'LES with fixed SST=',les_config(jg)%sst

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
    
    IF(les_config(jg)%set_geowind .AND. les_config(jg)%ugeo(1)==0._wp  &
                                  .AND. les_config(jg)%vgeo(1)==0._wp) &

      CALL message('mo_les_nml:WARNING:','Input values for Geostrophic wind are 0!')
   
    IF(les_config(jg)%set_geowind .AND. .NOT.is_plane_torus) &
      CALL finish(TRIM(routine),'set_geowind is only applicable for torus grid!')
 

  END SUBROUTINE configure_les

END MODULE mo_les_config
