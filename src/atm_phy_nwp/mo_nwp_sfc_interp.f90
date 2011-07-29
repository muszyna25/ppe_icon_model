!>
!! This module contains routines for the vertical interpolation of
!! surface/soil fields provided by external analyses to the ICON grid
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-29)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_nwp_sfc_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_config,          ONLY: msg_level
  USE mo_parallel_config,     ONLY: nproma 
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_prepicon_nml,        ONLY: nlevsoil_in, nlev_in
  USE mo_prepicon_utils,      ONLY: t_prepicon_state
  USE mo_lnd_nwp_config,      ONLY: nlev_soil

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: process_sfcfields

CONTAINS


  !-------------
  !>
  !! SUBROUTINE process_sfcfields
  !! Routine to convert surface fields interpolated horizontally by IFS2ICON
  !! to the ICON prognostic variables. Important ingredients are 
  !! - height adjustment of temperatures (partly done)
  !! - vertical interpolation of soil temperature and moisture (not yet done)
  !! - conversion of soil moisture information (not yet done)
  !! - height adjustment of snow cover information (not yet done)
  !!
  !! Other open items - just to document them somewhere
  !! - IFS2ICON needs to take into account land-sea-mask information for horizontal
  !!   interpolation of surface fields
  !! - Proper conversion of soil moisture may require information about soil types
  !!   and field capacity or similar things
  !! - And, the most complicated problem, lakes/islands not present at all in the
  !!   source data, is not yet addressed at all here!
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-20)
  !!
  !!
  SUBROUTINE process_sfcfields(p_patch, prepicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon


    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1
    INTEGER :: nlen, nlev

!-------------------------------------------------------------------------

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen)

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! 2D fields that do not require height adjustment
      ! these fields are simply copied
      DO jc = 1, nlen
        prepicon%sfc%skinres(jc,jb) = prepicon%sfc_in%skinres(jc,jb) 
        prepicon%sfc%ls_mask(jc,jb) = prepicon%sfc_in%ls_mask(jc,jb)
        IF (prepicon%sfc_in%seaice(jc,jb) >= 0._wp) THEN
          prepicon%sfc%seaice(jc,jb)  = prepicon%sfc_in%seaice(jc,jb) 
        ELSE
          prepicon%sfc%seaice(jc,jb)  = -0.1_wp 
        ENDIF
      ENDDO

      ! 2D fields that require height adjustment
      ! unfortunately, skin temperature is the only variable for which it is
      ! intuitively clear what to do ...
      DO jc = 1, nlen
        ! Adjust skin temperature with the difference between the atmospheric
        ! temperatures at the lowest model level
        prepicon%sfc%tskin(jc,jb)    = prepicon%sfc_in%tskin(jc,jb) +         &
          (prepicon%atm%temp(jc,nlev,jb) - prepicon%atm_in%temp(jc,nlev_in,jb))

        prepicon%sfc%tsnow(jc,jb)    = prepicon%sfc_in%tsnow(jc,jb) 
        prepicon%sfc%snowweq(jc,jb)  = prepicon%sfc_in%snowweq(jc,jb) 
        prepicon%sfc%snowdens(jc,jb) = prepicon%sfc_in%snowdens(jc,jb) 
      ENDDO

      ! 3D fields: soil temperature and moisture
      ! What is implemented here so far is definitely nonsense!
      ! Just for testing technical functionality...
      DO jk = 1, MIN(nlevsoil_in, nlev_soil)
        DO jc = 1, nlen
          prepicon%sfc%tsoil(jc,jb,jk) = prepicon%sfc_in%tsoil(jc,jb,jk)
          prepicon%sfc%wsoil(jc,jb,jk) = prepicon%sfc_in%wsoil(jc,jb,jk)
        ENDDO
      ENDDO

      jk1 = MIN(nlevsoil_in, nlev_soil)
      DO jk = jk1+1, nlev_soil+2
        DO jc = 1, nlen
          prepicon%sfc%tsoil(jc,jb,jk) = prepicon%sfc_in%tsoil(jc,jb,jk1)
          prepicon%sfc%wsoil(jc,jb,jk) = prepicon%sfc_in%wsoil(jc,jb,jk1)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO 
!$OMP END PARALLEL

  END SUBROUTINE process_sfcfields

END MODULE mo_nwp_sfc_interp
