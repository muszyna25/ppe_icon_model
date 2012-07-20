!>
!! @brief configuration setup for z/i/p-level output
!!
!! configuration setup for z/i/p-level output
!!
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-05)
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
MODULE mo_nh_pzlev_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for z/p-level output
  !!--------------------------------------------------------------------------
  TYPE :: t_nh_pzlev_config

    ! namelist variables
    !
    INTEGER :: nzlev                 !< number of z-levels

    INTEGER :: nplev                 !< number of p-levels

    INTEGER :: nilev                 !< number of isentropes

    REAL(wp):: zlevels(100)          !< zlevel heights [m] 

    REAL(wp):: plevels(100)          !< plevel heights [Pa] 

    REAL(wp):: ilevels(100)          !< isentropes [K] 

    ! derived variables
    !
    REAL(wp), POINTER ::      &
      &  p3d(:,:,:),          & !< 3D pressure level target field for output on p-levels
      &  z3d(:,:,:),          & !< 3D height level target field for output on z-levels
      &  i3d(:,:,:)             !< 3D theta level target field for output on isentropes

  END TYPE t_nh_pzlev_config

  !>
  !!
  TYPE(t_nh_pzlev_config), TARGET :: nh_pzlev_config(0:max_dom)


CONTAINS

  !>
  !! setup components for output on pressure/height levels and isentropes
  !!
  !! Setup of additional control variables for output on pressure/height levels 
  !! and isentropes.  
  !! These may depend on the nh_pzlev-namelist and potentially other namelists. 
  !! This routine is called, after all namelists have been read and a synoptic 
  !! consistency check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-09-07)
  !!
  SUBROUTINE configure_nh_pzlev( jg, nproma, npromz_c, nblks_c )
  !
    INTEGER, INTENT(IN) :: jg           !< patch 
    INTEGER, INTENT(IN) :: nproma
    INTEGER, INTENT(IN) :: npromz_c
    INTEGER, INTENT(IN) :: nblks_c      !< number of blocks

    ! Local variables
    INTEGER :: ist
    INTEGER :: nlen
    INTEGER :: z_nplev, z_nzlev, z_nilev
    INTEGER :: jb, jk           ! loop indices
    !-----------------------------------------------------------------------

    z_nplev = nh_pzlev_config(jg)%nplev
    z_nzlev = nh_pzlev_config(jg)%nzlev
    z_nilev = nh_pzlev_config(jg)%nilev


    ! allocate 3D pressure and z-level fields
    ALLOCATE(nh_pzlev_config(jg)%p3d(nproma,z_nplev,nblks_c),          &
      &      nh_pzlev_config(jg)%z3d(nproma,z_nzlev,nblks_c),          &
      &      nh_pzlev_config(jg)%i3d(nproma,z_nzlev,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( 'mo_nh_pzlev_nml: configure_nh_pzlev',       &
        &      'allocation of p3d, z3d, i3d failed' )
    ENDIF


    ! Fill z3d field of pressure-level data and pressure field of 
    ! height-level data
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1,nblks_c

      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

      DO jk = 1, nh_pzlev_config(jg)%nplev
        nh_pzlev_config(jg)%p3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%plevels(jk)
      ENDDO

      DO jk = 1, nh_pzlev_config(jg)%nzlev
        nh_pzlev_config(jg)%z3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%zlevels(jk)
      ENDDO

      DO jk = 1, nh_pzlev_config(jg)%nilev
        nh_pzlev_config(jg)%i3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%ilevels(jk)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    
  END SUBROUTINE configure_nh_pzlev

END MODULE mo_nh_pzlev_config
