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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nh_pzlev_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom
  USE mo_math_utilities,     ONLY: t_value_set
  USE mo_exception,          ONLY: finish
  USE mo_util_sort,          ONLY: quicksort
  USE mo_mpi,                ONLY: my_process_is_stdio

  IMPLICIT NONE
  PUBLIC



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for z/p-level output
  !!--------------------------------------------------------------------------
  TYPE :: t_nh_pzlev_config

    ! namelist variables
    !
    TYPE (t_value_set) :: zlevels    !< zlevel heights [m] 
    TYPE (t_value_set) :: plevels    !< plevel heights [Pa] 
    TYPE (t_value_set) :: ilevels    !< isentropes [K]

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

    z_nplev = nh_pzlev_config(jg)%plevels%nvalues
    z_nzlev = nh_pzlev_config(jg)%zlevels%nvalues
    z_nilev = nh_pzlev_config(jg)%ilevels%nvalues

    ! do status output
    !
    IF (((z_nplev > 0) .OR. (z_nzlev > 0) .OR. (z_nilev > 0)) .AND. (my_process_is_stdio())) THEN
      WRITE (0,'(a)')      " "
      WRITE (0,'(a,i0)') " Output on pressure/height levels and/or isentropes: domain ", jg
      IF (z_nplev > 0) THEN
        WRITE (0,'(a)')      " selected pressure levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%plevels%values(1:z_nplev)
        WRITE (0,'(a)')      " "
      END IF
      IF (z_nzlev > 0) THEN
        WRITE (0,'(a)')      " selected height levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%zlevels%values(1:z_nzlev)
        WRITE (0,'(a)')      " "
      END IF
      IF (z_nilev > 0) THEN
        WRITE (0,'(a)')      " selected isentropic levels: "
        WRITE (0,'(5(f10.2,","))')  nh_pzlev_config(jg)%ilevels%values(1:z_nilev)
        WRITE (0,'(a)')      " "
      END IF
    END IF

    ! allocate 3D pressure and z-level fields
    ALLOCATE(nh_pzlev_config(jg)%p3d(nproma,z_nplev,nblks_c),          &
      &      nh_pzlev_config(jg)%z3d(nproma,z_nzlev,nblks_c),          &
      &      nh_pzlev_config(jg)%i3d(nproma,z_nilev,nblks_c), STAT=ist )
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

      DO jk = 1, z_nplev
        nh_pzlev_config(jg)%p3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%plevels%values(jk)
      ENDDO

      DO jk = 1, z_nzlev
        nh_pzlev_config(jg)%z3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%zlevels%values(jk)
      ENDDO

      DO jk = 1, z_nilev
        nh_pzlev_config(jg)%i3d(1:nlen,jk,jb) = nh_pzlev_config(jg)%ilevels%values(jk)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    
  END SUBROUTINE configure_nh_pzlev

END MODULE mo_nh_pzlev_config
