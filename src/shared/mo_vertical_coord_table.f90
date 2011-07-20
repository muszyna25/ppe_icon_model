!>
!! module *mo_vertical_coord_table* - *loop indices and surface-pressure independent
!! variables associated with the vertical finite-difference scheme.
!!
!! @par Revision History
!!  A.J. Simmons, ECMWF (1981-11-16)
!!  H. Wan, MPI-Met (2006-02) adapted from ECHAM5.3.01
!!  H. Wan, MPI-Met (2007-07-19)
!!   - calling of *message* removed when something goes wrong.
!!   - no longer initialize the parameter arrays with infinity.
!!   - changed the name of this module from mo_hyb to mo_hyb_params.
!!  H. Wan, MPI-Met (2007-08)
!!   - parameters used only for the semi-implicit correction were moved to
!!     module mo_si_correction.
!!   - inihyb renamed init_hyb_params.
!!  A. Gassmann, MPI-Met (2008-04)
!!   - read hyb_file according to level number given
!!  M.A. Giorgetta, MPI-Met (2009-02-08)
!!   - change input and output formats of the vertical coordinate table
!!  Almut Gassmann, MPI-M (2009-03-19)
!!   - make read_hyb_params public for using it in nonhydrostatic version
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_vertical_coord_table

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !

  ! USE mo_parameters

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max, find_next_free_unit
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_impl_constants,     ONLY: success, max_char_length, ishallow_water, &
                                   ihs_atm_temp, ihs_atm_theta, inh_atmosphere
  USE mo_physical_constants, ONLY: grav, rcpd, rd

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: alloc_vct, init_vct

  INTEGER :: nlevm1         ! (number of levels)-1.
  INTEGER :: nplev          ! *number of pressure levels.
  INTEGER :: nplvp1         ! *nplev+1.*
  INTEGER :: nplvp2         ! *nplev+2.*
  INTEGER :: nplvpa         ! *nplvp1,* or 2 if *nplev=0.*
  INTEGER :: nlmsgl         ! *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         ! *nlmsgl+1.*
  INTEGER :: nlmsla         ! *nlmslp,* or 2 if *nlmslp=1.*

  REAL(wp) :: apzero        ! *reference pressure for computation of the
  !                         !  hybrid vertical levels.
  REAL(wp) :: t0icao        ! *surface temperatur of reference atmosphere
  REAL(wp) :: tsticao       ! *stratospheric temperature of reference atmosphere
  REAL(wp) :: rdtstic       ! *rd*tsticao
  REAL(wp) :: rdlnp0i       ! *rd*ln(surface pressure) of reference atmosphere
  REAL(wp) :: alrrdic       ! *lapse-rate parameter of reference atmosphere
  REAL(wp) :: rdt0ral       ! *rd*t0icao/alphaic
  REAL(wp) :: ptricao       ! *tropopause pressure of reference atmosphere
  REAL(wp) :: rdlnpti       ! *rd*ln(ptricao)
  REAL(wp) :: gsticao       ! *constant used in geopotential calculation

  REAL(wp), ALLOCATABLE :: ralpha(:) ! rd*alpha at pressure and sigma levels.
  REAL(wp), ALLOCATABLE :: rlnpr(:)  ! rd*ln(p(k+.5)/p(k-.5))
  REAL(wp), ALLOCATABLE :: dela(:)   ! a(k+.5)-a(k-.5).
  REAL(wp), ALLOCATABLE :: delb(:)   ! b(k+.5)-b(k-.5).
  REAL(wp), ALLOCATABLE :: rddelb(:) ! rd*delb.
  REAL(wp), ALLOCATABLE :: cpg(:)    ! a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL(wp), ALLOCATABLE :: delpr(:)  ! p(k+.5)-p(k-.5) of
                                     ! the refrence surface pressure.
  REAL(wp), ALLOCATABLE :: rdelpr(:) ! *reciprocal of *delpr.*
  REAL(wp), ALLOCATABLE :: alpham(:) ! *constant array for use by dyn.
  REAL(wp), ALLOCATABLE :: ardprc(:) ! *constant array for use by dyn.
  REAL(wp), ALLOCATABLE :: ceta(:)   ! *full hybrid vertical levels.
  REAL(wp), ALLOCATABLE :: cetah(:)  ! *half hybrid vertical levels.

  REAL(wp), ALLOCATABLE :: vct_a(:) ! param. A of the vertical coordinte
  REAL(wp), ALLOCATABLE :: vct_b(:) ! param. B of the vertical coordinate
  REAL(wp), ALLOCATABLE :: vct  (:) ! param. A and B of the vertical coordinate


  ! calling tree
  !
  ! (hydro_atmos) -> *init_vertical_coord* |- *read_hyb_params*
  !                                        |- *alloc_hyb*
  !                                        |- *init_hyb_params*


CONTAINS

  !-------------------------------------------------------------------------
  !BOC

  !EOC
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  init_vertical_coord
  !
  ! !SUBROUTINE INTERFACE:

  SUBROUTINE init_vertical_coord_table(iequations,klev)

    ! !DESCRIPTION:
    !  Initialization of the hybrid vertical coordinate
    !
    ! !REVISION HISTORY:
    !  Original version by Hui Wan, MPI-M, 2006-02-09
    !
    INTEGER, INTENT(IN) :: klev        !< number of ful levels
    INTEGER, INTENT(IN) :: iequations
    INTEGER :: jk
    !EOP
    !-----------------------------------------------------------------------
    !BOC

    ! read the A and B parameters of the vertical coordinate

    CALL read_vct(iequations,klev)

    CALL message('vertical_coord_table:init_vertical_coord', '')
    CALL message('', 'Vertical coordinate table')
    CALL message('', '   k     vct_a(k) [Pa]   vct_b(k) []')
    DO jk = 1, SIZE(vct_a)
      WRITE(message_text,'(i4,f18.10,f14.10)')  jk, vct_a(jk), vct_b(jk)
      CALL message('', TRIM(message_text))
    ENDDO

    ! allocate memory for the auxiliary parameters and arrays

    CALL alloc_vct(klev)

    ! assign values to the the auxiliary parameters and arrays

    CALL init_vct(klev)

  END SUBROUTINE init_vertical_coord_table

  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !! Read the A and B parameters of the hybrid vertical grid,.
  !!
  !! Read the A and B parameters of the hybrid vertical grid,
  !! which define the half level pressure: ph=A+B*ps [Pa]
  !!
  SUBROUTINE  read_vct (iequations,klev)

    INTEGER, INTENT(IN) :: klev
    INTEGER, INTENT(IN) :: iequations

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
         &   'mo_vertical_coord_table:read_vct'
    CHARACTER(len=filename_max)              :: vct_file, line

    INTEGER :: ist, iunit

    INTEGER :: ik, jk

    !-------------------------------------------------------------------------
    !BOC

    ! allocate memory

    ALLOCATE(vct_a(klev+1), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of vct_a failed')
    ENDIF

    ALLOCATE(vct_b(klev+1), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of vct_b failed')
    ENDIF

    ALLOCATE(vct((klev+1)*2), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of vct failed')
    ENDIF

    ! Open file
    WRITE(line,FMT='(i4)') klev
    !
    SELECT CASE(iequations)
    CASE(ishallow_water,ihs_atm_temp,ihs_atm_theta)
      ! use hybrid sigma pressure tables
      vct_file = 'atm_hyb_sp_'//TRIM(ADJUSTL(line))
    CASE(inh_atmosphere)
      ! use hybrid sigma height tables
      vct_file = 'atm_hyb_sz_'//TRIM(ADJUSTL(line))
    CASE DEFAULT
      CALL finish (TRIM(routine), 'no vct file name defined for specified iequations')
    END SELECT

    iunit = find_next_free_unit(10,20)
    OPEN (unit=iunit,file=TRIM(vct_file),access='SEQUENTIAL', &
      &  form='FORMATTED', IOSTAT=ist)

    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'open vertical coordinate table file failed')
    ENDIF

    ! Skip header line
    READ (iunit,*,IOSTAT=ist) line
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'reading header line failed')
    ENDIF

    ! Read A and B
    DO jk=1,klev+1
       READ (iunit,*,IOSTAT=ist) ik, vct_a(jk), vct_b(jk)
       IF(ist/=success)THEN
          CALL finish (TRIM(routine), 'reading vct_a and vct_b failed')
       ENDIF
    END DO

    CLOSE(iunit)

    vct(     1: klev+1   ) = vct_a(:)
    vct(klev+2:(klev+1)*2) = vct_b(:)

  END SUBROUTINE  read_vct

  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !!
  SUBROUTINE alloc_vct(klev)

    INTEGER, INTENT(IN) :: klev
    INTEGER :: ist
    INTEGER :: klevp1

    CHARACTER(len=max_char_length),PARAMETER :: routine = &
         & 'mo_vertical_coord_table:alloc_vct'

    !-----------------------------------------------------------------------
    !BOC

    klevp1 = klev+1

    ALLOCATE (ralpha(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ralpha failed')
    ENDIF

    ALLOCATE (rlnpr(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rlnpr failed')
    ENDIF

    ALLOCATE (dela(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of dela failed')
    ENDIF

    ALLOCATE (delb(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of delb failed')
    ENDIF

    ALLOCATE (rddelb(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rddelb failed')
    ENDIF

    ALLOCATE (cpg(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of cpg failed')
    ENDIF

    ALLOCATE (delpr(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (rdelpr(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (alpham(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of alpham failed')
    ENDIF

    ALLOCATE (ardprc(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ardprc failed')
    ENDIF

    ALLOCATE (ceta(klev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ceta failed')
    ENDIF

    ALLOCATE (cetah(klevp1), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of cetah failed')
    ENDIF

  END SUBROUTINE alloc_vct

  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !!  Initializes constants for vertical coordinate calculations.
  !!
  !!  Method:
  !!    Compute loop indices and surface-pressure independent
  !!    variables associated with the vertical finite-difference scheme.
  !!    Output is in module *mo_hyb*
  !!
  !! @par Revision History
  !!  A. J. Simmons, ECMWF, November 1981, original source
  !!  L. Kornblueh, MPI, May 1998, f90 rewrite
  !!  U. Schulzweida, MPI, May 1998, f90 rewrite
  !!  A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
  !!  H. Wan, MPI-M, Feb 2006, rewrite: "goto" removed
  !!  H. Wan, MPI-M, Aug 2007, new name: init_hyb_params
  !!  A. Gassmann, MPI-M, (2008-04-23), change apzero to 10^5Pa
  !! @par
  !!  for more details see file AUTHORS
  !!
  SUBROUTINE init_vct(klev)

    INTEGER, INTENT(IN) :: klev

    !  Local scalars:
    REAL(wp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER  :: ilev, klevp1, ilevp1, iplev, iplvp1, is, ism, ist, &
      &         jk, jlev

    !  Intrinsic functions
    INTRINSIC EXP, LOG

    !-----------------------------------------------------------------------
    !BOC

    !  Executable statements

    !-- 1. Initialize variables

!ag    apzero    = 101325._wp ! changed for NCAR summer colloquium!
    apzero    = 100000._wp
    zrd       = rd
    ralpha(1) = zrd*LOG(2._wp)
    rlnpr(1)  = 2._wp*ralpha(1)
    ilev      = klev
    ilevp1    = ilev + 1
    klevp1    = ilevp1
    nlevm1    = ilev - 1
    iplev     = 0
    iplvp1    = 1
    is        = klevp1 + ilevp1
    ism       = is - 1
    zpp       = vct(1)
    zsm       = vct(is)

    t0icao  = 288._wp
    tsticao = 216.5_wp
    zp0icao = 101320._wp
    rdlnp0i = rd*LOG(zp0icao)
    rdtstic = rd*tsticao
    alrrdic = 0.0065_wp/grav
    rdt0ral = t0icao/alrrdic
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic
    ptricao = EXP(rdlnpti/rd)
    gsticao = tsticao*(rdlnpti-1._wp/alrrdic)

    zb      = vct(klevp1+iplvp1+1)

    !-- 2. Calculate pressure-level values

    DO WHILE ( zb == 0._wp )

      iplev  = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) EXIT    ! if all levels are pressure levels

      zp            = zpp
      zpp           = vct(iplvp1)
      delpr(iplev)  = zpp - zp
      rdelpr(iplev) = 1._wp/delpr(iplev)

      IF ( iplev>1 ) THEN
        rlnpr(iplev)  = zrd*LOG(zpp/zp)
        ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
      END IF

      alpham(iplev) = ralpha(iplev)*rcpd
      ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
      zb            = vct(klevp1+iplvp1+1)

    ENDDO


    IF (iplvp1/=ilevp1) THEN   ! All levels are not pressure-levels

      nplev  = iplev
      nplvp1 = iplvp1
      nplvp2 = iplvp1 + 1
      IF (iplev==0) THEN
        nplvpa = 2
      ELSE
        nplvpa = iplvp1
      END IF

      !-- 3. Calculate sigma-level values

      za = vct(ism-klevp1)

      DO WHILE ( za == 0._wp )

        is  = ism
        ism = is - 1
        ist = is - klevp1
        zs  = zsm
        zsm = vct(is)
        IF (ist==1) THEN
          nlmsgl = 0
          nlmslp = 1
          nlmsla = 2
          EXIT
        ELSE
          rlnpr(ist)  = zrd*LOG(zs/zsm)
          ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
        END IF
        za = vct(ism-klevp1)

      END DO

      IF (za>0._wp) THEN
        nlmsgl = ism - klevp1
        nlmslp = nlmsgl + 1
        nlmsla = nlmslp
      END IF

      !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

      DO jk = 1, klev
        dela(jk)   = vct(jk+1) - vct(jk)
        delb(jk)   = vct(klevp1+jk+1) - vct(klevp1+jk)
        rddelb(jk) = rd*delb(jk)
        cpg(jk)    = vct(klevp1+jk)*vct(jk+1) - vct(klevp1+jk+1)*vct(jk)
      END DO

      DO jk = nlmslp, klev
        alpham(jk) = ralpha(jk)*delb(jk)
      END DO

    ENDIF  ! If all levels are not pressure-levels

    !-- 5. Compute full level values of the hybrid coordinate

    zetam    = vct(1)/apzero + vct(klevp1+1)
    cetah(1) = zetam

    DO jlev = 1, klev
      zetap         = vct(jlev+1)/apzero + vct(klevp1+1+jlev)
      ceta(jlev)    = (zetam+zetap)*.5_wp
      cetah(jlev+1) = zetap
      zetam = zetap
    END DO

  END SUBROUTINE init_vct

END MODULE mo_vertical_coord_table

