! ************************************************************************
MODULE messy_main_tools
! ************************************************************************

  ! MESSy-SMCL tools
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, 2003-...
  !         Hella Riede    , MPICH, Mainz, 2009-...
  !         Rolf Sander,     MPICH, Mainz, 2006-...

  USE messy_main_constants_mem, ONLY: SP, DP, STRLEN_LONG, STRLEN_VLONG, &
                                      STRLEN_MEDIUM

  IMPLICIT NONE
  PRIVATE

  ! mz_rs_20081222+
  ! 1D array also for integer:
  TYPE PTR_1D_ARRAY_INT
     INTEGER, DIMENSION(:), POINTER :: PTR
  END TYPE PTR_1D_ARRAY_INT
  ! mz_rs_20081222-

  TYPE PTR_1D_ARRAY
     REAL(DP), DIMENSION(:), POINTER :: PTR
  END TYPE PTR_1D_ARRAY

  TYPE PTR_2D_ARRAY
     REAL(DP), DIMENSION(:,:), POINTER :: PTR
  END TYPE PTR_2D_ARRAY

  TYPE PTR_3D_ARRAY
     REAL(DP), DIMENSION(:,:,:), POINTER :: PTR
  END TYPE PTR_3D_ARRAY

  TYPE PTR_4D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:), POINTER :: PTR
  END TYPE PTR_4D_ARRAY

  TYPE PTR_5D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: PTR
  END TYPE PTR_5D_ARRAY

  PUBLIC :: PTR_1D_ARRAY, PTR_2D_ARRAY, PTR_3D_ARRAY &
          , PTR_4D_ARRAY, PTR_5D_ARRAY, PTR_1D_ARRAY_INT

  ! mz_ht_20040414+
  ! FOR LOOKUP TABLE INITIALIZATION (CONVECTION at al.)
  ! lower and upper bound in K * 1000
  INTEGER,  PARAMETER, PUBLIC :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER,  PARAMETER, PUBLIC :: jptlucu2 = 400000  ! lookup table upper bound
  ! table - e_s*Rd/Rv
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucua(jptlucu1:jptlucu2)
  ! table - for derivative calculation
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucub(jptlucu1:jptlucu2)
  ! table - l/cp
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuc(jptlucu1:jptlucu2)
  ! table
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuaw(jptlucu1:jptlucu2)
  ! mz_ht_20040414-

  ! mz_rs_20060110+
  INTERFACE str
     MODULE PROCEDURE str_logical
     MODULE PROCEDURE str_integer
     MODULE PROCEDURE str_real_sp
     MODULE PROCEDURE str_real_dp
  END INTERFACE
  PUBLIC :: str
  ! mz_rs_20060110-

  ! mz_hr_20100704+ 
  ! psat and psatf replaced by psat_mk, converging to zero,
  ! always > 0
  PUBLIC :: psat_mk  ! Murphy-Koop 2005 sat. vap. press ice or liquid [Pa]
  PUBLIC :: cair_trad ! cair based on trad. definition of relhum
  PUBLIC :: cair_wmo ! cair based on WMO definition of relhum
  ! mz_hr_20100704- 

  INTERFACE iso2ind
     MODULE PROCEDURE iso2ind_1d
     MODULE PROCEDURE iso2ind_2d
  END INTERFACE

  INTERFACE ind2val
     MODULE PROCEDURE ind2val_1d
     MODULE PROCEDURE ind2val_2d
  END INTERFACE

! mz_ab_20100610+
  INTERFACE remap_bounds
     MODULE PROCEDURE remap_bounds1
     MODULE PROCEDURE remap_bounds2
     MODULE PROCEDURE remap_bounds3
     MODULE PROCEDURE remap_bounds4
  END INTERFACE
! mz_ab_20100610-

  ! mz_bk_20110707+
  INTERFACE str2num
     MODULE PROCEDURE str2num_real_dp
     MODULE PROCEDURE str2num_real_sp
     MODULE PROCEDURE str2num_integer
  END INTERFACE
  ! mz_bk_20110707-

  ! SUBROUTINES
  PUBLIC :: read_nml_open       ! Utilities ...
  PUBLIC :: read_nml_check      ! ... to simplify ...
  PUBLIC :: read_nml_close      ! ... namelist input

  PUBLIC :: iso2ind             ! find index
  PUBLIC :: ind2val             ! find value at index level
  PUBLIC :: int2str             ! convert integer to string
  PUBLIC :: strcrack            ! cracking strings into parts
  PUBLIC :: nn_index            ! look for nearest neighbour(s) in ordered list

  ! mz_ab_20100624+
  PUBLIC :: fliparray           ! reverses the order of 1-D array
  ! mz_ab_20100624-

  ! mz_ap_20070913+
  PUBLIC :: ns_index            ! look for sourrounding neighbour(s) in list
  ! mz_ap_20070913-

  ! mz_ht_20040414+
  PUBLIC :: init_convect_tables ! lookup table for convection et al. 
  ! mz_ht_20040414-

  PUBLIC :: match_wild          ! compare strings with wildcards
  PUBLIC :: str2chob            ! convert string to channel/object list
  PUBLIC :: bilin_weight        ! weights for bilinear interpolation

  ! mz_hr_20100704+ 
  PUBLIC :: ucase               ! turn a string into all uppercase
  PUBLIC :: spec2relhum         ! conversion from specific to relative
  PUBLIC :: spec2relhumwmo      ! same as spec2relhum with WMO def. of relhum
  PUBLIC :: rh2mr               ! rel. hum. -> mix. ratio (calls relhum*2mr)
  PUBLIC :: relhum2mr           ! relative humidity -> mol(H2O)/mol(dryair)
  PUBLIC :: relhumwmo2mr        ! WMO relative humidity -> mol(H2O)/mol(dryair)
  PUBLIC :: rel2spechum         ! conversion from relative to specific 
  PUBLIC :: rel2spechumwmo      ! conversion from relative (WMO) to specific 
  ! mz_hr_20100704- 

  PUBLIC :: find_next_free_unit

  ! mz_ab_20100610+
  PUBLIC :: remap_bounds  ! for obtaining pointers(-1:) etc, used for bounds
  ! mz_ab_20100610-

  ! mz_rs_20101118+
  PUBLIC :: density
  PUBLIC :: layerthickness
  ! mz_rs_20101118-

  ! mz_ab_20100831+
  PUBLIC :: full2half           ! converts variable on full level pressures
                                ! to half level pressures
  PUBLIC :: spline1D            ! spline interpolation
  PUBLIC :: splint1D
  ! mz_ab_20100831-
  ! mz_ab_20111125+
  PUBLIC :: CalcTimeAngles
  ! mz_ab_20111125+

  ! mz_bk_20110707+
  PUBLIC :: str2num             ! convert string to numerical value
  ! mz_bk_20110707-

CONTAINS

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_open(lex, substr, iou, nmlstr, modstr)

    USE messy_main_blather, ONLY: start_message
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    LOGICAL, INTENT(OUT)                   :: lex       ! file exists ?
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name


    CALL start_message(TRIM(modstr), 'INITIALISATION', substr)

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) '*** WARNING: FILE '''//TRIM(modstr)//'.nml'&
            &//'''  NOT FOUND !'
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(modstr)//'.nml')
    WRITE(*,*) 'Reading namelist '''//TRIM(nmlstr)//''''//&
         &' from '''//TRIM(modstr)//'.nml',''' (unit ',iou,') ...'

  END SUBROUTINE read_nml_open
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_check(fstat, substr, iou, nmlstr, modstr)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(IN)                    :: fstat     ! file status
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name

    IF (fstat /= 0) THEN
       WRITE(*,*) '*** ERROR: READ ERROR in NAMELIST '''//TRIM(nmlstr)//''''&
            &//' in FILE '''//TRIM(modstr)//'.nml'//''' !'
       CLOSE(iou)
    ELSE
       WRITE(*,*) ' ... OK !'
    END IF

  END SUBROUTINE read_nml_check
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_close(substr, iou, modstr)

    USE messy_main_blather, ONLY: end_message
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name

    CLOSE(iou)

    CALL end_message(TRIM(modstr), 'INITIALISATION', substr)

  END SUBROUTINE read_nml_close
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE iso2ind_1d(field, iso, k, f, lrev)
    
    ! ISOSURFACE TO INDEX
    ! OUTPUT IS THE LEVEL INDEX FOR A GIVEN ISOSURFACE
    ! NOTE:
    !   THIS ROUTINE IS WORKING PROPERLY ONLY IF THE 'FIELD'
    !   IS MONOTONIC
    !   (e.g., POTENTIAL TEMPERATURE, PRESSURE, POTENTIAL VORTICITY, ...)
    ! METHOD:
    !  lrev = .false. (default) -> SEARCH FROM LOW INDEX TO HIGH INDEX
    !  lrev = .true.            -> SEARCH FROM HIGH INDEX TO LOW INDEX
    !  ADJUST INDEX
    !  OPTIONAL: FRACTION OF BOX 'BELOW' ISO-SURFACE
    
    IMPLICIT NONE
    INTRINSIC :: NINT, ABS, MAX, MIN, PRESENT, REAL, SIZE, TINY

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: field  ! input field
    REAL(DP),               INTENT(IN)            :: iso    ! isosurface value
    INTEGER ,               INTENT(OUT)           :: k      ! isosurface level
    REAL(DP),               INTENT(OUT), OPTIONAL :: f      ! fraction in layer
    LOGICAL,                INTENT(IN),  OPTIONAL :: lrev   ! reverse order ?

    ! LOCAL
    INTEGER :: nk, jk, dk
    ! mz_rs_20080224: dk changed to integer
    REAL(DP) :: zf
    LOGICAL  :: llrev
    INTEGER  :: jstart, jstop, jstep

    k = 0
    nk = SIZE(field)    

    IF (PRESENT(lrev)) THEN
       llrev = lrev
    ELSE
       llrev = .FALSE. ! default
    END IF

    IF (llrev) THEN
       jstart = nk
       jstop  = 2
       jstep = -1
    ELSE
       jstart = 2
       jstop  = nk
       jstep = 1
    END IF

    DO jk = jstart, jstop, jstep
       IF ( ( (iso >= field(jk-1)) .AND. (iso <= field(jk)) ) .OR. &
            ( (iso <= field(jk-1)) .AND. (iso >= field(jk)) ) ) THEN
          k=jk
          EXIT
       END IF
    END DO

    IF ( k == 0 ) THEN   ! NOT FOUND
       IF (llrev) THEN
          k = 2
          IF (PRESENT(f)) f = 1.0
       ELSE
          k = nk
          IF (PRESENT(f)) f = 0.0
       END IF
       RETURN
    END IF

    ! ADJUST INDEX
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    !
    ! METHOD: LINEAR INTERPOLATION
    !
    ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
    ! SINCE THE FIRST LEVEL WITH 
    !    FIELD(k-1) <= ISO <= FLIELD(k)
    ! OR
    !    FIELD(k-1) >= ISO >= FLIELD(k)
    ! IS SEARCHED
    !
    IF ( ABS( field(k) - field(k-1) ) > TINY(0.0) ) THEN
       zf = ABS( (iso-field(k-1)) / (field(k)-field(k-1)) )    ! e [0,1)
    ELSE
       zf = 0.5_dp  ! SHOULD BE NEVER REACHED !!!
    END IF
    
    zf = MIN(1._dp,zf)
    zf = MAX(0._dp,zf)
    
    ! dk = INT(zf+0.5)
    dk = NINT(zf)
    ! zf e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
    ! zf e (0.5,1) -> dk = 1 -> KEEP LEVEL

    k = k - 1 + dk
    
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    ! ONE LEVEL ABOVE (dk = 0) -> zf e [0, 0.5]
    !                          -> f  e [0.5, 0]
    !       EXAMPLE: zf  = 0   -> ISO AT BOX MID         -> FRACT. = 0.5
    !                zf  = 0.5 -> ISO AT LOWER INTERFACE -> FRACT. = 0.0
    ! KEEP LEVEL      (dk = 1) -> zf e (0.5, 1)
    !                          -> f  e (1, 0.5)
    !       EXAMPLE: zf  = 0.5 -> ISO AT UPPER INTERFACE -> FRACT. = 1.0
    !                zf  = 1   -> ISO AT BOX MID         -> FRACT. = 0.5
    
    !                            FOR dk=1            FOR dk=0
    IF (PRESENT(f)) f = (1.5_dp-zf)*REAL(dk,DP) + (0.5_dp-zf)*REAL(1-dk,DP)
     
  END SUBROUTINE iso2ind_1d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE iso2ind_2d(kproma, field, iso, k, f, lrev)
    
    ! as iso2ind_1d, however for vectors ...
    
    IMPLICIT NONE
    INTRINSIC :: NINT, ABS, MAX, MIN, PRESENT, REAL, SIZE, TINY

    ! I/O
    INTEGER,                  INTENT(IN)            :: kproma
    REAL(DP), DIMENSION(:,:), INTENT(IN)            :: field ! input field
    REAL(DP), DIMENSION(:),   INTENT(IN)            :: iso   ! isosurface value
    INTEGER,  DIMENSION(:),   INTENT(OUT)           :: k     ! isosurface level
    REAL(DP), DIMENSION(:),   INTENT(OUT), OPTIONAL :: f     ! fraction in layer
    LOGICAL,                  INTENT(IN),  OPTIONAL :: lrev  ! reverse order ?

    ! LOCAL
    INTEGER  :: nk, jk, jl, dk
    ! mz_rs_20080224: dk changed to integer
    REAL(DP) :: zf
    LOGICAL  :: llrev
    INTEGER  :: jstart, jstop, jstep
    INTEGER, DIMENSION(kproma) :: ilfound

    k = 0
    ilfound(:) = 0
    nk = SIZE(field,2)    

    IF (PRESENT(lrev)) THEN
       llrev = lrev
    ELSE
       llrev = .FALSE. ! default
    END IF

    IF (llrev) THEN
       jstart = nk
       jstop  = 2
       jstep = -1
    ELSE
       jstart = 2
       jstop  = nk
       jstep = 1
    END IF

    DO jk = jstart, jstop, jstep
       DO jl = 1, kproma
          IF ( ilfound(jl) == 1) CYCLE
          IF ( ( (iso(jl) >= field(jl,jk-1)) .AND. &
               (iso(jl) <= field(jl,jk)) ) .OR. &
               ( (iso(jl) <= field(jl,jk-1)) .AND. &
               (iso(jl) >= field(jl,jk)) ) ) THEN
             k(jl) = jk
             ilfound(jl) = 1 
          END IF
       END DO
    END DO

    DO jl = 1, kproma

       IF ( k(jl) == 0 ) THEN   ! NOT FOUND
          IF (llrev) THEN
             k(jl) = 2
             IF (PRESENT(f)) f(jl) = 1.0_dp
          ELSE
             k(jl) = nk
             IF (PRESENT(f)) f(jl) = 0.0_dp
          END IF

       ELSE
          
          ! ADJUST INDEX
          ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
          !
          ! METHOD: LINEAR INTERPOLATION
          !
          ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
          ! SINCE THE FIRST LEVEL WITH 
          !    FIELD(k-1) <= ISO <= FLIELD(k)
          ! OR
          !    FIELD(k-1) >= ISO >= FLIELD(k)
          ! IS SEARCHED
          !
          IF ( ABS( field(jl,k(jl)) - field(jl,k(jl)-1) ) > TINY(0.0_dp) ) THEN
             zf = ABS( (iso(jl)-field(jl,k(jl)-1)) / &
                  (field(jl,k(jl))-field(jl,k(jl)-1)) )    ! e [0,1)
          ELSE
             zf = 0.5_dp  ! SHOULD BE NEVER REACHED !!!
          END IF
    
          zf = MIN(1._dp,zf)
          zf = MAX(0._dp,zf)
    
          ! dk = INT(zf+0.5)
          dk = NINT(zf)
          ! zf e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
          ! zf e (0.5,1) -> dk = 1 -> KEEP LEVEL
    
          k(jl) = k(jl) - 1 + dk
    
          ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
          ! ONE LEVEL ABOVE (dk = 0) -> zf e [0, 0.5]
          !                          -> f  e [0.5, 0]
          !       EXAMPLE: zf  = 0   -> ISO AT BOX MID         -> FRACT. = 0.5
          !                zf  = 0.5 -> ISO AT LOWER INTERFACE -> FRACT. = 0.0
          ! KEEP LEVEL      (dk = 1) -> zf e (0.5, 1)
          !                          -> f  e (1, 0.5)
          !       EXAMPLE: zf  = 0.5 -> ISO AT UPPER INTERFACE -> FRACT. = 1.0
          !                zf  = 1   -> ISO AT BOX MID         -> FRACT. = 0.5

          !                            FOR dk=1            FOR dk=0
          IF (PRESENT(f)) f(jl) = (1.5_dp-zf)*REAL(dk,dp) + &
               (0.5_dp-zf)*REAL(1-dk,dp)
          
       END IF

    END DO
     
  END SUBROUTINE iso2ind_2d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE ind2val_1d(val, field, k, f)

    ! CONVERT INDEX (AND FRACTION 'BELOW') IN A MONOTONOUS FIELD
    ! TO THE VALUE
    ! METHOD:
    !   - PICK OUT INDEX
    !   - LINEAR INTERPOLATION BETWEEN NEIGHBOURS, IF f IS PRESENT
    
    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    REAL(DP),               INTENT(OUT)          :: val   ! value
    REAL(DP), DIMENSION(:), INTENT(IN)           :: field ! field
    INTEGER,                INTENT(IN)           :: k     ! level
    REAL(DP),               INTENT(IN), OPTIONAL :: f     ! fraction
    
    ! LOCAL
    INTEGER  :: nk
    REAL(DP) :: ri, gf
    
    nk = SIZE(field)
    
    IF (PRESENT(f)) THEN
       ri  = 0.5_dp - f           ! e (-0.5,0.5) -> (top, bottom) of box
       IF (ri >= 0.0_dp) THEN
          IF (k == nk) THEN
             gf  = (field(k)-field(k-1))
          ELSE
             gf  = (field(k+1)-field(k))
          END IF
       ELSE
          IF (k == 1) THEN
             gf  = (field(k+1)-field(k))
          ELSE
             gf  = (field(k)-field(k-1))
          END IF
       END IF
       !
       val = field(k) + ri * gf
    ELSE
       val = field(k)
    END IF
    
  END SUBROUTINE ind2val_1d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE ind2val_2d(kproma, val, field, k, f)

    ! as ind2val_1d, however for vectors ...
    
    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    INTEGER,                  INTENT(IN)           :: kproma
    REAL(DP), DIMENSION(:),   INTENT(OUT)          :: val   ! value
    REAL(DP), DIMENSION(:,:), INTENT(IN)           :: field ! field
    INTEGER,  DIMENSION(:),   INTENT(IN)           :: k     ! level
    REAL(DP), DIMENSION(:),   INTENT(IN), OPTIONAL :: f     ! fraction
    
    ! LOCAL
    INTEGER  :: jl, nk
    REAL(DP) :: ri, gf
    
    nk = SIZE(field,2)
    
    IF (PRESENT(f)) THEN
       DO jl = 1, kproma
          ri  = 0.5_dp - f(jl)       ! e (-0.5,0.5) -> (top, bottom) of box
          IF (ri >= 0.0_dp) THEN
             IF (k(jl) == nk) THEN
                gf  = (field(jl,k(jl))-field(jl,k(jl)-1))
             ELSE
                gf  = (field(jl,k(jl)+1)-field(jl,k(jl)))
             END IF
          ELSE
             IF (k(jl) == 1) THEN
                gf  = (field(jl,k(jl)+1)-field(jl,k(jl)))
             ELSE
                gf  = (field(jl,k(jl))-field(jl,k(jl)-1))
             END IF
          END IF
          !
          val(jl) = field(jl,k(jl)) + ri * gf
       END DO
    ELSE
       DO jl = 1, kproma
          val(jl) = field(jl,k(jl))
       END DO
    END IF
    
  END SUBROUTINE ind2val_2d
! -----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE int2str(str, ii, cpad, cerr)
    
    IMPLICIT NONE
    INTRINSIC :: MOD, LEN, PRESENT
    
    ! I/O
    CHARACTER(LEN=*), INTENT(OUT)          :: str   ! STRING
    INTEGER,          INTENT(IN)           :: ii    ! INTEGER
    CHARACTER,        INTENT(IN), OPTIONAL :: cpad  ! CHAR FOR PADDING
    CHARACTER,        INTENT(IN), OPTIONAL :: cerr  ! CHAR FOR ERROR
    
    ! LOCAL
    INTEGER              :: n, zi, zn, k
    INTEGER              :: rest
    CHARACTER            :: zcpad
    
    IF (PRESENT(cpad)) THEN
       zcpad = cpad
    ELSE
       zcpad = '0'      ! DEFAULT PADDING
    END IF
    
    n  = LEN(str)
    zi = ii
    zn = n
    
    DO
       rest = MOD(zi, 10)
       zi   = zi/10
       WRITE(str(zn:zn),'(i1)') rest
       zn = zn - 1
       IF (zi == 0) EXIT
       IF (zn == 0) EXIT
    END DO
    
    IF (PRESENT(cerr)) THEN
       IF ( (zn == 0) .AND. (zi /= 0) ) THEN
          DO k = 1, n
             WRITE(str(k:k),'(a1)') cerr
          END DO
       END IF
    END IF
    
    DO k = zn, 1, -1
       WRITE(str(k:k),'(a1)') zcpad
    END DO
    
  END SUBROUTINE int2str
! ---------------------------------------------------------------------

! mz_rs_20060110+
! ---------------------------------------------------------------------
  ! various types of the function str:

  CHARACTER(LEN=5) FUNCTION str_logical(zlogical, fmt)
    ! create string from logical
    IMPLICIT NONE
    INTRINSIC :: PRESENT
    LOGICAL, INTENT(in) :: zlogical
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_logical,fmt) zlogical
    ELSE
      IF (zlogical) THEN
        str_logical = 'TRUE '
      ELSE
        str_logical = 'FALSE'
      ENDIF
    ENDIF
  END FUNCTION str_logical

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_integer(zinteger, fmt)
    ! create string from integer
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    INTEGER, INTENT(in) :: zinteger
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_integer,fmt) zinteger
    ELSE
      WRITE(str_integer,*) zinteger
      str_integer = ADJUSTL(str_integer) ! remove leading spaces
    ENDIF
  END FUNCTION str_integer

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_real_sp(zreal_sp, fmt)
    ! create string from real_sp
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    REAL(sp), INTENT(in) :: zreal_sp
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_real_sp,fmt) zreal_sp
    ELSE
      WRITE(str_real_sp,*) zreal_sp
      str_real_sp = ADJUSTL(str_real_sp) ! remove leading spaces
    ENDIF
  END FUNCTION str_real_sp

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_real_dp(zreal_dp, fmt)
    ! create string from real_dp
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    REAL(dp), INTENT(in) :: zreal_dp
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_real_dp,fmt) zreal_dp
    ELSE
      WRITE(str_real_dp,*) zreal_dp
      str_real_dp = ADJUSTL(str_real_dp) ! remove leading spaces
    ENDIF
  END FUNCTION str_real_dp
! ---------------------------------------------------------------------
! mz_rs_20060110-

! ---------------------------------------------------------------------
  SUBROUTINE strcrack(str, ch, el, n)

    ! strcrack = string crack

    ! Split the string <str> into small pieces which are separated by
    ! the character <ch>. Delete trailing spaces from the resulting <n>
    ! pieces, then put them into the array <el>.

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, ASSOCIATED, INDEX, LEN_TRIM, TRIM

    ! I/O
    CHARACTER(LEN=*),               INTENT(IN)  :: str ! string 2b cracked
    CHARACTER,                      INTENT(IN)  :: ch  ! separating char
    CHARACTER(LEN=*), DIMENSION(:), POINTER     :: el  ! field of substrings
    INTEGER,                        INTENT(OUT) :: n   ! # of substrings
    
    ! LOCAL
    INTEGER :: idx1, idx2, i
    
    ! INIT
    IF (ASSOCIATED(el)) DEALLOCATE(el)
    NULLIFY(el)
    n = 0
    
    ! EMPTY STRING
    IF ( (TRIM(str) == '') .OR. (TRIM(str) == ch) ) RETURN
    
    idx1 = 0
    idx2 = 0
    DO 
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
    END DO
    
    ! ALLOCATE SPACE
    ALLOCATE(el(n))
    DO i=1, n
       el(i) = ''
    END DO
    
    n = 0
    idx1 = 0
    idx2 = 0
    DO     
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
!       el(n) = str(idx1:idx2-1)
       el(n) = ADJUSTL(str(idx1:idx2-1))
       
    END DO

  END SUBROUTINE strcrack
! ---------------------------------------------------------------------  

! mz_ap_20070913+
! ---------------------------------------------------------------------
  SUBROUTINE ns_index(list, value, idx, idx2, l_periodic)

    IMPLICIT NONE
    INTRINSIC :: ABS, PRESENT, SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: list
    REAL(DP),               INTENT(IN)            :: value
    INTEGER,                INTENT(OUT)           :: idx
    INTEGER,                INTENT(OUT), OPTIONAL :: idx2
    LOGICAL,                INTENT(IN),  OPTIONAL :: l_periodic 
    ! list(i) < value < list(i+1)
    ! l_periodic : periodic boundaries

    ! LOCAL
    INTEGER     :: n, i

    ! INIT
    n = SIZE(list)
    idx = 0

    ! we suppouse that both lists are ordered
    IF (value<list(1)) THEN
      DO i=1, n
         IF (value <= list(i)) THEN
            idx = i
         END IF
      END DO
    ELSE
      DO i=1, n
         IF (value >= list(i)) THEN
            idx = i
         END IF
      END DO
    END IF
    IF (PRESENT(idx2)) THEN
       IF (idx == n) THEN
          IF (PRESENT(l_periodic)) THEN
            idx2 = 1
          ELSE
            idx2 = n
          ENDIF
       ELSE
          idx2 = idx+1
       ENDIF
    END IF
    IF (idx == 0 ) THEN
      IF (PRESENT(l_periodic)) THEN 
        idx=n
        IF (PRESENT(idx2)) idx2=1
      ELSE 
        idx=1
        IF (PRESENT(idx2)) idx2=2
      ENDIF
    ENDIF

  END SUBROUTINE ns_index
! mz_ap_20070913-
! --------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE nn_index(list, value, idx, idx2, fac, status)

    IMPLICIT NONE
    INTRINSIC :: ABS, PRESENT, SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: list
    REAL(DP),               INTENT(IN)            :: value
    INTEGER,                INTENT(OUT)           :: idx
    INTEGER,                INTENT(OUT), OPTIONAL :: idx2
    REAL(DP),               INTENT(IN),  OPTIONAL :: fac
    INTEGER,                INTENT(OUT), OPTIONAL :: status

    ! LOCAL
    REAL(DP)    :: zfac
    INTEGER     :: n, i
    REAL(DP)    :: dmin

    ! INIT
    n = SIZE(list)
    idx = 0
    IF (PRESENT(idx2)) idx2 = 0
    IF (PRESENT(fac)) THEN
       zfac = fac
    ELSE
       zfac = 2.0_dp     ! default
    END IF
    IF (PRESENT(status)) status = 0

    ! LIST MUST CONTAIN AT LEAST TWO ELEMENTS
    IF (n < 2) THEN
       idx = 1
       IF (PRESENT(idx2))   idx2 = 1
       IF (PRESENT(status)) status = 1
       RETURN
    END IF

    ! ALLOW zfac TIMES THE FIRST / LAST INTERVAL FOR EXTRAPOLATION ...
    IF (list(1) <= list(n)) THEN
       ! increasing ...
       IF (value < (list(1) - (list(2)-list(1))*zfac) ) THEN
          idx = 1
          IF (PRESENT(idx2))   idx2 = 1
          IF (PRESENT(status)) status = 2
          RETURN
       END IF
       IF (value > (list(n) + (list(n)-list(n-1))*zfac) ) THEN
          idx = n
          IF (PRESENT(idx2))   idx2 = n
          IF (PRESENT(status)) status = 3
          RETURN
       END IF
    ELSE
       ! decreasing ...
       IF (value > (list(1) - (list(2)-list(1))*zfac) ) THEN
          idx = 1
          IF (PRESENT(idx2))   idx2 = 1
          IF (PRESENT(status)) status = 4
          RETURN
       END IF
       IF (value < (list(n) + (list(n)-list(n-1))*zfac) ) THEN
          idx = n
          IF (PRESENT(idx2))   idx2 = n
          IF (PRESENT(status)) status = 5
          RETURN
       ENDIF
    END IF

    ! SEARCH MINIMUM
    dmin = ABS(list(1) - list(n))
    DO i=1, n
       IF (ABS(list(i)-value) <= dmin) THEN
          dmin = ABS(list(i)-value)
          idx = i
       END IF
    END DO

    ! NOTE: FOR IDX2 IT IS ASSUMED THAT THE LIST IS ORDERED
    IF (PRESENT(idx2)) THEN
       IF (idx == 1) THEN
          idx2 = 2
          RETURN
       END IF
       IF (idx == n) THEN
          idx2 = n-1
          RETURN
       END IF
       IF ( ABS(list(idx+1)-value) <= ABS(list(idx-1)-value) ) THEN
          idx2 = idx+1
       ELSE
          idx2 = idx-1
       END IF
    END IF

  END SUBROUTINE nn_index
! --------------------------------------------------------------------

! mz_ab_20100624+
! --------------------------------------------------------------------
  FUNCTION fliparray(list) RESULT(buffer)
    
    ! flips (reverses) array:
    ! (1,2,3,4) --> (4,3,2,1)

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN) :: list
    ! LOCAL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: buffer
    INTEGER :: i, n

    n = SIZE(list,1)

    ALLOCATE(buffer(n))

    DO i=1,n
       buffer(i) = list(n-i+1)
    END DO

  END FUNCTION fliparray
! --------------------------------------------------------------------
! mz_ab_20100624-

! mz_ht_20042510+
! ---------------------------------------------------------------------  
  SUBROUTINE init_convect_tables

    ! Lookup tables for convective adjustment code
    !
    ! D. Salmond, CRAY (UK), August 1991, original code

    USE messy_main_constants_mem, ONLY: rd, rv, tmelt, cpd => cp_air &
                                      , alv, als ! op_pj_20100209

    IMPLICIT NONE
    INTRINSIC :: EXP, LOG
                       
    REAL(dp), PARAMETER :: zavl1 = -6096.9385_dp
    REAL(dp), PARAMETER :: zavl2 =    21.2409642_dp
    REAL(dp), PARAMETER :: zavl3 =    -2.711193_dp
    REAL(dp), PARAMETER :: zavl4 =     1.673952_dp
    REAL(dp), PARAMETER :: zavl5 =     2.433502_dp 

    REAL(dp), PARAMETER :: zavi1 = -6024.5282_dp
    REAL(dp), PARAMETER :: zavi2 =    29.32707_dp
    REAL(dp), PARAMETER :: zavi3 =     1.0613868_dp
    REAL(dp), PARAMETER :: zavi4 =    -1.3198825_dp
    REAL(dp), PARAMETER :: zavi5 =    -0.49382577_dp   
    
! op_pj_20100209+
!!$    REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for
!!$    !                                          ! vaporisation in J/kg
!!$    REAL(dp), PARAMETER :: als   = 2.8345e6_dp ! latent heat for
!!$    !                                          ! sublimation in J/kg
! op_pj_20100209-
    ! Constants used for computation of saturation mixing ratio
    ! over liquid water (*c_les*) or ice(*c_ies*)
    REAL(dp),PARAMETER :: c3les = 17.269_dp           ! 
    REAL(dp),PARAMETER :: c3ies = 21.875_dp           ! 
    REAL(dp),PARAMETER :: c4les = 35.86_dp            ! 
    REAL(dp),PARAMETER :: c4ies = 7.66_dp             ! 
    REAL(dp),PARAMETER :: c5les = c3les*(tmelt-c4les) ! 
    REAL(dp),PARAMETER :: c5ies = c3ies*(tmelt-c4ies) ! 
   
    REAL(dp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL(dp) :: ztt, zldcp
!!$    REAL(dp) :: zcvm3, zcvm4, zcvm5
    REAL(dp) :: zcvm4, zcvm5
    REAL(dp) :: zavm1, zavm2, zavm3, zavm4, zavm5

    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    DO it = jptlucu1, jptlucu2
      ztt = 0.001_dp*it
      IF ((ztt-tmelt) > 0.0_dp) THEN
!!$        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
!!$        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF
      tlucuc(it)  = zldcp
      tlucua(it)  = EXP((zavm1/ztt+zavm2+zavm3*0.01_dp* &
           ztt+zavm4*ztt*ztt*1.e-5_dp+zavm5*LOG(ztt)))*rd/rv
      tlucub(it)  = zcvm5*(1.0_dp/(ztt-zcvm4))**2
      tlucuaw(it) = EXP((zavl1/ztt+zavl2+zavl3*0.01_dp* &
           ztt+zavl4*ztt*ztt*1.e-5_dp+zavl5*LOG(ztt)))*rd/rv
    END DO
    
  END SUBROUTINE init_convect_tables
! ---------------------------------------------------------------------  
! mz_ht_20042510-

! ---------------------------------------------------------------------  
LOGICAL FUNCTION match_wild(pattern, string)

   ! compare given string for match to pattern which may
   ! contain wildcard characters:
   ! "?" matching any one character, and
   ! "*" matching any zero or more characters.
   ! Both strings may have trailing spaces which are ignored.
   ! Authors: Clive Page, userid: cgp  domain: le.ac.uk, 2003 (original code)
   !          Rolf Sander, 2005 (bug fixes and pattern preprocessing)
   ! Minor bug fixed by Clive Page, 2005 Nov 29.

   ! This program is free software; you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as published by
   ! the Free Software Foundation; either version 2 of the License, or
   ! (at your option) any later version.
   !
   ! This program is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with this program; if not, write to the Free Software
   ! Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 
   ! 02110-1301  USA

   IMPLICIT NONE
   INTRINSIC :: INDEX, LEN, LEN_TRIM, REPEAT

   CHARACTER(LEN=*), INTENT(IN) :: pattern ! pattern may contain * and ?
   CHARACTER(LEN=*), INTENT(IN) :: string  ! string to be compared
   INTEGER :: lenp, lenp2, lens, n, p2, p, s
   INTEGER :: n_question, n_asterisk

   CHARACTER(LEN=LEN(pattern)) :: pattern2

   lens = LEN_TRIM(string)
   lenp = LEN_TRIM(pattern)

   ! If the pattern is empty, always return true
   IF (lenp == 0) THEN
     match_wild = .TRUE.
     RETURN
   ENDIF

   ! The pattern must be preprocessed. All consecutive occurences of
   ! one or more question marks ('?') and asterisks ('*') are sorted and
   ! compressed. The result is stored in pattern2.

   pattern2(:)=''
   p  = 1 ! current position in pattern
   p2 = 1 ! current position in pattern2
   DO
     IF ((pattern(p:p) == '?').OR.(pattern(p:p) == '*')) THEN
       ! a special character was found in the pattern
       n_question = 0
       n_asterisk = 0
       DO WHILE (p <= lenp)
         ! count the consecutive question marks and asterisks
         IF ((pattern(p:p) /= '?').AND.(pattern(p:p) /= '*')) EXIT
         IF (pattern(p:p) == '?') n_question = n_question + 1
         IF (pattern(p:p) == '*') n_asterisk = n_asterisk + 1
         p = p + 1
       ENDDO
       IF (n_question>0) THEN ! first, all the question marks
         pattern2(p2:p2+n_question-1) = REPEAT('?',n_question)
         p2 = p2 + n_question
       ENDIF
       IF (n_asterisk>0) THEN ! next, the asterisk (only one!)
         pattern2(p2:p2) = '*'
         p2 = p2 + 1
       ENDIF
     ELSE
       ! just a normal character
       pattern2(p2:p2) = pattern(p:p)
       p2 = p2 + 1
       p = p + 1
     ENDIF
     IF (p > lenp) EXIT
   ENDDO
   lenp2 = LEN_TRIM(pattern2)

   ! The modified wildcard in pattern2 is compared to the string:

   p2 = 1
   s = 1
   match_wild = .FALSE.
   DO
     IF (pattern2(p2:p2) == '?') THEN
       ! accept any char in string
       p2 = p2 + 1
       s = s + 1
     ELSEIF (pattern2(p2:p2) == "*") THEN
       p2 = p2 + 1
       IF (p2 > lenp2) THEN
         ! anything goes in rest of string
         match_wild = .TRUE.
         EXIT ! .TRUE.
       ELSE
         ! search string for char at p2
         n = INDEX(string(s:), pattern2(p2:p2))
         IF (n == 0) EXIT  ! .FALSE.
         s = n + s - 1
       ENDIF
     ELSEIF (pattern2(p2:p2) == string(s:s)) THEN
       ! single char match
       p2 = p2 + 1
       s = s + 1
     ELSE
       ! non-match
       EXIT ! .FALSE.
     ENDIF
     IF (p2 > lenp2 .AND. s > lens) THEN
       ! end of both pattern2 and string
       match_wild = .TRUE.
       EXIT ! .TRUE.
     ENDIF

     IF (s > lens .AND. p2 == lenp) THEN
       IF(pattern2(p2:p2) == "*") THEN
         ! "*" at end of pattern2 represents an empty string
         match_wild = .TRUE.
         EXIT ! .TRUE.
       END IF
     ENDIF
     IF (p2 > lenp2 .OR. s > lens) THEN
       ! end of either pattern2 or string
       EXIT ! .FALSE.
     ENDIF
   ENDDO

END FUNCTION match_wild
! ---------------------------------------------------------------------  

! ----------------------------------------------------------------------
  SUBROUTINE str2chob(status, str, n, c, o)

! op_pj_20100803+
!!$    USE messy_main_constants_mem, ONLY: STRLEN_ULONG
! op_pj_20100803-

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,           INTENT(OUT)               :: status
    CHARACTER(LEN=*),  INTENT(IN)                :: str
    INTEGER,           INTENT(OUT)               :: n
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: c
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: o

    ! LOCAL
    INTEGER                                            :: nc, no, two, j, i
! op_pj_20100803+
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp1 => NULL()
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp2 => NULL()
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp3 => NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp1 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp2 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp3 !=> NULL()
! op_pj_20100803-

    ! INIT
    IF (ASSOCIATED(c)) DEALLOCATE(c)
    IF (ASSOCIATED(o)) DEALLOCATE(o)
    NULLIFY(c)
    NULLIFY(o)
    n = 0
    ! op_pj_20100811+
    NULLIFY(tmp1)
    NULLIFY(tmp2)
    NULLIFY(tmp3)
    ! op_pj_20100811-

    status = 1 ! ERROR

    ! COUNT OBJECTS
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       IF (two /= 2) RETURN      ! ERROR: 0 or more than 1 ':' in string
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       n = n + no
    END DO

    ! ALLOCATE SPACE
    ALLOCATE(c(n))
    ALLOCATE(o(n))
    ! INIT
    DO i=1,n
       c(i) = ''
       o(i) = ''
    END DO

    ! PARSE STRING
    n = 0
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       DO j=1, no
          n = n + 1
          c(n) = TRIM(tmp2(1))
          o(n) = TRIM(tmp3(j))
       END DO
    END DO   

    IF (ASSOCIATED(tmp1)) THEN
      DEALLOCATE(tmp1) ; NULLIFY(tmp1)
    ENDIF
    IF (ASSOCIATED(tmp2)) THEN
      DEALLOCATE(tmp2) ; NULLIFY(tmp2)
    ENDIF
    IF (ASSOCIATED(tmp3)) THEN
      DEALLOCATE(tmp3) ; NULLIFY(tmp3)
    ENDIF

    status  = 0  ! NO ERROR

  END SUBROUTINE str2chob
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE bilin_weight(vn, v, w)
    
    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(4, 2), INTENT(IN)  :: vn ! neighbours
    REAL(dp), DIMENSION(2),    INTENT(IN)  :: v  ! position
    REAL(dp), DIMENSION(4),    INTENT(OUT) :: w  ! weight

    ! LOCAL
    REAL(DP) :: t, u

    t = (v(1) - vn(1,1)) / (vn(2,1) - vn(1,1))
    u = (v(2) - vn(1,2)) / (vn(3,2) - vn(1,2))

    w(1) = (1._dp-t)*(1._dp-u)
    w(2) = t*(1._dp-u)
    w(3) = t*u
    w(4) = (1._dp-t)*u

  END SUBROUTINE bilin_weight
! ----------------------------------------------------------------------

! mz_hr_20100704+ 
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION psat_mk(p_temp, p_l_liquid)

    ! psat_mk replaces psat and psatf

    ! saturation water vapor pressure [Pa] over liquid or ice as
    ! recommended by Murphy and Koop, Q. J. R. Meteorol. Soc. (2005),
    ! 131, pp. 1539-1565
    ! temperature in Kelvin [K]
    ! switch from ice to liquid above 273.15 K (= tmelt) in analogy
    !   to convect_init_tables
    ! valid temperature range liquid formula: 123 K < T < 332 K
    ! valid temperature range ice formula:    110 K < T < 273.16 K

    ! always use l_liquid=T for WMO relative humidity and for
    ! supercooling effects and supersaturation

    USE messy_main_constants_mem, ONLY: tmelt

    IMPLICIT NONE

    INTRINSIC LOG, EXP, TANH

    ! I/O
    REAL(dp), INTENT(in)           :: p_temp ! temperature [K]
    LOGICAL,  INTENT(IN), OPTIONAL :: p_l_liquid
    ! local
    REAL(dp), PARAMETER      :: lowlim_ice = 110._dp ! lower temp limit [K]
    REAL(dp), PARAMETER      :: lowlim_liq = 123._dp ! lower temp limit [K]
    REAL(dp), PARAMETER      :: uplim      = 332._dp ! upper temp limit [K]
    REAL(dp) :: lowlim
    LOGICAL :: l_liquid

    IF (PRESENT(p_l_liquid)) THEN
      l_liquid = p_l_liquid
    ELSE
      l_liquid = .FALSE.
    ENDIF

    IF (l_liquid) THEN
      lowlim = lowlim_liq
    ELSE
      lowlim = lowlim_ice
    ENDIF

    ! test temperature range and l_liquid:
    IF ( (p_temp < uplim) .AND. &
      ( (p_temp > tmelt) .OR. (l_liquid.AND.(p_temp > lowlim)) ) ) THEN
      ! water vapor saturation pressure over liquid:
      psat_mk = EXP( 54.842763_dp - 6763.22_dp/p_temp &
        - 4.21_dp*LOG(p_temp) + 0.000367_dp*p_temp    &
        + TANH(0.0415_dp*(p_temp-218.8_dp))           &
        * (53.878_dp - 1331.22_dp/p_temp              &
        - 9.44523_dp*LOG(p_temp) + 0.014025_dp*p_temp) )
    ELSEIF ((p_temp <= tmelt).AND.(.NOT.l_liquid).AND.(p_temp > lowlim)) THEN
      ! water vapor saturation pressure over ice:
      psat_mk = EXP( 9.550426_dp - 5723.265_dp/p_temp &
        + 3.53068_dp*LOG(p_temp) - 0.00728332_dp*p_temp )
    ELSE
      WRITE(*,'(A,F10.3,A,F10.3,A,F10.3,A)') '  WARNING psat_mk: &
        & temperature (', p_temp, ') out of valid temperature range [', &
        lowlim, ',', uplim, ']'
    ENDIF

  END FUNCTION psat_mk
! ---------------------------------------------------------------------
! mz_hr_20100704-

! mz_hr_20080221+
! ---------------------------------------------------------------------
  SUBROUTINE ucase(string)
    ! from book Stephen Chapman "F90/95"

    IMPLICIT NONE

    INTRINSIC LGE, LLE, ACHAR, IACHAR

    ! I/O
    CHARACTER(LEN=*),        INTENT(INOUT)  :: string

    ! LOCAL
    INTEGER :: i 
    INTEGER :: length

    ! get len of str
    length = LEN(string)

    ! shift lower case to upper case
    DO i=1, length
      IF (LGE(string(i:i),'a') .AND. LLE(string(i:i),'z')) THEN
        string(i:i) = ACHAR(IACHAR(string(i:i)) - 32)
      ENDIF
    ENDDO
    
  END SUBROUTINE ucase
! ---------------------------------------------------------------------
! mz_hr_20080221-

! mz_hr_20080226+
! ----------------------------------------------------------------------------
  REAL(dp) FUNCTION spec2relhumwmo(status, z_spechum, z_temp, z_press, &
                                   z_l_psat_emac) ! mz_hr_20100704
    ! mz_hr_20100704               z_l_psatf) !mz_hr_20100608

    ! calculates relative humidity from specific humidity
    ! relative humidity as defined by WMO:
    ! mass mixing ratio water vapor/saturation mass mixing ratio water vapor
    ! [Jacobson, Fundamentals of Atmospheric Modeling, CamUnivPress, 1999]

    ! molar mass of water (vapour) / molar mass of dry air
    ! mz_hr_20100704+ 
    !USE messy_main_constants_mem, ONLY: MM_eps, TINY_DP, FLAGGED_BAD
    USE messy_main_constants_mem, ONLY: MM_eps, TINY_DP, FLAGGED_BAD, Rd, Rv
    ! mz_hr_20100704- 

    IMPLICIT NONE

    ! e funct, natural log, absolute value
    INTRINSIC ABS, INT

    !I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_spechum    ! kg/kg
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    !mz_hr_20100608+
    !LOGICAL, INTENT(IN), OPTIONAL :: z_l_psatf
    !mz_hr_20100608-
    ! mz_hr_20100704+ 
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac
    ! mz_hr_20100704- 

    !LOCAL
    REAL(dp) :: omega_v, omega_vs, psat

    status = 0

    IF (ABS(1._dp - z_spechum) < TINY_DP) THEN
      WRITE(*,*) 'ERROR spec2relhumwmo: spechum = 1, division by 0'
      spec2relhumwmo = FLAGGED_BAD
      status = 1
      RETURN
    !mz_hr_20081125+
    ELSEIF (z_spechum > 1._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhumwmo: spechum > 1 : ', z_spechum
      spec2relhumwmo = FLAGGED_BAD
      status = 1
      RETURN
    ELSEIF (z_spechum < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhumwmo: spechum < 0 : ', z_spechum
      spec2relhumwmo = FLAGGED_BAD
      status = 1
      RETURN
    !mz_hr_20081125-
    ENDIF
    
    ! water vapor saturation mass mixing ratio
    ! always use psat over liquid water
    ! mz_hr_20100705+ 
    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucuaw(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
         ! l_relhum_wmo=.TRUE. since here WMO function
        psat = psat_mk(z_temp,.TRUE.)
      ENDIF
    ELSE
       ! l_relhum_wmo=.TRUE. since here WMO function
      psat = psat_mk(z_temp,.TRUE.)
    ENDIF 

    ! convert spechum to mass mixing ratio of water in dry air
    omega_v = z_spechum/(1._dp-z_spechum)
    ! saturation mass mixing ratio
    omega_vs = MM_eps * psat / (z_press - psat)
    ! mz_hr_20100705- 

    ! calc relhum, def by World Meteorological Organization WMO
    spec2relhumwmo = omega_v / omega_vs

    IF ((spec2relhumwmo >= 1._dp) .AND. (spec2relhumwmo <= 1.1_dp)) THEN
      WRITE(*,'(A, F10.3)') '  WARNING spec2relhumwmo: relhum >= 1'// &
        ' (< 1.1) : ', spec2relhumwmo
    ELSEIF ((spec2relhumwmo > 1.1_dp)) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhumwmo: relhum > 1.1 : ', &
        spec2relhumwmo
      spec2relhumwmo = FLAGGED_BAD
      status = 1
      RETURN
    ELSEIF (spec2relhumwmo < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhumwmo: relhum < 0 : ', &
        spec2relhumwmo
      spec2relhumwmo = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION spec2relhumwmo
! ----------------------------------------------------------------------------
! mz_hr_20080226-

! mz_hr_20100607+
! ---------------------------------------------------------------------
  REAL(dp) FUNCTION spec2relhum(status, z_spechum, z_temp, z_press, &
                                z_l_psat_emac)
    !                           z_l_psatf) ! mz_hr_20100704 

    ! calculates relative humidity from specific humidity
    ! traditional definition:
    ! partial press H20/saturation partial pressure H2O
    ! relhum = p(H2O)/psat(H2O)

    ! mz_hr_20100704+ 
    !USE messy_main_constants_mem, ONLY: MM_eps, FLAGGED_BAD
    USE messy_main_constants_mem, ONLY: MM_eps, FLAGGED_BAD, tmelt, Rv, Rd
    ! mz_hr_20100704- 

    IMPLICIT NONE

    !I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_spechum    ! kg/kg
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    !mz_hr_20100608+
    !LOGICAL, INTENT(IN), OPTIONAL :: z_l_psatf
    !mz_hr_20100608-
    ! mz_hr_20100704+ 
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac
    ! mz_hr_20100704- 

    !LOCAL
    REAL(dp) :: p_v, psat, omega_v

    status = 0

    !mz_hr_20081125+
    IF (z_spechum > 1._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhum: spechum > 1 : ', z_spechum
      spec2relhum = FLAGGED_BAD
      status = 1
      RETURN
    ELSEIF (z_spechum < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhum: spechum < 0 : ', z_spechum
      spec2relhum = FLAGGED_BAD
      status = 1
      RETURN
    !mz_hr_20081125-
    ENDIF
    
    ! mz_hr_20100705+ 
    ! for supersaturation use psat_mk with l_liquid=.TRUE. at all temp
    !IF (PRESENT(z_l_psatf)) THEN ! mz_hr_20100608+
    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucua(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
        psat = psat_mk(z_temp)
      ENDIF
    ELSE
      psat = psat_mk(z_temp)
    ENDIF

    ! water vapor partial pressure
    ! taking into account that press = press of humid air
    ! water vapor mass mixing ratio in dry air
    omega_v = z_spechum/(1._dp-z_spechum)
    p_v = z_press * omega_v / (MM_eps + omega_v)
    ! same but complicated:
    ! p_v = z_spechum/MM_eps * z_press / (1-z_spechum+z_spechum/MM_eps)
    ! wrong since press NOT dry air press:
    !p_v = z_spechum/MM_eps * z_press

    spec2relhum = p_v / psat
    !spec2relhum = p_v / psatf(z_temp)
    ! mz_hr_20100705- 

    ! mz_hr_20100706+ 
    !IF ((spec2relhum >= 1._dp) .AND. (spec2relhum <= 1.1_dp)) THEN
    IF (spec2relhum > 1._dp) THEN
      WRITE(*,'(A, F10.3)') '  WARNING spec2relhum: relhum > 1 : ', &
        spec2relhum
    ! mz_hr_20100706- 
    ! mz_hr_20100706+ error only if wmo, with this trad supersat ok
    !ELSEIF ((spec2relhum > 1.1_dp)) THEN
    !  WRITE(*,*) 'ERROR spec2relhum: calculated relhum > 1.1 : ', &
    !    spec2relhum
    !  spec2relhum = FLAGGED_BAD
    !  status = 1
    !  RETURN
    ! mz_hr_20100706- 
    ELSEIF (spec2relhum < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR spec2relhum: relhum < 0 : ', &
        spec2relhum
      spec2relhum = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION spec2relhum
! --------------------------------------------------------------------- 
! mz_hr_20100607-

! ---------------------------------------------------------------------
  REAL(dp) FUNCTION rel2spechumwmo(status, z_relhum, z_temp, z_press, &
                                  z_l_psat_emac) ! mz_hr_20100704
    ! mz_hr_20100704              z_l_psatf) !mz_hr_20100608

    ! calculates specific humidity from relative humidity
    ! for WMO definition, see, e.g., 
    ! [Jacobson, Fundamentals of Atmospheric Modeling, CamUnivPress, 1999]

    ! molar mass of water (vapour) / molar mass of dry air
    ! mz_hr_20100705+ 
    !USE messy_main_constants_mem, ONLY: MM_eps ! M_H20/M_air
    USE messy_main_constants_mem, ONLY: MM_eps, Rd, Rv
    ! mz_hr_20100705- 

    IMPLICIT NONE

    INTRINSIC INT

    !I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum     ! relhum [0-1] 
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    ! mz_hr_20100704+ 
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac
    !LOGICAL, INTENT(IN), OPTIONAL :: z_l_psatf ! mz_hr_20100608-
    ! mz_hr_20100704- 
    

    !LOCAL
    REAL(dp) ::  omega_vs, psat

    status = 0

    ! mz_hr_20100704+ 
    !IF (PRESENT(z_l_psatf)) THEN !mz_hr_20100608+
    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucuaw(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
         ! l_relhum_wmo=.TRUE. since here WMO function
        psat = psat_mk(z_temp,.TRUE.)
      ENDIF
    ELSE
       ! l_relhum_wmo=.TRUE. since here WMO function
      psat = psat_mk(z_temp,.TRUE.)
    ENDIF

    omega_vs = MM_eps * psat / (z_press - psat)

    ! WMO definition of relhum
    !rel2spechumwmo = 1._dp/(1._dp/(z_relhum * omega_vs) + 1._dp)
    ! same but clearer:
    rel2spechumwmo = z_relhum * omega_vs / (1._dp + z_relhum * omega_vs)
    ! mz_hr_20100704- 
 
  END FUNCTION rel2spechumwmo
! --------------------------------------------------------------------- 

! mz_hr_20100608+
! ---------------------------------------------------------------------
  REAL(dp) FUNCTION rel2spechum(status, z_relhum, z_temp, z_press, &
                                   z_l_psat_emac) ! mz_hr_20100704
    !                                 z_l_psatf) !mz_hr_20100608

    ! calculates specific humidity from relative humidity
    ! standard definition relhum = p_v/p_vs

    ! molar mass of water (vapour) / molar mass of dry air
    ! mz_hr_20100705+ 
    !USE messy_main_constants_mem, ONLY: MM_eps
    USE messy_main_constants_mem, ONLY: MM_eps, Rd, Rv
    ! mz_hr_20100705- 

    IMPLICIT NONE

    INTRINSIC INT

    !I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum    ! 
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    ! mz_hr_20100704+ 
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac
    ! mz_hr_20100704- 

    !LOCAL
    REAL(dp) ::  omega_v, psat

    status = 0

    ! mz_hr_20100705+ 
    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucua(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
        psat = psat_mk(z_temp)
      ENDIF
    ELSE
      psat = psat_mk(z_temp)
    ENDIF
    !mz_hr_20100608-

    ! accounting for press = press(humid atmosphere)
    omega_v = MM_eps * z_relhum * psat / (z_press - z_relhum * psat)
    rel2spechum = omega_v / (1._dp + omega_v)
    ! same in one step
    !rel2spechum = &
    !    MM_eps * z_relhum * psat / (z_press + z_relhum*psat*(MM_eps-1._dp))
    ! wrong since press = press(humid air):
    !rel2spechum = z_relhum * MM_eps * psat / z_press
    ! mz_hr_20100705- 
 
  END FUNCTION rel2spechum
! --------------------------------------------------------------------- 
! mz_hr_20100608-

! mz_rs_20100902+
! --------------------------------------------------------------------- 
  REAL(dp) FUNCTION rh2mr(status, p_relhum, p_temp, p_press, &
                          p_l_psat_emac, p_l_relhum_wmo)
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: p_relhum     ! [0-1]
    REAL(dp), INTENT(IN) :: p_temp       ! K
    REAL(dp), INTENT(IN) :: p_press      ! Pa
    LOGICAL, INTENT(IN), OPTIONAL :: p_l_psat_emac
    LOGICAL, INTENT(IN), OPTIONAL :: p_l_relhum_wmo
    ! local:
    LOGICAL :: z_l_psat_emac, z_l_relhum_wmo

    IF (PRESENT(p_l_psat_emac)) THEN
      z_l_psat_emac = p_l_psat_emac
    ELSE
      z_l_psat_emac  = .FALSE.
    ENDIF

    IF (PRESENT(p_l_relhum_wmo)) THEN
      z_l_relhum_wmo = p_l_relhum_wmo
    ELSE
      z_l_relhum_wmo = .FALSE.
    ENDIF

    IF (z_l_relhum_wmo) THEN
      ! using relhum = (mixing ratio H2O)/(sat. mixing ratio H2O): WMO def.
      rh2mr = relhumwmo2mr(status,p_relhum,p_temp,p_press,z_l_psat_emac)
    ELSE
      ! using relhum = p(H2O)/ps(H2O)
      rh2mr = relhum2mr(status,p_relhum,p_temp,p_press,z_l_psat_emac)
    ENDIF
  END FUNCTION rh2mr
! ---------------------------------------------------------------------
! mz_rs_20100902-

! mz_hr_20100704+
! ---------------------------------------------------------------------
  REAL(dp) FUNCTION relhum2mr(status, z_relhum, z_temp, z_press, &
                                z_l_psat_emac)

    ! calculates water mixing ratio (mol H2O)/(mol dryair) as function
    ! of traditional relative humidity:
    ! relhum = p(H2O)/psat(H2O)

    ! molar mass of water (vapour) / molar mass of dry air
    ! mz_hr_20100704+ 
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, Rv, Rd
    ! mz_hr_20100704- 

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum     ! [0-1]
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    LOGICAL,  INTENT(IN) :: z_l_psat_emac

    ! local
    REAL(dp) :: psat ! H2O sat pressure [Pa]

    status = 0

    ! mz_hr_20100706+ 
    !IF ((z_relhum >= 1._dp) .AND. (z_relhum <= 1.1_dp)) THEN
    IF (z_relhum > 1._dp) THEN
    ! mz_hr_20100706- 
      WRITE(*,'(A, F10.3)') '  WARNING relhum2mr: relhum > 1 : ', &
        z_relhum
    ! mz_hr_20100706+ not so strict for traditional relhum
    !ELSEIF ((z_relhum > 1.1_dp)) THEN
    !  WRITE(*,*) 'ERROR relhum2mr: rel. humidity > 1.1 : ', &
    !    z_relhum
    !  status = 1
    !  RETURN
    ! mz_hr_20100706- 
    ELSEIF (z_relhum < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhum2mr: relhum < 0 : ', &
        z_relhum
      status = 1
      RETURN
    ENDIF

    IF (z_l_psat_emac) THEN
      psat = tlucua(INT(z_temp*1000._dp))*Rv/Rd
    ELSE
      psat = psat_mk(z_temp)
    ENDIF
  
    relhum2mr = z_relhum * psat / (z_press - z_relhum*psat)    


    IF ((relhum2mr >= 1._dp)) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhum2mr: mixing ratio > 1.0 : ', &
        relhum2mr
      relhum2mr = FLAGGED_BAD
      status = 1
      RETURN
    ELSEIF (relhum2mr < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhum2mr: mixing ratio < 0 : ', &
        relhum2mr
      relhum2mr = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION relhum2mr
! --------------------------------------------------------------------- 
! mz_hr_20100704-

! mz_hr_20100705+
! ---------------------------------------------------------------------
  REAL(dp) FUNCTION relhumwmo2mr(status, z_relhum, z_temp, z_press, &
                                z_l_psat_emac)

    ! calculates water mixing ratio (mol H2O)/(mol dryair) as function
    ! of traditional relative humidity:
    ! relhum = w(H2O)/wsat(H2O)
    ! w: (kg H2O)/(kg dryair)

    ! molar mass of water (vapour) / molar mass of dry air
    ! mz_hr_20100704+ 
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, Rv, Rd
    ! mz_hr_20100704- 

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum     ! [0-1]
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    LOGICAL,  INTENT(IN) :: z_l_psat_emac

    ! local
    REAL(dp) :: psat ! H2O sat pressure [Pa]

    status = 0

    IF ((z_relhum > 1._dp) .AND. (z_relhum <= 1.1_dp)) THEN
      WRITE(*,'(A, F10.3)') '  WARNING relhumwmo2mr: relhum > 1 (<= 1.1) : ', &
        z_relhum
    ELSEIF ((z_relhum > 1.1_dp)) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhumwmo2mr: relhum > 1.1 : ', &
        z_relhum
      status = 1
      RETURN
    ELSEIF (z_relhum < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhumwmo2mr: relhum < 0 : ', &
        z_relhum
      status = 1
      RETURN
    ENDIF

    IF (z_l_psat_emac) THEN
      ! WMO definition of relative humidity always above liquid (w, l)
      psat = tlucuaw(INT(z_temp*1000._dp))*Rv/Rd
    ELSE
       ! l_relhum_wmo=.TRUE. since here WMO function
      psat = psat_mk(z_temp,.TRUE.)
    ENDIF
  
    relhumwmo2mr = z_relhum * psat / (z_press - psat)    

    IF ((relhumwmo2mr >= 1._dp)) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhumwmo2mr: mixing ratio > 1.0 : ', &
        relhumwmo2mr
      relhumwmo2mr = FLAGGED_BAD
      status = 1
      RETURN
    ELSEIF (relhumwmo2mr < 0._dp) THEN
      WRITE(*,'(A, F10.3)') 'ERROR relhumwmo2mr: mixing ratio < 0 : ', &
        relhumwmo2mr
      relhumwmo2mr = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION relhumwmo2mr
! --------------------------------------------------------------------- 
! mz_hr_20100705-

! mz_hr_20100705+
! --------------------------------------------------------------------- 
  REAL(dp) FUNCTION cair_trad(status, z_relhum, z_temp, z_press, z_l_psat_emac)

    ! based on traditional definition of relative humidity
    ! relhum = p_v / psat

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, Rv, Rd, R_gas, N_A

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum     ! [0-1]
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac

    ! local
    REAL(dp) :: psat

    status = 0

    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucua(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
        psat = psat_mk(z_temp)
      ENDIF
    ELSE
      psat = psat_mk(z_temp)
    ENDIF
    !print *, 'HERE cair_trad: psat_mk   = ', psat_mk(z_temp)
    !print *, 'HERE cair_trad: psat_emac = ', tlucua(INT(z_temp*1000._dp))*Rv/Rd

    cair_trad = (N_A/1.E6_dp) * z_press / (R_gas*z_temp) * (1._dp-z_relhum*psat/z_press)
    !print *, 'HERE cair_trad: cair_trad = ', cair_trad
    !print *, 'HERE cair_trad: cair old  = ', (N_A/1.E6_dp) * z_press / (R_gas*z_temp)

    IF (cair_trad <= 0) THEN
      WRITE(*,'(A, F10.3)') 'ERROR cair_trad: cair_trad < 0, cair_trad = ', cair_trad
      cair_trad = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION cair_trad
! --------------------------------------------------------------------- 
! mz_hr_20100705-

! mz_hr_20100705+
! --------------------------------------------------------------------- 
  REAL(dp) FUNCTION cair_wmo(status, z_relhum, z_temp, z_press, z_l_psat_emac)

    ! based on WMO definition of relative humidity
    ! relhum = omega_v/omega_vs
    ! omega_v: water vapor mass mixing ratio (kg H2O)/(kg dry air)
    ! omega_vs: saturation mass mixing ratio
    ! psat always above water

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, Rv, Rd, R_gas, N_A

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: z_relhum     ! [0-1]
    REAL(dp), INTENT(IN) :: z_temp       ! K
    REAL(dp), INTENT(IN) :: z_press      ! Pa
    LOGICAL, INTENT(IN), OPTIONAL :: z_l_psat_emac

    ! local
    REAL(dp) :: psat

    status = 0

    IF (PRESENT(z_l_psat_emac)) THEN
      IF (z_l_psat_emac) THEN
        psat = tlucuaw(INT(z_temp*1000._dp))*Rv/Rd
      ELSE
         ! l_relhum_wmo=.TRUE. since here WMO function
        psat = psat_mk(z_temp,.TRUE.)
      ENDIF
    ELSE
       ! l_relhum_wmo=.TRUE. since here WMO function
      psat = psat_mk(z_temp,.TRUE.)
    ENDIF

    cair_wmo = (N_A/1.E6_dp) * z_press / (R_gas*z_temp) &
      * (z_press-psat)/(z_press + psat*(z_relhum-1._dp))

    IF (cair_wmo <= 0) THEN
      WRITE(*,'(A, F10.3)') 'ERROR cair_wmo: cair_wmo < 0, cair_wmo = ', &
        cair_wmo
      cair_wmo = FLAGGED_BAD
      status = 1
      RETURN
    ENDIF

  END FUNCTION cair_wmo
! --------------------------------------------------------------------- 
! mz_hr_20100705-

! mz_ab_20100517+
! --------------------------------------------------------------------- 
  FUNCTION remap_bounds1(lb1,array) RESULT(ptr)
    INTEGER,                   INTENT(IN)          :: lb1
    REAL(dp), DIMENSION(lb1:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:),                POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds1
  
  FUNCTION remap_bounds2(lb1,lb2,array) RESULT(ptr)
    INTEGER,                        INTENT(IN)          :: lb1,lb2
    REAL(dp), DIMENSION(lb1:,lb2:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:),                   POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds2
  
  FUNCTION remap_bounds3(lb1,lb2,lb3,array) RESULT(ptr)
    INTEGER,                             INTENT(IN)          :: lb1,lb2,lb3
    REAL(dp), DIMENSION(lb1:,lb2:,lb3:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:,:),                      POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds3
  
  FUNCTION remap_bounds4(lb1,lb2,lb3,lb4,array) RESULT(ptr)
    INTEGER,                                  INTENT(IN)  :: lb1,lb2,lb3,lb4
    REAL(dp), DIMENSION(lb1:,lb2:,lb3:,lb4:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:,:,:),                         POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds4
! --------------------------------------------------------------------- 
! mz_ab_20100517-

! mz_ab_20100831+
! --------------------------------------------------------------------- 
  SUBROUTINE full2half(full,half,press,pressi)

    ! This subroutine returns a variable on half level pressures
    ! given a 3-D variable (but single jrow) on full level pressures.
    ! 
    ! half(:,k) lies between full(:,k) and full(:,k+1), so
    ! it is the lower boundary. Note that pressi is defined for upper
    ! boundary, therefore we use pressi(k+1) below.

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: full &
                                            , press, pressi
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: half 

    INTEGER :: nlev, k

    nlev = SIZE(press,2)

    DO k=1, nlev-1
       half(:,k)=(pressi(:,k+1)-press (:,k  ))/&
                   (press (:,k+1)-press (:,k  )) * full(:,k)   + &
                   (press (:,k+1)-pressi(:,k+1))/&
                   (press (:,k+1)-press (:,k  )) * full(:,k+1)
    ENDDO

    ! It is not possible to calculate this for lowest level, so
    ! this is set to the lowest-but-1 level
    half(:,nlev)=half(:,nlev-1)

  END SUBROUTINE full2half
! --------------------------------------------------------------------- 

! mz_ab_20111125+
! --------------------------------------------------------------------- 
    SUBROUTINE CalcTimeAngles(utsec, sda, ssa, rlt, csza, &
                              sza, daylen_sec, rlat, phi)

      USE messy_main_constants_mem, ONLY: DTR

      REAL(dp), INTENT(IN)  :: utsec, daylen_sec
      REAL(dp), INTENT(IN)  :: sda
      REAL(dp), INTENT(IN)  :: rlat(:), phi(:)

      REAL(dp), INTENT(OUT) :: rlt(:), ssa(:)       &
                             , csza(:,:), sza(:,:)

      INTEGER :: lat_dim, lon_dim, l, m

      lat_dim = SIZE(rlat,1)
      lon_dim = SIZE(phi,1)

      LONLOOP: DO l = 1, lon_dim

         ! Solar siderial angle

         ssa(l)    = phi(l) + 360.*(utsec/daylen_sec) + 180.
         IF ( ssa(l) >= 360.0 ) ssa(l) = ssa(l) - 360.0

         ! Local time

         rlt(l) = 180.0 + ssa(l)
         IF ( rlt(l) > 360.0 ) rlt(l) = rlt(l) - 360.0
         rlt(l) = rlt(l)*DTR

         LATLOOP: DO m = 1 , lat_dim

            ! Calculate solar zenith angle

            csza(m,l)  = -COS(rlat(m))*COS(sda)*COS(rlt(l))+ &
                 SIN(rlat(m))*SIN(sda)
            sza(m,l)   =  ACOS(csza(m,l))

         END DO LATLOOP

      END DO LONLOOP


    END SUBROUTINE CalcTimeAngles
! --------------------------------------------------------------------- 
! mz_ab_20111125-

! --------------------------------------------------------------------- 
  SUBROUTINE Spline1D(X,Y,N,YP1,YPN,Y2,natspline)
    ! Numerical recipies Spline1d routine, see also partner split1d
    !
    ! GIVEN ARRAYS X AND Y OF LENGTH N CONTAINING A TABULATED FUNCTION
    ! IE YI=F(XI), WITH X1 < X2 < ... < XN, AND GIVEN VALUES YP1 AND
    ! YPN FOR THE FIRST DERIVATIVE OF THE INTERPOLATIONG FUNCTIONS AT
    ! POINTS 1 AND N RESPECTIVELY, THIS ROUTINE RETURNS AN ARRAY Y2 OF
    ! LENGTH N WHICH CONTAINS THE SECOND DERIVATIVES OF THE INTERPOLATING
    ! FUNCTION AT THE TABULATED POINTS XI.  IF YP2 AND/OR YPN ARE EQUAL
    ! TO 1E30 OR LARGER, THE ROUTINE IS SIGNALLED TO SET THE CORRESPONDING
    ! BOUNDARY CONDITIONS FOR A NATURAL SPLINE, WITH ZERO SECOND DERIVATIVES
    ! ON THAT BOUNDARY.

    REAL(dp), INTENT(IN) :: X(:), Y(:)
    INTEGER,  INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: yp1,ypn
    REAL(dp), INTENT(OUT) :: Y2(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: natspline

    ! Local
    INTEGER, PARAMETER :: NMAX=400
    REAL(dp) :: U(NMAX)
    REAL(dp) :: xa, y2a, ya
    REAL(dp) :: sig,p,qn,un
    INTEGER  :: i,k
    INTEGER  :: klo,khi
    REAL(dp) :: h,a,b
    LOGICAL  :: lnatspline = .FALSE.

    IF (PRESENT(natspline)) THEN
       lnatspline=natspline
    ELSE
       lnatspline=.FALSE.
    ENDIF

    IF (lnatspline) THEN
       Y2(1)=0._dp
       U(1)=0._dp
    ELSE
       Y2(1)=-0.5_dp
       U(1)=(3.0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    ENDIF

    DO I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2._dp
       Y2(I)=(SIG-1.0)/P
       U(I)=(6._dp*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))   &
            /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    ENDDO

    IF (lnatspline) THEN
       QN=0._dp
       UN=0._dp
    ELSE
       QN=0.5_dp
       UN=(3._dp/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    ENDIF

    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1._dp)
    DO K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    ENDDO

    RETURN

  END SUBROUTINE Spline1D
! --------------------------------------------------------------------- 

! --------------------------------------------------------------------- 
  SUBROUTINE Splint1D(XA,YA,Y2A,N,X,Y,status)

    ! Numerical recipies spline routine. See also Spline1d
      
      ! GIVEN THE ARRAYS XA AND YA OF LENGTH N, WHICH TABULATE A FUNCTION
      ! (WITH THE XAI'S IN ORDER), AND GIVEN THE ARRAY Y2A, WHICH IS THE
      ! OUTPUT FROM SPLINE ABOVE, AND GIVEN A VALUE OF X, THIS ROUTINE
      ! RETURNS A CUBIC-SPLINE INTERPOLATED VALUE Y.
            
    REAL(dp), INTENT(IN)  :: xa(:), ya(:), y2a(:)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: X
    REAL(dp), INTENT(OUT) :: Y
    INTEGER,  INTENT(OUT) :: status

    INTEGER  :: nmax
    REAL(dp) :: Y2,U
    REAL(dp) :: yp1,ypn,sig,p,qn,un
    INTEGER  :: i,k
    INTEGER  :: klo,khi
    REAL(dp) :: h,a,b

    status = 1 ! ERROR
      
    KLO=1
    KHI=N
    DO WHILE (KHI-KLO > 1) 
       K=(KHI+KLO)/2
       IF (XA(K) > X) THEN
          KHI=K
       ELSE
          KLO=K
       ENDIF
    ENDDO
    
    H=XA(KHI)-XA(KLO)
    IF (H == 0 ) RETURN
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI) +                                  &
         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0
    
    status = 0 ! NO ERROR
      
    RETURN
      
  END SUBROUTINE Splint1D
!---------------------------------------------------------------------
! mz_ab_20100831-

!---------------------------------------------------------------------
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)

    ! I/O
    INTEGER :: istart, istop, unit

    ! LOCAL
    LOGICAL        :: found, opened
    INTEGER        :: i
    CHARACTER(256) :: info

    found = .FALSE.
    DO i=istart,istop
       INQUIRE(unit=i,opened=opened)
       IF (.NOT.opened) THEN
          unit = i
          found = .TRUE.
          EXIT
       END IF
    END DO

    IF (.NOT. found) THEN
       WRITE(info,'(a,i2.2,a,i2.2,a)') &
         'No unit in range <',istart,':',istop,'> free.'
       !CALL error_bi(info,'find_next_free_unit')
    END IF

  END FUNCTION find_next_free_unit
! --------------------------------------------------------------------------

! mz_rs_20101118+
! --------------------------------------------------------------------------
  ELEMENTAL REAL(dp) FUNCTION density(press, temp, sphum)

    ! calculate air density in kg/m3

    USE messy_main_constants_mem, ONLY: R_gas, M_air, M_H2O

    REAL(dp), INTENT(in) :: press
    REAL(dp), INTENT(in) :: temp
    REAL(dp), INTENT(in) :: sphum

    density = press /( R_gas * &
      temp * (1._dp +(M_air / M_H2O -1._dp) *sphum))

  END FUNCTION density
!---------------------------------------------------------------------------

! --------------------------------------------------------------------------
  ELEMENTAL REAL(dp) FUNCTION layerthickness(geopot_u, geopot_l)

    USE messy_main_constants_mem,  ONLY: g

    REAL(dp), INTENT(in) :: geopot_u
    REAL(dp), INTENT(in) :: geopot_l
    !LOCAL

    layerthickness  = (geopot_u-geopot_l)/ g
  END FUNCTION layerthickness
! --------------------------------------------------------------------------
! mz_rs_20101118-

! mz_bk_20110707+
  ! --------------------------------------------------------------------------
  ! str2num : converts string to number of kinds real(dp), real(sp) or integer
  ! procedures : str2num(character(len=*) str, real(dp) out, integer err)
  !              str2num(character(len=*) str, real(sp) out, integer err)
  !              str2num(character(len=*) str, integer out, integer err)

  SUBROUTINE str2num_real_dp(str, out, err)
    CHARACTER(LEN=*), INTENT(in)   :: str
    REAL(dp), INTENT(out)          :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_real_dp

  SUBROUTINE str2num_real_sp(str, out, err)
    CHARACTER(LEN=*), INTENT(in)   :: str
    REAL(sp), INTENT(out)          :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_real_sp

  SUBROUTINE str2num_integer(str, out, err)
    CHARACTER(LEN=*), INTENT(in)   :: str
    INTEGER, INTENT(out)           :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_integer
  ! --------------------------------------------------------------------------
! mz_bk_20110707-

! ************************************************************************
END MODULE messy_main_tools
! ************************************************************************
