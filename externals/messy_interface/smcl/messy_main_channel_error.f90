! **********************************************************************
MODULE messy_main_channel_error
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

  ! op_bk_20130816+
  PUBLIC
  ! op_bk_20130816-

CONTAINS

  ! -------------------------------------------------------------------
  FUNCTION channel_error_str(status)

    USE messy_main_tools,         ONLY: int2str
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=STRLEN_VLONG)           :: channel_error_str
    INTEGER,                   INTENT(IN) :: status

    ! LOCAL
    CHARACTER(LEN=4) :: echar = '    '

    channel_error_str = ''

    CALL int2str(echar, status, cpad='0', cerr='*')

    SELECT CASE(echar)
       ! NO ERROR
    CASE('0000')
       channel_error_str = 'E'//echar//': NO ERROR'

       ! ATTRIBUTE ERRORS
    CASE('0801')
       channel_error_str = 'E'//echar//': ATTRIBUTE TYPE IS AMBIGUOUS'
    CASE('0802')
       channel_error_str = 'E'//echar//': ATTRIBUTE NAME TOO LONG'
    CASE('0803')
       channel_error_str = 'E'//echar//': CHARACTER ATTRIBUTE TOO LONG'
    CASE('0804')
       channel_error_str = 'E'//echar//': ATTRIBUTE EXISTS ALREADY'
    CASE('0805')
       channel_error_str = 'E'//echar//': ATTRIBUTE DOES NOT EXIST'
    CASE('0806')
       channel_error_str = 'E'//echar//': ATTRIBUTE TYPE IS UNKNOWN'
    CASE('0807')
       channel_error_str = 'E'//echar//': ATTRIBUTE LIST IS ALREADY ASSOCIATED'

       ! DIMENSION ERRORS
    CASE('0901')
       channel_error_str = 'E'//echar//': DIMENSION LENGTH INVALID'
    CASE('0902')
       channel_error_str = 'E'//echar//': DIMENSION NAME TOO LONG'
    CASE('0904')
       channel_error_str = 'E'//echar//': DIMENSION EXISTS ALREADY'
    CASE('0905')
       channel_error_str = 'E'//echar//': DIMENSION DOES NOT EXIST'
    CASE('0907')
       channel_error_str = 'E'//echar//': DIMENSION LIST IS EMPTY'
    CASE('0909')
       channel_error_str = 'E'//echar//': DIMENSION ELEMENT NOT ASSOCIATED'
    CASE('0910')
       channel_error_str = 'E'//echar//': DIMENSION ID INVALID'

       ! DIMENSION VARIABLE ERRORS
    CASE('0951')
       channel_error_str = 'E'//echar//': DIMENSION VARIABLE NAME TOO LONG'
    CASE('0952')
       channel_error_str = 'E'//echar//': DIMENSION VARIABLE EXISTS ALREADY'
    CASE('0953')
       channel_error_str = 'E'//echar//': DIMENSION VARIABLE DOES NOT EXIST'
    CASE('0954')
       channel_error_str = 'E'//echar//': VECTOR LENGTH NOT CONFORM WITH DIMENSION LENGTH'
    CASE('0955')
       channel_error_str = 'E'//echar//': DIMENSION VARIABLE UPDATE ONLY FOR TIME DIMENSION'

       ! MEMORY ERRORS
    CASE('1000')
       channel_error_str = 'E'//echar//': MEMORY ALLOCATION FAILED'
    CASE('1001')
       channel_error_str = 'E'//echar//': MEMORY DEALLOCATION FAILED'
    CASE('1002')
       channel_error_str = 'E'//echar//': INTERNAL DIMENSION COUNT ERROR'
    CASE('1003')
       channel_error_str = 'E'//echar//': INTERNAL REPRESENTATION COUNT ERROR'
    CASE('1004')
       channel_error_str = 'E'//echar//': INTERNAL CHANNEL COUNT ERROR'
    CASE('1005')
       channel_error_str = 'E'//echar//': POINTER (TRACER) IS ASSOCIATED'
    CASE('1006')

       ! REPRESENTATION ERRORS
    CASE('2001')
       channel_error_str = 'E'//echar//': REPRESENTATION NAME TOO LONG'
    CASE('2002')
       channel_error_str = 'E'//echar//': REPRESENTATION EXISTS ALREADY'
    CASE('2003')
       channel_error_str = 'E'//echar//': REPRESENTATION DOES NOT EXIST'
    CASE('2004')
       channel_error_str = 'E'//echar//': REPRESENTATION LIST EMPTY'
    CASE('2005')
       channel_error_str = 'E'//echar//': DIMENSION NAME IS EMPTY'
    CASE('2006')
       channel_error_str = 'E'//echar//': LOCAL LENGTH OF NON-DECOMPOSED DIMENSION DOES NOT EQUAL GLOBAL LENGTH'
    CASE('2007')
       channel_error_str = 'E'//echar//': LOCAL DIMENSION IS LARGER THAN GLOBAL DIMENSION'
    CASE('2008')
       channel_error_str = 'E'//echar//': DEFAULT DIMENSION IS NOT ASSOCIATED'
    CASE('2009')
       channel_error_str = 'E'//echar//': DIMENSION ID IS <= ZERO'
    CASE('2010')
       channel_error_str = 'E'//echar//': MAXIMUM RANK LIMIT EXCEEDED'
    CASE('2011')
       channel_error_str = 'E'//echar//': SYNTAX ERROR IN LINK STRING'
    CASE('2012')
       channel_error_str = 'E'//echar//': LINK STRING NOT CONFORM WITH RANK'
    CASE('2013')
       channel_error_str = 'E'//echar//': INVALID REPRESENTATION ID'
    CASE('2014')
       channel_error_str = 'E'//echar//': INVALID ORDER OF OUTPUT DIMENSIONS'
    CASE('2015')
       channel_error_str = 'E'//echar//': INVALID REPR_REORDER FLAG'
    CASE('2016')
       channel_error_str = 'E'//echar//': RANK < 0'
    CASE('2017')
       channel_error_str = 'E'//echar//': POINTER p0 COULD NOT BE ASSOCIATED'
    CASE('2018')
       channel_error_str = 'E'//echar//': POINTER p1 COULD NOT BE ASSOCIATED'
    CASE('2019')
       channel_error_str = 'E'//echar//': POINTER p2 COULD NOT BE ASSOCIATED'
    CASE('2020')
       channel_error_str = 'E'//echar//': POINTER p3 COULD NOT BE ASSOCIATED'
    CASE('2021')
       channel_error_str = 'E'//echar//': POINTER p4 COULD NOT BE ASSOCIATED'
    CASE('2022')
       channel_error_str = 'E'//echar//': UNKNOWN REPRESENTATION DECOMPOSITION TYPE / RANK'
    CASE('2023')
       channel_error_str = 'E'//echar//': REPRESENTATION DECOMPOSITION TABLE EXISTS ALREADY'
    CASE('2024')
       channel_error_str = 'E'//echar//': REPRESENTATION SEGMENTATION MISMATCH'
    CASE('2025')
       channel_error_str = 'E'//echar//': REPRESENTATION DECOMPOSITION TABLE DOES NOT EXIST'
    CASE('2026')
       channel_error_str = 'E'//echar//': REPRESENTATION DECOMPOSITION MISMATCH'
    CASE('2100')
       channel_error_str = 'E'//echar//': NON-GENERIC FAST REORDER NOT IMPLEMENTED FOR SPECIFIED PERMUTATION'

       ! CHANNEL ERRORS
    CASE('3000')
       channel_error_str = 'E'//echar//': OPERATION NOT POSSIBLE AFTER CHANNEL FIXATION'
    CASE('3001')
       channel_error_str = 'E'//echar//': CHANNEL NAME TOO LONG'
    CASE('3002')
       channel_error_str = 'E'//echar//': CHANNEL EXISTS ALREADY'
    CASE('3003')
       channel_error_str = 'E'//echar//': CHANNEL (NAME) DOES NOT EXIST'
    CASE('3004')
       channel_error_str = 'E'//echar//': CHANNEL (ID) DOES NOT EXIST'
    CASE('3005')
       channel_error_str = 'E'//echar//': NUMBER OF CHANNELS MISMATCH'
    CASE('3006')
       channel_error_str = 'E'//echar//': CHANNEL POINTER NOT ASSOCIATED'
    CASE('3007')
       channel_error_str = 'E'//echar//': INVALID CHANNEL ID'
    CASE('3008')
       channel_error_str = 'E'//echar//': CHANNEL PRIMARY MEMORY ERROR'
    CASE('3009')
       channel_error_str = 'E'//echar//': CHANNEL SECONDARY MEMORY ERROR'
    CASE('3010')
       channel_error_str = 'E'//echar//': CHANNEL NAME IS EMPTY'


       ! CHANNEL OBJECT ERRORS
    CASE('3101')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT NAME TOO LONG'
    CASE('3102')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT EXISTS ALREADY'
    CASE('3103')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT DOES NOT EXIST'
    CASE('3104')
       channel_error_str = 'E'//echar//': CHANNEL DEFAULT REPRESENTATION IS UNDEFINED'
    CASE('3105')
       channel_error_str = 'E'//echar//': SHAPE OF MEMORY NOT CONFORM WITH REPRESENTATION'
    CASE('3106')
       channel_error_str = 'E'//echar//': UNKNOW FLAG'
    CASE('3107')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT POINTER NOT ASSOCIATED'
    CASE('3108')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT PRIMARY MEMORY ERROR'
    CASE('3109')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT SECONDARY MEMORY ERROR'
    CASE('3110')
       channel_error_str = 'E'//echar//': CHANNEL OBJECT NAME IS EMPTY'

       ! CHANNEL I/O ERRORS
    CASE('3200')
       channel_error_str = 'E'//echar//': CHANNEL OUTPUT FILE TYPE UNDEFINED'
    CASE('3201')
       channel_error_str = 'E'//echar//': CHANNEL OUTPUT FILE TYPE NOT YET IMPLEMENTED'
    CASE('3202')
       channel_error_str = 'E'//echar//': CHANNEL OUTPUT FILE TYPE UNKNOWN'
    CASE('3203')
       channel_error_str = 'E'//echar//': SECONDARY DATA TYPE UNKNOWN'
    CASE('3204')
       channel_error_str = 'E'//echar//': RESTART TIME STRING TOO LONG'
    CASE('3205')
       channel_error_str = 'E'//echar//': NO INPUT OF OUTPUT FILES'
    CASE('3206')
       channel_error_str = 'E'//echar//': RESTART FILE REQUIRED BUT NOT PRESENT'
    CASE('3207')
       channel_error_str = 'E'//echar//': MISSING ATTRIBUTE IN RESTART FILE'
    CASE('3208')
       channel_error_str = 'E'//echar//': RESTART ATTRIBUTE HAS NON-SCALAR RANK'
    CASE('3209')
       channel_error_str = 'E'//echar//': RESTART CHARACTER ATTRIBUTE TOO LONG'
    CASE('3210')
       channel_error_str = 'E'//echar//': RESTART ATTRIBUTE MISMATCH'
    CASE('3211')
       channel_error_str = 'E'//echar//': RESTART VARIABLE REQUIRED BUT NOT PRESENT'
    CASE('3212')
       channel_error_str = 'E'//echar//': CHANNEL OUTPUT FILE TYPE NOT AVAILABLE'

       ! FORMAT SPECIFIC I/O ERRORS
    CASE('4000')
       channel_error_str = 'E'//echar//': NETCDF ERROR'
    CASE('4001')
       channel_error_str = 'E'//echar//': UNDEFINED FILE-ID'
    CASE('4002')
       channel_error_str = 'E'//echar//': PARALLEL NETCDF ERROR'
    CASE('4003')
       channel_error_str = 'E'//echar//': UNKNOWN PRECISION FLAG FOR (P)NETCDF'

    CASE DEFAULT
       channel_error_str = 'E'//echar//': UNKONW ERROR STATUS'
    END SELECT

  END FUNCTION channel_error_str
  ! -------------------------------------------------------------------


! **********************************************************************
END MODULE messy_main_channel_error
! **********************************************************************
