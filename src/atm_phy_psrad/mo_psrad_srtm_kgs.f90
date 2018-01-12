!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_srtm_kgs

  USE mo_psrad_general, ONLY : wp, nbndsw
  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: &
    ng(nbndsw) = (/6, 12, 8, 8, 10, 10, 2, 10, 8, 6, 6, 8, 6, 12/), &
    no(nbndsw) = &
      (/16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16/)

  REAL(wp) :: kao16(9,5,13,no(1))
  REAL(wp) :: kbo16(5,47,no(1))
  REAL(wp) :: selfrefo16(10,no(1)), forrefo16(3,no(1))
  REAL(wp) :: sfluxrefo16(no(1))
  REAL(wp) :: rayl16

  REAL(wp) :: ka16(9,5,13,ng(1)) , absa16(585,ng(1))
  REAL(wp) :: kb16(5,47,ng(1)), absb16(235,ng(1))
  REAL(wp) :: selfref16(10,ng(1)), forref16(3,ng(1))
  REAL(wp) :: sfluxref16(ng(1))

  EQUIVALENCE (ka16(1,1,1,1),absa16(1,1)), (kb16(1,1,1),absb16(1,1))

  REAL(wp) :: kao17(9,5,13,no(2))
  REAL(wp) :: kbo17(5,5,47,no(2))
  REAL(wp) :: selfrefo17(10,no(2)), forrefo17(4,no(2))
  REAL(wp) :: sfluxrefo17(no(2),5)
  REAL(wp) :: rayl17

  REAL(wp) :: ka17(9,5,13,ng(2)) , absa17(585,ng(2))
  REAL(wp) :: kb17(5,5,47,ng(2)), absb17(1175,ng(2))
  REAL(wp) :: selfref17(10,ng(2)), forref17(4,ng(2))
  REAL(wp) :: sfluxref17(ng(2),5)

  EQUIVALENCE (ka17(1,1,1,1),absa17(1,1)), (kb17(1,1,1,1),absb17(1,1))

  REAL(wp) :: kao18(9,5,13,no(3))
  REAL(wp) :: kbo18(5,47,no(3))
  REAL(wp) :: selfrefo18(10,no(3)), forrefo18(3,no(3))
  REAL(wp) :: sfluxrefo18(no(3),9)
  REAL(wp) :: rayl18

  REAL(wp) :: ka18(9,5,13,ng(3)), absa18(585,ng(3))
  REAL(wp) :: kb18(5,47,ng(3)), absb18(235,ng(3))
  REAL(wp) :: selfref18(10,ng(3)), forref18(3,ng(3))
  REAL(wp) :: sfluxref18(ng(3),9)

  EQUIVALENCE (ka18(1,1,1,1),absa18(1,1)), (kb18(1,1,1),absb18(1,1))

  REAL(wp) :: kao19(9,5,13,no(4))
  REAL(wp) :: kbo19(5,47,no(4))
  REAL(wp) :: selfrefo19(10,no(4)), forrefo19(3,no(4))
  REAL(wp) :: sfluxrefo19(no(4),9)
  REAL(wp) :: rayl19

  REAL(wp) :: ka19(9,5,13,ng(4)), absa19(585,ng(4))
  REAL(wp) :: kb19(5,47,ng(4)), absb19(235,ng(4))
  REAL(wp) :: selfref19(10,ng(4)), forref19(3,ng(4))
  REAL(wp) :: sfluxref19(ng(4),9)

  EQUIVALENCE (ka19(1,1,1,1),absa19(1,1)), (kb19(1,1,1),absb19(1,1))

  REAL(wp) :: kao20(5,13,no(5))
  REAL(wp) :: kbo20(5,47,no(5))
  REAL(wp) :: selfrefo20(10,no(5)), forrefo20(4,no(5))
  REAL(wp) :: sfluxrefo20(no(5))
  REAL(wp) :: absch4o20(no(5))

  REAL(wp) :: rayl20

  REAL(wp) :: ka20(5,13,ng(5)), absa20(65,ng(5))
  REAL(wp) :: kb20(5,47,ng(5)), absb20(235,ng(5))
  REAL(wp) :: selfref20(10,ng(5)), forref20(4,ng(5))
  REAL(wp) :: sfluxref20(ng(5))
  REAL(wp) :: absch420(ng(5))

  EQUIVALENCE (ka20(1,1,1),absa20(1,1)), (kb20(1,1,1),absb20(1,1))

  REAL(wp) :: kao21(9,5,13,no(6))
  REAL(wp) :: kbo21(5,5,47,no(6))
  REAL(wp) :: selfrefo21(10,no(6)), forrefo21(4,no(6))
  REAL(wp) :: sfluxrefo21(no(6),9)
  REAL(wp) :: rayl21

  REAL(wp) :: ka21(9,5,13,ng(6)), absa21(585,ng(6))
  REAL(wp) :: kb21(5,5,47,ng(6)), absb21(1175,ng(6))
  REAL(wp) :: selfref21(10,ng(6)), forref21(4,ng(6))
  REAL(wp) :: sfluxref21(ng(6),9)

  EQUIVALENCE (ka21(1,1,1,1),absa21(1,1)), (kb21(1,1,1,1),absb21(1,1))

  REAL(wp) :: kao22(9,5,13,no(7))
  REAL(wp) :: kbo22(5,47,no(7))
  REAL(wp) :: selfrefo22(10,no(7)), forrefo22(3,no(7))
  REAL(wp) :: sfluxrefo22(no(7),9)

  REAL(wp) :: rayl22

  REAL(wp) :: ka22(9,5,13,ng(7)), absa22(585,ng(7))
  REAL(wp) :: kb22(5,47,ng(7)), absb22(235,ng(7))
  REAL(wp) :: selfref22(10,ng(7)), forref22(3,ng(7))
  REAL(wp) :: sfluxref22(ng(7),9)

  EQUIVALENCE (ka22(1,1,1,1),absa22(1,1)), (kb22(1,1,1),absb22(1,1))

  REAL(wp) :: kao23(5,13,no(8))
  REAL(wp) :: selfrefo23(10,no(8)), forrefo23(3,no(8))
  REAL(wp) :: sfluxrefo23(no(8))
  REAL(wp) :: raylo23(no(8))

  REAL(wp) :: ka23(5,13,ng(8)), absa23(65,ng(8))
  REAL(wp) :: selfref23(10,ng(8)), forref23(3,ng(8))
  REAL(wp) :: sfluxref23(ng(8)), rayl23(ng(8))

  EQUIVALENCE (ka23(1,1,1),absa23(1,1))

  REAL(wp) :: kao24(9,5,13,no(9))
  REAL(wp) :: kbo24(5,47,no(9))
  REAL(wp) :: selfrefo24(10,no(9)), forrefo24(3,no(9))
  REAL(wp) :: sfluxrefo24(no(9),9)
  REAL(wp) :: abso3ao24(no(9)), abso3bo24(no(9))
  REAL(wp) :: raylao24(no(9),9), raylbo24(no(9))

  REAL(wp) :: ka24(9,5,13,ng(9)), absa24(585,ng(9))
  REAL(wp) :: kb24(5,47,ng(9)), absb24(235,ng(9))
  REAL(wp) :: selfref24(10,ng(9)), forref24(3,ng(9))
  REAL(wp) :: sfluxref24(ng(9),9)
  REAL(wp) :: abso3a24(ng(9)), abso3b24(ng(9))
  REAL(wp) :: rayla24(ng(9),9), raylb24(ng(9))

  EQUIVALENCE (ka24(1,1,1,1),absa24(1,1)), (kb24(1,1,1),absb24(1,1))

  REAL(wp) :: kao25(5,13,no(10))
  REAL(wp) :: sfluxrefo25(no(10))
  REAL(wp) :: abso3ao25(no(10)), abso3bo25(no(10))
  REAL(wp) :: raylo25(no(10))

  REAL(wp) :: ka25(5,13,ng(10)), absa25(65,ng(10))
  REAL(wp) :: sfluxref25(ng(10))
  REAL(wp) :: abso3a25(ng(10)), abso3b25(ng(10))
  REAL(wp) :: rayl25(ng(10))

  EQUIVALENCE (ka25(1,1,1),absa25(1,1))

  REAL(wp) :: sfluxrefo26(no(11))
  REAL(wp) :: raylo26(no(11))

  REAL(wp) :: sfluxref26(ng(11))
  REAL(wp) :: rayl26(ng(11))

  REAL(wp) :: kao27(5,13,no(12))
  REAL(wp) :: kbo27(5,47,no(12))
  REAL(wp) :: sfluxrefo27(no(12))
  REAL(wp) :: raylo27(no(12))

  REAL(wp) :: ka27(5,13,ng(12)), absa27(65,ng(12))
  REAL(wp) :: kb27(5,47,ng(12)), absb27(235,ng(12))
  REAL(wp) :: sfluxref27(ng(12))
  REAL(wp) :: rayl27(ng(12))

  EQUIVALENCE (ka27(1,1,1),absa27(1,1)), (kb27(1,1,1),absb27(1,1))

  REAL(wp) :: kao28(9,5,13,no(13))
  REAL(wp) :: kbo28(5,5,47,no(13))
  REAL(wp) :: sfluxrefo28(no(13),5)
  REAL(wp) :: rayl28

  REAL(wp) :: ka28(9,5,13,ng(13)), absa28(585,ng(13))
  REAL(wp) :: kb28(5,5,47,ng(13)), absb28(1175,ng(13))
  REAL(wp) :: sfluxref28(ng(13),5)

  EQUIVALENCE (ka28(1,1,1,1),absa28(1,1)), (kb28(1,1,1,1),absb28(1,1))

  REAL(wp) :: kao29(5,13,no(14))
  REAL(wp) :: kbo29(5,47,no(14))
  REAL(wp) :: selfrefo29(10,no(14)), forrefo29(4,no(14))
  REAL(wp) :: sfluxrefo29(no(14))
  REAL(wp) :: absh2oo29(no(14)), absco2o29(no(14))
  REAL(wp) :: rayl29

  REAL(wp) :: ka29(5,13,ng(14)), absa29(65,ng(14))
  REAL(wp) :: kb29(5,47,ng(14)), absb29(235,ng(14))
  REAL(wp) :: selfref29(10,ng(14)), forref29(4,ng(14))
  REAL(wp) :: sfluxref29(ng(14))
  REAL(wp) :: absh2o29(ng(14)), absco229(ng(14))

  EQUIVALENCE (ka29(1,1,1),absa29(1,1)), (kb29(1,1,1),absb29(1,1))

END MODULE mo_psrad_srtm_kgs

