!>
!! @brief bgc tracer parameters
!!
!! Definition of parameters
!!
!! @author N.N.
!!
!! See the LICENSE and the WARRANTY conditions.
!!
MODULE mo_param1_bgc
  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC


  ! advected tracers
  INTEGER, PARAMETER :: i_base_adv = 15,              &
       &                isco212    = 1,               &
       &                ialkali    = 2,               &
       &                iphosph    = 3,               &
       &                iano3      = 4,               &
       &                igasnit    = 5,               &
       &                iphy       = 6,               &
       &                izoo       = 7,               &
       &                icya       = 8,               &
       &                ioxygen    = 9,               &
       &                isilica    = 10,              &
       &                idoc       = 11,              &
       &                ian2o      = 12,              &
       &                idet       = 13,              &         
       &                iiron      = 14,              &               
       &                idms       = 15               


  ! total number of advected tracers
  INTEGER, PARAMETER :: ntraad=i_base_adv

  ! non-advected (fast sinking) tracers
  INTEGER, PARAMETER ::  icalc    = ntraad+1,         &
       &                 iopal    = ntraad+2,         &
       &                 idust    = ntraad+3,         &
       &                 i_base   = 3  


  INTEGER, PARAMETER :: n_bgctra = ntraad + i_base 



  ! atmosphere
  INTEGER, PARAMETER :: iatmco2    = 1,               &
       &                iatmo2     = 2,               &
       &                iatmn2     = 3,               &
       &                i_base_atm = 3

  INTEGER, PARAMETER :: natm = i_base_atm 

  ! sediment
  INTEGER, PARAMETER :: issso12   = 1,                &
       &                isssc12   = 2,                &
       &                issssil   = 3,                &
       &                issster   = 4,                &
       &                nsss_base = 4
  INTEGER, PARAMETER :: nsedtra = nsss_base 

  ! pore water tracers, index must be the same as for ocetra otherwise problems in dipowa.f90!
  INTEGER, PARAMETER :: ipowaic    = 1,               &
       &                ipowaal    = 2,               &
       &                ipowaph    = 3,               &
       &                ipowaox    = 4,               &
       &                ipown2     = 5,               &
       &                ipowno3    = 6,               &
       &                ipowasi    = 7,               &
       &                ipowafe    = 8,               &
       &                npowa_base = 8
  
  INTEGER, PARAMETER :: npowtra=npowa_base

 ! Diagnostics

  INTEGER, PARAMETER :: kphosy    = 1,               &
      &                 ksred     = 2,               &
      &                 kdenit    = 3,               &
      &                 kremin    = 4,               &
      &                 kgraz     = 5,               &
      &                 knfix     = 6,               &
      &                 kbacfra   = 7,               &
      &                 kpho_cya  = 8,               &
      &                 kcyaloss  = 9,               &
      &                 kn2b      = 10,               &
      &                 kh2ob     = 11,               &
      &                 kdelsil   = 12,               &
      &                 kdelcar   = 13,               &
      &                 kbacfrac  = 14,               &
      &                 kdmsprod  = 15,               &
      &                 kdmsbac   = 16,               &
      &                 kdmsuv    = 17,               &
      &                 keuexp    = 18,               &
      &                 kplim     = 19,               &
      &                 kflim     = 20,               &
      &                 knlim     = 21,               &
      &                 kgraton   = 22,               &
      &                 kexudp    = 23,               &
      &                 kexudz    = 24,               &
      &                 kzdy      = 25,               &
      &                 kpdy      = 26,               &
      &                 kaou      = 27,               &
      &                 kcLlim    = 28,               &
      &                 kcTlim    = 29,               &
      &                 kcPlim    = 30,               &
      &                 kcFlim    = 31,               &
                        nbgctend  = 31 
 
  INTEGER, PARAMETER :: kcflux     = 1,               &
      &                 koflux     = 2,               &
      &                 knflux     = 3,               &
      &                 knfixd     = 4,               &
      &                 kprcaca    = 5,               &
      &                 kprorca    = 6,               &
      &                 ksilpro    = 7,               &
      &                 kcoex90    = 8,               &
      &                 kcalex90   = 9,               &
      &                 kopex90    = 10,               &
      &                 kn2oflux   = 11,               &
      &                 korginp    = 12,              &
      &                 ksilinp    = 13,              &
      &                 kcalinp    = 14,              &
      &                 kprodus    = 15,               &
      &                 kdmsflux   = 16,               &
      &                 kcoex1000  = 17,               &
      &                 kcalex1000 = 18,               &
      &                 kopex1000  = 19,               &
      &                 kcoex2000  = 20,               &
      &                 kcalex2000 = 21,               &
      &                 kopex2000  = 22,               &
                        nbgcflux   = 22  
END MODULE mo_param1_bgc
