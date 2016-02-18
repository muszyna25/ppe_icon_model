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
  INTEGER, PARAMETER :: i_base_adv = 16,              &
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
       &                idoccya    = 14,              &               
       &                iiron      = 15,              &               
       &                idms       = 16               


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
                        nbgctend  = 18 
 
  INTEGER, PARAMETER :: kcflux     = 1,               &
      &                 koflux     = 2,               &
      &                 knflux     = 3,               &
      &                 knfixd     = 4,               &
      &                 kprcaca    = 5,               &
      &                 kprorca    = 6,               &
      &                 ksilpro    = 7,               &
      &                 kcoex90    = 8,               &
      &                 kn2oflux   = 9,               &
      &                 korginp   = 10,              &
      &                 ksilinp   = 11,              &
      &                 kcalinp   = 12,              &
      &                 kprodus   = 13,               &
      &                 kdmsflux  = 14,               &
                        nbgcflux   = 14  
END MODULE mo_param1_bgc
