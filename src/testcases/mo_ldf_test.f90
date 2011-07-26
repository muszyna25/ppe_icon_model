!>
!! This module contains a subroutine used in the local
!! diabatic forcing test
!!
!! mo_local_forcing_test is called from mo_hierarchy_management
!!
!! @author Constantin Junk, MPI-M
!!
!! @par Revision History
!! Original implementation by Constantin Junk, MPI-M (2011-01-06)
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
MODULE mo_ldf_test

  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: pi
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,  ONLY: nproma
  USE mo_run_config,         ONLY: iforcing
  USE mo_impl_constants,     ONLY: ildf_echam
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_ha_testcases,       ONLY: ldf_symm, ildf_init_type


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldf_temp

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !>
  !! Local diabatic forcing subroutine -- if the variable
  !! iforcing=4 in the namelist, the dry dynamical
  !! core is externally forced by a local diabatic forcing.
  !! For iforcing=4, one has the choice to initialise the ICOHDC
  !! with a resting, isothermal atmosphere (ildf_init_type=0) or
  !! or with an zonal wind field similar to JWs (ildf_init_type=1).
  !!
  !! If iforcing=5, the dry dynamical core is forced by ECHAM
  !! physics and the local forcing! The model is initialized with
  !! with a zonal wind field and a user specified moisture.
  !!
  !! The forcing is centered either around the equator 
  !! (ldf_symm=.TRUE.) or around 30N (ldf_symm=.FALSE.)
  !! and is strongest in the middle troposphere (~500hPa level).
  !! for both cases. 
  !!
  !! The horizontal part of the forcing function is mainly taken
  !! taken from Webster (1972, eq.34). The strength of the diabatic
  !! forcing is specified by the variable 'zTempA'
  !!
  !! @par See also
  !! Gill,A.E.,1980:    Some simple solutions for heat-induced
  !!                    tropical circulations. Quart.J.R.Met.Soc
  !! Webster,P.J.,1972: Response of the Tropical Atmosphere
  !!                    to Local, Steady Forcing.
  !!
  !!
  !!
  !!
  !! @par Revision History
  !! Original implementation by Constantin Junk, MPI-M (2010-01-06)
  !! Modification by Constantin Junk, MPI-M (2010-04-04)
  !! - to use the same implementation as for the physics part,
  !!   the forcing is added to the dynamics tendency.
  !!

  SUBROUTINE ldf_temp(  p_patch,      &!in
                      & p_prog,       &!in
                      & dyn_diag,     &!in
                      & acc_tend )     !inout

    ! Arguments

    TYPE(t_patch),TARGET,INTENT(IN)           :: p_patch
    TYPE(t_hydro_atm_prog),INTENT(IN)         :: p_prog
    TYPE(t_hydro_atm_diag),INTENT(IN)         :: dyn_diag
    TYPE(t_hydro_atm_prog),INTENT(INOUT)      :: acc_tend

    ! Local scalar

    INTEGER :: jk, is, ie, jc
    INTEGER :: jb, jbs, nblks_c

    ! Local arrays

    REAL(wp) :: zsigma (nproma)              ! sigma = pres/pres_sfc
    REAL(wp) :: zlat   (nproma)              ! latitude
    REAL(wp) :: zlon   (nproma)              ! longitude

    REAL(wp) :: zsinsig2(nproma), zsinsig4(nproma)
    REAL(wp) :: zTempHrz(nproma),zTempLon(nproma),zTempLat(nproma)


    REAL(wp) :: zTempA    !Amplitude of the local diabatic forcing (in Kelvin/s)

    nblks_c = p_patch%nblks_int_c                   
    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,zlat,zlon,zsigma,zsinsig2,zsinsig4,zTempLon,zTempLat,zTempHrz)
   DO jb = jbs,nblks_c

       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       DO jk = 1,p_patch%nlev

         DO jc = is,ie

           zlat(jc) = p_patch%cells%center(jc,jb)%lat
           zlon(jc) = p_patch%cells%center(jc,jb)%lon

           !Vertical part of the forcing function
           !
           zsigma(jc) = dyn_diag%pres_mc(jc,jk,jb) &
                           /p_prog%pres_sfc(jc,jb)
           zsinsig2(jc) = SIN(zsigma(jc)*pi)**2
           zsinsig4(jc) = zsinsig2(jc)**2

           ! horizontal part of the forcing function
    
           IF (ldf_symm) THEN  
              !forcing symmetric about the equator
              zTempLat(jc) = -25._wp * (zlat(jc))**2

           ELSE                             
              !antisymmetric about the equator, placed at 30 N
              zTempLat(jc) = -25._wp * (zlat(jc) - pi/9._wp)**2 

           END IF

           zTempLon(jc) = (-81._wp * (zlon(jc))**2) / (4._wp*(pi**2))
           zTempHrz(jc) = EXP( zTempLon(jc) + zTempLat(jc) )
 
           !Amplitude of the forcing depending on testcase setup

           IF (iforcing == ildf_echam) THEN
                 zTempA  = 5.E-5_wp
           ELSE
              IF  (ildf_init_type == 0) THEN
                 zTempA  = 1.E-6_wp
              ELSE  !(ildf_init_type == 1)
                 zTempA  = 5.E-5_wp  
              END IF
           END IF

           ! accumulate the tendencies by adding the local diabatic forcing

           acc_tend%temp(jc,jk,jb) = acc_tend%temp(jc,jk,jb)+ zTempA * zsinsig4(jc) * zTempHrz(jc)

         ENDDO !jc-loop

       ENDDO !vertical layer loop

   ENDDO !jb-loop
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE ldf_temp

END MODULE mo_ldf_test



