!! @file mo_sedmnt_diffusion.f90
!! @brief sediment pore water diffusion
!! routines: dipowa, powadi

MODULE mo_sedmnt_diffusion

 USE mo_kind, ONLY           : wp
 USE mo_hamocc_nml, ONLY: ks,porwat
 
 IMPLICIT NONE

 PRIVATE

 PUBLIC:: dipowa, powadi


CONTAINS

SUBROUTINE DIPOWA (start_idx,end_idx)
!! @brief diffusion of pore water
!!
!! vertical diffusion of sediment pore water tracers
!! calculate vertical diffusion of sediment pore water properties
!! and diffusive flux through the ocean/sediment interface.
!! integration.
!!
!! implicit formulation;
!! constant diffusion coefficient : 1.e-9 set in mo_sedment.
!! diffusion coefficient : zcoefsu/zcoeflo for upper/lower
!! sediment layer boundary.



  USE mo_memory_bgc, ONLY         : bgctra, sedfluxo, satoxy, bgctend, &
  &                                 kbo, bolay

  USE mo_sedmnt, ONLY         : sedict, seddzi, seddw, &
       &                        porwah, powtra

  USE mo_control_bgc, ONLY    : dtbgc

  USE mo_param1_bgc, ONLY     : npowtra, ipowaox,  ioxygen,  &
  &                             ipowno3, ipowasi, iphosph, iano3, &
  &                             isilica, ipowafe, iiron, ialkali, &  
  &                             isco212, igasnit, ipowaph, ipowaal, &
  &                             ipown2, ipowaic, kaou 

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx    !< start index for j loop (ICON cells, MPIOM lat dir)          
  INTEGER, INTENT(in)  :: end_idx      !< end index  for j loop  (ICON cells, MPIOM lat dir) 
         

  !! Local variables

  INTEGER :: j,k,l,iv
  INTEGER :: iv_oc                         !< index of bgctra in powtra loop

  REAL(wp) :: sedb1(0:ks,npowtra)          !< 
  REAL(wp) :: zcoefsu(0:ks),zcoeflo(0:ks)  !< diffusion coefficients (upper/lower)
  REAL(wp) :: tredsy(0:ks,3)               !< redsy for 'reduced system'

  REAL(wp) :: aprior                       !< start value of oceanic tracer in bottom layer


  !
  ! --------------------------------------------------------------------
  !
  zcoefsu(0) = 0.0_wp

  DO  k = 1, ks
     ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
     zcoefsu(k  ) = -sedict * seddzi(k) * porwah(k)
     ! the lowerdiffusive flux of layer k is identical to
     ! the upper diff. flux of layer k+1
     zcoeflo(k-1) = -sedict * seddzi(k) * porwah(k)
  END DO

  ! diffusion coefficient for bottom sediment layer
  zcoeflo(ks) = 0.0_wp

  DO j = start_idx, end_idx

        k = 0

        tredsy(k,1) = zcoefsu(k)
        tredsy(k,3) = zcoeflo(k)
        ! dz(kbo) - diff upper - diff lower
        tredsy(k,2) =  bolay(j) - tredsy(k,1) - tredsy(k,3)


        DO iv = 1, npowtra      ! loop over pore water tracers

          iv_oc = iv

          if(iv == ipowaox) iv_oc = ioxygen
          if(iv == ipowno3) iv_oc = iano3
          if(iv == ipowasi) iv_oc = isilica
          if(iv == ipowafe) iv_oc = iiron
          if(iv == ipowaal) iv_oc = ialkali
          if(iv == ipowaph) iv_oc = iphosph
          if(iv == ipown2) iv_oc = igasnit
          if(iv == ipowaic) iv_oc = isco212

          sedb1( k, iv) = 0._wp
           ! tracer_concentration(kbo) * dz(kbo)
          IF (bolay( j) > 0._wp)                                     &
                & sedb1(k,iv) = bgctra(j,kbo(j),iv_oc) * bolay(j) 
   
           
        END DO

        DO k = 1, ks
           tredsy(k,1) = zcoefsu(k)
           tredsy(k,3) = zcoeflo(k)
           tredsy(k,2) = seddw(k) * porwat(k) - tredsy(k,1) - tredsy(k,3)
       END DO

       DO iv= 1, npowtra
        DO k = 1, ks
              ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
              sedb1(k,iv) = powtra(j,k,iv) * porwat(k) * seddw(k)
        END DO
       END DO

       DO k = 1, ks

           IF (bolay( j) > 0._wp) THEN
              ! this overwrites tredsy(k=0) for k=1
              tredsy(k-1,1) = tredsy(k,1) / tredsy(k-1,2)
              !                 diff upper    / conc (k-1)
              tredsy(k,2)   = tredsy(k,2)                      &
                   &          - tredsy(k-1,3) * tredsy(k,1) / tredsy(k-1,2)
              !   concentration -diff lower     * diff upper    / conc(k-1)
           ENDIF

       END DO

     ! diffusion from above
     DO iv = 1, npowtra
        DO k = 1, ks
              sedb1(k,iv) = sedb1(k,iv)                        &
                   &        - tredsy(k-1,1) * sedb1(k-1,iv)
        END DO
     END DO

     ! sediment bottom layer
     k = ks
     DO iv = 1, npowtra

           IF (bolay( j) > 0._wp) THEN
              powtra(j,k,iv) = sedb1(k,iv) / tredsy(k,2)
           ENDIF

     END DO


     ! sediment column
     DO iv = 1, npowtra
        DO k = 1, ks-1
           l = ks-k
              IF (bolay( j) > 0._wp) THEN
                 powtra(j,l,iv) = ( sedb1(l,iv)            &
                      &             - tredsy(l,3) * powtra(j,l+1,iv) )    &
                      &             / tredsy(l,2)
              ENDIF
        END DO
     END DO

     ! sediment ocean interface
     DO iv = 1, npowtra
        ! 
        ! check mo_param1_bgc.f90 for consistency
        iv_oc = iv
        if(iv == ipowaox) iv_oc = ioxygen
        if(iv == ipowno3) iv_oc = iano3
        if(iv == ipowasi) iv_oc = isilica
        if(iv == ipowafe) iv_oc = iiron
        if(iv == ipowaal) iv_oc = ialkali
        if(iv == ipowaph) iv_oc = iphosph
        if(iv == ipown2) iv_oc = igasnit
        if(iv == ipowaic) iv_oc = isco212

           l = 0
           IF (bolay( j) > 0._wp) THEN

              aprior = bgctra(j,kbo(j),iv_oc)
              bgctra(j,kbo(j),iv_oc) =                                  &
                   &         ( sedb1(l,iv) - tredsy(l,3) * powtra(j,l+1,iv) ) &
                   &         / tredsy(l,2)

              sedfluxo(j,iv) = (bgctra(j,kbo(j),iv_oc)-aprior)*bolay(j)/dtbgc
           ENDIF       !----------------------------------------------------------------------



     END DO
     ! update AOU
     IF (bolay( j) > 0._wp)bgctend(j,kbo(j),kaou) = satoxy(j,kbo(j)) - bgctra(j,kbo(j),ioxygen)
     


  END DO ! j loop

END SUBROUTINE DIPOWA

SUBROUTINE powadi ( j,  solrat, sedb1, sediso, bolven)
!! @file powadi.f90
!! @brief vertical diffusion with simultaneous dissolution,
!! implicit discretisation.
!!
!!


  USE mo_sedmnt, ONLY     : sedict, seddzi, seddw, porwah

  USE mo_memory_bgc, ONLY     : bolay


  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)     :: j             !< zonal grid index

  REAL(wp), INTENT(in)    :: solrat(ks)    !< dissolution rate
  REAL(wp), INTENT(inout) :: sedb1(0:ks)   !< tracer at entry
  REAL(wp), INTENT(inout) :: sediso(0:ks)  !<
  REAL(wp), INTENT(in)    :: bolven        !< bottom layer ventilation rate

  !! Local variables

  INTEGER,  SAVE  :: k,l   
  REAL(wp) :: tredsy(0:ks,3)
  REAL(wp), SAVE :: asu,alo
  !
  !----------------------------------------------------------------------
  !
  DO k = 1, ks

     asu=sedict*seddzi(k)*porwah(k)
     alo = 0._wp

     IF (k < ks) alo = sedict*seddzi(k+1)*porwah(k+1)

        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = seddw(k)*porwat(k) - tredsy(k,1)       &
             &        - tredsy(k,3) + solrat(k)*porwat(k)*seddw(k)

  END DO

  k = 0

  asu = 0._wp
  alo = sedict*seddzi(1)*porwah(1)

     IF(bolay(j) > 0._wp) THEN
        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = bolven*bolay(j) - tredsy(k,1) - tredsy(k,3)
     ELSE
        tredsy(k,1) = 0._wp
        tredsy(k,3) = 0._wp
        tredsy(k,2) = 0._wp
     ENDIF


  DO k = 1, ks
        IF (bolay(j) > 0._wp) THEN
           tredsy(k-1,1) = tredsy(k,1) / tredsy(k-1,2)
           tredsy(k,2)   = tredsy(k,2)                       &
                &          - tredsy(k-1,3) * tredsy(k,1) / tredsy(k-1,2)
        ENDIF
  END DO

  DO k = 1, ks
        sedb1(k) = sedb1(k) - tredsy(k-1,1) * sedb1(k-1)
  END DO

  k = ks

     IF (bolay(j) > 0._wp)sediso(k) = sedb1(k) / tredsy(k,2)

  DO k = 1, ks
     l = ks-k
        IF (bolay(j) > 0._wp) sediso(l) =                           &
             &           ( sedb1(l) - tredsy(l,3) * sediso(l+1) )       &
             &           / tredsy(l,2)
  END DO

END SUBROUTINE powadi

END MODULE mo_sedmnt_diffusion
