#include "hamocc_omp_definitions.inc"

  SUBROUTINE swr_absorption(start_idx,end_idx, klevs, pfswr, psicomo, dzw)

    USE mo_memory_bgc, ONLY : bgctra, swr_frac, meanswr, strahl
    USE mo_param1_bgc, ONLY : iphy, icya
    USE mo_hamocc_nml, ONLY : l_cyadyn
    USE mo_kind,    ONLY: wp                                               
    USE mo_control_bgc, ONLY: bgc_zlevs, bgc_nproma


    INTEGER, INTENT(in):: start_idx
    INTEGER, INTENT(IN):: end_idx
    INTEGER :: klevs(bgc_nproma)
    REAL(wp), INTENT(in):: pfswr(bgc_nproma)
    REAL(wp), INTENT(in):: psicomo(bgc_nproma)
    REAL(wp), INTENT(in):: dzw(bgc_nproma,bgc_zlevs)

    !! Analogue to Zielinski et al., Deep-Sea Research II 49 (2002), 3529-3542

    REAL(wp), PARAMETER :: redfrac=0.4_wp !< red fraction of the spectral domain (> 580nm)

    REAL(wp), PARAMETER :: c_to_chl=12.0_wp/60.0_wp   !< ration Carbon to Chlorophyll
    REAL(wp), PARAMETER :: r_car=122.0_wp   !< Redfield ratio
    REAL(wp), PARAMETER :: pho_to_chl=r_car*c_to_chl*1.e6_wp !< 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll]

    REAL(wp), PARAMETER :: atten_r=0.35_wp !< attenuation of red light [m-1]
    REAL(wp), PARAMETER :: atten_w=0.03_wp !< attenuation of blue/green light
                                           !! in clear water between 400nm and 580nm [m-1]
    REAL(wp), PARAMETER :: atten_c=0.04_wp !< attenuation of blue/green light
                                           !! by chlorophyll [m-1]

    REAL(wp) :: swr_r
    REAL(wp) :: swr_b
    REAL(wp) :: rcyano

    INTEGER :: k, kpke, j

    ! if prognostic cyanobacteria are calculated
    ! use them in absorption (rcyano=1)
    rcyano=merge(1._wp,0._wp,l_cyadyn)

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,swr_r,swr_b) HAMOCC_OMP_DEFAULT_SCHEDULE
    DO j = start_idx, end_idx

      strahl(j) = pfswr(j) * (1._wp - psicomo(j))

      swr_frac(j,1) = 1.0_wp
      
      
      kpke = klevs(j)
    
      IF(kpke > 0) then

      swr_r = redfrac
      swr_b = (1._wp-redfrac)

      DO k=2,kpke
 
           swr_r = swr_r * EXP(-dzw(j,k-1) *  atten_r)
           swr_b = swr_b * EXP(-dzw(j,k-1) * (atten_w +&
        &    atten_c*pho_to_chl*MAX(0.0_wp,(bgctra(j,k-1,iphy)+rcyano*bgctra(j,k-1,icya)))))
           swr_frac(j,k) = swr_r + swr_b

      END DO
      DO k=1,kpke-1
           meanswr(j,k) = (swr_frac(j,k) + swr_frac(j,k+1))/2._wp
      END DO
      meanswr(j,kpke) = swr_frac(j,k) 

      ENDIF
   ENDDO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
  END SUBROUTINE swr_absorption



