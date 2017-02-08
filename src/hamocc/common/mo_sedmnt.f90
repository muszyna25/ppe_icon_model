!>
!! @brief Variables for sediment modules
!!
!! declaration and memory allocation
#include "hamocc_omp_definitions.inc"
MODULE mo_sedmnt

  USE mo_kind, ONLY        : wp
  USE mo_param1_bgc, ONLY  : nsedtra, npowtra, n_bgctra
  USE mo_control_bgc, ONLY: dtbgc, bgc_nproma, bgc_zlevs 
  USE mo_hamocc_nml, ONLY : isac
  USE mo_memory_bgc, ONLY : kbo, bolay, bolaymin, wdust
  USE mo_bgc_constants, ONLY: g,rhoref_water

  IMPLICIT NONE

  PUBLIC
  INTEGER, PARAMETER::  &
  &        nsed_diag =3,&
  &        isremino =1,&
  &        isreminn =2,&
  &        isremins =3   

  !TODO: These are used IN currently 29 different source code file - the absolutely need more descriptive names
  INTEGER, PARAMETER :: ks=12, ksp=ks+1 ! sediment layering

  REAL(wp), ALLOCATABLE, TARGET :: sedlay (:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: powtra (:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: sedtend(:,:,:)

  REAL(wp), ALLOCATABLE,TARGET :: sedhpl(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: seddenit(:)

  REAL(wp), ALLOCATABLE :: seddw(:)
  REAL(wp), ALLOCATABLE :: porsol(:)
  REAL(wp), ALLOCATABLE :: porwah(:)
  REAL(wp), ALLOCATABLE :: porwat(:)
  REAL(wp), ALLOCATABLE :: pors2w(:)

  REAL(wp), ALLOCATABLE :: dzs(:)
  REAL(wp)              :: dzsed(100)
  REAL(wp), ALLOCATABLE :: seddzi(:)
  REAL(wp), ALLOCATABLE :: z_sed(:)

  REAL(wp), ALLOCATABLE :: silpro(:)
  REAL(wp), ALLOCATABLE :: prorca(:)
  REAL(wp), ALLOCATABLE :: prcaca(:)
  REAL(wp), ALLOCATABLE :: produs(:)

  REAL(wp), ALLOCATABLE, TARGET :: burial(:,:)

  ! pown2bud closes the mass balance for the alkalinity for biogenic induced changes in N2 in sediment
  REAL(wp), ALLOCATABLE :: pown2bud(:,:)
  ! powh2obud closes the mass balance for oxygen
  REAL(wp), ALLOCATABLE :: powh2obud(:,:)

  REAL(wp) :: sedict, calcon, rno3, o2ut, ansed, sedac, sedifti
  REAL(wp) :: calcwei, opalwei, orgwei
  REAL(wp) :: calcdens, opaldens, orgdens, claydens
  REAL(wp) :: calfa, oplfa, orgfa, clafa, solfu
  REAL(wp) :: sred_sed, disso_op,disso_cal
  REAL(wp) :: silsat, o2thresh
CONTAINS

SUBROUTINE sediment_bottom

 REAL(wp) :: sumsed, dustd1, dustd2
 INTEGER:: k
  !
  !---------------------------------------------------------------------
  !
  ! put into nml to allow for variable layer numbers (also ks, ksp)
 
  dzsed(:)=-1._wp

  dzs(1) = 0.001_wp
  dzs(2) = 0.003_wp
  dzs(3) = 0.005_wp
  dzs(4) = 0.007_wp
  dzs(5) = 0.009_wp
  dzs(6) = 0.011_wp
  dzs(7) = 0.013_wp
  dzs(8) = 0.015_wp
  dzs(9) = 0.017_wp
  dzs(10) = 0.019_wp
  dzs(11) = 0.021_wp
  dzs(12) = 0.023_wp
  dzs(13) = 0.025_wp

  dzsed(1:13) = dzs(:)

  porwat(1) = 0.85_wp
  porwat(2) = 0.83_wp
  porwat(3) = 0.8_wp
  porwat(4) = 0.79_wp
  porwat(5) = 0.77_wp
  porwat(6) = 0.75_wp
  porwat(7) = 0.73_wp
  porwat(8) = 0.7_wp
  porwat(9) = 0.68_wp
  porwat(10) = 0.66_wp
  porwat(11) = 0.64_wp
  porwat(12) = 0.62_wp


  seddzi(1) = 500._wp

  DO k = 1, ks
     porsol(k) = 1._wp - porwat(k)
     IF ( k == 1 ) porwah(k) = 0.5_wp * (1._wp + porwat(1))
     IF ( k >= 2 ) porwah(k) = 0.5_wp * (porwat(k) + porwat(k-1))
     ! inverse thickness of sediment layer
     seddzi(k+1) = 1._wp / dzs(k+1)
     ! define sediment layer thickness [m]
     seddw(k) = 0.5_wp * (dzs(k) + dzs(k+1))
  END DO

  ! ******************************************************************
  ! the following section is to include the SEDIMENT ACCELERATION
  ! mechanism to accelerate the sediment:

  sedac = 1._wp / REAL(isac, wp)

  ! determine total solid sediment thickness sumsed
  ! and reduced sediment volume

  sumsed = 0._wp

  DO k = 1, ks
     porwat(k) = porwat(k) * sedac
     porsol(k) = porsol(k) * sedac
     pors2w(k) = 1._wp/porwat(k) -1._wp ! ratio of porsol/porwat
     sumsed = sumsed + seddw(k)
  ENDDO

  ! depth of sediment layers for output

  z_sed(1) = dzs(1)

  DO k = 2, ks
     z_sed(k) = z_sed(k-1) + dzs(k)
  ENDDO

  ! determine reduced diffusion rate sedict,
  ! and scaling for bottom layer ventilation, sedifti

  sedict=1.e-9_wp*dtbgc                        ! units? m/s2 ?
  sedifti = sedict / (sumsed**2)               ! not used
  sedict=sedict*sedac

  ! ******************************************************************
  !
  ! densities etc. for SEDIMENT SHIFTING
  !
  ! define wei(ght) of calcium carbonate, opal, and poc [kg/kmol]
  !
  calcwei = 100._wp ! 40+12+3*16 kg/kmol C
  opalwei = 60._wp  ! 28 + 2*16  kg/kmol Si
  orgwei  = 30._wp  ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                    ! after Alldredge, 1998: POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3
  !
  ! define dens(ities) of caco3, opal, caco3 [kg/m3]
  !
  calcdens = 2600._wp
  opaldens = 2200._wp
  orgdens  = 1000._wp
  claydens = 2600._wp !quartz
  !
  ! define volumes occupied by solid constituents [m3/kmol]
  !
  calfa = calcwei/calcdens         ! 3.85e-2
  oplfa = opalwei/opaldens         ! 2.73e-2
  orgfa = orgwei/orgdens           ! 3.0e-2
  clafa = 1._wp/claydens           ! 3.85e-4 (clay is calculated in kg/m3)
  !
  ! determine total solid column integrated sediment volume  (1-d)
  !
  solfu = 0._wp
  !
  DO k=1,ks
     solfu = solfu + seddw(k)*porsol(k)
  ENDDO

  dustd1 = 0.0001_wp !cm = 1 um, boundary between clay and silt
  dustd2 = dustd1*dustd1
  wdust = (g * 86400._wp / 18._wp                          &  ! g * sec per day / 18. | js: Stoke's law for small particles
         &   * (claydens - rhoref_water) / 1.567_wp * 1000._wp  &  ! excess density / dyn. visc. | -> cm/s to m/day
         &   * dustd2 * 1.e-4_wp)                              ! *diameter**2 |*1000 *1.e-4?


END SUBROUTINE





SUBROUTINE  ini_bottom(start_idx,end_idx,klevs,pddpo)

 REAL(wp), INTENT(IN):: pddpo(bgc_nproma,bgc_zlevs)

 INTEGER, INTENT(IN) :: start_idx, end_idx
 INTEGER,  TARGET::klevs(bgc_nproma)
 INTEGER ::  j, k, kpke


  DO j = start_idx, end_idx
   k=klevs(j)
        kbo(j) = 1
        bolay(j) = 0._wp
        IF(k.gt.0)THEN
        IF ( pddpo(j,k) > 0.5_wp ) THEN
           bolay(j) = pddpo(j,k)
           kbo(j) = k
        ENDIF
        ENDIF
  END DO


 ! find depth of last wet layer
  ! and evaluate global min depth of last layer


  bolaymin=8000._wp
!HAMOCC_OMP_PARALLEL	
!HAMOCC_OMP_DO PRIVATE(j,kpke,k,bolaymin) HAMOCC_OMP_DEFAULT_SCHEDULE
  DO j = start_idx, end_idx
   kpke=klevs(j)
   IF(kpke>0)THEN ! always valid for MPIOM
   DO k = kpke-1, 1, -1
      IF ( pddpo(j,k) > 0.5_wp .AND. pddpo(j,k+1) <= 0.5_wp ) THEN
         bolay(j) = pddpo(j,k) ! local thickness of bottom layer
         kbo(j) = k
         bolaymin = MIN(bolaymin,bolay(j))
      ENDIF
    END DO
   ENDIF
  END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

 END SUBROUTINE
 SUBROUTINE ALLOC_MEM_SEDMNT


    ALLOCATE (sedlay(bgc_nproma,ks,nsedtra))
    sedlay(:,:,:) = 0._wp

    ALLOCATE (sedhpl(bgc_nproma,ks))
    ALLOCATE (burial(bgc_nproma,nsedtra))
    burial(:,:) = 0._wp
    ALLOCATE (powtra(bgc_nproma,ks,npowtra))
    ALLOCATE (silpro(bgc_nproma))
    silpro(:) = 0._wp
    ALLOCATE (prorca(bgc_nproma))
    prorca(:) = 0._wp

    ALLOCATE (prcaca(bgc_nproma))
    prcaca(:) = 0._wp
    ALLOCATE (produs(bgc_nproma))
    produs(:) = 0._wp
    ALLOCATE (dzs(ksp))
    ALLOCATE (seddzi(ksp))
    ALLOCATE (z_sed(ks))
    ALLOCATE (seddw(ks))
    ALLOCATE (porsol(ks))
    ALLOCATE (porwah(ks))
    ALLOCATE (porwat(ks))
    ALLOCATE (pors2w(ks))
    ALLOCATE (pown2bud(bgc_nproma,ks))
    pown2bud(:,:) = 0._wp
    ALLOCATE (powh2obud(bgc_nproma,ks))
    powh2obud(:,:) = 0._wp
    ALLOCATE (sedtend(bgc_nproma,ks,nsed_diag))
    sedtend(:,:,:)=0._wp
    ALLOCATE (seddenit(bgc_nproma))
    seddenit(:)=0._wp

  END SUBROUTINE ALLOC_MEM_SEDMNT

END MODULE mo_sedmnt
