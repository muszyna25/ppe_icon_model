!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-09-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aerosol_util

  USE mo_impl_constants,       ONLY: min_rlcell, iss, iorg, ibc, iso4, idu, dzsoil, nclass_aero
  USE mo_math_constants,       ONLY: rad2deg
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_lrtm_par,             ONLY: jpband => nbndlw
  USE mo_model_domain,         ONLY: t_patch
  USE mo_radiation_rg_par,     ONLY: jpspec
  USE mo_srtm_config,          ONLY: jpsw
  USE sfc_terra_data,          ONLY: cadp, cfcap
  USE mo_lnd_nwp_config,       ONLY: ntiles_lnd
  USE mo_nwp_tuning_config,    ONLY: tune_dust_abs

  IMPLICIT NONE

  PRIVATE


 
  !RRTM
  REAL  (wp)              ::           &
  zaea_rrtm(jpsw+jpband,5), &  ! ratio of optical thickness for the absorption in spectral
                               ! interval jpspec  and total optical thickness at 0.55m*1.E-06 
                               ! for an aerosoltyp specified by second array index
  zaes_rrtm(jpsw+jpband,5), &  ! analog for the optical thickness of scattering 
  zaeg_rrtm(jpsw+jpband,5)!, zaef_rrtm(jpsw+jpband,5)

  !RG
  REAL  (wp)              ::           &
  zaea_rg(jpspec,5), &  ! ratio of optical thickness for the absorption in spectral
                     ! interval jpspec  and total optical thickness at 0.55m*1.E-06 
                     ! for an aerosoltyp specified by second array index
  zaes_rg(jpspec,5), &  ! analog for the optical thickness of scattering 
  zaeg_rg(jpspec,5), zaef_rg(jpspec,5)


  PUBLIC :: zaea_rrtm, zaes_rrtm, zaeg_rrtm, zaea_rg, zaes_rg, zaeg_rg, zaef_rg, &
    &       init_aerosol_dstrb_tanre,init_aerosol_props_tanre_rrtm, &
    &       init_aerosol_props_tanre_rg, &
    &       init_aerosol_props_tegen_rrtm, init_aerosol_props_tegen_rg, prog_aerosol_2D, tune_dust
  
CONTAINS

  !>
  !! Initializes aerosol (climatologically).
  !!
  !! Taken from COSMO model (version 4.16), adaptated for better vectorization.
  !!
  !! @par Revision History
  !! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-02-28)
  !!      
  SUBROUTINE init_aerosol_dstrb_tanre (kbdim,pt_patch,aersea,aerlan,aerurb,aerdes)
    
    INTEGER,            INTENT(in)    :: kbdim

    TYPE(t_patch),      INTENT(in)    :: pt_patch    ! Patch

    REAL(wp),           INTENT(out)   :: &
      & aersea(kbdim,pt_patch%nblks_c),  &
      & aerlan(kbdim,pt_patch%nblks_c),  &
      & aerurb(kbdim,pt_patch%nblks_c),  &
      & aerdes(kbdim,pt_patch%nblks_c)

    ! Local arrays and scalars:
    ! -------------------------
    INTEGER                   ::  &
      jc, jb,                     & ! loop indices
      jzj, jzm1, jzm2, jzm, jzn,  & ! indices for Legendre coefficients
      imn, imnc, imns, jmm, jnn

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL (wp)              ::  &
                                ! arrays for the T10 distribution of
      zaesc(66) , zaess (55) , & ! sea    type aerosols                     
      zaelc(66) , zaels (55) , & ! land   type aerosols
      zaeuc(66) , zaeus (55) , & ! urban  type aerosols
      zaedc(66) , zaeds (55) , & ! desert type aerosols
      zfaes(kbdim,21) , zfael (kbdim,21) , & ! coefficients for spectral
      zfaeu(kbdim,21) , zfaed (kbdim,21) , & ! expansion
      zalp (kbdim,66) ,              & !
      zsinphi(kbdim)   , zcosphi(kbdim)    , & !
      zm, z2m, zre1, ze1, ze2, & !
      zf1m(kbdim), zf2m(kbdim), zn, zn2,     & !
      zsin1, zsin2, zsin3, zsin4, zsin5,  & ! 
      zsin6, zsin7, zsin8, zsin9, zsin10, & !
      zcos1, zcos2, zcos3, zcos4, zcos5,  & ! 
      zcos6, zcos7, zcos8, zcos9, zcos10    !

    !------------------------------------------------------------------------------
    ! Section 0: Data for the Fourier coefficients of the four aerosol types          
    !------------------------------------------------------------------------------

    DATA zaesc/ &
      +.6688E+00_wp,-.1172E+00_wp,-.1013E+00_wp,+.1636E-01_wp,-.3699E-01_wp,+.1775E-01_wp,   &
      -.9635E-02_wp,+.1290E-02_wp,+.4681E-04_wp,-.9106E-04_wp,+.9355E-04_wp,                 &
      -.7076E-01_wp,-.1782E-01_wp,+.1856E-01_wp,+.1372E-01_wp,+.8210E-04_wp,+.2149E-02_wp,   &
      +.4856E-03_wp,+.2231E-03_wp,+.1824E-03_wp,+.1960E-05_wp,                               &
      +.2057E-01_wp,+.2703E-01_wp,+.2424E-01_wp,+.9716E-02_wp,+.1312E-02_wp,-.8846E-03_wp,   &
      -.3347E-03_wp,+.6231E-04_wp,+.6397E-04_wp,                                             &
      -.3341E-02_wp,-.1295E-01_wp,-.4598E-02_wp,+.3242E-03_wp,+.8122E-03_wp,-.2975E-03_wp,   &
      -.7757E-04_wp,+.7793E-04_wp,                                                           &
      +.4455E-02_wp,-.1584E-01_wp,-.2551E-02_wp,+.1174E-02_wp,+.1335E-04_wp,+.5112E-04_wp,   &
      +.5605E-04_wp,                                                                         &
      +.7412E-04_wp,+.1857E-02_wp,-.1917E-03_wp,+.4460E-03_wp,+.1767E-04_wp,-.5281E-04_wp,   &
      -.5043E-03_wp,+.2467E-03_wp,-.2497E-03_wp,-.2377E-04_wp,-.3954E-04_wp,                 &
      +.2666E-03_wp,-.8186E-03_wp,-.1441E-03_wp,-.1904E-04_wp,                               &
      +.3337E-03_wp,-.1696E-03_wp,-.2503E-04_wp,                                             &
      +.1239E-03_wp,-.9983E-04_wp,                                                           &
      -.5283E-04_wp  /

    DATA zaess/                                                                              &
      -.3374E-01_wp,-.3247E-01_wp,-.1012E-01_wp,+.6002E-02_wp,+.5190E-02_wp,+.7784E-03_wp,   &
      -.1090E-02_wp,+.3294E-03_wp,+.1719E-03_wp,-.5866E-05_wp,                               &
      -.4124E-03_wp,-.3742E-01_wp,-.5054E-02_wp,+.3430E-02_wp,+.5513E-03_wp,-.6235E-03_wp,   &
      +.2892E-03_wp,-.9730E-04_wp,+.7078E-04_wp,                                             &
      -.3300E-01_wp,+.5104E-03_wp,-.2156E-02_wp,-.3194E-02_wp,-.5079E-03_wp,-.5517E-03_wp,   &
      +.4632E-04_wp,+.5369E-04_wp,                                                           &
      -.2731E-01_wp,+.5126E-02_wp,+.2241E-02_wp,-.5789E-03_wp,-.3048E-03_wp,-.1774E-03_wp,   &
      +.1946E-05_wp,                                                                         &
      -.8247E-02_wp,+.2338E-02_wp,+.1021E-02_wp,+.1575E-04_wp,+.2612E-05_wp,+.1995E-04_wp,   &
      -.1319E-02_wp,+.1384E-02_wp,-.4159E-03_wp,-.2337E-03_wp,+.5764E-04_wp,                 &
      +.1495E-02_wp,-.3727E-03_wp,+.6075E-04_wp,-.4642E-04_wp,                               &
      +.5368E-03_wp,-.7619E-04_wp,+.3774E-04_wp,                                             &
      +.1206E-03_wp,-.4104E-06_wp,                                                           &
      +.2158E-04_wp  /

    DATA zaelc/                                                                              &
      +.1542E+00_wp,+.8245E-01_wp,-.1879E-03_wp,+.4864E-02_wp,-.5527E-02_wp,-.7966E-02_wp,   &
      -.2683E-02_wp,-.2011E-02_wp,-.8889E-03_wp,-.1058E-03_wp,-.1614E-04_wp,                 &
      +.4206E-01_wp,+.1912E-01_wp,-.9476E-02_wp,-.6780E-02_wp,+.1767E-03_wp,-.5422E-03_wp,   &
      -.7753E-03_wp,-.2106E-03_wp,-.9870E-04_wp,-.1721E-04_wp,                               &
      -.9536E-02_wp,-.9580E-02_wp,-.1050E-01_wp,-.5747E-02_wp,-.1282E-02_wp,+.2248E-03_wp,   &
      +.1694E-03_wp,-.4782E-04_wp,-.2441E-04_wp,                                             &
      +.5781E-03_wp,+.6212E-02_wp,+.1921E-02_wp,-.1102E-02_wp,-.8145E-03_wp,+.2497E-03_wp,   &
      +.1539E-03_wp,-.2538E-04_wp,                                                           &
      -.3993E-02_wp,+.9777E-02_wp,+.4837E-03_wp,-.1304E-02_wp,+.2417E-04_wp,-.1370E-04_wp,   &
      -.3731E-05_wp,                                                                         &
      +.1922E-02_wp,-.5167E-03_wp,+.4295E-03_wp,-.1888E-03_wp,+.2427E-04_wp,+.4012E-04_wp,   &
      +.1529E-02_wp,-.2120E-03_wp,+.8166E-04_wp,+.2579E-04_wp,+.3488E-04_wp,                 &
      +.2140E-03_wp,+.2274E-03_wp,-.3447E-05_wp,-.1075E-04_wp,                               &
      -.1018E-03_wp,+.2864E-04_wp,+.3442E-04_wp,                                             &
      -.1002E-03_wp,+.7117E-04_wp,                                                           &
      +.2045E-04_wp  /

    DATA zaels/                                                                              &
      +.1637E-01_wp,+.1935E-01_wp,+.1080E-01_wp,+.2784E-02_wp,+.1606E-03_wp,+.1860E-02_wp,   &
      +.1263E-02_wp,-.2707E-03_wp,-.2290E-03_wp,-.9761E-05_wp,                               &
      -.7317E-02_wp,+.2465E-01_wp,+.6799E-02_wp,-.1913E-02_wp,+.1382E-02_wp,+.6691E-03_wp,   &
      +.1414E-03_wp,+.3527E-04_wp,-.5210E-04_wp,                                             &
      +.1873E-01_wp,+.2977E-02_wp,+.4650E-02_wp,+.2509E-02_wp,+.3680E-03_wp,+.1481E-03_wp,   &
      -.6594E-04_wp,-.5634E-04_wp,                                                           &
      +.1592E-01_wp,-.1875E-02_wp,-.1093E-02_wp,+.3022E-03_wp,+.2625E-03_wp,+.3252E-04_wp,   &
      -.3803E-04_wp,                                                                         &
      +.4218E-02_wp,-.1843E-02_wp,-.1351E-02_wp,-.2952E-03_wp,-.8171E-05_wp,-.1473E-04_wp,   &
      +.9076E-03_wp,-.1057E-02_wp,+.2676E-03_wp,+.1307E-03_wp,-.3628E-04_wp,                 &
      -.9158E-03_wp,+.4335E-03_wp,+.2927E-04_wp,+.6602E-04_wp,                               &
      -.3570E-03_wp,+.5760E-04_wp,-.3465E-04_wp,                                             &
      -.8535E-04_wp,-.2011E-04_wp,                                                           &
      +.6612E-06_wp  /  

    DATA zaeuc/                                                                              &
      +.8005E-01_wp,+.7095E-01_wp,+.2014E-01_wp,-.1412E-01_wp,-.2425E-01_wp,-.1332E-01_wp,   &
      -.2904E-02_wp,+.5068E-03_wp,+.9369E-03_wp,+.4114E-03_wp,+.7549E-04_wp,                 &
      +.1922E-01_wp,+.2534E-01_wp,+.2088E-01_wp,+.1064E-01_wp,+.1063E-02_wp,-.2526E-02_wp,   &
      -.2091E-02_wp,-.9660E-03_wp,-.2030E-03_wp,+.3865E-04_wp,                               &
      -.9900E-02_wp,-.5964E-02_wp,+.2223E-02_wp,+.4941E-02_wp,+.3277E-02_wp,+.1038E-02_wp,   &
      -.1480E-03_wp,-.2844E-03_wp,-.1208E-03_wp,                                             &
      +.3999E-02_wp,+.6282E-02_wp,+.2813E-02_wp,+.1475E-02_wp,+.4571E-03_wp,-.1349E-03_wp,   &
      -.9011E-04_wp,-.1936E-04_wp,                                                           &
      +.1994E-02_wp,+.3540E-02_wp,+.8837E-03_wp,+.1992E-03_wp,+.3092E-04_wp,-.7979E-04_wp,   &
      -.2664E-04_wp,                                                                         &
      -.5006E-04_wp,+.6447E-03_wp,+.5550E-03_wp,+.1197E-03_wp,+.6657E-04_wp,+.1488E-04_wp,   &
      -.9141E-04_wp,-.2896E-03_wp,-.1561E-03_wp,-.6524E-04_wp,-.1559E-04_wp,                 &
      -.1082E-03_wp,-.4126E-03_wp,-.1732E-03_wp,-.8286E-04_wp,                               &
      -.1993E-04_wp,+.3850E-04_wp,+.2870E-04_wp,                                             &
      +.4493E-04_wp,+.4721E-04_wp,                                                           &
      +.1338E-04_wp  /

    DATA zaeus/                                                                              &
      +.6646E-02_wp,+.8373E-02_wp,+.5463E-02_wp,+.4554E-02_wp,+.3301E-02_wp,+.5725E-03_wp,   &
      -.7482E-03_wp,-.6222E-03_wp,-.2603E-03_wp,-.5127E-04_wp,                               &
      -.3849E-04_wp,+.9741E-02_wp,+.8190E-02_wp,+.5712E-02_wp,+.3039E-02_wp,+.5290E-03_wp,   &
      -.2044E-03_wp,-.2309E-03_wp,-.1160E-03_wp,                                             &
      +.9160E-02_wp,+.1286E-01_wp,+.1170E-01_wp,+.5491E-02_wp,+.1393E-02_wp,-.6288E-04_wp,   &
      -.2715E-03_wp,-.1047E-03_wp,                                                           &
      +.4873E-02_wp,+.3545E-02_wp,+.3069E-02_wp,+.1819E-02_wp,+.6947E-03_wp,+.1416E-03_wp,   &
      -.1538E-04_wp,                                                                         &
      -.4351E-03_wp,-.1907E-02_wp,-.5774E-03_wp,-.2247E-03_wp,+.5345E-04_wp,+.9052E-04_wp,   &
      -.3972E-04_wp,-.9665E-04_wp,+.7912E-04_wp,-.1094E-04_wp,-.6776E-05_wp,                 &
      +.2724E-03_wp,+.1973E-03_wp,+.6837E-04_wp,+.4313E-04_wp,                               &
      -.7174E-05_wp,+.8527E-05_wp,-.2160E-05_wp,                                             &
      -.7852E-04_wp,+.3453E-06_wp,                                                           &
      -.2402E-05_wp  /

    DATA zaedc/                                                                              &
      +.2840E-01_wp,+.1775E-01_wp,-.1069E-01_wp,-.1553E-01_wp,-.3299E-02_wp,+.3583E-02_wp,   &
      +.2274E-02_wp,+.5767E-04_wp,-.3678E-03_wp,-.1050E-03_wp,+.2133E-04_wp,                 &
      +.2326E-01_wp,+.1566E-01_wp,-.3130E-02_wp,-.8253E-02_wp,-.2615E-02_wp,+.1247E-02_wp,   &
      +.1059E-02_wp,+.1196E-03_wp,-.1303E-03_wp,-.5094E-04_wp,                               &
      +.1185E-01_wp,+.7238E-02_wp,-.1562E-02_wp,-.3665E-02_wp,-.1182E-02_wp,+.4678E-03_wp,   &
      +.4448E-03_wp,+.8307E-04_wp,-.3468E-04_wp,                                             &
      +.5273E-02_wp,+.3037E-02_wp,-.4014E-03_wp,-.1202E-02_wp,-.4647E-03_wp,+.5148E-04_wp,   &
      +.1014E-03_wp,+.2996E-04_wp,                                                           &
      +.2505E-02_wp,+.1495E-02_wp,+.2438E-03_wp,-.1223E-03_wp,-.7669E-04_wp,-.1638E-04_wp,   &
      +.1869E-05_wp,                                                                         &
      +.1094E-02_wp,+.6131E-03_wp,+.1508E-03_wp,+.1765E-04_wp,+.1360E-05_wp,-.7998E-06_wp,   &
      +.4475E-03_wp,+.2737E-03_wp,+.6430E-04_wp,-.6759E-05_wp,-.6761E-05_wp,                 &
      +.1992E-03_wp,+.1531E-03_wp,+.4828E-04_wp,+.5103E-06_wp,                               &
      +.7454E-04_wp,+.5917E-04_wp,+.2152E-04_wp,                                             &
      +.9300E-05_wp,+.9790E-05_wp,                                                           &
      -.8853E-05_wp /

    DATA zaeds/                                                                              &
      +.9815E-02_wp,+.8436E-02_wp,+.1087E-02_wp,-.2717E-02_wp,-.1755E-02_wp,-.1559E-03_wp,   &
      +.2367E-03_wp,+.8808E-04_wp,+.2001E-05_wp,-.1244E-05_wp,                               &
      +.1041E-01_wp,+.8039E-02_wp,+.1005E-02_wp,-.1981E-02_wp,-.1090E-02_wp,+.1595E-05_wp,   &
      +.1787E-03_wp,+.4644E-04_wp,-.1052E-04_wp,                                             &
      +.6593E-02_wp,+.3983E-02_wp,-.1527E-03_wp,-.1235E-02_wp,-.5078E-03_wp,+.3649E-04_wp,   &
      +.1005E-03_wp,+.3182E-04_wp,                                                           &
      +.3225E-02_wp,+.1672E-02_wp,-.7752E-04_wp,-.4312E-03_wp,-.1872E-03_wp,-.1666E-04_wp,   &
      +.1872E-04_wp,                                                                         &
      +.1133E-02_wp,+.5643E-03_wp,+.7747E-04_wp,-.2980E-04_wp,-.2092E-04_wp,-.8590E-05_wp,   &
      +.2988E-03_wp,+.6714E-04_wp,-.6249E-05_wp,+.1052E-04_wp,+.8790E-05_wp,                 &
      +.1569E-03_wp,-.1175E-04_wp,-.3033E-04_wp,-.9777E-06_wp,                               &
      +.1101E-03_wp,+.6827E-05_wp,-.1023E-04_wp,                                             &
      +.4231E-04_wp,+.4905E-05_wp,                                                           &
      +.6229E-05_wp  /


    !------------------------------------------------------------------------------
    ! Section 2: Calculation of the inverse Legendre and Fourier transformation  
    !------------------------------------------------------------------------------
 
    i_nchdom  = MAX(1,pt_patch%n_childdom)

    ! nest boudaries have to be included for reduced-grid option
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    zalp (:,1) = 1.0_wp

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      i_startidx=1 ! because start index is missing in RRTM
      
      DO jc=i_startidx,i_endidx

        ! Calculation of the values zalp for the sine of latitude (zsinphi) of the
        ! normalized Legendre associated functions. The limit wave number is 10.

        zsinphi(jc)  = SIN(pt_patch%cells%center(jc,jb)%lat )
        zcosphi(jc)  = SQRT(1._wp-zsinphi(jc)**2)

        zf1m(jc)     = SQRT(3.0_wp)
        !        zalp (1) = 1.0_wp
        zalp (jc,2) = zf1m(jc)*zsinphi(jc)

      ENDDO

      jzj      = 2

      wave_number_loop : DO jzm1 = 1, 11

        jzm  = jzm1-1
        zm   = REAL(jzm,wp)
        z2m  = zm + zm
        zre1 = SQRT(z2m+3.0_wp)
        ze1  = 1.0_wp/zre1
        IF (jzm.NE.0) THEN
          jzj       = jzj + 1
          DO jc=i_startidx,i_endidx
            zf2m(jc)      = zf1m(jc)*zcosphi(jc)/SQRT(z2m)
            zf1m(jc)      = zf2m(jc)*zre1
            zalp(jc,jzj) = zf2m(jc)
          ENDDO
          IF(jzm ==10) CYCLE wave_number_loop
          jzj       = jzj + 1
          DO jc=i_startidx,i_endidx
            zalp(jc,jzj) = zf1m(jc)*zsinphi(jc)
          ENDDO
          IF(jzm1==10) CYCLE wave_number_loop
        ENDIF
        jzm2 = jzm+2
        DO jzn = jzm2, 10
          zn        = REAL(jzn,wp)
          zn2       = zn**2
          ze2       = SQRT( (4.0_wp*zn2-1.0_wp)/(zn2-zm**2) )
          jzj       = jzj+1
          DO jc=i_startidx,i_endidx
            zalp(jc,jzj) = ze2*(zsinphi(jc)*zalp(jc,jzj-1)-ze1*zalp(jc,jzj-2))
          ENDDO
          ze1       = 1.0_wp/ze2
        ENDDO
      ENDDO wave_number_loop


      ! Legendre transform of aerosols

      DO jc=i_startidx,i_endidx
        zfaes(jc,:) = 0.0_wp
        zfael(jc,:) = 0.0_wp
        zfaeu(jc,:) = 0.0_wp
        zfaed(jc,:) = 0.0_wp
      ENDDO !jc

      imn  = 0
      imnc = 0
      imns = 0

      DO jmm = 1, 11
        imn  = imn  + 1
        DO jnn = jmm, 11
          imnc       = imnc + 1
          DO jc=i_startidx,i_endidx
            zfaes(jc,imn) = zfaes(jc,imn)+zalp(jc,imnc)*zaesc(imnc)
            zfael(jc,imn) = zfael(jc,imn)+zalp(jc,imnc)*zaelc(imnc)
            zfaeu(jc,imn) = zfaeu(jc,imn)+zalp(jc,imnc)*zaeuc(imnc)
            zfaed(jc,imn) = zfaed(jc,imn)+zalp(jc,imnc)*zaedc(imnc)
          ENDDO
        ENDDO
        IF(jmm.NE.1) THEN
          imn  = imn+1
          DO jnn = jmm, 11
            imns       = imns + 1
            DO jc=i_startidx,i_endidx
              zfaes(jc,imn) = zfaes(jc,imn)+zalp(jc,imns+11)*zaess(imns)
              zfael(jc,imn) = zfael(jc,imn)+zalp(jc,imns+11)*zaels(imns)
              zfaeu(jc,imn) = zfaeu(jc,imn)+zalp(jc,imns+11)*zaeus(imns)
              zfaed(jc,imn) = zfaed(jc,imn)+zalp(jc,imns+11)*zaeds(imns)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      DO jc=i_startidx,i_endidx

        ! Inverse Fourier transformation

        zcos1   = COS(pt_patch%cells%center(jc,jb)%lon )
        zsin1   = SIN(pt_patch%cells%center(jc,jb)%lon )
        zcos2   = zcos1*zcos1 - zsin1*zsin1
        zsin2   = zsin1*zcos1 + zcos1*zsin1
        zcos3   = zcos2*zcos1 - zsin2*zsin1
        zsin3   = zsin2*zcos1 + zcos2*zsin1
        zcos4   = zcos3*zcos1 - zsin3*zsin1
        zsin4   = zsin3*zcos1 + zcos3*zsin1
        zcos5   = zcos4*zcos1 - zsin4*zsin1
        zsin5   = zsin4*zcos1 + zcos4*zsin1
        zcos6   = zcos5*zcos1 - zsin5*zsin1
        zsin6   = zsin5*zcos1 + zcos5*zsin1
        zcos7   = zcos6*zcos1 - zsin6*zsin1
        zsin7   = zsin6*zcos1 + zcos6*zsin1
        zcos8   = zcos7*zcos1 - zsin7*zsin1
        zsin8   = zsin7*zcos1 + zcos7*zsin1
        zcos9   = zcos8*zcos1 - zsin8*zsin1
        zsin9   = zsin8*zcos1 + zcos8*zsin1
        zcos10  = zcos9*zcos1 - zsin9*zsin1
        zsin10  = zsin9*zcos1 + zcos9*zsin1

        aersea(jc,jb) = zfaes(jc,1) + 2._wp* ( zfaes(jc,2 ) * zcos1 + zfaes(jc,3 ) * zsin1  &
          + zfaes(jc,4 ) * zcos2 + zfaes(jc,5 ) * zsin2    &
          + zfaes(jc,6 ) * zcos3 + zfaes(jc,7 ) * zsin3    &
          + zfaes(jc,8 ) * zcos4 + zfaes(jc,9 ) * zsin4    &
          + zfaes(jc,10) * zcos5 + zfaes(jc,11) * zsin5    &
          + zfaes(jc,12) * zcos6 + zfaes(jc,13) * zsin6    &
          + zfaes(jc,14) * zcos7 + zfaes(jc,15) * zsin7    &
          + zfaes(jc,16) * zcos8 + zfaes(jc,17) * zsin8    &
          + zfaes(jc,18) * zcos9 + zfaes(jc,19) * zsin9    &
          + zfaes(jc,20) * zcos10+ zfaes(jc,21) * zsin10 )

        aerlan(jc,jb) = zfael(jc,1) + 2._wp* ( zfael(jc,2 ) * zcos1 + zfael(jc,3 ) * zsin1  &
          + zfael(jc,4 ) * zcos2 + zfael(jc,5 ) * zsin2    &
          + zfael(jc,6 ) * zcos3 + zfael(jc,7 ) * zsin3    &
          + zfael(jc,8 ) * zcos4 + zfael(jc,9 ) * zsin4    &
          + zfael(jc,10) * zcos5 + zfael(jc,11) * zsin5    &
          + zfael(jc,12) * zcos6 + zfael(jc,13) * zsin6    &
          + zfael(jc,14) * zcos7 + zfael(jc,15) * zsin7    &
          + zfael(jc,16) * zcos8 + zfael(jc,17) * zsin8    &
          + zfael(jc,18) * zcos9 + zfael(jc,19) * zsin9    &
          + zfael(jc,20) * zcos10+ zfael(jc,21) * zsin10 )

        aerurb(jc,jb) = zfaeu(jc,1) + 2._wp* ( zfaeu(jc,2 ) * zcos1 + zfaeu(jc,3 ) * zsin1  &
          + zfaeu(jc,4 ) * zcos2 + zfaeu(jc,5 ) * zsin2    &
          + zfaeu(jc,6 ) * zcos3 + zfaeu(jc,7 ) * zsin3    &
          + zfaeu(jc,8 ) * zcos4 + zfaeu(jc,9 ) * zsin4    &
          + zfaeu(jc,10) * zcos5 + zfaeu(jc,11) * zsin5    &
          + zfaeu(jc,12) * zcos6 + zfaeu(jc,13) * zsin6    &
          + zfaeu(jc,14) * zcos7 + zfaeu(jc,15) * zsin7    &
          + zfaeu(jc,16) * zcos8 + zfaeu(jc,17) * zsin8    &
          + zfaeu(jc,18) * zcos9 + zfaeu(jc,19) * zsin9    &
          + zfaeu(jc,20) * zcos10+ zfaeu(jc,21) * zsin10 )

        aerdes(jc,jb) = zfaed(jc,1) + 2._wp* ( zfaed(jc,2 ) * zcos1 + zfaed(jc,3 ) * zsin1  &
          + zfaed(jc,4 ) * zcos2 + zfaed(jc,5 ) * zsin2    &
          + zfaed(jc,6 ) * zcos3 + zfaed(jc,7 ) * zsin3    &
          + zfaed(jc,8 ) * zcos4 + zfaed(jc,9 ) * zsin4    &
          + zfaed(jc,10) * zcos5 + zfaed(jc,11) * zsin5    &
          + zfaed(jc,12) * zcos6 + zfaed(jc,13) * zsin6    &
          + zfaed(jc,14) * zcos7 + zfaed(jc,15) * zsin7    &
          + zfaed(jc,16) * zcos8 + zfaed(jc,17) * zsin8    &
          + zfaed(jc,18) * zcos9 + zfaed(jc,19) * zsin9    &
          + zfaed(jc,20) * zcos10+ zfaed(jc,21) * zsin10 )

        aersea(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aersea(jc,jb) ) )
        aerlan(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerlan(jc,jb) ) )
        aerurb(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerurb(jc,jb) ) )
        aerdes(jc,jb) = MAX( 0.0_wp, MIN( 1.0_wp, aerdes(jc,jb) ) )

      ENDDO !jc

    ENDDO !jb 

  END SUBROUTINE init_aerosol_dstrb_tanre

  SUBROUTINE init_aerosol_props_tanre_rrtm

   ! the following aerosol types (second array index) are considered:
   ! 1 : continental
   ! 2 : maritime
   ! 3 : urban
   ! 4 : vulcano ashes
   ! 5 : stratospheric background aerosol (SB)

   !absorption
   zaea_rrtm=RESHAPE( (/ &
     &0.0469_wp,0.0447_wp,0.0397_wp,0.0432_wp,0.0336_wp,0.0384_wp,0.0560_wp,0.0734_wp,&
     &0.0484_wp,0.0421_wp,0.0314_wp,0.0249_wp,0.0247_wp,0.0240_wp,0.0224_wp,0.0260_wp,&
     &0.0260_wp,0.0505_wp,0.0356_wp,0.0375_wp,0.0545_wp,0.0730_wp,0.0792_wp,0.0895_wp,&
     &0.0936_wp,0.1111_wp,0.1401_wp,0.2029_wp,0.4395_wp,0.0458_wp,                    &
     &0.1229_wp,0.1375_wp,0.1778_wp,0.1841_wp,0.1618_wp,0.0932_wp,0.0634_wp,0.0723_wp,&
     &0.0620_wp,0.0636_wp,0.1372_wp,0.0383_wp,0.0405_wp,0.0345_wp,0.0219_wp,0.1185_wp,&
     &0.1185_wp,0.1457_wp,0.0146_wp,0.0131_wp,0.0119_wp,0.0143_wp,0.0143_wp,0.0132_wp,&
     &0.0100_wp,0.0106_wp,0.0139_wp,0.0269_wp,0.0925_wp,0.0781_wp,                    &
     &0.0164_wp,0.0175_wp,0.0197_wp,0.0240_wp,0.0189_wp,0.0217_wp,0.0319_wp,0.0471_wp,&
     &0.0386_wp,0.0361_wp,0.0339_wp,0.0368_wp,0.0418_wp,0.0435_wp,0.0453_wp,0.0543_wp,&
     &0.0543_wp,0.0804_wp,0.0802_wp,0.0894_wp,0.1116_wp,0.1133_wp,0.0901_wp,0.2015_wp,&
     &0.3239_wp,0.3727_wp,0.5078_wp,0.6559_wp,0.8995_wp,0.0322_wp,                    &
     &0.0125_wp,0.0110_wp,0.0166_wp,0.0159_wp,0.0154_wp,0.0242_wp,0.0292_wp,0.0248_wp,&
     &0.0174_wp,0.0097_wp,0.0058_wp,0.0037_wp,0.0032_wp,0.0033_wp,0.0037_wp,0.0074_wp,&
     &0.0074_wp,0.0088_wp,0.0031_wp,0.0010_wp,0.0012_wp,0.0087_wp,0.0194_wp,0.0280_wp,&
     &0.0410_wp,0.0555_wp,0.0733_wp,0.1103_wp,0.2602_wp,0.0200_wp,                    &
     &0.0060_wp,0.0117_wp,0.0269_wp,0.0222_wp,0.0195_wp,0.0398_wp,0.0733_wp,0.1091_wp,&     ! SB
     &0.1124_wp,0.0415_wp,0.0424_wp,0.0495_wp,0.0451_wp,0.0484_wp,0.0540_wp,0.0735_wp,&     ! SB
     &0.0735_wp,0.0188_wp,0.0021_wp,0.0014_wp,0.0007_wp,0.0002_wp,0.0000_wp,0.0000_wp,&     ! SB
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0628_wp/),(/jpsw+jpband,5/))      ! SB

   !scattering
   zaes_rrtm=RESHAPE( (/ &
     &0.0235_wp,0.0321_wp,0.0381_wp,0.0401_wp,0.0402_wp,0.0455_wp,0.0480_wp,0.0481_wp,&
     &0.0308_wp,0.0430_wp,0.0353_wp,0.0583_wp,0.0685_wp,0.0718_wp,0.0761_wp,0.0827_wp,&
     &0.0827_wp,0.0872_wp,0.1118_wp,0.1308_wp,0.1685_wp,0.2320_wp,0.2774_wp,0.4229_wp,&
     &0.6696_wp,0.9410_wp,1.2715_wp,1.5532_wp,1.6360_wp,0.0448_wp,                    &
     &0.0302_wp,0.0577_wp,0.0715_wp,0.0592_wp,0.0500_wp,0.0719_wp,0.1451_wp,0.1901_wp,&
     &0.2028_wp,0.2536_wp,0.2315_wp,0.3699_wp,0.4437_wp,0.4762_wp,0.5228_wp,0.5163_wp,&
     &0.5163_wp,0.4484_wp,0.6568_wp,0.7082_wp,0.7559_wp,0.8051_wp,0.8305_wp,0.8788_wp,&
     &0.9381_wp,0.9985_wp,1.0708_wp,1.1366_wp,1.1403_wp,0.1669_wp,                    &
     &0.0019_wp,0.0027_wp,0.0034_wp,0.0034_wp,0.0037_wp,0.0048_wp,0.0063_wp,0.0081_wp,&
     &0.0030_wp,0.0052_wp,0.0056_wp,0.0094_wp,0.0126_wp,0.0140_wp,0.0160_wp,0.0206_wp,&
     &0.0206_wp,0.0266_wp,0.0407_wp,0.0525_wp,0.0784_wp,0.1682_wp,0.2616_wp,0.2775_wp,&
     &0.3969_wp,0.6920_wp,0.9869_wp,1.2506_wp,1.3541_wp,0.0058_wp,                    &
     &0.0000_wp,0.0002_wp,0.0003_wp,0.0003_wp,0.0007_wp,0.0024_wp,0.0048_wp,0.0050_wp,&
     &0.0017_wp,0.0030_wp,0.0048_wp,0.0104_wp,0.0173_wp,0.0200_wp,0.0241_wp,0.0376_wp,&
     &0.0376_wp,0.0628_wp,0.0318_wp,0.0132_wp,0.0196_wp,0.1666_wp,0.3869_wp,0.5767_wp,&
     &0.8068_wp,0.9621_wp,1.0634_wp,1.0735_wp,0.9208_wp,0.0060_wp,                    &
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0001_wp,0.0003_wp,0.0006_wp,0.0008_wp,&     ! SB
     &0.0005_wp,0.0003_wp,0.0008_wp,0.0013_wp,0.0024_wp,0.0030_wp,0.0040_wp,0.0059_wp,&     ! SB
     &0.0059_wp,0.0123_wp,0.0236_wp,0.0384_wp,0.0651_wp,0.1246_wp,0.1801_wp,0.3807_wp,&     ! SB
     &0.7105_wp,1.0514_wp,1.3754_wp,1.5334_wp,1.5495_wp,0.0009_wp/),(/jpsw+jpband,5/))      ! SB

   !asymmetry factor
   zaeg_rrtm=RESHAPE( (/ &
     &0.5856_wp,0.6513_wp,0.7064_wp,0.7354_wp,0.7608_wp,0.7356_wp,0.7118_wp,0.7233_wp,&
     &0.8772_wp,0.8244_wp,0.8601_wp,0.8315_wp,0.8064_wp,0.7988_wp,0.7896_wp,0.7796_wp,&
     &0.7796_wp,0.7841_wp,0.7459_wp,0.7204_wp,0.6758_wp,0.6426_wp,0.6362_wp,0.6324_wp,&
     &0.6319_wp,0.6382_wp,0.6463_wp,0.6571_wp,0.6845_wp,0.7636_wp,                    &
     &0.4380_wp,0.5305_wp,0.5771_wp,0.6318_wp,0.6832_wp,0.7491_wp,0.7445_wp,0.7366_wp,&
     &0.7646_wp,0.7618_wp,0.7419_wp,0.7689_wp,0.7620_wp,0.7583_wp,0.7522_wp,0.7556_wp,&
     &0.7556_wp,0.8332_wp,0.7968_wp,0.7853_wp,0.7778_wp,0.7684_wp,0.7635_wp,0.7571_wp,&
     &0.7499_wp,0.7453_wp,0.7439_wp,0.7444_wp,0.7564_wp,0.7517_wp,                    &
     &0.5757_wp,0.6206_wp,0.6523_wp,0.7001_wp,0.6923_wp,0.6310_wp,0.5574_wp,0.5094_wp,&
     &0.7863_wp,0.6673_wp,0.6354_wp,0.6324_wp,0.6077_wp,0.6005_wp,0.5932_wp,0.5872_wp,&
     &0.5872_wp,0.5883_wp,0.5853_wp,0.5844_wp,0.5729_wp,0.6021_wp,0.6293_wp,0.5850_wp,&
     &0.5784_wp,0.5920_wp,0.6004_wp,0.6123_wp,0.6413_wp,0.6070_wp,                    &
     &0.0158_wp,0.0374_wp,0.0451_wp,0.0520_wp,0.0733_wp,0.1215_wp,0.1578_wp,0.1668_wp,&
     &0.1406_wp,0.1711_wp,0.2084_wp,0.2638_wp,0.3108_wp,0.3265_wp,0.3492_wp,0.4042_wp,&
     &0.4042_wp,0.4702_wp,0.5106_wp,0.5516_wp,0.5821_wp,0.6213_wp,0.6334_wp,0.6658_wp,&
     &0.6886_wp,0.6989_wp,0.7050_wp,0.7143_wp,0.7553_wp,0.2521_wp,                    &
     &0.0021_wp,0.0039_wp,0.0061_wp,0.0078_wp,0.0109_wp,0.0161_wp,0.0201_wp,0.0206_wp,&     ! SB
     &0.0217_wp,0.0320_wp,0.0428_wp,0.0583_wp,0.0773_wp,0.0856_wp,0.0985_wp,0.1310_wp,&     ! SB
     &0.1310_wp,0.1906_wp,0.2625_wp,0.3154_wp,0.3869_wp,0.4787_wp,0.5279_wp,0.6272_wp,&     ! SB
     &0.6941_wp,0.7286_wp,0.7358_wp,0.7177_wp,0.6955_wp,0.0616_wp/),(/jpsw+jpband,5/))      ! SB

  END SUBROUTINE init_aerosol_props_tanre_rrtm

  SUBROUTINE init_aerosol_props_tanre_rg

   ! the following aerosol types (second array index) are considered:
   ! 1 : continental
   ! 2 : maritime
   ! 3 : urban
   ! 4 : vulcano ashes
   ! 5 : stratospheric background aerosol (SB)

    !absorption
    zaea_rg=RESHAPE((/0.0477_wp, 0.0875_wp,  0.1198_wp, 0.0458_wp, &
      0.0387_wp, 0.0439_wp,  0.0599_wp, 0.0396_wp, &
      0.0381_wp, 0.0129_wp,  0.0130_wp, 0.1304_wp, &
      0.1757_wp, 0.0949_wp,  0.0653_wp, 0.0795_wp, &
      0.0962_wp, 0.2046_wp,  0.4116_wp, 0.0169_wp, &
      0.0204_wp, 0.0263_wp,  0.0348_wp, 0.0361_wp, &
      0.0030_wp, 0.0271_wp,  0.0613_wp, 0.0118_wp, &
      0.0160_wp, 0.0231_wp,  0.0287_wp, 0.0127_wp, &
      0.0103_wp, 0.000016_wp,0.0000_wp, 0.0087_wp, &
      0.0238_wp, 0.0511_wp,  0.0734_wp, 0.0809_wp/),(/jpspec,5/))

    !scattering
    zaes_rg=RESHAPE((/0.1407_wp, 0.4256_wp,  1.0066_wp, 0.0279_wp, &
      0.0391_wp, 0.0445_wp,  0.0485_wp, 0.0362_wp, &
      0.6746_wp, 0.8761_wp,  1.0139_wp, 0.0443_wp, &
      0.0624_wp, 0.0921_wp,  0.1491_wp, 0.2327_wp, &
      0.0605_wp, 0.2761_wp,  0.7449_wp, 0.0023_wp, &
      0.0034_wp, 0.0051_wp,  0.0065_wp, 0.0045_wp, &
      0.0284_wp, 0.5524_wp,  0.9683_wp, 0.0001_wp, &
      0.0004_wp, 0.0024_wp,  0.0049_wp, 0.0030_wp, &
      0.0467_wp, 0.3854_wp,  1.1008_wp, 0.0000_wp, &
      0.00005_wp,0.0004_wp,  0.0006_wp, 0.0006_wp/),(/jpspec,5/))

   !asymmetry factor   
    zaeg_rg=RESHAPE((/0.6989_wp, 0.6329_wp,  0.6418_wp, 0.6243_wp, &
      0.7299_wp, 0.7430_wp,  0.7086_wp, 0.8569_wp, &
      0.7833_wp, 0.7575_wp,  0.7456_wp, 0.4997_wp, &
      0.6130_wp, 0.7440_wp,  0.7426_wp, 0.7590_wp, &
      0.5753_wp, 0.5867_wp,  0.5957_wp, 0.6027_wp, &
      0.6766_wp, 0.6117_wp,  0.5439_wp, 0.6905_wp, &
      0.5170_wp, 0.6674_wp,  0.7004_wp, 0.0340_wp, &
      0.0570_wp, 0.1289_wp,  0.1597_wp, 0.1906_wp, &
      0.3751_wp, 0.6353_wp,  0.7259_wp, 0.0037_wp, &
      0.0083_wp, 0.0177_wp,  0.0201_wp, 0.0332_wp/),(/jpspec,5/))
    
  END SUBROUTINE init_aerosol_props_tanre_rg

  SUBROUTINE init_aerosol_props_tegen_rg

  ! the following aerosol types (second array index) are considered:
  ! 1. continental, 2. maritime, 3. desert, 4. urban, 5. stratospheric background

    zaea_rg=RESHAPE((/0.0345_wp,0.0511_wp,0.0847_wp,0.0336_wp,&
      0.0499_wp,0.0364_wp,0.0382_wp,0.0260_wp,&
      0.0457_wp,0.0018_wp,0.0015_wp,0.1361_wp,&
      0.2346_wp,0.1177_wp,0.0684_wp,0.0808_wp,&
      0.0707_wp,0.0689_wp,0.1557_wp,0.1258_wp,&
      0.1588_wp,0.1973_wp,0.2766_wp,0.1134_wp,&
      0.0597_wp,0.1077_wp,0.2095_wp,0.0299_wp,&
      0.0456_wp,0.0358_wp,0.0377_wp,0.0304_wp,&
      0.0103_wp, 0.000016_wp,0.0000_wp, 0.0087_wp, &
      0.0238_wp, 0.0511_wp,  0.0734_wp, 0.0809_wp/),(/jpspec,5/))

     
    zaes_rg=RESHAPE((/0.1030_wp,0.3977_wp,1.0680_wp,0.0084_wp,&
      0.0142_wp,0.0191_wp,0.0234_wp,0.0140_wp,&
      0.7894_wp,0.9734_wp,1.0110_wp,0.0307_wp,&
      0.0531_wp,0.0546_wp,0.0839_wp,0.2142_wp,&
      0.7157_wp,0.8698_wp,0.8604_wp,0.0645_wp,&
      0.0781_wp,0.1256_wp,0.2317_wp,0.1409_wp,&
      0.0859_wp,0.3442_wp,0.9496_wp,0.0067_wp,&
      0.0113_wp,0.0153_wp,0.0187_wp,0.0113_wp,&
      0.0467_wp, 0.3854_wp,  1.1008_wp, 0.0000_wp, &
      0.00005_wp,0.0004_wp,  0.0006_wp, 0.0006_wp/),(/jpspec,5/))
     
     
    zaeg_rg=RESHAPE((/0.6562_wp,0.6614_wp,0.7109_wp,0.5043_wp,&
      0.6486_wp,0.6814_wp,0.6489_wp,0.7799_wp,&
      0.8105_wp,0.7906_wp,0.7947_wp,0.4374_wp,&
      0.5203_wp,0.7076_wp,0.7246_wp,0.7535_wp,&
      0.6932_wp,0.6962_wp,0.7402_wp,0.4029_wp,&
      0.5587_wp,0.5618_wp,0.4520_wp,0.7120_wp,&
      0.6462_wp,0.6510_wp,0.6955_wp,0.5041_wp,&
      0.6482_wp,0.6805_wp,0.6477_wp,0.7753_wp,&
      0.3751_wp, 0.6353_wp,  0.7259_wp, 0.0037_wp, &
      0.0083_wp, 0.0177_wp,  0.0201_wp, 0.0332_wp/),(/jpspec,5/))
    
  END SUBROUTINE init_aerosol_props_tegen_rg
  
  SUBROUTINE init_aerosol_props_tegen_rrtm

  ! the following aerosol types (second array index) are considered:
  ! 1. continental, 2. maritime, 3. desert, 4. urban, 5. stratospheric background (SB)

   !absorption
   zaea_rrtm=RESHAPE( (/ &
     &0.0304_wp,0.0367_wp,0.0462_wp,0.0566_wp,0.0496_wp,0.0336_wp,0.0355_wp,0.0456_wp,&
     &0.0272_wp,0.0264_wp,0.0290_wp,0.0156_wp,0.0165_wp,0.0157_wp,0.0138_wp,0.0401_wp,&
     &0.0401_wp,0.0760_wp,0.0214_wp,0.0227_wp,0.0295_wp,0.0394_wp,0.0431_wp,0.0519_wp,&
     &0.0611_wp,0.0774_wp,0.1012_wp,0.1412_wp,0.2632_wp,0.0324_wp,                    &
     &0.1096_wp,0.1614_wp,0.2294_wp,0.2506_wp,0.2242_wp,0.1190_wp,0.0680_wp,0.0664_wp,&
     &0.0656_wp,0.0749_wp,0.1250_wp,0.0425_wp,0.0498_wp,0.0425_wp,0.0259_wp,0.1619_wp,&
     &0.1619_wp,0.2152_wp,0.0139_wp,0.0119_wp,0.0046_wp,0.0036_wp,0.0020_wp,0.0016_wp,&
     &0.0012_wp,0.0013_wp,0.0016_wp,0.0035_wp,0.0147_wp,0.0882_wp,                    &
     &0.0974_wp,0.1529_wp,0.1643_wp,0.1373_wp,0.1753_wp,0.1923_wp,0.2804_wp,0.2426_wp,&
     &0.1263_wp,0.1321_wp,0.0979_wp,0.0664_wp,0.0360_wp,0.0311_wp,0.0325_wp,0.0833_wp,&
     &0.0833_wp,0.1170_wp,0.0739_wp,0.0631_wp,0.0604_wp,0.0628_wp,0.0645_wp,0.0677_wp,&
     &0.0843_wp,0.1328_wp,0.2224_wp,0.3022_wp,0.3579_wp,0.1820_wp,                    &
     &0.0267_wp,0.0329_wp,0.0420_wp,0.0515_wp,0.0461_wp,0.0332_wp,0.0354_wp,0.0447_wp,&
     &0.0303_wp,0.0306_wp,0.0342_wp,0.0248_wp,0.0274_wp,0.0276_wp,0.0271_wp,0.0526_wp,&
     &0.0526_wp,0.0903_wp,0.0450_wp,0.0492_wp,0.0596_wp,0.0754_wp,0.0842_wp,0.1082_wp,&
     &0.1429_wp,0.1926_wp,0.2595_wp,0.3379_wp,0.4761_wp,0.0340_wp,                    &
     &0.0060_wp,0.0117_wp,0.0269_wp,0.0222_wp,0.0195_wp,0.0398_wp,0.0733_wp,0.1091_wp,&     ! SB
     &0.1124_wp,0.0415_wp,0.0424_wp,0.0495_wp,0.0451_wp,0.0484_wp,0.0540_wp,0.0735_wp,&     ! SB
     &0.0735_wp,0.0188_wp,0.0021_wp,0.0014_wp,0.0007_wp,0.0002_wp,0.0000_wp,0.0000_wp,&     ! SB
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0628_wp/),(/jpsw+jpband,5/))      ! SB

   !scattering
   zaes_rrtm=RESHAPE( (/ &
     &0.0060_wp,0.0107_wp,0.0134_wp,0.0150_wp,0.0152_wp,0.0200_wp,0.0232_wp,0.0211_wp,&
     &0.0112_wp,0.0186_wp,0.0128_wp,0.0260_wp,0.0339_wp,0.0368_wp,0.0409_wp,0.0527_wp,&
     &0.0527_wp,0.0621_wp,0.0715_wp,0.0929_wp,0.1276_wp,0.1895_wp,0.2350_wp,0.3930_wp,&
     &0.6641_wp,0.9834_wp,1.3737_wp,1.7160_wp,1.9115_wp,0.0198_wp,                    &
     &0.0188_wp,0.0421_wp,0.0576_wp,0.0547_wp,0.0430_wp,0.0367_wp,0.0806_wp,0.1209_wp,&
     &0.1681_wp,0.2257_wp,0.2440_wp,0.3622_wp,0.4540_wp,0.5026_wp,0.5765_wp,0.5986_wp,&
     &0.5986_wp,0.5225_wp,0.7420_wp,0.8311_wp,0.8970_wp,0.9444_wp,0.9637_wp,0.9763_wp,&
     &0.9855_wp,1.0034_wp,1.0337_wp,1.0640_wp,1.0795_wp,0.1312_wp,                    &
     &0.0458_wp,0.0823_wp,0.0667_wp,0.0642_wp,0.1080_wp,0.1471_wp,0.2422_wp,0.1216_wp,&
     &0.0717_wp,0.1616_wp,0.2027_wp,0.3042_wp,0.4045_wp,0.4369_wp,0.4685_wp,0.5043_wp,&
     &0.5043_wp,0.5782_wp,0.6898_wp,0.7477_wp,0.7926_wp,0.8320_wp,0.8503_wp,0.8736_wp,&
     &0.8874_wp,0.8737_wp,0.8278_wp,0.7857_wp,0.7571_wp,0.1714_wp,                    &
     &0.0048_wp,0.0085_wp,0.0107_wp,0.0119_wp,0.0121_wp,0.0160_wp,0.0185_wp,0.0170_wp,&
     &0.0090_wp,0.0150_wp,0.0103_wp,0.0210_wp,0.0274_wp,0.0298_wp,0.0332_wp,0.0430_wp,&
     &0.0430_wp,0.0485_wp,0.0593_wp,0.0776_wp,0.1073_wp,0.1610_wp,0.2008_wp,0.3398_wp,&
     &0.5809_wp,0.8701_wp,1.2309_wp,1.5535_wp,1.7368_wp,0.0159_wp,                    &
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0001_wp,0.0003_wp,0.0006_wp,0.0008_wp,&     ! SB
     &0.0005_wp,0.0003_wp,0.0008_wp,0.0013_wp,0.0024_wp,0.0030_wp,0.0040_wp,0.0059_wp,&     ! SB
     &0.0059_wp,0.0123_wp,0.0236_wp,0.0384_wp,0.0651_wp,0.1246_wp,0.1801_wp,0.3807_wp,&     ! SB
     &0.7105_wp,1.0514_wp,1.3754_wp,1.5334_wp,1.5495_wp,0.0009_wp/),(/jpsw+jpband,5/))      ! SB

   !asymmetry factor
   zaeg_rrtm=RESHAPE( (/ &
     &0.4388_wp,0.5396_wp,0.6191_wp,0.6535_wp,0.6876_wp,0.6718_wp,0.6493_wp,0.6782_wp,&
     &0.7958_wp,0.7537_wp,0.7757_wp,0.7821_wp,0.7583_wp,0.7487_wp,0.7351_wp,0.6917_wp,&
     &0.6917_wp,0.6989_wp,0.6982_wp,0.6726_wp,0.6426_wp,0.6294_wp,0.6337_wp,0.6582_wp,&
     &0.6850_wp,0.7061_wp,0.7212_wp,0.7306_wp,0.7417_wp,0.6978_wp,                    &
     &0.4062_wp,0.4507_wp,0.4878_wp,0.5302_wp,0.5850_wp,0.6962_wp,0.7242_wp,0.7293_wp,&
     &0.7414_wp,0.7484_wp,0.7607_wp,0.7785_wp,0.7805_wp,0.7785_wp,0.7724_wp,0.7690_wp,&
     &0.7690_wp,0.8348_wp,0.8316_wp,0.8170_wp,0.8074_wp,0.7990_wp,0.7954_wp,0.7897_wp,&
     &0.7884_wp,0.7927_wp,0.8001_wp,0.8057_wp,0.8076_wp,0.7462_wp,                    &
     &0.4219_wp,0.3928_wp,0.5306_wp,0.6229_wp,0.5544_wp,0.5454_wp,0.4353_wp,0.5736_wp,&
     &0.7502_wp,0.6957_wp,0.7038_wp,0.6881_wp,0.6740_wp,0.6739_wp,0.6784_wp,0.6969_wp,&
     &0.6969_wp,0.7068_wp,0.6965_wp,0.6918_wp,0.6904_wp,0.6911_wp,0.6915_wp,0.6952_wp,&
     &0.7080_wp,0.7326_wp,0.7689_wp,0.8000_wp,0.8206_wp,0.5788_wp,                    &
     &0.4387_wp,0.5394_wp,0.6187_wp,0.6531_wp,0.6871_wp,0.6712_wp,0.6482_wp,0.6756_wp,&
     &0.7930_wp,0.7498_wp,0.7685_wp,0.7766_wp,0.7520_wp,0.7419_wp,0.7277_wp,0.6828_wp,&
     &0.6828_wp,0.6875_wp,0.6872_wp,0.6622_wp,0.6333_wp,0.6209_wp,0.6250_wp,0.6479_wp,&
     &0.6725_wp,0.6912_wp,0.7043_wp,0.7129_wp,0.7254_wp,0.6956_wp,                    &
     &0.0021_wp,0.0039_wp,0.0061_wp,0.0078_wp,0.0109_wp,0.0161_wp,0.0201_wp,0.0206_wp,&     ! SB
     &0.0217_wp,0.0320_wp,0.0428_wp,0.0583_wp,0.0773_wp,0.0856_wp,0.0985_wp,0.1310_wp,&     ! SB
     &0.1310_wp,0.1906_wp,0.2625_wp,0.3154_wp,0.3869_wp,0.4787_wp,0.5279_wp,0.6272_wp,&     ! SB
     &0.6941_wp,0.7286_wp,0.7358_wp,0.7177_wp,0.6955_wp,0.0616_wp/),(/jpsw+jpband,5/))      ! SB

  END SUBROUTINE init_aerosol_props_tegen_rrtm

  ! Very simple parameterization of source and sink terms for prognostic 2D aerosol fields
  !
  SUBROUTINE prog_aerosol_2D (nproma,jcs,jce,dtime,iprog_aero,aerosol,aercl_ss,aercl_or,aercl_bc,aercl_su,aercl_du,       &
                              rr_gsp,sr_gsp,rr_con,sr_con,gust_dyn,gust_con,soiltype,plcov_t,frac_t,w_so_t,t_so_t,h_snow_t)
                              
    
    INTEGER,  INTENT(in)    :: nproma, jcs, jce, iprog_aero
    REAL(wp), INTENT(in)    :: dtime

    REAL(wp), INTENT(inout) :: aerosol(:,:)
    REAL(wp), INTENT(in)    :: aercl_ss(:),aercl_or(:),aercl_bc(:),aercl_su(:),aercl_du(:),  &
                               rr_gsp(:),sr_gsp(:),rr_con(:),sr_con(:),                      &
                               gust_dyn(:),gust_con(:),plcov_t(:,:),frac_t(:,:),w_so_t(:,:), &
                               t_so_t(:,:),h_snow_t(:,:)
    INTEGER,  INTENT(in)    :: soiltype(:)

    INTEGER :: jc, js, jt

    REAL(wp) :: relax_scale(2), relax_fac(2), od_clim(nclass_aero), vfac, wsofac, tsofac, minfrac, &
                plcfac, snwfac, dustsrc(ntiles_lnd), ts_dustsrc, ts_saltsrc, ts_orgsrc, ts_bcsrc, ts_susrc, &
                washout, washout_scale

    ! Fractions of climatology used as target values for relaxation 
    ! (tuned to approximately balance the source terms for non-dust aerosol classes)
    relax_scale(1)  = 0.7_wp
    relax_scale(2)  = 0.4_wp

    ! Relaxation time scales
    relax_fac(1) = 1._wp/(3._wp*86400._wp)  ! 3 days for aerosols with shallow vertical extent
    relax_fac(2) = 1._wp/(10._wp*86400._wp) ! 10 days for aerosols with deep vertical extent (dust)

    ! Time scales for sources
    ts_dustsrc = 1._wp/(12.5_wp*86400._wp) ! 12.5 days for dust source terms
    ts_saltsrc = 1._wp/(2.5_wp*86400._wp)  ! 2.5  days for sea salt
    ts_orgsrc  = 1._wp/(2.5_wp*86400._wp)  ! 2.5  days for organic aerosol
    ts_bcsrc   = 1._wp/(2.5_wp*86400._wp)  ! 2.5  days for black carbon
    ts_susrc   = 1._wp/(2.5_wp*86400._wp)  ! 2.5  days for sulfate aerosol

    ! Washout scale for dust
    washout_scale = 1._wp/7.5_wp  ! e-folding scale 7.5 mm WE precipitation

    minfrac = 0.025_wp ! minimum allowed fraction of climatological aerosol optical depth

    ! optical depth offsets for climatology-based source terms
    od_clim(iss)  = 0.005_wp
    od_clim(iorg) = 0.015_wp
    od_clim(ibc)  = 0.002_wp
    od_clim(iso4) = 0.015_wp
    od_clim(idu)  = 0.075_wp

    ! Prediction of mineral dust; other aerosol classes are treated prognostically only if iprog_aero=2

    ! Relaxation to scaled climatology
    DO jc = jcs, jce
      aerosol(jc,idu)  = aerosol(jc,idu)  + dtime*relax_fac(2)*(relax_scale(2)*aercl_du(jc)-aerosol(jc,idu))
    ENDDO

    ! Sources and sinks for mineral dust, including weak climatology-based source term
    DO jc = jcs, jce
      IF (soiltype(jc) >= 3 .AND. soiltype(jc) <= 8) THEN ! land excluding glaciers and rock
        js = soiltype(jc)
        DO jt = 1, ntiles_lnd ! snow tiles are excluded a priori; the snow depth criterion is relevant without snow tiles only
          wsofac = MAX(0._wp,w_so_t(jc,jt)/dzsoil(1)-1.5_wp*cadp(js)) / (cfcap(js)-1.5_wp*cadp(js)) ! soil moisture < field cap.
          plcfac = MAX(0._wp, (0.75_wp-plcov_t(jc,jt))/0.75_wp)             ! plant cover < 75%
          snwfac = MAX(0._wp, (0.05_wp-h_snow_t(jc,jt))*20._wp)             ! snow depth < 5 cm
          tsofac = MAX(0._wp, MIN(1._wp,0.4_wp*(t_so_t(jc,jt)-270.65_wp)))  ! top soil layer warmer than -2.5 deg C
          vfac   = MAX(0._wp, gust_dyn(jc)+gust_con(jc) - (7.5_wp+17.5_wp*wsofac)) ! gusts > soil-moisture dependent threshold
          dustsrc(jt) = dtime*ts_dustsrc*( vfac*plcfac*snwfac*tsofac +        &
                   MAX(0._wp,aercl_du(jc)-MAX(od_clim(idu),aerosol(jc,idu)) ) )
        ENDDO
        ! washout of mineral dust by precipitation: convective precip is counted only by 50% because it
        ! is assumed not to cover the whole grid box
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+0.5_wp*(rr_con(jc)+sr_con(jc)))*aerosol(jc,idu)
        aerosol(jc,idu)  = aerosol(jc,idu) - washout + SUM(dustsrc(1:ntiles_lnd)*frac_t(jc,1:ntiles_lnd))
      ENDIF
    ENDDO

    ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
    DO jc = jcs, jce
      aerosol(jc,idu)  = MAX(aerosol(jc,idu),  minfrac*aercl_du(jc))
    ENDDO

    IF (iprog_aero == 2) THEN

      ! Relaxation to scaled climatology
      DO jc = jcs, jce
        aerosol(jc,iss)  = aerosol(jc,iss)  + dtime*relax_fac(1)*(relax_scale(1)*aercl_ss(jc)-aerosol(jc,iss))
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime*relax_fac(1)*(relax_scale(1)*aercl_or(jc)-aerosol(jc,iorg))
        aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime*relax_fac(1)*(relax_scale(1)*aercl_bc(jc)-aerosol(jc,ibc))
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime*relax_fac(1)*(relax_scale(1)*aercl_su(jc)-aerosol(jc,iso4))
      ENDDO

      ! Sources are specified where either the climatological value exceeds the global average
      ! sources of sea salt and black carbon are restricted to water and land surfaces, respectively
      DO jc = jcs, jce
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime*ts_susrc*MAX(0._wp,aercl_su(jc)-MAX(od_clim(iso4),0.3_wp*aerosol(jc,iso4)))
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime*ts_orgsrc*MAX(0._wp,aercl_or(jc)-MAX(od_clim(iorg),0.3_wp*aerosol(jc,iorg)))
        IF (soiltype(jc) == 9) THEN ! water
          aerosol(jc,iss)  = aerosol(jc,iss)  + dtime*ts_saltsrc*MAX(0._wp,aercl_ss(jc)-od_clim(iss))
        ELSE                        ! land
          aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime*ts_bcsrc* MAX(0._wp,aercl_bc(jc)-od_clim(ibc))
        ENDIF
      ENDDO

      ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
      DO jc = jcs, jce
        aerosol(jc,iss)  = MAX(aerosol(jc,iss),  minfrac*aercl_ss(jc))
        aerosol(jc,iorg) = MAX(aerosol(jc,iorg), minfrac*aercl_or(jc))
        aerosol(jc,ibc)  = MAX(aerosol(jc,ibc),  minfrac*aercl_bc(jc))
        aerosol(jc,iso4) = MAX(aerosol(jc,iso4), minfrac*aercl_su(jc))
      ENDDO

    ENDIF

  END SUBROUTINE prog_aerosol_2D


  ! Tuning of longwave absorption coefficient of mineral dust in order to reduce cold bias in the Saharan region
  !
  SUBROUTINE tune_dust (lat,lon,iend,tunefac)

    REAL(wp), INTENT(in) :: lat(:), lon(:)
    INTEGER,  INTENT(in) :: iend

    REAL(wp), INTENT(out) :: tunefac(:,:)

    INTEGER :: jc, jb
    REAL(wp) :: maxfac

    DO jb = 1, jpband
      maxfac = tune_dust_abs*5._wp*(jpband-MAX(8,jb))/REAL(jpband-8,wp)
      DO jc = 1, iend
        tunefac(jc,jb) = 1._wp + maxfac*(1._wp - MIN(1._wp,((rad2deg*lat(jc)-15._wp)/20._wp)**4)) * &
         (1._wp - MIN(1._wp,((rad2deg*lon(jc)-20._wp)/50._wp)**4))
      ENDDO
    ENDDO


  END SUBROUTINE tune_dust

END MODULE mo_aerosol_util

