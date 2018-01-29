!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Prepares and provides cloud optical properties
!!
!! @remarks
!!   This code makes use tabulated Mie Calculations for the RRTM band structure
!!   prepared by Stefan Kinne.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Rewrite and synthesis of ECHAM5 code, particularly rad_int.f90 whose
!!   contributers included:  M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,
!!   MPI-M (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas
!!   MPI-M (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); 
!!   S.J. Lorenz, MPI-M (2007-11).
!! 
!
MODULE mo_psrad_cloud_optics

  USE mo_psrad_general, ONLY: wp, finish, pi, rhoh2o
#ifndef PSRAD_ONLY
  USE mo_echam_cld_config, ONLY: echam_cld_config
#endif

  USE mo_psrad_io, ONLY: psrad_io_open, psrad_io_close, &
    psrad_io_copy_double

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_cloud_optics, cloud_optics

  INTEGER, PARAMETER :: &
    n_mdl_bnds = 30, & ! n of bands in effective radius table
    n_sizes = 61 ! n of bands in effective radius table

  REAL (wp), PARAMETER :: &
    ccwmin = 1.e-7_wp, & ! min condensate for lw cloud opacity
    zkap_cont = 1.143_wp, & ! continental (Martin et al. ) breadth param
    zkap_mrtm = 1.077_wp ! maritime (Martin et al.) breadth parameter

  REAL (wp), PARAMETER, DIMENSION(16) ::rebcug = (/ &
    0.718_wp, 0.726_wp, 1.136_wp, 1.320_wp, &
    1.505_wp, 1.290_wp, 0.911_wp, 0.949_wp, &
    1.021_wp, 1.193_wp, 1.279_wp, 0.626_wp, &
    0.647_wp, 0.668_wp, 0.690_wp, 0.690_wp/)

  REAL (wp), PARAMETER, DIMENSION(16) ::rebcuh = (/ &
    0.0069_wp, 0.0060_wp, 0.0024_wp, 0.0004_wp, &
    -0.0016_wp, 0.0003_wp, 0.0043_wp, 0.0038_wp, &
    0.0030_wp, 0.0013_wp, 0.0005_wp, 0.0054_wp, &
    0.0052_wp, 0.0050_wp, 0.0048_wp, 0.0048_wp/)

  LOGICAL, SAVE   :: l_variable_inhoml = .TRUE.

  !LOOKS LIKE A HACK...
  real(wp), public, save :: & ! values for liquid-cloud inhomogeneity factor:
    zinhomi, & ! ice-cloud inhomogeneity factor
    zinhoml1, & ! without convection
    zinhoml2, & ! with shallow convection  
    zinhoml3, & ! deep/mid-level convection  
    zinpar ! exponent for variable lquid-cloud inhomogeneity
  REAL (wp), SAVE :: &
    wavenumber(n_mdl_bnds), & ! effective wavenumber for table band
    wavelength(n_mdl_bnds), & ! effective wavelength for table band
    re_droplet(n_sizes), & ! effective drop radius for table sizes
    re_crystal(n_sizes), & ! effective ice radius for table sizes
    reimin, & ! min ice effective radius (microns)
    reimax, & ! max ice effective radius (microns)
    relmin, & ! min liquid effective radius (microns)
    relmax, & ! max liquid effective radius (microns)
    del_rei, & ! ice effective radius (microns) increment
    del_rel, & ! liq effective radius (microns) increment
    z_ext_l(n_sizes, n_mdl_bnds), & ! tabulated extinction liquid
    z_ext_i(n_sizes, n_mdl_bnds), & ! tabulated extinction ice
    z_coa_l(n_sizes, n_mdl_bnds), & ! tabulated co-albedo liquid
    z_coa_i(n_sizes, n_mdl_bnds), & ! tabulated co-albedo ice
    z_asy_l(n_sizes, n_mdl_bnds), & ! tabulated asymmetry liquid
    z_asy_i(n_sizes, n_mdl_bnds) ! tabulated asymmetry ice

!!$  TYPE(file_info) :: optical_tbl

  INTEGER :: fileid     !< id number of netcdf file

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief sets resolution dependent parameters for cloud optics
  !
#ifdef DUMP_COEFFS
  SUBROUTINE dump
    CALL dump_source("cloud_optics", "wavenumber", wavenumber)
    CALL dump_source("cloud_optics", "wavelength", wavelength)
    CALL dump_source("cloud_optics", "re_droplet", re_droplet)
    CALL dump_source("cloud_optics", "re_crystal", re_crystal)
    CALL dump_source("cloud_optics", "z_ext_l", z_ext_l)
    CALL dump_source("cloud_optics", "z_coa_l", z_coa_l)
    CALL dump_source("cloud_optics", "z_asy_l", z_asy_l)
    CALL dump_source("cloud_optics", "z_ext_i", z_ext_i)
    CALL dump_source("cloud_optics", "z_coa_i", z_coa_i)
    CALL dump_source("cloud_optics", "z_asy_i", z_asy_i)
  END SUBROUTINE
#endif

#ifdef PSRAD_ONLY
  SUBROUTINE setup_cloud_optics

    zinhoml1 = 0.6 ! 0.7
    zinhoml2 = 0.6 ! 0.7 for nn=31 in mo_control...
    zinhomi  = 0.7
    zinpar   = 0.10_wp

#else
  SUBROUTINE setup_cloud_optics(jg)

    INTEGER, OPTIONAL, INTENT(in) :: jg

    IF (PRESENT(jg)) THEN
       zinhoml1      = echam_cld_config(jg)% cinhoml1
       zinhoml2      = echam_cld_config(jg)% cinhoml2
       zinhoml3      = echam_cld_config(jg)% cinhoml3
       zinhomi       = echam_cld_config(jg)% cinhomi
    ELSE
       zinhoml1      = 0.8_wp
       zinhoml2      = 0.4_wp
       zinhoml3      = 0.8_wp
       zinhomi       = 0.8_wp
    ENDIF
#endif

    ! Variable liquid cloud inhomogeneity is not used:
    l_variable_inhoml = .FALSE.


!!$    IF (p_parallel_io) THEN
!!$      CALL io_open ('ECHAM6_CldOptProps.nc', optical_tbl, io_read)
    CALL psrad_io_open('ECHAM6_CldOptProps.nc', fileid)
    IF (fileid == 0) THEN
      CALL finish('mo_psrad_cloud_optics/setup_cloud_optics', 'File ECHAM6_CldOptProps.nc cannot be opened')
    END IF

    CALL psrad_io_copy_double(fileid, "wavenumber", &
      (/1/), (/n_mdl_bnds/), wavenumber)
    CALL psrad_io_copy_double(fileid, "wavelength", &
      (/1/), (/n_mdl_bnds/), wavelength)
    CALL psrad_io_copy_double(fileid, "re_droplet", &
      (/1/), (/n_sizes/), re_droplet)
    CALL psrad_io_copy_double(fileid, "re_crystal", &
      (/1/), (/n_sizes/), re_crystal)

    CALL psrad_io_copy_double(fileid, "extinction_per_mass_droplet", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_ext_l)
    CALL psrad_io_copy_double(fileid, "co_albedo_droplet", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_coa_l)
    CALL psrad_io_copy_double(fileid, "asymmetry_factor_droplet", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_asy_l)

    CALL psrad_io_copy_double(fileid, "extinction_per_mass_crystal", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_ext_i)
    CALL psrad_io_copy_double(fileid, "co_albedo_crystal", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_coa_i)
    CALL psrad_io_copy_double(fileid, "asymmetry_factor_crystal", &
      (/1,1/), (/n_sizes,n_mdl_bnds/), z_asy_i)

    CALL psrad_io_close(fileid)
#ifdef DUMP_COEFFS
    call dump()
#endif

    reimin = MINVAL(re_crystal)
    reimax = MAXVAL(re_crystal)
    del_rei= (re_crystal(2) - re_crystal(1))
    relmin = MINVAL(re_droplet)
    relmax = MAXVAL(re_droplet)
    del_rel= (re_droplet(2) - re_droplet(1))

  END SUBROUTINE setup_cloud_optics

  !! @brief Calculates cloud optical from cloud physical properties
  !! @remarks  
  !!   Currently this model assumes four bands in the SW and maps these to the
  !!   six-band model of ECHAM5.
  SUBROUTINE cloud_optics(laglac, laland, kproma, kbdim, klev, ktype, &
    icldlyr, zlwp, ziwp, zlwc, ziwc, zcdnc, tau_lw, tau_sw, omg, asy, &
    re_droplets2d, re_crystals2d)

    USE mo_psrad_general, ONLY : nbndsw, nbndlw
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, &
      ktype(KBDIM), & ! type of convection
      icldlyr(KBDIM,klev)
    LOGICAL, DIMENSION(KBDIM), INTENT(IN) :: & 
      laglac, & ! logical for glacier points
      laland ! logical for land points
    REAL(wp), DIMENSION(KBDIM,klev), INTENT (IN)  :: &
      zlwp, & ! liquid water path
      ziwp, & ! ice water path
      zcdnc, & ! cloud drop number concentration
      zlwc, & ! liquid water content
      ziwc ! ice water content

    REAL(wp), INTENT (OUT) :: &
      tau_lw(KBDIM,klev,nbndlw), & ! LW optical depth
      tau_sw(KBDIM,klev,nbndsw), & ! SW optical depth
      omg(KBDIM,klev,nbndsw), & ! cloud single scattering albedo
      asy(KBDIM,klev,nbndsw), & ! cloud asymmetry factor
      re_droplets2d(KBDIM,klev), & ! effective radius of liquid distribution
      re_crystals2d(KBDIM,klev) ! effective radius of ice distribution

    INTEGER  :: iband, ii, jk, jl, ml1, ml2, mi1, mi2
    REAL(wp) :: ztol, ztoi, zol, zoi, zgl, zgi, wl1, wl2, wi1, wi2, &
      zfact, zmsald, zmsaid, &
      zkap(KBDIM), & ! breath parameter for scaling effective radius
      zlwpt(KBDIM), & ! liquid water path
      zinhoml(KBDIM), & ! cloud inhomogeneity factor (liquid)
      re_droplets, & ! effective radius of liquid distribution
      re_crystals, & ! effective radius of ice distribution
      zscratch, & ! effective radius of ice distribution
      ztau(KBDIM,klev,n_mdl_bnds), & ! SW optical depth
      zomg(KBDIM,klev,n_mdl_bnds), & ! cloud single scattering albedo
      zasy(KBDIM,klev,n_mdl_bnds) ! cloud asymmetry factor

    IF (relmin < 1.5 .OR. relmin > 2.5 ) THEN
      CALL finish('Apparently unsuccessful loading of optical tables')
    END IF

    ! Basic cloud properties
    IF (l_variable_inhoml) THEN
      zlwpt(1:kproma) = 0.0_wp
      DO jk = 1, klev
        DO jl = 1, kproma
          zlwpt(jl) = zlwpt(jl)+zlwp(jl,jk)
        END DO
      END DO
      WHERE (zlwpt(1:kproma) > 1.0_wp) 
        zinhoml(1:kproma) = zlwpt(1:kproma)**(-zinpar)
      ELSEWHERE
        zinhoml(1:kproma) = 1.0_wp
      END WHERE
    ELSE
      DO jl = 1, kproma
        IF(ktype(jl) .EQ. 0) THEN          !no convection; ktype=0
          zinhoml(jl) = zinhoml1
#ifndef PSRAD_ONLY
        ! BUG? Bjorn's version did not contain these lines!
        ELSE IF(ktype(jl) .NE. 4) THEN !convection; ktype=1,2,3
          zinhoml(jl) = zinhoml3
#endif
        ELSE !shallow convection and clwprat>0.; ktype=4
          zinhoml(jl) = zinhoml2
        END IF
      END DO
    END IF

    WHERE (laland(1:kproma).AND.(.NOT.laglac(1:kproma))) 
      zkap(1:kproma)=zkap_cont ! continental breadth factor
    ELSEWHERE
      zkap(1:kproma)=zkap_mrtm ! maritime breadth factor 
    END WHERE

    ! Cloud Optical Properties by interpolating tables in effective radius
    zfact = 1.0e6_wp*(3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp) 
    DO jk=1,klev
      DO jl=1,kproma
        IF (icldlyr(jl,jk)==1 .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN
          
          re_crystals = MAX(reimin,&
            MIN(reimax,83.8_wp*ziwc(jl,jk)**0.216_wp))
          re_droplets = MAX(relmin,&
            MIN(relmax,zfact*zkap(jl)*(zlwc(jl,jk) / &
              zcdnc(jl,jk))**(1.0_wp/3.0_wp)))
          re_crystals2d(jl,jk) = re_crystals
          re_droplets2d(jl,jk) = re_droplets

          ml1 = MAX(1,&
            MIN(n_sizes-1,FLOOR(1.0_wp+(re_droplets-relmin)/del_rel)))
          ml2 = ml1 + 1
          wl1 = 1.0_wp - (re_droplets - (relmin + del_rel* (ml1-1)) )/del_rel
          wl2 = 1.0_wp - wl1

          mi1 = MAX(1,MIN(n_sizes-1,FLOOR(1.0_wp+(re_crystals-reimin)/del_rei)))
          mi2 = mi1 + 1
          wi1 = 1.0_wp - (re_crystals - (reimin + del_rei * (mi1-1)) )/del_rei
          wi2 = 1.0_wp - wi1

          DO iband = 1,n_mdl_bnds
            ztol = zlwp(jl,jk) * &
              (wl1 * z_ext_l(ml1,iband) + wl2 * z_ext_l(ml2,iband))
            ztoi = ziwp(jl,jk) * &
              (wi1 * z_ext_i(mi1,iband) + wi2 * z_ext_i(mi2,iband))
            zol = 1.0_wp - &
              (wl1 * z_coa_l(ml1,iband) + wl2 * z_coa_l(ml2,iband))
            zoi = 1.0_wp - &
              (wi1 * z_coa_i(mi1,iband) + wi2 * z_coa_i(mi2,iband))
            zgl = wl1 * z_asy_l(ml1,iband) + wl2 * z_asy_l(ml2,iband)
            zgi = wi1 * z_asy_i(mi1,iband) + wi2 * z_asy_i(mi2,iband)

            zscratch = (ztol * zol + ztoi * zoi)
            ztau(jl,jk,iband) = ztol * zinhoml(jl) + ztoi * zinhomi
            zomg(jl,jk,iband) = zscratch / (ztol + ztoi)
            zasy(jl,jk,iband) = (ztol * zol * zgl + ztoi * zoi * zgi) / zscratch
          END DO

          ! overwrite Kinne Optics with old Cloud Optics for LW Only
          zmsald = 0.025520637_wp + &
            0.2854650784_wp * EXP(-0.088968393014_wp * re_droplets)
          DO ii = 1,16
            iband = ii
            IF (ii == 16) iband=30 
            zmsaid = (rebcuh(ii)+rebcug(ii)/re_crystals) ! WRONG PARENTHESIS?
            ztau(jl,jk,iband) = zmsald * zlwp(jl,jk) * zinhoml(jl) + &
              zmsaid * ziwp(jl,jk) * zinhomi
          END DO
        ELSE
          ztau(jl,jk,:) = 0.0_wp
          zomg(jl,jk,:) = 1.0_wp
          zasy(jl,jk,:) = 0.0_wp
          re_crystals2d(jl,jk) = 0.0_wp
          re_droplets2d(jl,jk) = 0.0_wp
        END IF
        tau_lw(jl,jk,1:nbndlw-1) = ztau(jl,jk,1:15)
        tau_lw(jl,jk,nbndlw) = ztau(jl,jk,30)
        tau_sw(jl,jk,1:nbndsw) = ztau(jl,jk,16:29)
        omg(jl,jk,1:nbndsw) = zomg(jl,jk,16:29)
        asy(jl,jk,1:nbndsw) = zasy(jl,jk,16:29)
      END DO
    END DO

  END SUBROUTINE cloud_optics

END MODULE mo_psrad_cloud_optics
