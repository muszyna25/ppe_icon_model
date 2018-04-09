!>
!! @brief Subroutine echam_phy_diag contains small diagnostic routines,
!!  which are executed on a single block of cells.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_echam_phy_diag

  USE mo_kind                ,ONLY: wp

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend,      &
    &                               cdimissval

  USE mo_physical_constants  ,ONLY: cpd, cpv, cvd, cvv, Tf, tmelt
  USE mo_run_config          ,ONLY: iqv
  USE mo_echam_cld_config    ,ONLY: echam_cld_config
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice, ilnd

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: surface_fractions, &
       &    droplet_number,    &
       &    cpair_cvair_qconv, &
       &    initialize,        &
       &    finalize

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE surface_fractions(jg,jb,jcs,jce ,&
       &                       nproma,nlev   )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field), POINTER    :: field
    REAL(wp)                            :: zfrw (nproma) !< cell area fraction of open water
    REAL(wp)                            :: zfri (nproma) !< cell area fraction of ice covered water
    REAL(wp)                            :: zfrl (nproma) !< cell area fraction of land
    INTEGER                             :: jc
 
    field => prm_field(jg)
 
    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box.
      ! lsmask: land-sea mask, depends on input data, either:
      ! fractional, including or excluding lakes in the land part or
      ! non-fractional, each grid cell is either land, sea, or sea-ice.
      ! See mo_echam_phy_init or input data set for details.

      zfrl(jc) = MAX(field% lsmask(jc,jb),0._wp)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      IF (iwtr.LE.nsfc_type) THEN
         zfrw(jc) = MAX(1._wp-zfrl(jc),0._wp)*MAX(1._wp-(field%seaice(jc,jb)+field%lake_ice_frc(jc,jb)),0._wp)
      ELSE
         zfrw(jc) = 0._wp
      END IF

      ! fraction of sea ice in the grid box
      zfri(jc) = MAX(1._wp-zfrl(jc)-zfrw(jc),0._wp)
      !
      IF (iice.LE.nsfc_type) THEN
         ! security for ice temperature with changing ice mask
         IF(zfri(jc) > 0._wp .AND. field%ts_tile(jc,jb,iice) == cdimissval ) THEN
            field% ts_tile(jc,jb,iice)  = tmelt + Tf    ! = 271.35 K
         END IF
      END IF

    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) field%frac_tile(jcs:jce,jb,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iice) = zfri(jcs:jce)

    NULLIFY(field)

  END SUBROUTINE surface_fractions
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE droplet_number(jg,jb,jcs,jce ,&
       &                    nproma,nlev   )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field), POINTER    :: field
    INTEGER                             :: jc, jk
    LOGICAL                             :: lland(nproma)
    LOGICAL                             :: lglac(nproma)
    REAL(wp)                            :: zprat, zn1, zn2, zcdnc

    ! Shortcuts to components of echam_cld_config
    !
    REAL(wp), POINTER :: cn1lnd, cn2lnd, cn1sea, cn2sea
    !
    cn1lnd => echam_cld_config(jg)% cn1lnd
    cn2lnd => echam_cld_config(jg)% cn2lnd
    cn1sea => echam_cld_config(jg)% cn1sea
    cn2sea => echam_cld_config(jg)% cn2sea

    field => prm_field(jg)
    
    DO jc=jcs,jce
      lland(jc) = field%sftlf (jc,jb) > 0._wp
      lglac(jc) = field%sftgif(jc,jb) > 0._wp
    END DO

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**2

        IF (lland(jc).AND.(.NOT.lglac(jc))) THEN
          zn1= cn1lnd
          zn2= cn2lnd
        ELSE
          zn1= cn1sea
          zn2= cn2sea
        END IF
        IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        END IF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO

    NULLIFY(cn1lnd, cn2lnd, cn1sea, cn2sea)
    NULLIFY(field)

  END SUBROUTINE droplet_number
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE cpair_cvair_qconv(jg,jb,jcs,jce ,&
       &                       nproma,nlev   )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field), POINTER    :: field
    INTEGER                             :: jc, jk

    field => prm_field(jg)
    
    DO jk = 1,nlev
       DO jc = jcs,jce
          !
          field%cpair(jc,jk,jb) = cpd+(cpv-cpd)*field%qtrc(jc,jk,jb,iqv)
          field%cvair(jc,jk,jb) = cvd+(cvv-cvd)*field%qtrc(jc,jk,jb,iqv)
          !
          field%qconv(jc,jk,jb) = 1._wp/(field%mair(jc,jk,jb)*field%cpair(jc,jk,jb))
          !
       END DO
    END DO
    
    NULLIFY(field)

  END SUBROUTINE cpair_cvair_qconv
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE initialize(jg,jb,jcs,jce ,&
       &                nproma,nlev   )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field), POINTER    :: field
    INTEGER                             :: jc, jk

    field => prm_field(jg)
    
    DO jk = 1,nlev
       DO jc = jcs,jce
          !
          ! heating
          field% q_phy(jc,jk,jb) = 0._wp
          !
       END DO
    END DO

    DO jc = jcs,jce
       !
       ! vertical integral of heating
       field% q_phy_vi(jc,jb) = 0._wp
       !
    END DO

    NULLIFY(field)

  END SUBROUTINE initialize
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE finalize(jg,jb,jcs,jce ,&
       &              nproma,nlev   )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field), POINTER    :: field
    TYPE(t_echam_phy_tend) , POINTER    :: tend
    INTEGER                             :: jc, jk

    field => prm_field(jg)
    tend  => prm_tend (jg)
    
    ! precipitation flux from all processes
    !
    DO jc = jcs,jce
       field% pr(jc,jb) =  field% rsfl(jc,jb) & ! rain large scale
            &             +field% ssfl(jc,jb) & ! snow large scale
            &             +field% rsfc(jc,jb) & ! rain convection
            &             +field% ssfc(jc,jb)   ! snow convection
    END DO
 
    DO jk = 1,nlev
       DO jc = jcs,jce
          !
          ! vertical integral of heating
          field% q_phy_vi(jc,jb) = field% q_phy_vi(jc,jb) + field% q_phy(jc,jk,jb)
          !
          ! now convert the temperature tendency from physics, as computed for constant pressure conditions,
          ! to constant volume conditions, as needed for the coupling to the dynamics
          tend% ta_phy(jc,jk,jb) = tend% ta_phy(jc,jk,jb) * field% cpair(jc,jk,jb) / field% cvair(jc,jk,jb)
          !
       END DO
    END DO

    NULLIFY(field)
    NULLIFY(tend )

  END SUBROUTINE finalize
  !---------------------------------------------------------------------

END MODULE mo_echam_phy_diag
