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

      ! fraction of solid land in the grid box, i.e. land without lakes if lakes are used
      ! see mo_echam_phy_init or input data set for details.
      !
      IF (ilnd.LE.nsfc_type) THEN
         zfrl(jc) = MAX(MIN(            &
              &     field%lsmask(jc,jb) &
              &     ,1._wp),0._wp)
      ELSE
         zfrl(jc) = 0._wp
      END IF

      ! fraction of open water in the grid box, for sea and lakes
      !
      IF (iwtr.LE.nsfc_type) THEN
         zfrw(jc) = MAX(MIN(                                                                          &
              &      (1._wp-field%lsmask(jc,jb)-field%alake(jc,jb))*(1._wp-field%seaice(jc,jb))       & ! ocean
              &     +field%alake(jc,jb)                            *(1._wp-field%lake_ice_frc(jc,jb)) & ! lakes
              &     ,1._wp),0._wp)
         !
         ! security for water temperature with changing ice mask
         ! (over lakes; over ocean this is not an issue since, over ocean,
         ! ts_tile(iwtr) is overwritten again with SST from ocean in 
         ! coupling interface after update_surface)
         IF (zfrw(jc) > 0._wp .AND. field%ts_tile(jc,jb,iwtr) == cdimissval) THEN
           ! lake was completely frozen in previous time step but only partially 
           ! frozen in current time step
           field%ts_tile(jc,jb,iwtr) = tmelt
         END IF

      ELSE
         zfrw(jc) = 0._wp
      END IF

      ! fraction of ice covered water in the grid box, for sea and lakes
      !
      IF (iice.LE.nsfc_type) THEN
         zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
         !
         ! security for ice temperature with changing ice mask
         IF(zfri(jc) > 0._wp .AND. field%ts_tile(jc,jb,iice) == cdimissval ) THEN
            field%ts_tile(jc,jb,iice)  = tmelt + Tf    ! = 271.35 K
         END IF
      ELSE
         zfri(jc) = 0._wp
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
