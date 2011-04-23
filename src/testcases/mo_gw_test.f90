!>
!!  This module contains parameters, initialization subroutines and.
!!
!!  This module contains parameters, initialization subroutines and
!!  functions to be used in the 3D gravity wave test case of the
!!  hydrostatic dynamical core.
!!
!! @par Revision History
!!  Original version by Hui Wan, MPI-M (2008-01-28)
!!  Modification by Almut Gassmann, MPI-M (2008-04-25)
!!   - changings according to NCAR Workshop experiments
!!  Modification by Almut Gassmann, MPI-M (2008-10-09)
!!   - blocking and other restructurings
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
!!
MODULE mo_gw_test
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2005
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data,            ONLY: t_external_data
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog
  USE mo_physical_constants,  ONLY: grav, rd, cpd, re, omega
  USE mo_math_constants,      ONLY: pi
  USE mo_vertical_coord_table,ONLY: vct_a,vct_b
  USE mo_run_nml,             ONLY: lcorio, nproma

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  REAL(wp) :: gw_brunt_vais  ! Brunt Vaisala frequency (1/s)
  REAL(wp) :: gw_u0          ! mean zonal wind (m/s)
  REAL(wp) :: gw_lon_deg, gw_lat_deg

  PUBLIC   :: gw_brunt_vais, gw_u0, gw_lon_deg, gw_lat_deg
  PUBLIC   :: init_hydro_state_prog_gwtest

 CONTAINS

!--------------------------------------------------------------------

!-------------------------------------------------------------------------
!

 !>
 !!               Initialization of prognostic state vector.
 !!
 !!
 !! @par Revision History
 !!  Original version  by Hui Wan, MPI-M (2008-01-28)
 !!  Renewed by Almut Gassmann, MPI-M (2008-04-25)
 !! @par
 !!  arguments
 !!
 SUBROUTINE init_hydro_state_prog_gwtest(pt_patch, pt_ext_data, pt_prog)


     IMPLICIT NONE

     TYPE(t_patch),            INTENT(INOUT) :: pt_patch
     TYPE(t_external_data),    INTENT(INOUT) :: pt_ext_data
     TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog

     REAL (wp) :: z_ps_0, z_t_0, z_sot, z_kappa, z_l_z, z_nsqog, z_pres, &
                  z_latc, z_lonc, z_lat, z_lon, z_rr, z_hshape, z_d_the, &
                  z_pk, z_pk1, z_height, z_thetab, z_exner, z_vshape, z_r, &
                  z_omega

     INTEGER   :: jb, je, jk, jc, nlen, npromz_e, nblks_e, npromz_c, nblks_c
     INTEGER   :: nlev              !< number of full levels

!--------------------------------------------------------------------
!
!---initialize the prognostic variables
!
  IF (lcorio) THEN
     z_omega = omega
  ELSE
     z_omega = 0.0_wp
  ENDIF
  z_ps_0  = 100000.0_wp
  z_t_0   = 300.0_wp
  z_l_z   = 20000.0_wp
  z_latc  = gw_lat_deg*pi/180._wp
  z_d_the = 10.0_wp
  z_lonc  = gw_lon_deg*pi/180._wp
  z_sot   = (grav/gw_brunt_vais)**2/cpd/z_t_0
  z_kappa = rd/cpd
  z_nsqog = gw_brunt_vais**2/grav
  z_rr    = re/3.0_wp

  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c
  nblks_e   = pt_patch%nblks_int_e
  npromz_e  = pt_patch%npromz_int_e

  ! number of vertical levels
  nlev = pt_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,z_lat)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jc = 1, nlen

           ! topography
           pt_ext_data%atm%topography_c(jc,jb) = 0.0_wp

           z_lat   = pt_patch%cells%center(jc,jb)%lat

           ! surface pressure (ncells_all,nblks)
           pt_prog%pres_sfc(jc,jb) = z_ps_0 * EXP ( &
                  -re*(gw_brunt_vais**2)*gw_u0*0.5_wp/(grav**2)/z_kappa*&
                  ( gw_u0/re +2.0_wp*z_omega) * ((SIN(z_lat))**2) )
        ENDDO
     ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,je,nlen,z_lat)
     DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
           DO je = 1, nlen

              z_lat   = pt_patch%edges%center(je,jb)%lat

              ! normal wind (nedges_all,nvertical_layers,nblks)
              pt_prog%vn(je,jk,jb) = gw_u0 *  COS(z_lat) *&
                          pt_patch%edges%primal_normal(je,jb)%v1

           ENDDO
        ENDDO
     ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen,z_pk,z_pk1,z_pres,z_exner,z_height,z_thetab,&
!$OMP            z_vshape,z_lat,z_lon,z_r,z_hshape)

     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF
        DO jk = 1,nlev

           ! Full level pressure
           z_pk   = vct_a(jk  ) + vct_b(jk  )*z_ps_0
           z_pk1  = vct_a(jk+1) + vct_b(jk+1)*z_ps_0
           z_pres = EXP((z_pk1*LOG(z_pk1)-z_pk*LOG(z_pk))/(z_pk1-z_pk)-1.0_wp)

           z_exner  = (z_pres/z_ps_0)**z_kappa
           z_height = -LOG((z_exner-1.0_wp)/z_sot+1.0_wp)/z_nsqog
           z_thetab = z_t_0/((z_exner-1.0_wp)/z_sot+1.0_wp)
           z_vshape = SIN(2.0_wp*pi*z_height/z_l_z)

           DO jc = 1,nlen

              z_lat   = pt_patch%cells%center(jc,jb)%lat
              z_lon   = pt_patch%cells%center(jc,jb)%lon

              ! Horizontal shape function
              z_r     = re*ACOS(SIN(z_latc)*SIN(z_lat)+&
                                COS(z_latc)*COS(z_lat)*COS(z_lon-z_lonc))
              IF ( z_r < z_rr ) THEN
                 z_hshape = 0.5_wp*(1.0_wp+COS(pi*z_r/z_rr))
              ELSE
                 z_hshape = 0.0_wp
              ENDIF

              ! temperature (ncells_all,nvertical_layers,nblks)
              pt_prog%temp(jc,jk,jb) = &
                          z_exner*(z_thetab+z_d_the*z_hshape*z_vshape)

           ENDDO
        ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE init_hydro_state_prog_gwtest

END MODULE mo_gw_test
