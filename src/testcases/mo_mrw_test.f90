!>
!!  This module contains parameters, initialization subroutines and.
!!
!!  This module contains parameters, initialization subroutines and
!!  functions to be used in the mountain induced Rossby wave train
!!  case of the hydrostatic dynamical core.
!!
!! @par Revision History
!!  Original version by Almut Gassmann, MPI-M (2008-04-24)
!!  Modified by Almut Gassmann, MPI-M (2008-10-09)
!!  - blocking
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_mrw_test
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
  USE mo_physical_constants,  ONLY: cpd, grav, rd, re, omega
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_configuration,  ONLY: nproma

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  REAL (wp), PUBLIC :: mountctr_lon_deg, mountctr_lat_deg, mountctr_height, &
                       mount_half_width, mount_u0

  PUBLIC :: init_hydro_state_prog_mrossby, init_hydro_state_prog_mrossby2


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
 !! @par
 !!  arguments
 !!
 SUBROUTINE init_hydro_state_prog_mrossby(pt_patch, pt_prog, pt_ext_data)


  TYPE(t_patch), INTENT(INOUT)    :: pt_patch
  TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
  TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

  REAL (wp) :: bruntvaissq, temp_0, pres_sp, u_0, kappa, &
               zcoslat, zhelp1, zlon_mc, zlat_mc, zlon, zlat, zr, zexp

  INTEGER   :: je, jb, jc, nlen, nblks_e, nblks_c, npromz_e, npromz_c
!--------------------------------------------------------------------
!
!---initialize the prognostic variables
!
  temp_0      = 288.0_wp
  bruntvaissq = grav*grav/cpd/temp_0
  pres_sp     = 93000.0_wp
  u_0         = mount_u0
  kappa       = rd/cpd
  zhelp1      = bruntvaissq/grav/grav/kappa
  zlon_mc     = mountctr_lon_deg*pi/180._wp
  zlat_mc     = mountctr_lat_deg*pi/180._wp

  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c
  nblks_e   = pt_patch%nblks_int_e
  npromz_e  = pt_patch%npromz_int_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,zlat,zlon,zr,zexp,zcoslat)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! topography
        DO jc = 1, nlen

           zlat = pt_patch%cells%center(jc,jb)%lat
           zlon = pt_patch%cells%center(jc,jb)%lon

           zr = SIN(zlat_mc)*SIN(zlat)+COS(zlat_mc)*COS(zlat)*COS(zlon-zlon_mc)
           zexp = re*ACOS(zr)/mount_half_width

           pt_ext_data%atm%topography_c(jc,jb) = &
                        mountctr_height*EXP( - zexp*zexp )

           ! surface pressure
           zcoslat  = COS(zlat)
           pt_prog%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1 * ( u_0 *&
                              ( 0.5_wp*u_0 + re*omega) * zcoslat*zcoslat - &
                              pt_ext_data%atm%topography_c(jc,jb)*grav))
        ENDDO
     ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,je,nlen,zcoslat)
     ! normal wind
     DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO je = 1, nlen
           zcoslat  = COS(pt_patch%edges%center(je,jb)%lat)
           pt_prog%vn(je,:,jb) = u_0*zcoslat*&
                                     pt_patch%edges%primal_normal(je,jb)%v1
        ENDDO
     ENDDO
!$OMP END DO

     ! temperature
!$OMP WORKSHARE
     pt_prog%temp(:,:,:) = temp_0
!$OMP END WORKSHARE
!$OMP END PARALLEL


  END SUBROUTINE init_hydro_state_prog_mrossby


!-------------------------------------------------------------------------
!

 !>
 !!               Same as init_hydro_state_prog_mrossby,.
 !!
 !!               Same as init_hydro_state_prog_mrossby,
 !! but using mountain with small-scale corrugations
 !!
 !! @par Revision History
 !!  Created by Guenther Zaengl (2008-12-08)
 !! @par
 !!  arguments
 !!
 SUBROUTINE init_hydro_state_prog_mrossby2(pt_patch, pt_prog, pt_ext_data)


  TYPE(t_patch), INTENT(INOUT)    :: pt_patch
  TYPE(t_hydro_atm_prog), INTENT(INOUT) :: pt_prog
  TYPE(t_external_data), INTENT(INOUT) :: pt_ext_data !< external data

  REAL (wp) :: bruntvaissq, temp_0, pres_sp, u_0, kappa, &
               zcoslat, zhelp1, zlon_mc, zlat_mc, zlon, zlat, zr, zexp, &
               topo_aux(nproma,7)

  INTEGER   :: je, jb, jc, nlen, nblks_e, nblks_c, npromz_e, npromz_c,  &
               iiv, ibv, iv
  INTEGER   :: nlev              !< number of full levels
!--------------------------------------------------------------------
!
!---initialize the prognostic variables
!
  temp_0      = 288.0_wp
  bruntvaissq = grav*grav/cpd/temp_0
  pres_sp     = 93000.0_wp
  u_0         = mount_u0
  kappa       = rd/cpd
  zhelp1      = bruntvaissq/grav/grav/kappa
  zlon_mc     = mountctr_lon_deg*pi/180._wp
  zlat_mc     = mountctr_lat_deg*pi/180._wp

  nblks_c   = pt_patch%nblks_int_c
  npromz_c  = pt_patch%npromz_int_c
  nblks_e   = pt_patch%nblks_int_e
  npromz_e  = pt_patch%npromz_int_e

  nlev = pt_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,nlen,iv,zlat,zlon,zr,zexp,zcoslat,iiv,ibv,topo_aux)
     DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! topography
        DO jc = 1, nlen

           ! To obtain a topography height representative for the cell
           ! an average over the cell center and the vertices is taken

           zlat = pt_patch%cells%center(jc,jb)%lat
           zlon = pt_patch%cells%center(jc,jb)%lon

           zr = SIN(zlat_mc)*SIN(zlat)+COS(zlat_mc)*COS(zlat)*COS(zlon-zlon_mc)
           zexp = re*ACOS(zr)/mount_half_width

           topo_aux(jc,1) = mountctr_height*EXP( - zexp*zexp )*&
                            0.5_wp*(1._wp+COS(pi*zexp*2._wp))
        ENDDO


        DO iv = 1, pt_patch%cell_type

           DO jc = 1, nlen

             IF (iv > pt_patch%cells%num_edges(jc,jb)) THEN
               topo_aux(jc,pt_patch%cell_type+1) = topo_aux(jc,1)
               CYCLE
             ENDIF

             iiv = pt_patch%cells%vertex_idx(jc,jb,iv)
             ibv = pt_patch%cells%vertex_blk(jc,jb,iv)
             zlat = pt_patch%verts%vertex(iiv,ibv)%lat
             zlon = pt_patch%verts%vertex(iiv,ibv)%lon

             zr = SIN(zlat_mc)*SIN(zlat)+COS(zlat_mc)*COS(zlat)*COS(zlon-zlon_mc)
             zexp = re*ACOS(zr)/mount_half_width

             topo_aux(jc,iv+1) = mountctr_height*EXP( - zexp*zexp )*&
                                 0.5_wp*(1._wp+COS(pi*zexp*2._wp))

           ENDDO
        ENDDO

        DO jc = 1, nlen

           pt_ext_data%atm%topography_c(jc,jb) = 0.5_wp*topo_aux(jc,1) + &
             0.5_wp/REAL(pt_patch%cell_type,wp)*(SUM(topo_aux(jc,2:pt_patch%cell_type+1)))
           ! surface pressure

           zlat = pt_patch%cells%center(jc,jb)%lat
           zlon = pt_patch%cells%center(jc,jb)%lon

           zcoslat  = COS(zlat)
           pt_prog%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1 * ( u_0 *&
                              ( 0.5_wp*u_0 + re*omega) * zcoslat*zcoslat - &
                              pt_ext_data%atm%topography_c(jc,jb)*grav))
        ENDDO
     ENDDO

!$OMP END DO
!$OMP DO PRIVATE(jb,je,nlen,zcoslat)

     ! normal wind
     DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
           nlen = nproma
        ELSE
           nlen = npromz_e
        ENDIF
        DO je = 1, nlen
           zcoslat  = COS(pt_patch%edges%center(je,jb)%lat)
           pt_prog%vn(je,1:nlev,jb) = u_0*zcoslat*&
                                      pt_patch%edges%primal_normal(je,jb)%v1
        ENDDO
     ENDDO
!$OMP END DO

     ! temperature
!$OMP WORKSHARE
     pt_prog%temp(:,:,:) = temp_0
!$OMP END WORKSHARE
!$OMP END PARALLEL


  END SUBROUTINE init_hydro_state_prog_mrossby2

END MODULE mo_mrw_test
