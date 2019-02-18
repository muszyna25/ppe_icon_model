!>
!! Geometric computations which are specific to the transport algorithm.
!!
!! Module contains procedures for geometric computations which are
!! specific to the horizontal transport schemes.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2013-10-30)
!! Modification by Daniel Reinert, DWD (2013-10-30)
!! - moved divide_flux_area from mo_advection_traj to this module
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_geometry

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int, max_char_length
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_math_types,          ONLY: t_line, t_geographical_coordinates
  USE mo_math_utilities,      ONLY: lintersect, line_intersect 
  USE mo_advection_utils,     ONLY: t_list2D
  USE mo_fortran_tools,       ONLY: copy


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: divide_flux_area
  PUBLIC :: divide_flux_area_list

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Divide flux area
  !!
  !! Flux area (aka. departure region) is subdivided according to its overlap
  !! with the underlying grid.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-04-03)
  !!
  SUBROUTINE divide_flux_area(p_patch, p_int, p_vn, p_vt,            &
    &                         dreg_patch0, dreg_patch1, dreg_patch2, &
    &                         patch1_cell_idx, patch1_cell_blk,      &
    &                         patch2_cell_idx, patch2_cell_blk,      &
    &                         opt_rlstart, opt_rlend, opt_slev,      &
    &                         opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    REAL(vp), INTENT(INOUT) ::    &  !< patch 0 of subdivided departure region
      & dreg_patch0(:,:,:,:,:)       !< coordinates

    REAL(vp), INTENT(OUT) ::    &  !< patch 0,1,2 of subdivided departure region
      & dreg_patch1(:,:,:,:,:), &  !< coordinates
      & dreg_patch2(:,:,:,:,:)     !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of underlying cell
      & patch1_cell_idx(:,:,:), &  !< dim: (nproma,nlev,ptr_p%nblks_e)
      & patch1_cell_blk(:,:,:), &
      & patch2_cell_idx(:,:,:), &
      & patch2_cell_blk(:,:,:)


    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) ::       &   !< coordinates of arrival points. The origin
      &  arrival_pts(nproma,2,2,p_patch%nlev,p_patch%nblks_e)
                          !< of the coordinate system is at the circumcenter of
                          !< the upwind cell. Unit vectors point to local East
                          !< and North. (geographical coordinates)

    REAL(wp) ::    &   !< coordinates of departure points. The origin
      &  depart_pts(nproma,2,2,p_patch%nlev,p_patch%nblks_e)
                       !< of the coordinate system is at the circumcenter of
                       !< the upwind cell. Unit vectors point to local East
                       !< and North. (geographical coordinates)

    TYPE(t_line) ::                  & !< departure-line segment
      &  fl_line(nproma,p_patch%nlev)

    TYPE(t_line) ::                  & !< departure area edges
      &  fl_e1(nproma,p_patch%nlev), & !< edge 1
      &  fl_e2(nproma,p_patch%nlev)    !< edge 2

    TYPE(t_line) ::                  & !< triangle edge
      &  tri_line1(nproma,p_patch%nlev), &
      &  tri_line2(nproma,p_patch%nlev)

    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of vertex3
      &  ptr_v3(:,:,:)
    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of
      &  ptr_bfcc(:,:,:,:)                         !< of butterfly cell centers


    REAL(wp) :: ps1(2),              & !< coordinates of intersection
      &         ps2(2)                 !< points S1, S2

    REAL(wp) :: pi1(2),              & !< coordinates of intersection
      &         pi2(2)                 !< points I1, I2

    REAL(wp) :: bf_cc(2,2)             !< coordinates of butterfly cell centers

    INTEGER :: je, jk, jb, jl          !< loop index of edge, vert level, block, lists
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: slev, elev              !< vertical start and end level

    LOGICAL :: lintersect_line1, lintersect_line2
    LOGICAL :: lintersect_e2_line1, lintersect_e1_line2
    LOGICAL :: lvn_pos, lvn_sys_pos
    INTEGER :: icnt_c1, icnt_c2p, icnt_c3p, icnt_c2m, icnt_c3m
    INTEGER :: icnt_rem, icnt_err, icnt_vn0

    INTEGER ::           &         !< je index list
      &  idxlist_c1 (nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_c2p(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_c3p(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_c2m(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_c3m(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_rem(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_vn0(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  idxlist_err(nproma*p_patch%nlev,p_patch%nblks_e)


    INTEGER ::           &         !< jk index list
      &  levlist_c1 (nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_c2p(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_c3p(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_c2m(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_c3m(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_rem(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_vn0(nproma*p_patch%nlev,p_patch%nblks_e), &
      &  levlist_err(nproma*p_patch%nlev,p_patch%nblks_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_traj: divide_flux_area'

  !-------------------------------------------------------------------------


    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF


    ! number of child domains
    i_nchdom   = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    ! pointer to coordinates of vertex3 (i.e. vertex which belongs to
    ! the upwind cell but does not belong to the current edge.)
    !
    ptr_v3   => p_int%pos_on_tplane_c_edge(:,:,:,3)

    ! pointer to coordinates of butterfly cell centers
    !
    ptr_bfcc => p_int%pos_on_tplane_c_edge(:,:,:,4:5)

!$OMP PARALLEL
    ! get arrival and departure points. Note that the indices of the departure
    ! points have to be switched so that departure point 1 belongs to arrival
    ! point one and departure point 2 to arrival point 2.
    CALL copy(dreg_patch0(:,1:2,1:2,:,:),    arrival_pts(:,1:2,1:2,:,:))
    CALL copy(dreg_patch0(:,4:3:-1,1:2,:,:), depart_pts (:,1:2,1:2,:,:))
!$OMP BARRIER

!$OMP DO PRIVATE(jb,jk,je,jl,i_startidx,i_endidx,lvn_pos,fl_line,tri_line1, &
!$OMP            tri_line2,fl_e1,fl_e2,lintersect_line1,lintersect_line2,   &
!$OMP            lintersect_e2_line1,lintersect_e1_line2,icnt_c1,icnt_c2p,  &
!$OMP            icnt_c2m,icnt_rem,icnt_c3p,icnt_c3m,icnt_vn0,icnt_err,     &
!$OMP            lvn_sys_pos,ps1,ps2,pi1,pi2,bf_cc)
    DO jb = i_startblk, i_endblk

     CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

      ! Reset counter
      icnt_c1  = 0
      icnt_c2p = 0
      icnt_c2m = 0
      icnt_c3p = 0
      icnt_c3m = 0
      icnt_rem = 0
      icnt_err = 0
      icnt_vn0 = 0

      DO jk = slev, elev
         DO je = i_startidx, i_endidx

          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          !
          ! get flux area departure-line segment
          !
          fl_line(je,jk)%p1%lon = depart_pts(je,1,1,jk,jb)
          fl_line(je,jk)%p1%lat = depart_pts(je,1,2,jk,jb)
          fl_line(je,jk)%p2%lon = depart_pts(je,2,1,jk,jb)
          fl_line(je,jk)%p2%lat = depart_pts(je,2,2,jk,jb)

          ! get triangle edge 1 (A1V3)
          !
          tri_line1(je,jk)%p1%lon = arrival_pts(je,1,1,jk,jb)
          tri_line1(je,jk)%p1%lat = arrival_pts(je,1,2,jk,jb)
          tri_line1(je,jk)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
          tri_line1(je,jk)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)

          ! get triangle edge 2 (A2V3)
          !
          tri_line2(je,jk)%p1%lon = arrival_pts(je,2,1,jk,jb)
          tri_line2(je,jk)%p1%lat = arrival_pts(je,2,2,jk,jb)
          tri_line2(je,jk)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
          tri_line2(je,jk)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)



          !
          ! Compute index lists
          !


          ! Check whether departure line intersects with edge A1V3 and/or A2V3

          ! does departure-line segment intersect with A1V3?
          !
!CDIR NEXPAND(lintersect)
          lintersect_line1 = lintersect(fl_line(je,jk), tri_line1(je,jk))

          ! does departure-line segment intersect with A2V3?
          !
!CDIR NEXPAND(lintersect)
          lintersect_line2 = lintersect(fl_line(je,jk), tri_line2(je,jk))


          IF ( lintersect_line1 .AND. lintersect_line2 ) THEN

            !CASE I
            icnt_c1 = icnt_c1 + 1
            idxlist_c1(icnt_c1,jb) = je
            levlist_c1(icnt_c1,jb) = jk

          ELSE IF ( lintersect_line1 .AND. (.NOT. lintersect_line2) ) THEN

            !CASE IIa
            icnt_c2p = icnt_c2p + 1
            idxlist_c2p(icnt_c2p,jb) = je
            levlist_c2p(icnt_c2p,jb) = jk

          ELSE IF ( lintersect_line2 .AND. (.NOT. lintersect_line1) ) THEN

            !CASE IIb
            icnt_c2m = icnt_c2m + 1
            idxlist_c2m(icnt_c2m,jb) = je
            levlist_c2m(icnt_c2m,jb) = jk

          ELSE

            ! remaining cases needing further processing
            icnt_rem = icnt_rem + 1
            idxlist_rem(icnt_rem,jb) = je
            levlist_rem(icnt_rem,jb) = jk
          ENDIF

        ENDDO ! loop over edges

      ENDDO  ! loop over vertical levels




      ! Second step of index list computation
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_rem
        je = idxlist_rem(jl,jb)
        jk = levlist_rem(jl,jb)

        ! get flux area edge 1
        !
        fl_e1(je,jk)%p1%lon = arrival_pts(je,1,1,jk,jb)
        fl_e1(je,jk)%p1%lat = arrival_pts(je,1,2,jk,jb)
        fl_e1(je,jk)%p2%lon = depart_pts (je,1,1,jk,jb)
        fl_e1(je,jk)%p2%lat = depart_pts (je,1,2,jk,jb)

        ! get flux area edge 2
        !
        fl_e2(je,jk)%p1%lon = arrival_pts(je,2,1,jk,jb)
        fl_e2(je,jk)%p1%lat = arrival_pts(je,2,2,jk,jb)
        fl_e2(je,jk)%p2%lon = depart_pts (je,2,1,jk,jb)
        fl_e2(je,jk)%p2%lat = depart_pts (je,2,2,jk,jb)


        ! Check whether flux area edge 2 intersects with triangle edge 1
        !
!CDIR NEXPAND(lintersect)
        lintersect_e2_line1 = lintersect(fl_e2(je,jk), tri_line1(je,jk))


        ! Check whether flux area edge 1 intersects with triangle edge 2
        !
!CDIR NEXPAND(lintersect)
        lintersect_e1_line2 = lintersect(fl_e1(je,jk), tri_line2(je,jk))


        IF ( lintersect_e2_line1 ) THEN

          !CASE IIIa
          icnt_c3p = icnt_c3p + 1
          idxlist_c3p(icnt_c3p,jb) = je
          levlist_c3p(icnt_c3p,jb) = jk

        ELSE IF ( lintersect_e1_line2 ) THEN

          !CASE IIIb
          icnt_c3m = icnt_c3m + 1
          idxlist_c3m(icnt_c3m,jb) = je
          levlist_c3m(icnt_c3m,jb) = jk

        ELSE IF ( ABS(p_vn(je,jk,jb)) < 0.1_wp ) THEN

          ! CASE IV
          ! special case of very small normal velocity
          icnt_vn0 = icnt_vn0 + 1
          idxlist_vn0(icnt_vn0,jb) = je
          levlist_vn0(icnt_vn0,jb) = jk
        ELSE     ! error index list

          ! ERROR
          icnt_err = icnt_err + 1
          idxlist_err(icnt_err,jb) = je
          levlist_err(icnt_err,jb) = jk

          ! adding the error points to the weak-vn list is done in order to ensure
          ! reproducible (though bad) results in cases of too high wind speed
          icnt_vn0 = icnt_vn0 + 1
          idxlist_vn0(icnt_vn0,jb) = je
          levlist_vn0(icnt_vn0,jb) = jk
        ENDIF

      ENDDO  !jl


      IF ( icnt_err>0 ) THEN
        ! Check for unassigned grid points (i.e. collected in list_Err) because of CFL violation
        DO jl = 1, icnt_err

          je = idxlist_err(jl,jb)
          jk = levlist_err(jl,jb)

          ! Note: direct write to standard error output is used here by intention
          ! because warnings are otherwise suppressed for all PEs but PE0
          WRITE(0,'(a,a,i5,a,i5,a,i5,a,f8.2,a,f8.2,a,f10.2,a,f10.2)') &
               & 'horizontal CFL number exceeded at:',                &
               & ' je =',je,' jk =',jk,' jb =',jb,                    &
               & ' lon(deg)=',p_patch%edges%center(je,jb)%lon*rad2deg,&
               & ' lat(deg)=',p_patch%edges%center(je,jb)%lat*rad2deg,&
               & ' vn(m/s)=',p_vn(je,jk,jb),' vt(m/s)=',p_vt(je,jk,jb)
        ENDDO
      ENDIF


      ! Get corners of flux area patches
      !

      !
      ! CASE 1
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c1
        je = idxlist_c1(jl,jb)
        jk = levlist_c1(jl,jb)


        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line1
        ! Compute intersection point of fl_line with tri_line2
        !
        ps1(1:2) = line_intersect(fl_line(je,jk), tri_line1(je,jk))
        ps2(1:2) = line_intersect(fl_line(je,jk), tri_line2(je,jk))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 S2 S1
        ! vn < 0: A1 S1 S2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           ps1(2),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = ps2(1:2)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(1),arrival_pts(je,2,1,jk,jb), &
          &                           lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(2),arrival_pts(je,2,2,jk,jb), &
          &                           lvn_sys_pos)

        ! patch 1
        ! vn > 0: A1 S1 D1 A1 (degenerated)
        ! vn < 0: A1 D1 S1 A1 (degenerated)
        !
        dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch1(je,4,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch1(je,2,1,jk,jb)   = MERGE(ps1(1),                    &
          &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
        dreg_patch1(je,2,2,jk,jb)   = MERGE(ps1(2),                    &
          &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
        dreg_patch1(je,3,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch1(je,3,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
          &                           ps1(2), lvn_sys_pos)


        ! patch 2
        ! vn > 0: A2 D2 S2 A2 (degenerated)
        ! vn < 0: A2 S2 D2 A2 (degenerated)
        !
        dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
        dreg_patch2(je,4,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
        dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
          &                           ps2(1), lvn_sys_pos)
        dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
          &                           ps2(2), lvn_sys_pos)
        dreg_patch2(je,3,1,jk,jb)   = MERGE(ps2(1),                    &
          &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
        dreg_patch2(je,3,2,jk,jb)   = MERGE(ps2(2),                    &
          &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 2a
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c2p
        je = idxlist_c2p(jl,jb)
        jk = levlist_c2p(jl,jb)


        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line1
        !
        ps1(1:2) = line_intersect(fl_line(je,jk), tri_line1(je,jk))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 D2 S1
        ! vn < 0: A1 S1 D2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           ps1(2),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = depart_pts(je,2,1:2,jk,jb)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(1),arrival_pts(je,2,1,jk,jb), &
          &                           lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(2),arrival_pts(je,2,2,jk,jb), &
          &                           lvn_sys_pos)

          ! patch 1
          ! vn > 0: A1 S1 D1 A1 (degenerated)
          ! vn < 0: A1 D1 S1 A1 (degenerated)
          !
        dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch1(je,4,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch1(je,2,1,jk,jb)   = MERGE(ps1(1),                    &
          &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
        dreg_patch1(je,2,2,jk,jb)   = MERGE(ps1(2),                    &
          &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
        dreg_patch1(je,3,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch1(je,3,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
          &                           ps1(2), lvn_sys_pos)


        ! patch 2 (non-existing)
        !
        dreg_patch2(je,1,1:2,jk,jb) = 0._wp
        dreg_patch2(je,2,1:2,jk,jb) = 0._wp
        dreg_patch2(je,3,1:2,jk,jb) = 0._wp
        dreg_patch2(je,4,1:2,jk,jb) = 0._wp

      ENDDO  ! jl



      !
      ! CASE 2b
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c2m
        je = idxlist_c2m(jl,jb)
        jk = levlist_c2m(jl,jb)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line2
        !
        ps2(1:2) = line_intersect(fl_line(je,jk), tri_line2(je,jk))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 S2 D1
        ! vn < 0: A1 D1 S2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           depart_pts(je,1,1,jk,jb),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           depart_pts(je,1,2,jk,jb),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = ps2(1:2)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
          &                           arrival_pts(je,2,1,jk,jb), lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
          &                           arrival_pts(je,2,2,jk,jb), lvn_sys_pos)


        ! patch 1 (non-existing)
        !
        dreg_patch1(je,1,1:2,jk,jb) = 0._wp
        dreg_patch1(je,2,1:2,jk,jb) = 0._wp
        dreg_patch1(je,3,1:2,jk,jb) = 0._wp
        dreg_patch1(je,4,1:2,jk,jb) = 0._wp


        ! patch 2
        ! vn > 0: A2 D2 S2 A2 (degenerated)
        ! vn < 0: A2 S2 D2 A2 (degenerated)
        !
        dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
        dreg_patch2(je,4,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
        dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
          &                           ps2(1), lvn_sys_pos)
        dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
          &                           ps2(2), lvn_sys_pos)
        dreg_patch2(je,3,1,jk,jb)   = MERGE(ps2(1),                    &
          &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
        dreg_patch2(je,3,2,jk,jb)   = MERGE(ps2(2),                    &
          &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 3a
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c3p
        je = idxlist_c3p(jl,jb)
        jk = levlist_c3p(jl,jb)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_e2 with tri_line1
        !
        pi1(1:2) = line_intersect(fl_e2(je,jk), tri_line1(je,jk))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 I1 A1 (degenerated)
        ! vn < 0: A1 I1 A2 A1 (degenerated)
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,4,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           pi1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           pi1(2),lvn_sys_pos)
        dreg_patch0(je,3,1,jk,jb)   = MERGE(pi1(1),                    &
          &                           arrival_pts(je,2,1,jk,jb),lvn_sys_pos)
        dreg_patch0(je,3,2,jk,jb)   = MERGE(pi1(2),                    &
          &                           arrival_pts(je,2,2,jk,jb),lvn_sys_pos)


        ! patch 1
        ! vn > 0: A1 I1 D2 D1
        ! vn < 0: A1 D1 D2 I1
        !
        dreg_patch1(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch1(je,2,1,jk,jb)   = MERGE(pi1(1),                    &
          &                           depart_pts(je,1,1,jk,jb), lvn_sys_pos)
        dreg_patch1(je,2,2,jk,jb)   = MERGE(pi1(2),                    &
          &                           depart_pts(je,1,2,jk,jb), lvn_sys_pos)
        dreg_patch1(je,3,1:2,jk,jb) = depart_pts(je,2,1:2,jk,jb)
        dreg_patch1(je,4,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb),  &
          &                           pi1(1),lvn_sys_pos)
        dreg_patch1(je,4,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb),  &
          &                           pi1(2), lvn_sys_pos)


        ! patch 2 (non-existing)
        !
        dreg_patch2(je,1,1:2,jk,jb) = 0._wp
        dreg_patch2(je,2,1:2,jk,jb) = 0._wp
        dreg_patch2(je,3,1:2,jk,jb) = 0._wp
        dreg_patch2(je,4,1:2,jk,jb) = 0._wp

      ENDDO  ! jl


      !
      ! CASE 3b
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c3m
        je = idxlist_c3m(jl,jb)
        jk = levlist_c3m(jl,jb)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_e1 with tri_line2
        pi2(1:2) = line_intersect(fl_e1(je,jk), tri_line2(je,jk))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 I2 A1 (degenerated)
        ! vn < 0: A1 I2 A2 A1 (degenerated)
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,4,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           pi2(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           pi2(2),lvn_sys_pos)
        dreg_patch0(je,3,1,jk,jb)   = MERGE(pi2(1),                    &
          &                           arrival_pts(je,2,1,jk,jb),lvn_sys_pos)
        dreg_patch0(je,3,2,jk,jb)   = MERGE(pi2(2),                    &
          &                           arrival_pts(je,2,2,jk,jb),lvn_sys_pos)


        ! patch 1 (non-existing)
        !
        dreg_patch1(je,1,1:2,jk,jb) = 0._wp
        dreg_patch1(je,2,1:2,jk,jb) = 0._wp
        dreg_patch1(je,3,1:2,jk,jb) = 0._wp
        dreg_patch1(je,4,1:2,jk,jb) = 0._wp



        ! patch 2
        ! vn > 0: A2 D2 D1 I2
        ! vn < 0: A2 I2 D1 D2
        !
        dreg_patch2(je,1,1:2,jk,jb) = arrival_pts(je,2,1:2,jk,jb)
        dreg_patch2(je,2,1,jk,jb)   = MERGE(depart_pts(je,2,1,jk,jb),  &
          &                           pi2(1), lvn_sys_pos)
        dreg_patch2(je,2,2,jk,jb)   = MERGE(depart_pts(je,2,2,jk,jb),  &
          &                           pi2(2), lvn_sys_pos)
        dreg_patch2(je,3,1:2,jk,jb) = depart_pts(je,1,1:2,jk,jb)
        dreg_patch2(je,4,1,jk,jb)   = MERGE(pi2(1),                    &
          &                           depart_pts(je,2,1,jk,jb), lvn_sys_pos)
        dreg_patch2(je,4,2,jk,jb)   = MERGE(pi2(2),                    &
          &                           depart_pts(je,2,2,jk,jb), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 4  (very small normal velocity)
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_vn0
        je = idxlist_vn0(jl,jb)
        jk = levlist_vn0(jl,jb)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! patch 0 (non-existing)
        !
        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 D1 D2 A2
        ! vn < 0: A1 A2 D2 D1
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(je,1,1:2,jk,jb)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(depart_pts(je,1,1,jk,jb), &
          &                           arrival_pts(je,2,1,jk,jb),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(depart_pts(je,1,2,jk,jb), &
          &                           arrival_pts(je,2,2,jk,jb),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = depart_pts(je,2,1:2,jk,jb)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(arrival_pts(je,2,1,jk,jb), &
          &                           depart_pts(je,1,1,jk,jb),lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(arrival_pts(je,2,2,jk,jb), &
          &                           depart_pts(je,1,2,jk,jb),lvn_sys_pos)


        ! patch 1 (non-existing)
        !
        dreg_patch1(je,1,1:2,jk,jb) = 0._wp
        dreg_patch1(je,2,1:2,jk,jb) = 0._wp
        dreg_patch1(je,3,1:2,jk,jb) = 0._wp
        dreg_patch1(je,4,1:2,jk,jb) = 0._wp

        ! patch 2 (non-existing)
        !
        dreg_patch2(je,1,1:2,jk,jb) = 0._wp
        dreg_patch2(je,2,1:2,jk,jb) = 0._wp
        dreg_patch2(je,3,1:2,jk,jb) = 0._wp
        dreg_patch2(je,4,1:2,jk,jb) = 0._wp

      ENDDO  ! jl

     ! end of index list stuff


      DO jk = slev, elev

         DO je = i_startidx, i_endidx

           lvn_pos = p_vn(je,jk,jb) >= 0._wp

           ! get coordinates of butterfly cell centers (from upwind cell)
           bf_cc(1,1) = MERGE(ptr_bfcc(je,jb,1,1)%lon,         &
             &                ptr_bfcc(je,jb,2,1)%lon, lvn_pos )
           bf_cc(1,2) = MERGE(ptr_bfcc(je,jb,1,1)%lat,         &
             &                ptr_bfcc(je,jb,2,1)%lat, lvn_pos )
           bf_cc(2,1) = MERGE(ptr_bfcc(je,jb,1,2)%lon,         &
             &                ptr_bfcc(je,jb,2,2)%lon, lvn_pos )
           bf_cc(2,2) = MERGE(ptr_bfcc(je,jb,1,2)%lat,         &
             &                ptr_bfcc(je,jb,2,2)%lat, lvn_pos )


           ! patch 1 in translated system
           !
           dreg_patch1(je,1,1:2,jk,jb) = dreg_patch1(je,1,1:2,jk,jb) - bf_cc(1,1:2)
           dreg_patch1(je,2,1:2,jk,jb) = dreg_patch1(je,2,1:2,jk,jb) - bf_cc(1,1:2)
           dreg_patch1(je,3,1:2,jk,jb) = dreg_patch1(je,3,1:2,jk,jb) - bf_cc(1,1:2)
           dreg_patch1(je,4,1:2,jk,jb) = dreg_patch1(je,4,1:2,jk,jb) - bf_cc(1,1:2)


           ! patch 2 in translated system
           !
           dreg_patch2(je,1,1:2,jk,jb) = dreg_patch2(je,1,1:2,jk,jb) - bf_cc(2,1:2)
           dreg_patch2(je,2,1:2,jk,jb) = dreg_patch2(je,2,1:2,jk,jb) - bf_cc(2,1:2)
           dreg_patch2(je,3,1:2,jk,jb) = dreg_patch2(je,3,1:2,jk,jb) - bf_cc(2,1:2)
           dreg_patch2(je,4,1:2,jk,jb) = dreg_patch2(je,4,1:2,jk,jb) - bf_cc(2,1:2)



           ! store global index of the underlying grid cell
           !
           patch1_cell_idx(je,jk,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,1), &
             &                               p_patch%edges%butterfly_idx(je,jb,2,1), &
             &                               lvn_pos)
           patch2_cell_idx(je,jk,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,2), &
             &                               p_patch%edges%butterfly_idx(je,jb,2,2), &
             &                               lvn_pos)

           patch1_cell_blk(je,jk,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,1), &
             &                               p_patch%edges%butterfly_blk(je,jb,2,1), &
             &                               lvn_pos)
           patch2_cell_blk(je,jk,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,2), &
             &                               p_patch%edges%butterfly_blk(je,jb,2,2), &
             &                               lvn_pos)


         ENDDO ! loop over edges

      ENDDO   ! loop over vertical levels
    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE divide_flux_area




  !-------------------------------------------------------------------------
  !>
  !! List based version of divide_flux_area
  !!
  !! Flux area (aka. departure region) is subdivided according to its overlap
  !! with the underlying grid.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-11-01)
  !!
  SUBROUTINE divide_flux_area_list(p_patch, p_int, p_vn, p_vt, falist, &
    &                         dreg_patch0, dreg_patch1, dreg_patch2,   &
    &                         patch1_cell_idx, patch1_cell_blk,        &
    &                         patch2_cell_idx, patch2_cell_blk,        &
    &                         opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
      &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
      &  p_vt(:,:,:)               !< edge midpoints

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard
      &  falist                     !< Miura-type treatment of flux areas is
                                    !< insufficient

    REAL(vp), INTENT(INOUT) ::    &  !< patch 0 of subdivided departure region
      & dreg_patch0(:,:,:,:,:)       !< coordinates

    REAL(vp), INTENT(OUT) ::  &  !< patch 0,1,2 of subdivided departure region
      & dreg_patch1(:,:,:,:), &  !< coordinates
      & dreg_patch2(:,:,:,:)     !< dim: (npoints,4,2,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of underlying cell
      & patch1_cell_idx(:,:), &  !< dim: (npoints,ptr_p%nblks_e)
      & patch1_cell_blk(:,:), &
      & patch2_cell_idx(:,:), &
      & patch2_cell_blk(:,:)


    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)


    REAL(wp) ::       &   !< coordinates of arrival points. The origin
      &  arrival_pts(falist%npoints,2,2)
                          !< of the coordinate system is at the circumcenter of
                          !< the upwind cell. Unit vectors point to local East
                          !< and North. (geographical coordinates)

    REAL(wp) ::    &   !< coordinates of departure points. The origin
      &  depart_pts(falist%npoints,2,2)
                       !< of the coordinate system is at the circumcenter of
                       !< the upwind cell. Unit vectors point to local East
                       !< and North. (geographical coordinates)

    TYPE(t_line) ::                  & !< departure-line segment
      &  fl_line(falist%npoints)

    TYPE(t_line) ::             & !< departure area edges
      &  fl_e1(falist%npoints), & !< edge 1
      &  fl_e2(falist%npoints)    !< edge 2

    TYPE(t_line) ::                 & !< triangle edge
      &  tri_line1(falist%npoints), &
      &  tri_line2(falist%npoints)

    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of vertex3
      &  ptr_v3(:,:,:)
    TYPE(t_geographical_coordinates), POINTER :: & !< pointer to coordinates of
      &  ptr_bfcc(:,:,:,:)                         !< of butterfly cell centers


    REAL(wp) :: ps1(2),              & !< coordinates of intersection
      &         ps2(2)                 !< points S1, S2

    REAL(wp) :: pi1(2),              & !< coordinates of intersection
      &         pi2(2)                 !< points I1, I2

    REAL(wp) :: bf_cc(2,2)             !< coordinates of butterfly cell centers

    INTEGER :: je, jk, jb, jl          !< loop index of edge, vert level, block, lists
    INTEGER :: ie                      !< index list loop counter
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_rlstart, i_rlend, i_nchdom

    LOGICAL :: lintersect_line1, lintersect_line2
    LOGICAL :: lintersect_e2_line1, lintersect_e1_line2
    LOGICAL :: lvn_pos, lvn_sys_pos
    INTEGER :: icnt_c1, icnt_c2p, icnt_c3p, icnt_c2m, icnt_c3m
    INTEGER :: icnt_rem, icnt_err, icnt_vn0

    INTEGER ::           &         !< ie index list
      &  ielist_c1 (falist%npoints), &
      &  ielist_c2p(falist%npoints), &
      &  ielist_c3p(falist%npoints), &
      &  ielist_c2m(falist%npoints), &
      &  ielist_c3m(falist%npoints), &
      &  ielist_rem(falist%npoints), &
      &  ielist_vn0(falist%npoints), &
      &  ielist_err(falist%npoints)

    INTEGER ::           &         !< je index list
      &  idxlist_c1 (falist%npoints), &
      &  idxlist_c2p(falist%npoints), &
      &  idxlist_c3p(falist%npoints), &
      &  idxlist_c2m(falist%npoints), &
      &  idxlist_c3m(falist%npoints), &
      &  idxlist_rem(falist%npoints), &
      &  idxlist_vn0(falist%npoints), &
      &  idxlist_err(falist%npoints)


    INTEGER ::           &         !< jk index list
      &  levlist_c1 (falist%npoints), &
      &  levlist_c2p(falist%npoints), &
      &  levlist_c3p(falist%npoints), &
      &  levlist_c2m(falist%npoints), &
      &  levlist_c3m(falist%npoints), &
      &  levlist_rem(falist%npoints), &
      &  levlist_vn0(falist%npoints), &
      &  levlist_err(falist%npoints)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_traj: divide_flux_area'

  !-------------------------------------------------------------------------


    ! Check for optional arguments
    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF


    ! number of child domains
    i_nchdom   = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    ! pointer to coordinates of vertex3 (i.e. vertex which belongs to
    ! the upwind cell but does not belong to the current edge.)
    !
    ptr_v3   => p_int%pos_on_tplane_c_edge(:,:,:,3)

    ! pointer to coordinates of butterfly cell centers
    !
    ptr_bfcc => p_int%pos_on_tplane_c_edge(:,:,:,4:5)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,jl,ie,lvn_pos,fl_line,tri_line1,                  &
!$OMP            tri_line2,fl_e1,fl_e2,lintersect_line1,lintersect_line2,   &
!$OMP            lintersect_e2_line1,lintersect_e1_line2,icnt_c1,icnt_c2p,  &
!$OMP            icnt_c2m,icnt_rem,icnt_c3p,icnt_c3m,icnt_vn0,icnt_err,     &
!$OMP            lvn_sys_pos,ps1,ps2,pi1,pi2,bf_cc,arrival_pts,depart_pts,  &
!$OMP            idxlist_c1,levlist_c1,idxlist_c2p,levlist_c2p,             &
!$OMP            idxlist_c3p,levlist_c3p,idxlist_c2m,levlist_c2m,           &
!$OMP            idxlist_c3m,levlist_c3m,idxlist_rem,levlist_rem,           &
!$OMP            idxlist_vn0,levlist_vn0,idxlist_err,levlist_err,           &
!$OMP            ielist_c1,ielist_c2p,ielist_c3p,ielist_c2m,                &
!$OMP            ielist_c3m,ielist_rem,ielist_vn0,ielist_err), SCHEDULE(guided)

    DO jb = i_startblk, i_endblk


      ! get arrival and departure points. Note that the indices of the departure
      ! points have to be switched so that departure point 1 belongs to arrival
      ! point one and departure point 2 to arrival point 2.
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        arrival_pts(ie,1:2,1:2) = dreg_patch0(je,1:2,1:2,jk,jb)
        depart_pts (ie,1:2,1:2) = dreg_patch0(je,4:3:-1,1:2,jk,jb)

      ENDDO

      ! Reset counter
      icnt_c1  = 0
      icnt_c2p = 0
      icnt_c2m = 0
      icnt_c3p = 0
      icnt_c3m = 0
      icnt_rem = 0
      icnt_err = 0
      icnt_vn0 = 0


      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)


        lvn_pos = p_vn(je,jk,jb) >= 0._wp

        !
        ! get flux area departure-line segment
        !
        fl_line(ie)%p1%lon = depart_pts(ie,1,1)
        fl_line(ie)%p1%lat = depart_pts(ie,1,2)
        fl_line(ie)%p2%lon = depart_pts(ie,2,1)
        fl_line(ie)%p2%lat = depart_pts(ie,2,2)

        ! get triangle edge 1 (A1V3)
        !
        tri_line1(ie)%p1%lon = arrival_pts(ie,1,1)
        tri_line1(ie)%p1%lat = arrival_pts(ie,1,2)
        tri_line1(ie)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
        tri_line1(ie)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)

        ! get triangle edge 2 (A2V3)
        !
        tri_line2(ie)%p1%lon = arrival_pts(ie,2,1)
        tri_line2(ie)%p1%lat = arrival_pts(ie,2,2)
        tri_line2(ie)%p2%lon = MERGE(ptr_v3(je,jb,1)%lon,ptr_v3(je,jb,2)%lon,lvn_pos)
        tri_line2(ie)%p2%lat = MERGE(ptr_v3(je,jb,1)%lat,ptr_v3(je,jb,2)%lat,lvn_pos)



        !
        ! Compute index lists
        !


        ! Check whether departure line intersects with edge A1V3 and/or A2V3

        ! does departure-line segment intersect with A1V3?
        !
!CDIR NEXPAND(lintersect)
        lintersect_line1 = lintersect(fl_line(ie), tri_line1(ie))

        ! does departure-line segment intersect with A2V3?
        !
!CDIR NEXPAND(lintersect)
        lintersect_line2 = lintersect(fl_line(ie), tri_line2(ie))


        IF ( lintersect_line1 .AND. lintersect_line2 ) THEN

          !CASE I
          icnt_c1 = icnt_c1 + 1
          ielist_c1(icnt_c1)  = ie
          idxlist_c1(icnt_c1) = je
          levlist_c1(icnt_c1) = jk

        ELSE IF ( lintersect_line1 .AND. (.NOT. lintersect_line2) ) THEN

          !CASE IIa
          icnt_c2p = icnt_c2p + 1
          ielist_c2p(icnt_c2p)  = ie
          idxlist_c2p(icnt_c2p) = je
          levlist_c2p(icnt_c2p) = jk

        ELSE IF ( lintersect_line2 .AND. (.NOT. lintersect_line1) ) THEN

          !CASE IIb
          icnt_c2m = icnt_c2m + 1
          ielist_c2m(icnt_c2m)  = ie
          idxlist_c2m(icnt_c2m) = je
          levlist_c2m(icnt_c2m) = jk

        ELSE

          ! remaining cases needing further processing
          icnt_rem = icnt_rem + 1
          ielist_rem(icnt_rem)  = ie
          idxlist_rem(icnt_rem) = je
          levlist_rem(icnt_rem) = jk
        ENDIF

      ENDDO ! loop over index list



      ! Second step of index list computation
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_rem
        ie = ielist_rem(jl)
        je = idxlist_rem(jl)
        jk = levlist_rem(jl)

        ! get flux area edge 1
        !
        fl_e1(ie)%p1%lon = arrival_pts(ie,1,1)
        fl_e1(ie)%p1%lat = arrival_pts(ie,1,2)
        fl_e1(ie)%p2%lon = depart_pts (ie,1,1)
        fl_e1(ie)%p2%lat = depart_pts (ie,1,2)

        ! get flux area edge 2
        !
        fl_e2(ie)%p1%lon = arrival_pts(ie,2,1)
        fl_e2(ie)%p1%lat = arrival_pts(ie,2,2)
        fl_e2(ie)%p2%lon = depart_pts (ie,2,1)
        fl_e2(ie)%p2%lat = depart_pts (ie,2,2)


        ! Check whether flux area edge 2 intersects with triangle edge 1
        !
!CDIR NEXPAND(lintersect)
        lintersect_e2_line1 = lintersect(fl_e2(ie), tri_line1(ie))


        ! Check whether flux area edge 1 intersects with triangle edge 2
        !
!CDIR NEXPAND(lintersect)
        lintersect_e1_line2 = lintersect(fl_e1(ie), tri_line2(ie))


        IF ( lintersect_e2_line1 ) THEN

          !CASE IIIa
          icnt_c3p = icnt_c3p + 1
          ielist_c3p(icnt_c3p)  = ie
          idxlist_c3p(icnt_c3p) = je
          levlist_c3p(icnt_c3p) = jk

        ELSE IF ( lintersect_e1_line2 ) THEN

          !CASE IIIb
          icnt_c3m = icnt_c3m + 1
          ielist_c3m(icnt_c3m)  = ie
          idxlist_c3m(icnt_c3m) = je
          levlist_c3m(icnt_c3m) = jk

        ELSE IF ( ABS(p_vn(je,jk,jb)) < 0.1_wp ) THEN

          ! CASE IV
          ! special case of very small normal velocity
          icnt_vn0 = icnt_vn0 + 1
          ielist_vn0(icnt_vn0)  = ie
          idxlist_vn0(icnt_vn0) = je
          levlist_vn0(icnt_vn0) = jk

        ELSE     ! error index list

          ! ERROR
          icnt_err = icnt_err + 1
          ielist_err(icnt_err)  = ie
          idxlist_err(icnt_err) = je
          levlist_err(icnt_err) = jk

          ! adding the error points to the weak-vn list is done in order to ensure
          ! reproducible (though bad) results in cases of too high wind speed
          icnt_vn0 = icnt_vn0 + 1
          ielist_vn0(icnt_vn0)  = ie
          idxlist_vn0(icnt_vn0) = je
          levlist_vn0(icnt_vn0) = jk
        ENDIF

      ENDDO  !jl


      IF ( icnt_err>0 ) THEN
        ! Check for unassigned grid points (i.e. collected in list_Err) because of CFL violation
        DO jl = 1, icnt_err

          je = idxlist_err(jl)
          jk = levlist_err(jl)

          ! Note: direct write to standard error output is used here by intention
          ! because warnings are otherwise suppressed for all PEs but PE0
          WRITE(0,'(a,a,i5,a,i5,a,i5,a,f8.2,a,f8.2,a,f10.2,a,f10.2)') &
               & 'horizontal CFL number exceeded at:',                &
               & ' je =',je,' jk =',jk,' jb =',jb,                    &
               & ' lon(deg)=',p_patch%edges%center(je,jb)%lon*rad2deg,&
               & ' lat(deg)=',p_patch%edges%center(je,jb)%lat*rad2deg,&
               & ' vn(m/s)=',p_vn(je,jk,jb),' vt(m/s)=',p_vt(je,jk,jb)
        ENDDO
      ENDIF


      ! Get corners of flux area patches
      !

      !
      ! CASE 1
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c1
        ie = ielist_c1(jl)
        je = idxlist_c1(jl)
        jk = levlist_c1(jl)


        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line1
        ! Compute intersection point of fl_line with tri_line2
        !
        ps1(1:2) = line_intersect(fl_line(ie), tri_line1(ie))
        ps2(1:2) = line_intersect(fl_line(ie), tri_line2(ie))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 S2 S1
        ! vn < 0: A1 S1 S2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(ie,2,1), &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(ie,2,2), &
          &                           ps1(2),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = ps2(1:2)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(1),arrival_pts(ie,2,1), &
          &                           lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(2),arrival_pts(ie,2,2), &
          &                           lvn_sys_pos)

        ! patch 1
        ! vn > 0: A1 S1 D1 A1 (degenerated)
        ! vn < 0: A1 D1 S1 A1 (degenerated)
        !
        dreg_patch1(ie,1,1:2,jb) = arrival_pts(ie,1,1:2)
        dreg_patch1(ie,4,1:2,jb) = arrival_pts(ie,1,1:2)
        dreg_patch1(ie,2,1,jb)   = MERGE(ps1(1),depart_pts(ie,1,1), lvn_sys_pos)
        dreg_patch1(ie,2,2,jb)   = MERGE(ps1(2),depart_pts(ie,1,2), lvn_sys_pos)
        dreg_patch1(ie,3,1,jb)   = MERGE(depart_pts(ie,1,1),ps1(1), lvn_sys_pos)
        dreg_patch1(ie,3,2,jb)   = MERGE(depart_pts(ie,1,2),ps1(2), lvn_sys_pos)


        ! patch 2
        ! vn > 0: A2 D2 S2 A2 (degenerated)
        ! vn < 0: A2 S2 D2 A2 (degenerated)
        !
        dreg_patch2(ie,1,1:2,jb) = arrival_pts(ie,2,1:2)
        dreg_patch2(ie,4,1:2,jb) = arrival_pts(ie,2,1:2)
        dreg_patch2(ie,2,1,jb)   = MERGE(depart_pts(ie,2,1),ps2(1), lvn_sys_pos)
        dreg_patch2(ie,2,2,jb)   = MERGE(depart_pts(ie,2,2),ps2(2), lvn_sys_pos)
        dreg_patch2(ie,3,1,jb)   = MERGE(ps2(1),depart_pts(ie,2,1), lvn_sys_pos)
        dreg_patch2(ie,3,2,jb)   = MERGE(ps2(2),depart_pts(ie,2,2), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 2a
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c2p
        ie = ielist_c2p(jl)
        je = idxlist_c2p(jl)
        jk = levlist_c2p(jl)


        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line1
        !
        ps1(1:2) = line_intersect(fl_line(ie), tri_line1(ie))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 D2 S1
        ! vn < 0: A1 S1 D2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(ie,2,1), &
          &                           ps1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(ie,2,2), &
          &                           ps1(2),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = depart_pts(ie,2,1:2)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(ps1(1),arrival_pts(ie,2,1), &
          &                           lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(ps1(2),arrival_pts(ie,2,2), &
          &                           lvn_sys_pos)

          ! patch 1
          ! vn > 0: A1 S1 D1 A1 (degenerated)
          ! vn < 0: A1 D1 S1 A1 (degenerated)
          !
        dreg_patch1(ie,1,1:2,jb) = arrival_pts(ie,1,1:2)
        dreg_patch1(ie,4,1:2,jb) = arrival_pts(ie,1,1:2)
        dreg_patch1(ie,2,1,jb)   = MERGE(ps1(1),depart_pts(ie,1,1), lvn_sys_pos)
        dreg_patch1(ie,2,2,jb)   = MERGE(ps1(2),depart_pts(ie,1,2), lvn_sys_pos)
        dreg_patch1(ie,3,1,jb)   = MERGE(depart_pts(ie,1,1),ps1(1), lvn_sys_pos)
        dreg_patch1(ie,3,2,jb)   = MERGE(depart_pts(ie,1,2),ps1(2), lvn_sys_pos)


        ! patch 2 (non-existing)
        !
        dreg_patch2(ie,1,1:2,jb) = 0._wp
        dreg_patch2(ie,2,1:2,jb) = 0._wp
        dreg_patch2(ie,3,1:2,jb) = 0._wp
        dreg_patch2(ie,4,1:2,jb) = 0._wp

      ENDDO  ! jl



      !
      ! CASE 2b
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c2m
        ie = ielist_c2m(jl)
        je = idxlist_c2m(jl)
        jk = levlist_c2m(jl)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_line with tri_line2
        !
        ps2(1:2) = line_intersect(fl_line(ie), tri_line2(ie))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 S2 D1
        ! vn < 0: A1 D1 S2 A2
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(ie,2,1), &
          &                           depart_pts(ie,1,1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(ie,2,2), &
          &                           depart_pts(ie,1,2),lvn_sys_pos)
        dreg_patch0(je,3,1:2,jk,jb) = ps2(1:2)
        dreg_patch0(je,4,1,jk,jb)   = MERGE(depart_pts(ie,1,1),  &
          &                           arrival_pts(ie,2,1), lvn_sys_pos)
        dreg_patch0(je,4,2,jk,jb)   = MERGE(depart_pts(ie,1,2),  &
          &                           arrival_pts(ie,2,2), lvn_sys_pos)


        ! patch 1 (non-existing)
        !
        dreg_patch1(ie,1,1:2,jb) = 0._wp
        dreg_patch1(ie,2,1:2,jb) = 0._wp
        dreg_patch1(ie,3,1:2,jb) = 0._wp
        dreg_patch1(ie,4,1:2,jb) = 0._wp


        ! patch 2
        ! vn > 0: A2 D2 S2 A2 (degenerated)
        ! vn < 0: A2 S2 D2 A2 (degenerated)
        !
        dreg_patch2(ie,1,1:2,jb) = arrival_pts(ie,2,1:2)
        dreg_patch2(ie,4,1:2,jb) = arrival_pts(ie,2,1:2)
        dreg_patch2(ie,2,1,jb)   = MERGE(depart_pts(ie,2,1),ps2(1), lvn_sys_pos)
        dreg_patch2(ie,2,2,jb)   = MERGE(depart_pts(ie,2,2),ps2(2), lvn_sys_pos)
        dreg_patch2(ie,3,1,jb)   = MERGE(ps2(1),depart_pts(ie,2,1), lvn_sys_pos)
        dreg_patch2(ie,3,2,jb)   = MERGE(ps2(2),depart_pts(ie,2,2), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 3a
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c3p
        ie = ielist_c3p(jl)
        je = idxlist_c3p(jl)
        jk = levlist_c3p(jl)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_e2 with tri_line1
        !
        pi1(1:2) = line_intersect(fl_e2(ie), tri_line1(ie))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 I1 A1 (degenerated)
        ! vn < 0: A1 I1 A2 A1 (degenerated)
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,4,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(ie,2,1),pi1(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(ie,2,2),pi1(2),lvn_sys_pos)
        dreg_patch0(je,3,1,jk,jb)   = MERGE(pi1(1),arrival_pts(ie,2,1),lvn_sys_pos)
        dreg_patch0(je,3,2,jk,jb)   = MERGE(pi1(2),arrival_pts(ie,2,2),lvn_sys_pos)


        ! patch 1
        ! vn > 0: A1 I1 D2 D1
        ! vn < 0: A1 D1 D2 I1
        !
        dreg_patch1(ie,1,1:2,jb) = arrival_pts(ie,1,1:2)
        dreg_patch1(ie,2,1,jb)   = MERGE(pi1(1),depart_pts(ie,1,1), lvn_sys_pos)
        dreg_patch1(ie,2,2,jb)   = MERGE(pi1(2),depart_pts(ie,1,2), lvn_sys_pos)
        dreg_patch1(ie,3,1:2,jb) = depart_pts(ie,2,1:2)
        dreg_patch1(ie,4,1,jb)   = MERGE(depart_pts(ie,1,1),pi1(1), lvn_sys_pos)
        dreg_patch1(ie,4,2,jb)   = MERGE(depart_pts(ie,1,2),pi1(2), lvn_sys_pos)


        ! patch 2 (non-existing)
        !
        dreg_patch2(ie,1,1:2,jb) = 0._wp
        dreg_patch2(ie,2,1:2,jb) = 0._wp
        dreg_patch2(ie,3,1:2,jb) = 0._wp
        dreg_patch2(ie,4,1:2,jb) = 0._wp

      ENDDO  ! jl


      !
      ! CASE 3b
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_c3m
        ie = ielist_c3m(jl)
        je = idxlist_c3m(jl)
        jk = levlist_c3m(jl)

        lvn_sys_pos = (p_vn(je,jk,jb) * p_patch%edges%tangent_orientation(je,jb)) >= 0._wp

        ! Compute intersection point of fl_e1 with tri_line2
        pi2(1:2) = line_intersect(fl_e1(ie), tri_line2(ie))

        ! store corners of flux area patches (counterclockwise)
        ! patch 0
        ! vn > 0: A1 A2 I2 A1 (degenerated)
        ! vn < 0: A1 I2 A2 A1 (degenerated)
        !
        dreg_patch0(je,1,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,4,1:2,jk,jb) = arrival_pts(ie,1,1:2)
        dreg_patch0(je,2,1,jk,jb)   = MERGE(arrival_pts(ie,2,1),pi2(1),lvn_sys_pos)
        dreg_patch0(je,2,2,jk,jb)   = MERGE(arrival_pts(ie,2,2),pi2(2),lvn_sys_pos)
        dreg_patch0(je,3,1,jk,jb)   = MERGE(pi2(1),arrival_pts(ie,2,1),lvn_sys_pos)
        dreg_patch0(je,3,2,jk,jb)   = MERGE(pi2(2),arrival_pts(ie,2,2),lvn_sys_pos)


        ! patch 1 (non-existing)
        !
        dreg_patch1(ie,1,1:2,jb) = 0._wp
        dreg_patch1(ie,2,1:2,jb) = 0._wp
        dreg_patch1(ie,3,1:2,jb) = 0._wp
        dreg_patch1(ie,4,1:2,jb) = 0._wp



        ! patch 2
        ! vn > 0: A2 D2 D1 I2
        ! vn < 0: A2 I2 D1 D2
        !
        dreg_patch2(ie,1,1:2,jb) = arrival_pts(ie,2,1:2)
        dreg_patch2(ie,2,1,jb)   = MERGE(depart_pts(ie,2,1),pi2(1), lvn_sys_pos)
        dreg_patch2(ie,2,2,jb)   = MERGE(depart_pts(ie,2,2),pi2(2), lvn_sys_pos)
        dreg_patch2(ie,3,1:2,jb) = depart_pts(ie,1,1:2)
        dreg_patch2(ie,4,1,jb)   = MERGE(pi2(1),depart_pts(ie,2,1), lvn_sys_pos)
        dreg_patch2(ie,4,2,jb)   = MERGE(pi2(2),depart_pts(ie,2,2), lvn_sys_pos)
      ENDDO  ! jl



      !
      ! CASE 4  (very small normal velocity)
      !
!CDIR NODEP,VOVERTAKE,VOB
      DO jl = 1, icnt_vn0
        ie = ielist_vn0(jl)

        ! patch 1 (non-existing)
        !
        dreg_patch1(ie,1,1:2,jb) = 0._wp
        dreg_patch1(ie,2,1:2,jb) = 0._wp
        dreg_patch1(ie,3,1:2,jb) = 0._wp
        dreg_patch1(ie,4,1:2,jb) = 0._wp

        ! patch 2 (non-existing)
        !
        dreg_patch2(ie,1,1:2,jb) = 0._wp
        dreg_patch2(ie,2,1:2,jb) = 0._wp
        dreg_patch2(ie,3,1:2,jb) = 0._wp
        dreg_patch2(ie,4,1:2,jb) = 0._wp

      ENDDO  ! jl

     ! end of index list stuff


!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        lvn_pos = p_vn(je,jk,jb) >= 0._wp

        ! get coordinates of butterfly cell centers (from upwind cell)
        bf_cc(1,1) = MERGE(ptr_bfcc(je,jb,1,1)%lon,         &
          &                ptr_bfcc(je,jb,2,1)%lon, lvn_pos )
        bf_cc(1,2) = MERGE(ptr_bfcc(je,jb,1,1)%lat,         &
          &                ptr_bfcc(je,jb,2,1)%lat, lvn_pos )
        bf_cc(2,1) = MERGE(ptr_bfcc(je,jb,1,2)%lon,         &
          &                ptr_bfcc(je,jb,2,2)%lon, lvn_pos )
        bf_cc(2,2) = MERGE(ptr_bfcc(je,jb,1,2)%lat,         &
          &                ptr_bfcc(je,jb,2,2)%lat, lvn_pos )


        ! patch 1 in translated system
        !
        dreg_patch1(ie,1,1:2,jb) = dreg_patch1(ie,1,1:2,jb) - bf_cc(1,1:2)
        dreg_patch1(ie,2,1:2,jb) = dreg_patch1(ie,2,1:2,jb) - bf_cc(1,1:2)
        dreg_patch1(ie,3,1:2,jb) = dreg_patch1(ie,3,1:2,jb) - bf_cc(1,1:2)
        dreg_patch1(ie,4,1:2,jb) = dreg_patch1(ie,4,1:2,jb) - bf_cc(1,1:2)


        ! patch 2 in translated system
        !
        dreg_patch2(ie,1,1:2,jb) = dreg_patch2(ie,1,1:2,jb) - bf_cc(2,1:2)
        dreg_patch2(ie,2,1:2,jb) = dreg_patch2(ie,2,1:2,jb) - bf_cc(2,1:2)
        dreg_patch2(ie,3,1:2,jb) = dreg_patch2(ie,3,1:2,jb) - bf_cc(2,1:2)
        dreg_patch2(ie,4,1:2,jb) = dreg_patch2(ie,4,1:2,jb) - bf_cc(2,1:2)



        ! store global index of the underlying grid cell
        !
        patch1_cell_idx(ie,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,1), &
          &                            p_patch%edges%butterfly_idx(je,jb,2,1), &
          &                            lvn_pos)
        patch2_cell_idx(ie,jb) = MERGE(p_patch%edges%butterfly_idx(je,jb,1,2), &
          &                            p_patch%edges%butterfly_idx(je,jb,2,2), &
          &                            lvn_pos)

        patch1_cell_blk(ie,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,1), &
          &                            p_patch%edges%butterfly_blk(je,jb,2,1), &
          &                            lvn_pos)
        patch2_cell_blk(ie,jb) = MERGE(p_patch%edges%butterfly_blk(je,jb,1,2), &
          &                            p_patch%edges%butterfly_blk(je,jb,2,2), &
          &                            lvn_pos)


      ENDDO ! loop over index list

    END DO    ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE divide_flux_area_list

END MODULE mo_advection_geometry
