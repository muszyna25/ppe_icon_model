!>
!! Contains routine for topography smoothing
!!
!! Contains routine for topography smoothing
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Guenther Zaengl, DWD
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-02-23)
!! Modification by Daniel Reinert, DWD (2012-02-23)
!! - Routine smooth_topography moved here from mo_ext_data
!!
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_smooth_topo

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_extpar_config,      ONLY: fac_smooth_topo, n_iter_smooth_topo,  &
    &                              heightdiff_threshold
  USE mo_math_laplace,       ONLY: nabla2_scalar, nabla4_scalar
  USE mo_intp,               ONLY: cells2verts_scalar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: smooth_topography 


  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !-------------------------------------------------------------------------

   SUBROUTINE smooth_topography (p_patch, p_int, topography_c, topography_v, sso_stdh)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(INOUT) :: topography_c(:,:), topography_v(:,:)
    REAL(wp), OPTIONAL , INTENT(INOUT) :: sso_stdh(:,:)

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il, npts, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c),    &
      &         z_topo_old(nproma,1,p_patch%nblks_c),z_nabla2_topo(nproma,1,p_patch%nblks_c),&
      &         z_topo_c_sv(nproma,p_patch%nblks_c),z_hdiffmax(nproma,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new,zdcoeff,z_heightdiff_threshold,zhdiff
    LOGICAL  :: lnabla2_mask(nproma,p_patch%nblks_c)


    jg = p_patch%id

    z_topo_v(:,1,:)  = topography_v(:,:)
    z_topo_c_sv(:,:) = topography_c(:,:)

    z_nabla4_topo(:,1,:) = 0._wp
    z_nabla2_topo(:,1,:) = 0._wp

    nblks_c    = p_patch%nblks_c

    zdcoeff = 0.05_wp ! diffusion coefficient for nabla2 diffusion
    niter   = 20      ! number of iterations for local nabla2 diffusion

    IF (n_iter_smooth_topo(jg) > 0) THEN
      WRITE(message_text,'(a,i3,a,i3)') 'number of topography smoothing steps in domain ', &
        jg, ': ', n_iter_smooth_topo(jg)
      CALL message('', TRIM(message_text))
    ENDIF

    ! Step 1: local nabla2 diffusion at grid points where the height difference to the
    ! neighbors exceeds a certain threshold value

    z_topo(:,1,:)   = topography_c(:,:)

    DO iter = 1, niter ! perform niter iterations
                       ! note: a variable number of iterations (with an exit condition) potentially 
                       ! causes trouble with MPI reproducibility

      z_heightdiff_threshold = heightdiff_threshold(jg)
      IF (iter >= niter-1 .AND. heightdiff_threshold(jg) < 2500._wp) &
        z_heightdiff_threshold = 0.75_wp*heightdiff_threshold(jg)

      z_hdiffmax(:,:)   = 0._wp
      lnabla2_mask(:,:) = .FALSE.

      CALL nabla2_scalar( z_topo, p_patch, p_int, z_nabla2_topo )

      npts = 0

      i_startblk = p_patch%cells%start_blk(2,1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO il=1,p_patch%cell_type
          DO jc = i_startidx, i_endidx
            z_hdiffmax(jc,jb) = MAX(z_hdiffmax(jc,jb), ABS(z_topo(jc,1,jb) - &
              &  z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1,              &
              &         p_patch%cells%neighbor_blk(jc,jb,il))) )

            IF (z_hdiffmax(jc,jb) > z_heightdiff_threshold) lnabla2_mask(jc,jb) = .TRUE.

          ENDDO
        ENDDO

      ENDDO

      ! set diffusion mask also true if one of the neighboring grid points 
      ! fulfills the height difference criterion
      i_startblk = p_patch%cells%start_blk(3,1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO il=1,p_patch%cell_type
          DO jc = i_startidx, i_endidx
            IF (z_hdiffmax(p_patch%cells%neighbor_idx(jc,jb,il),                   &
                p_patch%cells%neighbor_blk(jc,jb,il)) > z_heightdiff_threshold) THEN
              lnabla2_mask(jc,jb) = .TRUE.
            ENDIF
          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          IF (lnabla2_mask(jc,jb)) THEN
            npts = npts + 1
            topography_c(jc,jb) = z_topo(jc,1,jb) + zdcoeff *                      &
                                  p_patch%cells%area(jc,jb) * z_nabla2_topo(jc,1,jb)
          ENDIF
        ENDDO

      ENDDO

      z_topo(:,1,:) = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:) = z_topo(:,1,:)

    ENDDO ! iteration of local nabla2 smoothing

    i_startblk = p_patch%cells%start_blk(3,1)

    ! Step 2: local nabla4 diffusion with monotonous limiter
    DO iter = 1, n_iter_smooth_topo(jg)
      z_topo(:,1,:)   = topography_c(:,:)
      z_topo_old(:,1,:) = z_topo(:,1,:)

      CALL nabla4_scalar(z_topo, p_patch, p_int, z_nabla4_topo, &
        & opt_nabla2=z_nabla2_topo )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema

          ! compute maximum (zmaxtop) and minimum (zmintop) of neighbor cells topography
          !set zmaxtop and zmintop to first neighbor's value
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now max resp min of all neighbors
          ENDDO

          z_topo_new = z_topo (jc,1,jb) -  &
              fac_smooth_topo * z_nabla4_topo(jc,1,jb) * &
              p_patch%cells%area(jc,jb)*p_patch%cells%area(jc,jb)

          !If it was a local maximum in the old field, dont make it higher
          IF ( zmaxtop < z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) CYCLE
          ENDIF
          !If it was a local minimum in the old field, dont make it lower:
          IF ( zmintop > z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) > 0.0_wp ) CYCLE
          ENDIF

          !If it became a local maximum in the new field with regard to old neighbors, avoid it:
          IF (( zmaxtop < z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp )) THEN
            topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
            CYCLE
          ENDIF
          !If it became a local minimum in the new field with regard to old neighbors, avoid it:
          IF (( zmintop > z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
            topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
            CYCLE
          ENDIF

          topography_c(jc,jb)=z_topo_new

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)
        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now mx resp min of all neighbors
          ENDDO

           !If it became a local maximum in the new field, avoid it
           IF ( ( zmaxtop < z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) ) THEN
             topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
           !If it became a local minimum in the new field, avoid it
           ELSEIF (( zmintop > z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
             topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
           ENDIF

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

    ENDDO !iter


    ! Interpolate smooth topography from cells to vertices
    z_topo(:,1,:)   = topography_c(:,:)
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

    topography_v(:,:) = z_topo_v(:,1,:)

    ! Apply correction to SSO standard deviation
    IF (PRESENT(sso_stdh)) THEN

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          zhdiff = topography_c(jc,jb) - z_topo_c_sv(jc,jb)
          sso_stdh(jc,jb) = SQRT(sso_stdh(jc,jb)**2 + zhdiff**2)

        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE smooth_topography


END MODULE mo_smooth_topo

