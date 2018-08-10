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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_smooth_topo

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, sync_patch_array
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_extpar_config,      ONLY: fac_smooth_topo, n_iter_smooth_topo,           &
    &                              heightdiff_threshold, hgtdiff_max_smooth_topo, &
    &                              lrevert_sea_height
  USE mo_math_laplace,       ONLY: nabla2_scalar, nabla4_scalar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: smooth_topo 
  PUBLIC :: smooth_topo_real_data

CONTAINS

  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Topography smoothing
  !!
  !! Topography smoothing by selectively applying nabla2 and nabla4 
  !! operators.
  !!
  !! @par Revision History
  !! Developed by G. Zaengl
  !!
  SUBROUTINE smooth_topo (p_patch, p_int, topography_c)

    TYPE(t_patch)      , INTENT(INOUT) :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(INOUT) :: topography_c(:,:) ! original topography on input
                                                            ! smoothed one on output

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il, npts, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c),    &
      &         z_topo_old(nproma,1,p_patch%nblks_c),z_nabla2_topo(nproma,1,p_patch%nblks_c),&
      &         z_hdiffmax(nproma,p_patch%nblks_c)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new,zdcoeff,z_heightdiff_threshold,rms_hdiff,smooth_fac
    LOGICAL  :: lnabla2_mask(nproma,p_patch%nblks_c)


    jg = p_patch%id

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

        DO il=1,p_patch%geometry_info%cell_type
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

        DO il=1,p_patch%geometry_info%cell_type
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
          rms_hdiff = (z_topo(jc,1,jb)-z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
                      p_patch%cells%neighbor_blk(jc,jb,1)))**2
          DO il=2,p_patch%geometry_info%cell_type

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

            rms_hdiff = rms_hdiff + (z_topo(jc,1,jb)-z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                      p_patch%cells%neighbor_blk(jc,jb,il)))**2
          ENDDO

          rms_hdiff = SQRT(rms_hdiff/REAL(p_patch%geometry_info%cell_type,wp))

          IF (hgtdiff_max_smooth_topo(jg) < 1._wp) THEN
            smooth_fac = fac_smooth_topo
          ELSE
            smooth_fac = fac_smooth_topo * MIN(1._wp,rms_hdiff/hgtdiff_max_smooth_topo(jg))
          ENDIF

          z_topo_new = z_topo (jc,1,jb) -  &
              smooth_fac * z_nabla4_topo(jc,1,jb) * &
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
          DO il=2,p_patch%geometry_info%cell_type

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

   END SUBROUTINE smooth_topo



  !-----------------------------------------------------------------------
  !>
  !! Topography smoothing for real-case runs
  !!
  !! Topography smoothing for real-case runs. Apart from smoothing the 
  !! topography field, 
  !! - the height of sea-points is reset to that of the raw topography 
  !!   data set
  !! - the SSO standard deviation field is updated based on the smoothed 
  !!   topography
  !!
  !! @par Revision History
  !! Developed by G. Zaengl  (2004).
  !! Modification by Daniel Reinert, DWD (2016-07-06)
  !! - bring sea points back to zero height after smoothing
  !!
  SUBROUTINE smooth_topo_real_data (p_patch, p_int, fr_land, fr_lake, topography_c, sso_stdh)

    TYPE(t_patch)      , INTENT(INOUT) :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(IN)    :: fr_land(:,:)
    REAL(wp)           , INTENT(IN)    :: fr_lake(:,:)
    REAL(wp)           , INTENT(INOUT) :: topography_c(:,:)
    REAL(wp)           , INTENT(INOUT) :: sso_stdh(:,:)

    ! local variables
    INTEGER  :: jb, jc
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo_c_sv(nproma,p_patch%nblks_c)
    REAL(wp) :: zhdiff


    nblks_c    = p_patch%nblks_c

    i_startblk = p_patch%cells%start_blk(3,1)

    ! save original raw topography
    z_topo_c_sv(:,:) = topography_c(:,:)


    ! topography smoothing
    !
    CALL smooth_topo(p_patch, p_int, topography_c)


    ! bring sea-points back to zero height (i.e. to original extpar values)
    !
    IF (lrevert_sea_height) THEN

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc=i_startidx, i_endidx
          !
          ! bring grid cell back to zero height, if it is entirely covered by sea. 
          IF ( (fr_land(jc,jb) + fr_lake(jc,jb)) == 0._wp ) THEN
            topography_c(jc,jb) = z_topo_c_sv(jc,jb)
          ENDIF
        ENDDO  !jc

      ENDDO  !jb
    ENDIF


    ! re-compute SSO_STDH based on smoothed topography
    !
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 3)

      DO jc = i_startidx, i_endidx

        zhdiff = topography_c(jc,jb) - z_topo_c_sv(jc,jb)
        sso_stdh(jc,jb) = SQRT(sso_stdh(jc,jb)**2 + zhdiff**2)

      ENDDO
    ENDDO


  END SUBROUTINE smooth_topo_real_data



END MODULE mo_smooth_topo

