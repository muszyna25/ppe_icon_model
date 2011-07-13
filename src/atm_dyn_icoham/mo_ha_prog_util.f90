!>
!! This module contains subroutines for manipulating the state vector of the
!! hydrostatic atmospheric dynamical core.
!!
!! @author Hui Wan, MPI-M
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
MODULE mo_ha_prog_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog
  USE mo_run_config,          ONLY: iqc, iqi, iqr, iqs, msg_level
  USE mo_exception,           ONLY: message, finish, message_text

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: copy_prog_state
  PUBLIC :: update_prog_state
  PUBLIC :: init_hydro_state_prog_isoRest
  PUBLIC :: hydro_state_prog_add_random
  PUBLIC :: add_tend, reset_tend, diag_tend

CONTAINS

  !>
  !! Copy values of the prognostic variables.
  !!
  !! This is used in
  !! <ol>
  !! <li> initializing 3-time-level integration schemes,
  !! <li> the Runge-Kutta method,
  !! <li> and preparing model output.
  !! </ol>
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2008-11-04)
  !!
  SUBROUTINE copy_prog_state( pt_prog_source, pt_prog_dest, &
                              lcopy_theta, lcopy_tracer    )

    TYPE(t_hydro_atm_prog) :: pt_prog_source
    TYPE(t_hydro_atm_prog) :: pt_prog_dest
    LOGICAL,INTENT(in)     :: lcopy_theta
    LOGICAL,INTENT(in)     :: lcopy_tracer

!$OMP PARALLEL
!$OMP WORKSHARE
    pt_prog_dest%vn       = pt_prog_source%vn
    pt_prog_dest%pres_sfc = pt_prog_source%pres_sfc
    pt_prog_dest%temp     = pt_prog_source%temp
!$OMP END WORKSHARE

    IF (lcopy_theta) THEN
!$OMP WORKSHARE
      pt_prog_dest%theta  = pt_prog_source%theta
!$OMP END WORKSHARE
    ENDIF

    IF (lcopy_tracer) THEN
!$OMP WORKSHARE
      pt_prog_dest%tracer = pt_prog_source%tracer
!$OMP END WORKSHARE
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE copy_prog_state
  !-------------
  !>
  !! Update the prognostic variables using the given tendency
  !! This is used in physics-dynamics coupling
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2010-07-05)
  !!  Optional argument opt_prog added by Hui Wan, MPI-M (2011-01-12)
  !!
  SUBROUTINE update_prog_state( pdtime, p_patch, p_tend,                   &
                                lupdate_ps, lupdate_theta, lupdate_tracer, &
                                p_prog, opt_prog )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_ha_prog_util:update_prog_state'

    REAL(wp),INTENT(IN)      :: pdtime
    TYPE(t_patch),INTENT(IN) :: p_patch

    LOGICAL,INTENT(IN) :: lupdate_ps
    LOGICAL,INTENT(in) :: lupdate_theta
    LOGICAL,INTENT(in) :: lupdate_tracer

    TYPE(t_hydro_atm_prog),INTENT(IN)             :: p_tend
    TYPE(t_hydro_atm_prog),INTENT(INOUT),TARGET   :: p_prog
    TYPE(t_hydro_atm_prog),INTENT(INOUT),TARGET,OPTIONAL :: opt_prog

    TYPE(t_hydro_atm_prog),POINTER :: ptr_out=>NULL()

    INTEGER :: jb, jbs, is, ie
    INTEGER :: nblks_e, nblks_c
    !------
    ! Determine in which variable will the updated values be stored.

    IF (PRESENT(opt_prog)) THEN
      ptr_out => opt_prog
    ELSE
      ptr_out => p_prog
    END IF

    ! Dimension parameters related to refinement and MPI parallelisation

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c

    ! Update cell-based variables

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       IF (lupdate_ps)                                                   &
         ptr_out%pres_sfc(is:ie,jb) =   p_prog%pres_sfc(is:ie,jb)        &
                                    & + p_tend%pres_sfc(is:ie,jb)*pdtime

       IF (lupdate_theta) THEN
         ptr_out%theta(is:ie,:,jb) =   p_prog%theta(is:ie,:,jb)        &
                                   & + p_tend%theta(is:ie,:,jb)*pdtime
       ELSE
         ptr_out%temp(is:ie,:,jb)  =   p_prog%temp(is:ie,:,jb)         &
                                   & + p_tend%temp(is:ie,:,jb)*pdtime
       ENDIF

       IF (lupdate_tracer)                                                  &
         ptr_out%tracer(is:ie,:,jb,:) =   p_prog%tracer(is:ie,:,jb,:)       &
                                      & + p_tend%tracer(is:ie,:,jb,:)*pdtime
     ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (msg_level > 10) THEN

      ! check positive definitness

      WRITE(message_text,'(a,2(E10.3,2x))') ' min tracer qc,qi = ',&
           MINVAL(ptr_out%tracer(:,:,:,iqc)),MINVAL(ptr_out%tracer(:,:,:,iqi))
      CALL message(TRIM(routine), TRIM(message_text))
      WRITE(message_text,'(a,2(E10.3,2x))') ' min tracer qr,qs = ',&
           MINVAL(ptr_out%tracer(:,:,:,iqr)),MINVAL(ptr_out%tracer(:,:,:,iqs))
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    ! Update edge-based variable

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )
       ptr_out%vn(is:ie,:,jb) =   p_prog%vn(is:ie,:,jb)       &
                              & + p_tend%vn(is:ie,:,jb)*pdtime
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    NULLIFY(ptr_out)

  END SUBROUTINE update_prog_state
  !---------
  !>
  !! Set the prognostic state vector to a resting isothermal state.
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2008-06-03)
  !!
  SUBROUTINE init_hydro_state_prog_isoRest( ptemp, ppres_sfc, pt_prog )

    REAL(wp),INTENT(IN) :: ptemp
    REAL(wp),INTENT(IN) :: ppres_sfc

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_prog

    pt_prog%vn       = 0._wp
    pt_prog%temp     = ptemp
    pt_prog%pres_sfc = ppres_sfc

  END SUBROUTINE init_hydro_state_prog_isoRest
  !-------------
  !>
  !! Add random perturbation to the normal wind.
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2008-11-04)
  !!  Date and time functions IDATE and ITIME replaced by
  !!  DATE_AND_TIME by Hui Wan, MPI-M (2009-02-05)
  !!
  SUBROUTINE hydro_state_prog_add_random(pt_patch, & ! in
                                         pt_prog,  & ! inout
                                         pscale, nproma, nlev) ! input

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_ha_prog_util:hydro_state_prog_add_random'

    TYPE(t_patch)            :: pt_patch
    TYPE(t_hydro_atm_prog) :: pt_prog

    REAL(wp), INTENT(IN) :: pscale ! magnitude of the perturbation
    INTEGER,  INTENT(IN) :: nlev   ! number of vertical layers

    ! LOCAL VARIABLES

    INTEGER :: jk
    INTEGER :: nproma, npromz, nblks

    INTEGER :: DateTimeArray(8)    ! Holds the date and time

    INTEGER :: seed_size, js, seed_trigger, ist

    INTEGER, ALLOCATABLE :: seed_array(:)

    REAL(wp) :: zrand(nproma,pt_patch%nblks_e)

    CHARACTER(len=MAX_CHAR_LENGTH) :: string
    !-----
    CALL message(TRIM(routine),'=========== generating random number =============')

    !-----------------------------------------------------------
    ! 1. prepare memory for the seed
    !-----------------------------------------------------------

    CALL RANDOM_SEED(SIZE=seed_size)
    WRITE(string,*) 'The size of the intrinsic seed is', seed_size
    CALL message(TRIM(routine),TRIM(string))

    ALLOCATE( seed_array(seed_size), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish('random number:','allocation of seed_array failed')
    ENDIF

    !----------------------------------------------------------
    ! 2. get seed
    !----------------------------------------------------------
    ! First inquire the current date and time.
    ! The output of the DATE_AND_TIME routine contains 8 elements:
    !   1: year (e.g. 2008)
    !   2: month
    !   3: day
    !   4: time difference with UTC in minutes
    !   5: hour
    !   6: minute
    !   7: seconds
    !   8: milliseconds

    CALL DATE_AND_TIME( VALUES=DateTimeArray )

    seed_trigger =   DateTimeArray(2)*100000000 + DateTimeArray(3)*1000000 &
                 & + DateTimeArray(5)*10000     + DateTimeArray(6)*100     &
                 & + DateTimeArray(7)

    WRITE(string,*) 'the seed trigger is', seed_trigger
    CALL message(TRIM(routine),TRIM(string))

    DO js=1,seed_size
       seed_array(js)=ABS(seed_trigger)+(js-1)
    ENDDO

    CALL message(TRIM(routine),'Seed generated')

    !-----------------------------------------------------------
    ! 3. generate random numbers and perturb the normal wind
    !-----------------------------------------------------------

    CALL RANDOM_SEED( PUT=seed_array )
    CALL RANDOM_NUMBER( zrand )

    zrand = (zrand - 0.5_wp)*pscale

    nblks = pt_patch%nblks_e
    npromz = pt_patch%npromz_e

    ! add the same random noise to all vertical layers
    DO jk=1,nlev
       pt_prog%vn(:,jk,1:nblks-1)    = pt_prog%vn(:,jk,1:nblks-1)    + zrand(:,1:nblks-1)
       pt_prog%vn(1:npromz,jk,nblks) = pt_prog%vn(1:npromz,jk,nblks) + zrand(1:npromz,nblks)
    ENDDO

    !-----------------------------------------------------------
    ! 4. clean up
    !-----------------------------------------------------------

    DEALLOCATE( seed_array, STAT=ist)
    IF(ist/=SUCCESS) CALL finish('random number:','deallocation of seed_array failed')
    CALL message(TRIM(routine),'=========================================')

  END SUBROUTINE hydro_state_prog_add_random
  !-------------
  !>
  !! @brief Add diabatic forcing to the dynamics tendency state.
  !! This is used in physics-dynamics coupling
  !!
  !! @par
  !! Note that in the boundary zone of a nested domain tendencies are interpolated
  !! from the parent domain. These values must not be overwritten. The width of the
  !! boundary zone is grf_bdywidth_c and grf_bdywidth_e for cell- and edge-based
  !! variables, respectively.
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2010-07-13)
  !!
  SUBROUTINE add_tend( p_patch, ltheta, ltracer, &! in
                     & p_tend_in, p_tend_acc )

    TYPE(t_patch),INTENT(IN) :: p_patch
    LOGICAL,INTENT(IN) :: ltheta
    LOGICAL,INTENT(IN) :: ltracer

    TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_tend_in
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend_acc

    INTEGER :: jb, jbs, is, ie
    INTEGER :: nblks_e, nblks_c
    !------
    ! Dimension parameters related to refinement and MPI parallelisation

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c

    ! Accumulate cell-based tendencies

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       IF (ltheta) THEN
         p_tend_acc% theta(is:ie,:,jb) =   p_tend_acc% theta(is:ie,:,jb) &
                                       & + p_tend_in%  theta(is:ie,:,jb)
       ELSE
         p_tend_acc% temp (is:ie,:,jb) =   p_tend_acc% temp(is:ie,:,jb)  &
                                       & + p_tend_in%  temp(is:ie,:,jb)
       ENDIF

       IF (ltracer) THEN
         p_tend_acc% tracer(is:ie,:,jb,:) =   p_tend_acc% tracer(is:ie,:,jb,:)  &
                                          & + p_tend_in%  tracer(is:ie,:,jb,:)
       ENDIF
     ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Accumulate wind tendency

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )
       p_tend_acc% vn(is:ie,:,jb) =   p_tend_acc% vn(is:ie,:,jb) &
                                  & + p_tend_in%  vn(is:ie,:,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE add_tend
  !-------------
  !>
  !! @brief Reset tendencies to zero.
  !!
  !! @par
  !! Note that in the boundary zone of a nested domain tendencies are interpolated
  !! from the parent domain. These values must not be overwritten. The width of the
  !! boundary zone is grf_bdywidth_c and grf_bdywidth_e for cell- and edge-based
  !! variables, respectively.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07-26)
  !!
  SUBROUTINE reset_tend( p_patch, lreset_theta, lreset_tracer, &! in
                       & p_tend                               ) ! inout

    TYPE(t_patch),INTENT(IN) :: p_patch
    LOGICAL,INTENT(in) :: lreset_theta
    LOGICAL,INTENT(in) :: lreset_tracer

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend

    INTEGER :: jb,jbs,jbe
    INTEGER :: jcs,jce

!$OMP PARALLEL

    ! Cell-based variables

    jbs = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
    jbe = p_patch%nblks_int_c

!$OMP DO PRIVATE(jb,jcs,jce)
    DO jb = jbs,jbe
      CALL get_indices_c( p_patch, jb,jbs,jbe, jcs,jce,grf_bdywidth_c+1 )

      IF (lreset_theta)  THEN
        p_tend% theta(jcs:jce,:,jb) = 0._wp
      ELSE
        p_tend% temp (jcs:jce,:,jb) = 0._wp
      ENDIF
      IF (lreset_tracer) p_tend% tracer(jcs:jce,:,jb,:) = 0._wp
    ENDDO
!$OMP END DO

    ! Edge-based variable

    jbs = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
    jbe = p_patch%nblks_int_e

!$OMP DO PRIVATE(jb,jcs,jce)
    DO jb = jbs,jbe
      CALL get_indices_e( p_patch, jb,jbs,jbe, jcs,jce,grf_bdywidth_e+1 )
      p_tend% vn(jcs:jce,:,jb) = 0._wp
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE reset_tend
  !-------------
  !>
  !! @brief Diagnose tendencies from the prognosic variables.
  !! This is used in physics-dynamics coupling
  !!
  !! @par
  !! Note that in the boundary zone of a nested domain tendencies 
  !! are interpolated from the parent domain. These values must 
  !! not be overwritten. The width of the boundary zone is 
  !! grf_bdywidth_c and grf_bdywidth_e for cell- and edge-based
  !! variables, respectively.
  !!
  !! @par Revision History
  !!  Original version by Hui Wan, MPI-M (2011-01-13)
  !!
  SUBROUTINE diag_tend( p_patch, ltheta, ltracer, pdtime, &! in
                      & p_prog_now, p_prog_new,           &! in
                      & p_tend                            )! inout

    TYPE(t_patch),INTENT(IN) :: p_patch
    LOGICAL,INTENT(IN) :: ltheta
    LOGICAL,INTENT(IN) :: ltracer
    REAL(wp) :: pdtime

    TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_prog_now
    TYPE(t_hydro_atm_prog),INTENT(IN)    :: p_prog_new
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_tend

    INTEGER :: jb, jbs, is, ie
    INTEGER :: nblks_e, nblks_c
    REAL(wp) :: zrdtime
    !------
    ! Reciprocal of time step

    zrdtime = 1._wp/pdtime

    ! Dimension parameters related to refinement and MPI parallelisation

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c

    ! Diagnose cell-based tendencies

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       IF (ltheta) THEN
         p_tend% theta(is:ie,:,jb) = ( p_prog_new% theta(is:ie,:,jb)   &
                                   &  -p_prog_now% theta(is:ie,:,jb) ) &
                                   & *zrdtime
       ELSE
         p_tend% temp (is:ie,:,jb) = ( p_prog_new% temp(is:ie,:,jb)   &
                                   &  -p_prog_now% temp(is:ie,:,jb) ) &
                                   & *zrdtime
       ENDIF

       IF (ltracer) THEN
         p_tend% tracer(is:ie,:,jb,:) = ( p_prog_new% tracer(is:ie,:,jb,:)   &
                                      &  -p_prog_now% tracer(is:ie,:,jb,:) ) &
                                      & *zrdtime
       ENDIF
    ENDDO
!$OMP END DO

    ! Diagnose wind tendency

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )
       p_tend% vn(is:ie,:,jb) = ( p_prog_new% vn(is:ie,:,jb)   &
                              &  -p_prog_now% vn(is:ie,:,jb) ) &
                              & *zrdtime
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE diag_tend

END MODULE mo_ha_prog_util
