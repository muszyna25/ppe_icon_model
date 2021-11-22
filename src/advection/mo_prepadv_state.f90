!>
!! Construction of a data object which is used to store mass fluxes 
!! and velocities that are used to drive the tracer transport schemes. 
!! In contrast to the mass fluxes in the dynamical core, these mass fluxes 
!! are averaged over the dynamics substeps.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2021-08-11)
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
MODULE mo_prepadv_state

  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,               ONLY: message, finish
  USE mo_model_domain,            ONLY: t_patch
  USE mo_prepadv_types,           ONLY: t_prepare_adv, t_step_adv
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_run_config,              ONLY: ntracer
  USE mo_var_list_register,       ONLY: vlr_add, vlr_del
  USE mo_var_list,                ONLY: add_var, t_var_list_ptr
  USE mo_zaxis_type,              ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_SURFACE
  USE mo_cdi,                     ONLY: DATATYPE_PACK16, GRID_UNSTRUCTURED, TSTEP_INSTANT, &
    &                                   DATATYPE_FLT32
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_CELL, GRID_EDGE
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var, grib2_var

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prep_adv
  PUBLIC :: prep_adv_list
  PUBLIC :: jstep_adv

  PUBLIC :: construct_prepadv_state
  PUBLIC :: destruct_prepadv_state


  ! state variable
  TYPE(t_prepare_adv), ALLOCATABLE :: prep_adv(:)  ! n_dom
  TYPE(t_step_adv),    ALLOCATABLE :: jstep_adv(:) ! n_dom

  ! variable list
  TYPE (t_var_list_ptr), ALLOCATABLE :: prep_adv_list(:)  ! n_dom

CONTAINS

  !>
  !! Constructor for prepadv state
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2021-08-11)
  !!
  SUBROUTINE construct_prepadv_state (p_patch)

    TYPE(t_patch), INTENT(in) :: p_patch(:)

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:construct_prepadv_state'


    ! Allocate pointer arrays prep_adv, as well as the corresponding list arrays.
    !   
    ALLOCATE(prep_adv(n_dom), prep_adv_list(n_dom),STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'allocation of prep_adv array and list failed')
    ENDIF

    !$ACC ENTER DATA COPYIN(prep_adv)

    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'prepadv_of_domain_',jg
      CALL new_prep_adv_list( p_patch(jg), listname, prep_adv_list(jg), prep_adv(jg))
    ENDDO

    CALL message(routine, 'construction of prep_adv state finished')

  END SUBROUTINE construct_prepadv_state



  !>
  !! Destructor for prepadv state
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2021-08-11)
  !!
  SUBROUTINE destruct_prepadv_state ()

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:destruct_prepadv_state'

    !--------------------------------------------------------------

    ! delete prep_adv varlist
    DO jg = 1, n_dom
      CALL vlr_del(prep_adv_list(jg))
    ENDDO

    DEALLOCATE(prep_adv, prep_adv_list, STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'deallocation of prep_adv array and list failed')
    ENDIF

    !$ACC EXIT DATA DELETE(prep_adv)
 
    CALL message(routine, 'destruction of prep_adv state finished')

  END SUBROUTINE destruct_prepadv_state



  !>
  !! Constructor for prepadv state
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2021-08-11)
  !!
  SUBROUTINE new_prep_adv_list (p_patch, listname, prep_adv_list, prep_adv)

    TYPE(t_patch)         , INTENT(IN   ) :: p_patch
    CHARACTER(len=*)      , INTENT(IN   ) :: listname
    TYPE(t_var_list_ptr)  , INTENT(INOUT) :: prep_adv_list
    TYPE(t_prepare_adv)   , INTENT(INOUT) :: prep_adv

    ! local variables
    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_e(3), shape3d_chalf(3), shape3d_tracer(3)
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt

    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:new_prep_adv_list'

    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits        = DATATYPE_PACK16   ! "entropy" of horizontal slice
    datatype_flt = DATATYPE_FLT32 

    shape3d_e     = (/nproma, nlev   , nblks_e/)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c/)
    shape3d_tracer= (/nproma, nblks_c, MAX(1,ntracer)/)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(prep_adv%mass_flx_me, &
      &     prep_adv%mass_flx_ic, &
      &     prep_adv%vn_traj,     &
      &     prep_adv%topflx_tra   )

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(prep_adv_list, TRIM(listname), patch_id=p_patch%id, lrestart=.FALSE.)


    ! mass_flx_me      prep_adv%mass_flx_me(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('mass_flx_me', 'kg m-2 s-1', &
      &                   'horizontal mass flux (averaged over dynamics substeps)', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( prep_adv_list, 'mass_flx_me', prep_adv%mass_flx_me,              &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_e, loutput=.FALSE.,                                &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%mass_flx_me)


    ! mass_flx_ic      prep_adv%mass_flx_ic(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('mass_flx_ic', 'kg m-2 s-1', &
      &                   'vertical mass flux (averaged over dynamics substeps)',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'mass_flx_ic', prep_adv%mass_flx_ic,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, loutput=.FALSE.,                            &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%mass_flx_ic)


    ! vn_traj          prep_adv%vn_traj(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('vn_traj', 'm s-1', &
      &                   'velocity normal to edge (averaged over dynamics substeps)',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( prep_adv_list, 'vn_traj', prep_adv%vn_traj,                      &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_e, loutput=.FALSE.,                                &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%vn_traj)


    ! topflx_tra       prep_adv%topflx_tra(nproma,nblks_c,ntracer)
    cf_desc    = t_cf_var('topflx_tra', 'kg m-2 s-1', &
      &                   'tracer mass flux at upper boundary',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'topflx_tra', prep_adv%topflx_tra,                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
                & ldims=shape3d_tracer, loutput=.FALSE.,                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%topflx_tra)

  END SUBROUTINE new_prep_adv_list


END MODULE mo_prepadv_state

