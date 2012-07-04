!>
!! Implementation of physics utility routines.
!!
!! @par Revision History
!!  Initial revision  :  F. Prill, DWD (2012-07-03)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_util_phys

  USE mo_kind,      ONLY: wp
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 cpd, p0ref, rd
  USE mo_satad,                 ONLY: sat_pres_water
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag
  USE mo_run_config,            ONLY: iqv
  USE mo_loopindices,           ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rel_hum
  PUBLIC :: compute_field_rel_hum

  CHARACTER(len=*), PARAMETER :: version = &
    &  '$Id$'

CONTAINS

  !> POINTWISE computation of relative humidity as r=e/e_sat
  !!
  !! (domain independent and elemental)
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-03) 
  ELEMENTAL FUNCTION rel_hum(temp, qv, p_ex)
    REAL(wp) :: rel_hum
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    e_s = sat_pres_water(temp)
    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum = 100._wp * e/e_s

  END FUNCTION rel_hum


  !> computation of relative humidity as r=e/e_sat
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-04) 
  SUBROUTINE compute_field_rel_hum(ptr_patch, p_prog, p_diag, out_var, &
    &                              opt_slev, opt_elev, opt_rlstart, opt_rlend)

    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN) :: p_prog
    TYPE(t_nh_diag), INTENT(IN) :: p_diag
    ! output variable, dim: (nproma,nlev,nblks_c):
    REAL(wp),INTENT(INOUT) :: out_var(:,:,:)
    ! optional vertical start/end level:
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev, opt_elev
    ! start and end values of refin_ctrl flag:
    INTEGER, INTENT(in), OPTIONAL     :: opt_rlstart, opt_rlend
   
    ! local variables
    REAL(wp) :: temp, qv, p_ex
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, jk, jb

    ! default values
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(slev,     opt_slev)
    CALL assign_if_present(elev,     opt_elev)
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)
    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
          out_var(jc,jk,jb) = rel_hum(temp, qv, p_ex)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_rel_hum

END MODULE mo_util_phys
