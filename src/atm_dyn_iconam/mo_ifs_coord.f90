!>
!! Contains subroutines for calculating the pressure values and
!! some other auxiliary variables related to the pressure-sigma
!! hybrid vertical coordinate (the eta coordinate).
!!
!! @par Revision History
!!   Original version from ECHAM5.3.01
!!   Adapted for ICOHDC by Hui Wan, 2006-02
!!   Modifications include:
!!   - Calculation of half-level geopotential added to *geopot*
!!   - Calculation of logorithm of half-level pressure added to *auxhyb*
!!   - Subroutine <i>full_level_pressure</i> added
!!   Modifications by Hui Wan (MPI-M, 2010-01-29)
!!   - Renamed subroutines and variables.
!!   - Changed the sequence of items in the argument lists.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ifs_coord

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message_text, message
  USE mo_impl_constants,     ONLY: success
  USE mo_physical_constants, ONLY: rd
  USE mo_fortran_tools,      ONLY: DO_DEALLOCATE
  USE mo_mpi,                ONLY: p_bcast, p_comm_rank
  USE mo_read_interface,     ONLY: nf

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  ! Data type containing all data specifying a pressure-sigma hybrid
  ! vertical coordinate (the eta coordinate).
  !
  TYPE :: t_vct
    INTEGER :: nlevm1         ! (number of levels)-1.
    INTEGER :: nplev          ! *number of pressure levels.
    INTEGER :: nplvp1         ! *nplev+1.*
    INTEGER :: nplvp2         ! *nplev+2.*
    INTEGER :: nlmsgl         ! *nlev* - (number of sigma levels).
    INTEGER :: nlmslp         ! *nlmsgl+1.*
    INTEGER :: nlmsla         ! *nlmslp,* or 2 if *nlmslp=1.*
    
    REAL(wp), ALLOCATABLE :: ralpha(:) ! rd*alpha at pressure and sigma levels.
    REAL(wp), ALLOCATABLE :: rlnpr(:)  ! rd*ln(p(k+.5)/p(k-.5))
    REAL(wp), ALLOCATABLE :: delpr(:)  ! p(k+.5)-p(k-.5) of
                                       ! the reference surface pressure.
    REAL(wp), ALLOCATABLE :: rdelpr(:) ! *reciprocal of *delpr.*
    
    REAL(wp), ALLOCATABLE :: vct  (:) ! param. A and B of the vertical coordinate

  CONTAINS
    PROCEDURE :: t_vct_construct_with_arrays
    PROCEDURE :: t_vct_construct_ncfile
    GENERIC :: construct => t_vct_construct_with_arrays, t_vct_construct_ncfile

    PROCEDURE :: finalize            => t_vct_finalize
    PROCEDURE :: half_level_pressure => t_vct_half_level_pressure
    PROCEDURE :: full_level_pressure => t_vct_full_level_pressure
    PROCEDURE :: auxhyb              => t_vct_auxhyb
   
  END TYPE t_vct


  PUBLIC :: t_vct
  PUBLIC :: geopot


  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ifs_coord'

CONTAINS

  ! ----------------------------------------------------------------------
  !> Allocates and initializes vertical coordinate data structure
  !> based on a NetCDF file.
  !
  !  Note: This call is MPI-collective: only the first PE reads the
  !  data from file and broadcasts it to the other PEs.
  !
  SUBROUTINE t_vct_construct_ncfile(vct, ncid, p_io, mpi_comm)
    CLASS(t_vct), INTENT(INOUT) :: vct
    INTEGER, INTENT(IN) :: ncid           !< NetCDF file ID (file already opened)
    INTEGER, INTENT(IN) :: p_io           !< MPI rank of reading process
    INTEGER, INTENT(IN) :: mpi_comm       !< MPI communicator

    ! local variables
    CHARACTER(LEN=*),PARAMETER :: routine =  modname//'::t_vct_construct_ncfile'
    REAL(wp), ALLOCATABLE :: vct_ab(:,:) ! param. A and B of the vertical coordinte
    REAL(wp), ALLOCATABLE :: lev_ifs(:)
    INTEGER,  ALLOCATABLE :: lev_hyi(:)
    REAL(wp), ALLOCATABLE :: hyab(:)
    INTEGER  :: varid, nlev_in, dimid, nhyi, ierrstat, var_ndims, var_dimlen(NF_MAX_VAR_DIMS), &
      &         var_dimids(NF_MAX_VAR_DIMS), i
    LOGICAL  :: lread_process  !< .TRUE. on the reading PE

    lread_process = (p_comm_rank(mpi_comm) == p_io)

    IF (lread_process) THEN
      CALL nf(nf_inq_varid(ncid, 'lev', varid), routine)
      ! retrieve number of levels
      CALL nf(nf_inq_varndims(ncid, varid, var_ndims), routine)
      CALL nf(nf_inq_vardimid(ncid, varid, var_dimids), routine)
      DO i = 1, var_ndims
        CALL nf(nf_inq_dimlen (ncid, var_dimids(i), var_dimlen(i)), routine)
      END DO
      nlev_in = var_dimlen(1)
    END IF

    CALL p_bcast(nlev_in, p_io, mpi_comm)

    ALLOCATE( vct_ab(nlev_in+1, 2), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    IF (lread_process) THEN
      CALL nf(nf_inq_dimid(ncid, 'nhyi', dimid), routine)
      CALL nf(nf_inq_dimlen(ncid, dimid, nhyi), routine)

      ALLOCATE( lev_ifs(nlev_in), lev_hyi(nlev_in+1), hyab(nhyi), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      CALL nf(nf_get_var_double(ncid, varid, lev_ifs), routine)
      lev_hyi(1:nlev_in) = NINT( lev_ifs(:) )
      lev_hyi(nlev_in+1) = lev_hyi(nlev_in) + 1
      IF ( nlev_in+1 /= nhyi) THEN
        WRITE(message_text,*) 'Reading only IFS levels ', lev_hyi(1:nlev_in)
        CALL message(routine, TRIM(message_text))
      END IF

      CALL nf(nf_inq_varid(ncid, 'hyai', varid), routine)
      CALL nf(nf_get_var_double(ncid, varid, hyab), routine)
      vct_ab(:,1) = hyab( lev_hyi(:))

      CALL nf(nf_inq_varid(ncid, 'hybi', varid), routine)
      CALL nf(nf_get_var_double(ncid, varid, hyab), routine)
      vct_ab(:,2) = hyab( lev_hyi(:))
    ENDIF

    CALL p_bcast(vct_ab, p_io, mpi_comm)

    ! now call constructor based on arrays:
    CALL vct%construct(nlev_in, vct_ab(:,1), vct_ab(:,2))

    DEALLOCATE(vct_ab, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE t_vct_construct_ncfile


  ! ----------------------------------------------------------------------
  !> Allocates and initializes vertical coordinate data structure
  !> based on "vct_a", "vct_b" arrays.
  !
  SUBROUTINE t_vct_construct_with_arrays(vct, nlev, vct_a, vct_b)
    CLASS(t_vct), INTENT(INOUT) :: vct
    REAL(wp), INTENT(IN) :: vct_a(:) ! param. A of the vertical coordinte
    REAL(wp), INTENT(IN) :: vct_b(:) ! param. B of the vertical coordinate

    INTEGER, INTENT(IN) :: nlev

    INTEGER :: nlevp1

    CHARACTER(LEN=*),PARAMETER :: routine =  modname//'::t_vct_construct_with_arrays'

    !  Local scalars:
    REAL(wp) :: za, zb, zp, zpp, zrd, zs, zsm
    INTEGER  :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, &
      &         nvclev

    !  Intrinsic functions
    INTRINSIC LOG

    nlevp1 = nlev+1

    ALLOCATE( vct%vct(2*(nlevp1)), STAT=ist)
    IF (ist/=success)  CALL finish(routine, 'ALLOCATE failed!') 

    vct%vct(1:nlevp1)        = vct_a(:)
    vct%vct(nlev+2:2*nlevp1) = vct_b(:)

    ALLOCATE (vct%ralpha(nlev), vct%rlnpr(nlev),  &
      &       vct%delpr(nlev),   vct%rdelpr(nlev), STAT=ist)
    IF (ist/=success)  CALL finish(routine, 'ALLOCATE failed!') 


    ! --------------------------------------------------------------------------------

    !-- 1. Initialize variables
    nvclev = nlev+1
    zrd       = rd
    vct%ralpha(1) = zrd*LOG(2._wp)
    vct%rlnpr(1)  = 2._wp*vct%ralpha(1)
    ilev      = nlev
    ilevp1    = ilev + 1
    vct%nlevm1    = ilev - 1
    iplev     = 0
    iplvp1    = 1
    is        = nvclev + ilevp1
    ism       = is - 1
    zpp       = vct%vct(1)
    zsm       = vct%vct(is)

    zb      = vct%vct(nvclev+iplvp1+1)

    !-- 2. Calculate pressure-level values

    DO WHILE ( zb == 0._wp )

      iplev  = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) EXIT    ! if all levels are pressure levels

      zp            = zpp
      zpp           = vct%vct(iplvp1)
      vct%delpr(iplev)  = zpp - zp
      vct%rdelpr(iplev) = 1._wp/vct%delpr(iplev)

      IF ( iplev>1 ) THEN
        vct%rlnpr(iplev)  = zrd*LOG(zpp/zp)
        vct%ralpha(iplev) = zrd - zp*vct%rlnpr(iplev)/vct%delpr(iplev)
      END IF

      zb            = vct%vct(nvclev+iplvp1+1)

    ENDDO


    IF (iplvp1/=ilevp1) THEN   ! All levels are not pressure-levels

      vct%nplev  = iplev
      vct%nplvp1 = iplvp1
      vct%nplvp2 = iplvp1 + 1


      !-- 3. Calculate sigma-level values

      za = vct%vct(ism-nvclev)

      DO WHILE ( za == 0._wp )

        is  = ism
        ism = is - 1
        ist = is - nvclev
        zs  = zsm
        zsm = vct%vct(is)
        IF (ist==1) THEN
          vct%nlmsgl = 0
          vct%nlmslp = 1
          vct%nlmsla = 2
          EXIT
        ELSE
          vct%rlnpr(ist)  = zrd*LOG(zs/zsm)
          vct%ralpha(ist) = zrd - zsm*vct%rlnpr(ist)/(zs-zsm)
        END IF
        za = vct%vct(ism-nvclev)

      END DO

      IF (za>0._wp) THEN
        vct%nlmsgl = ism - nvclev
        vct%nlmslp = vct%nlmsgl + 1
        vct%nlmsla = vct%nlmslp
      END IF

    ENDIF  ! If all levels are not pressure-levels

  END SUBROUTINE t_vct_construct_with_arrays


  SUBROUTINE t_vct_finalize(vct)
    CLASS(t_vct), INTENT(INOUT) :: vct

    CHARACTER(LEN=*),PARAMETER :: routine =  modname//'::t_vct_finalize'

    CALL DO_DEALLOCATE(vct%vct)
    CALL DO_DEALLOCATE(vct%ralpha)
    CALL DO_DEALLOCATE(vct%rlnpr)
    CALL DO_DEALLOCATE(vct%delpr)
    CALL DO_DEALLOCATE(vct%rdelpr)

  END SUBROUTINE t_vct_finalize


  !-------------------------------------------------------------------------
  !>
  !! Calculate half-level pressures at all model levels
  !! for a given surface pressure.
  !!
  !! @par Method
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !! @par Arguments
  !!   *ps*        surface pressure.
  !!   *kdimp*     first dimension of 2-d array *ph.*
  !!   *klen*      number of points for which calculation is
  !!               performed.
  !!   *ph*        computed half-level pressures.
  !!
  !! @par Parameters
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialized by a call of
  !!  subroutine init_vertical_coord.
  !!
  !! @par Results
  !!  Results are computed for *klen* consecutive points at
  !!  each model half level.
  !!
  !! @see
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI, Feb 2006, adapted for ICOHDC
  !!    H. Wan, MPI, Jan 2010, renamed the interface
  !!
  SUBROUTINE t_vct_half_level_pressure(vct, ps,kdimp,klen,nlev_in, ph)
    CLASS(t_vct), INTENT(IN) :: vct

    INTEGER ,INTENT(in)  :: kdimp
    REAL(wp),INTENT(in)  :: ps(:)   !< surface pressure
    INTEGER ,INTENT(in)  :: klen, nlev_in

    REAL(wp),INTENT(inout) :: ph(kdimp,nlev_in+1) !< half-level pressure

    CHARACTER(LEN=*),PARAMETER :: routine =  modname//'::t_vct_half_level_pressure'
    REAL(wp):: zb, zp
    INTEGER :: jk, jl, nvclev, nlevp1

    ! consistency checks:

    IF (SIZE(ph,1) < klen) THEN
      WRITE (message_text,*) "Wrong dimension: SIZE(ph,1)=",SIZE(ph,1), ", klen=",klen
      CALL finish(routine, message_text)
    END IF
    IF (SIZE(ps) /= kdimp) THEN
      WRITE (message_text,*) "Wrong dimension: SIZE(ps)=",SIZE(ps), ", kdimp=",kdimp
      CALL finish(routine, message_text)
    END IF
    IF (.NOT. ALLOCATED(vct%vct)) THEN
      CALL finish(routine, "vct%vct not allocated!")
    END IF
    IF (SIZE(vct%vct) < vct%nlmsgl) THEN
      WRITE (message_text,*) "Wrong dimension: SIZE(vct%vct)=",SIZE(vct%vct),"; vct%nlmsgl=",vct%nlmsgl
      CALL finish(routine, message_text)
    END IF

    nvclev = nlev_in+1
    nlevp1 = nvclev

    ! Transfer pressure level values

    DO jk = 1, vct%nplvp1
      zp = vct%vct(jk)
      DO jl = 1, klen
        ph(jl,jk) = zp
      END DO
    END DO

    ! Compute hybrid level values

    DO jk = vct%nplvp2, vct%nlmsgl
      zp = vct%vct(jk)
      zb = vct%vct(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zp + zb*ps(jl)
      END DO
    END DO

    ! Compute sigma-level values

    DO jk = vct%nlmslp, nlevp1
      zb = vct%vct(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zb*ps(jl)
      END DO
    END DO

  END SUBROUTINE t_vct_half_level_pressure


  !>
  !! Calculate full-level pressures for all vertical layers
  !! Method: Simmons and Burridge (Mon.Wea.Rev.,1981,p761,Eqn.(3.18))
  !!
  !! @par Parameters
  !!    *pres_i*    half-level pressure values.
  !!    *pres_m*    computed full-level pressure values.
  !!    *kdimp*     first dimension of 2-d array *pres_i*
  !!    *klen*      number of points for which calculation is
  !!                performed.
  !!
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialiazed
  !!  by a call of subroutine *inihyb*.
  !!
  !! @par Revision History
  !!    H. Wan, MPI, 2006-08-17
  !!
  SUBROUTINE t_vct_full_level_pressure( vct, pres_i, kdimp, klen, nlev_in, pres_m)
    CLASS(t_vct), INTENT(IN) :: vct

    INTEGER ,INTENT(in) :: kdimp, klen, nlev_in    !< dimension parameters
    REAL(wp),INTENT(in) :: pres_i(kdimp,nlev_in+1) !< half-level pressure

    REAL(wp),INTENT(inout) :: pres_m(kdimp,nlev_in) !< full(/mid)-level pressure

    REAL(wp):: ztmp, zpres_i_top_min
    INTEGER :: jk, jl, ikp, ik_top, nlev

    !-----

    nlev = nlev_in
    zpres_i_top_min = vct%vct(1)
    IF ( zpres_i_top_min > 0._wp ) THEN
      ik_top = 1
    ELSE
      ik_top = 2
      pres_m(1:klen,1) = pres_i(1:klen,2)*0.5_wp
    END IF

    DO jk = ik_top, nlev
       ikp = jk+1
       DO jl = 1, klen
         ztmp = ( pres_i(jl,ikp)*LOG(pres_i(jl,ikp))   &
         &       -pres_i(jl,jk )*LOG(pres_i(jl,jk )) ) &
         &     /( pres_i(jl,ikp)-pres_i(jl,jk) )
         pres_m(jl,jk) = EXP(ztmp-1._wp)
       END DO
    END DO

  END SUBROUTINE t_vct_full_level_pressure


  !-------------------------------------------------------------------------
  !>
  !! Calculates auxiliary variables connected with
  !! the vertical finite-difference scheme.
  !!
  !! @par Arguments
  !!
  !!   *ph*          *specified half-level pressures.
  !!   *kdim*        *first dimension of 2-d arrays *pdelp,*
  !!                  *plnpr,* *palpha,* and *ph.*
  !!   *klen*        *number of points for which calculation is performed.
  !!   *pdelp*       *computed pressure difference across layers.
  !!   *prdelp*      *reciprocal of *pdelp.*
  !!   *plnpr*       *computed logarithm of ratio of pressures.
  !!   *palpha*      *computed alphas for integration of the
  !!                  hydrostatic equation and related terms.
  !!
  !!  Required constants are obtained from modules
  !!  *mo_constants* and *mo_hyb*. The latter should have been
  !!  initialized by a call of subroutine *inihyb*.
  !!
  !!  Results are computed for *klen* consecutive points for
  !!  all required levels.
  !!
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation can be found
  !!  in MPI-M technical report No. 349
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!
  SUBROUTINE t_vct_auxhyb( vct, ph,kdim,klen,nlev_in,                   &
    &                      pdelp,prdelp,plnph,plnpr,palpha )
    
    CLASS(t_vct), INTENT(IN) :: vct

    INTEGER ,INTENT(in)  :: kdim, klen, nlev_in
    REAL(wp),INTENT(in)  :: ph(kdim, nlev_in+1)

    REAL(wp), INTENT(inout) :: pdelp (kdim, nlev_in  ), prdelp(kdim, nlev_in)
    REAL(wp), INTENT(inout) :: plnph (kdim, nlev_in+1), plnpr (kdim, nlev_in)
    REAL(wp), INTENT(inout) :: palpha(kdim, nlev_in  )

    REAL(wp) :: za, zd, zl, zr
    INTEGER  :: jk, jl, nlev, nlevp1

    nlev = nlev_in
    nlevp1 = nlev+1

    !-----
    ! Set pressure-level values or other top-level values

    DO jk = 1, vct%nplev
      zd = vct%delpr(jk)
      zr = vct%rdelpr(jk)
      zl = vct%rlnpr(jk)
      za = vct%ralpha(jk)
      DO jl = 1, klen
        pdelp(jl,jk) = zd
        prdelp(jl,jk) = zr
        plnpr(jl,jk) = zl
        palpha(jl,jk) = za
      END DO
    END DO

    ! Calculate hybrid-level values

    DO jk = vct%nplvp1, vct%nlmsgl

      IF (jk == 1 .AND. vct%vct(1) == 0.0_wp ) THEN

        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          palpha(jl,jk) = rd*LOG(2._wp)
          plnpr(jl,jk)  = 2._wp*palpha(jl,jk)
        END DO

      ELSE
        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          plnpr(jl,jk)  = rd*LOG(ph(jl,jk+1)/ph(jl,jk))
          palpha(jl,jk) = rd - ph(jl,jk)*plnpr(jl,jk)*prdelp(jl,jk)
        END DO

      ENDIF
    END DO

    ! Set sigma-level values

    DO jk = vct%nlmsla, nlev
      zl = vct%rlnpr(jk)
      za = vct%ralpha(jk)
      DO jl = 1, klen
        pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
        prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
        plnpr(jl,jk)  = zl
        palpha(jl,jk) = za
      END DO
    END DO

    DO jk = 2,nlevp1
      DO jl = 1, klen
        plnph(jl,jk) = LOG(ph(jl,jk))
      END DO
    END DO

  END SUBROUTINE t_vct_auxhyb


  !-------------------------------------------------------------------------
  !>
  !! Calculates full- and half-level geopotential
  !!
  !! Method: Integrate the hydrostatic equation in the vertical
  !! to obtain full- or half-level values of geopotential.
  !!
  !! *geopot* is called during the calculation of adiabatic
  !! tendencies, prior to the calculation of the physical paramet-
  !! erizations, and in the post-processing.
  !!
  !! Parameters are
  !!  *pgeop_sfc*   *surface geopotential.
  !!  *pgeop_m*     *computed geopotential on full-levels
  !!  *pgeop_i*     *computed geopotential on half-levels
  !!  *ptv*         *virtual temperature.
  !!  *plnpr*       *logarithm of ratio of pressures, computed by *auxhyb*.
  !!  *palpha*      *for full-level values use *alpha* as computed by *auxhyb*.
  !!  *kdim*        *first dimension of 2-d arrays *phi,*
  !!                 *ptv,* *plnpr,* and *palpha.*
  !!  *kstart,kend* *start and end points for which calculation is performed.
  !!
  !! Required constants are obtained from module *mo_hyb*.
  !! The latter should have been initialized by a call of subroutine *inihyb*.
  !!
  !! Results are computed for *klen* consecutive points for *nlev* levels.
  !!
  !! The choice of full- or half-level values is determined
  !! by the specification of the input array *alpha.*
  !!
  !! External documentation of the model equations and the
  !! organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, January 1982, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI-M, 2006-08-10, included in m_vertical
  !!    H. Wan, MPI-M, 2006-08-15, half-level geopotential added to the output
  !!    H. Wan, MPI-M, 2007-07-19, calculation of the full-level geopotential
  !!                               simplified.
  !!    G. Zaengl, DWD, 2009-07-06, replace klen with kstart/kend for correct
  !!                                execution on nested domains
  !!    A. Seifert, DWD, 2010-06-21, add missing lower boundary by restructuring
  !!                                 the loop for half and full levels

  SUBROUTINE geopot( ptv,plnpr,palpha,pgeop_sfc,kdim,kstart,kend,nlev_in, &
                     pgeop_m, pgeop_i )

    INTEGER ,INTENT(in) :: kdim, kstart, kend, nlev_in

    REAL(wp) ,INTENT(in)    :: ptv   (kdim,nlev_in),     plnpr(kdim,nlev_in)
    REAL(wp) ,INTENT(in)    :: palpha(kdim,nlev_in), pgeop_sfc(kdim)
    REAL(wp) ,INTENT(inout) :: pgeop_m(kdim,nlev_in),   pgeop_i(kdim,nlev_in+1)

    INTEGER :: jk, jl, jkp, nlev, nlevp1, nlevm1

    nlev = nlev_in
    nlevp1 = nlev+1
    nlevm1 = nlev-1

    ! Integrate hydrostatic equation

    DO jl = kstart, kend
      pgeop_i(jl,nlevp1) = pgeop_sfc(jl)
      pgeop_m(jl,nlev)   = palpha(jl,nlev)*ptv(jl,nlev) + pgeop_sfc(jl)
    END DO

    ! half levels
    DO jk = nlev, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_i(jl,jk) = plnpr(jl,jk)*ptv(jl,jk) + pgeop_i(jl,jkp)
      END DO
    END DO

    ! full levels
    DO jk = nlevm1, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_m(jl,jk)  = pgeop_i(jl,jkp) + palpha(jl,jk)*ptv(jl,jk)
      END DO
    END DO

  END SUBROUTINE geopot


END MODULE mo_ifs_coord

