!>
!!                 Contains global variables used by the shallow water model.
!!
!!                 More comments on many of the variable defined here
!!                 to be found in the <i>user_introduction</i>.
!!
!! @par Revision History
!!  Initial version  by Luca Bonaventura (2005).
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Thomas Heinze (2006-05-04):
!!  - changed nerr to nerror to avoid conflict with nerr: error output unit
!!  Modification by Peter Korn, MPI-M (2007-02):
!!  - ocean parameters introduced.
!!  - topography_selection introduced.
!!  Modification by Hui Wan, MPI-M (2007-02-23):
!!  Reorganization:
!!  - declaration of namelists moved to here from <i>mo_io_utilities</i>.
!!  - changed the names of the namelists.
!!  - subroutine <i>setup_control</i> and <i>setup_dyn</i> introduced
!!    to parse namelist <i>run_ctl</i> and <i>dyn_ctl, ocean_ctl</i>.
!!    The contents were parts of subroutine <i>init_dyn</i> and
!!    <i>check_initialization</i> in <i>mo_io_utilities</i>.
!!  Modification by Hui Wan, MPI-M (2007-11)
!!  - added additional namelist variables for additional smoothing to
!!    divergence operator.
!!  - Comment: CFL number of horizontal diffusion needs to be checked
!!    for the three-time-level schemes.
!!  Modification by Hui Wan, MPI-M (2008-02-04)
!!  - added namelist variable <i>dt_hdiff</i> which is the user-specified
!!    time step for horizontal diffusion
!!  Modification by Hui Wan, MPI-M (2008-04-04)
!!  - loce renamed locean
!!  - topography_selection renamed itopo
!!  - default value of nproma set to 1
!!  - dt renamed dtime
!!  - modtype renamed itime_scheme
!!  - epsass renamed asselin_coeff
!!  - difftype renamed hdiff_order
!!  - modtype_hdiff renamed hdiff_time_scheme
!!  - efold renamed hdiff_efold
!!  Modification by Almut Gassman, MPI-M (2008-09-19)
!!  - switch lshallow_water introduced with intention to have only one
!!    executable for atmosphere and shallow water
!!  - switch lidealized for using artificial model generated initial data
!!  - switch for running options as triangular or hexagonal model
!!  Modification by A. Gassmann, MPI-M (2008-09-23)
!!  - clean namelists and removed not necessary variables
!!  - rename this subroutine in mo_global_variables
!!    (no reference to shallow water any longer)
!!  - rename theta in si_2tls
!!  Modification by A. Gassmann, MPI-M (2008-09-30)
!!  - inserted Namelists originally found at hydro_control
!!  - split out a Namelist for horizontal diffusion
!!  Modification by M. Giorgetta, MPI-M (2009-02-23)
!!  - lidealized renamed to ltestcase
!!  Modification by M. Giorgetta, MPI-M (2009-02-25)
!!  - revision of &run_ctl
!!    - rename ltracer to ltransport
!!    - rename tracer_ctl to transport_ctl
!!    - added ldynamics
!!    - rename dyn_ctl to dynamics_ctl
!!    - renamed setup_dyn to setup_dynamics
!!    - added lhydrostatic
!!    - renamed hydro_ctl to hydrostatic_ctl
!!    - renamed setup_hydro to setup_hydrostatic
!!    - moved nlev from namelist hydrostatic_ctl to namelist run_ctl
!!    - removed nlev_ocean from ocean_ctl
!!    - added iequations
!!    - added latmosphere
!!  Modification by Almut Gassmann, MPI-M (2009-03-05)
!!  - start preparation for NH Model
!!  Modification by M. Giorgetta, MPI-M (2009-03-22)
!!  - add parameters 'modelname' and 'modelversion'
!!  Modification by Almut Gassmann, MPI-M (2009-10-05)
!!  - remove luse_rbf_vec_int.. switches
!!  Modification by Stephan Lorenz, MPI-M (2010-02-23)
!!  - extension of namelist ocean_ctl and default vertical z-levels
!!  Modification by Daniel Reinert, DWD (2010-06-07)
!!  - added variables iqv, iqc, iqr, iqi, iqs which can be used to
!!    acess the desired fields in the 4D tracer array
!!  Modification by Kristina Froehlich, DWD (2010-09-22)
!!  - added variables itconv, itccov, itrad, itsso plus the respective field
!!    to  introduce the calling-time control array
!!  Modification by Constantin Junk, MPI-M (2011-02-22)
!!  - separated run_ctl and the subroutine setup_run
!!    into new module 'mo_run_nml'.
!!  Modification by Constantin Junk, MPI-M (2011-02-24)
!!  - separated io_ctl and the subroutine setup_io
!!    into new module 'mo_io_nml'.
!!  Modification by Constantin Junk, MPI-M (2011-03-18)
!!  - separated mpiom_phy_ctl, ocean_ctl und octst_ctl and
!!    the subroutine setup_ocean into new module 'mo_ocean_nml'
!!  Modification by Constantin Junk, MPI-M (2011-03-23)
!!  - removed hydrostatic_ctl and setup_hydrostatic. Included
!!    the former hydrostatic_ctl variables into mo_dynamics_nml and
!!    the consistency checks of setup_hydrostatic into setup_dynamics_nml
!!  Modification by Constantin Junk, MPI-M (2011-03-28)
!!  - separated nonhydrostatic_ctl and the subroutine setup_nonhydrostatic
!!    into new module 'mo_nonhydrostatic_nml'
!!  Modification by Constantin Junk, MPI-M (2011-03-29)
!!  - separated sleve_ctl and the subroutine setup_sleve
!!    into new module 'mo_sleve_nml'
!!
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
!!
MODULE mo_global_variables
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_run_nml,            ONLY: inwp,iecham, iforcing,lshallow_water,iequations,    &
    &                              inoforcing,iheldsuarez,impiom, ildf_dry,ildf_echam
  USE mo_ocean_nml,          ONLY: lmpiom_radiation, lmpiom_convection,           &
    &                              lmpiom_gentmcwill, mpiom_phy_ctl
  USE mo_dynamics_nml,       ONLY: ldry_dycore

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'




!--------------------------------------------------------------------

  CONTAINS


 SUBROUTINE setup_physics
!!-----------------------------------------------------------------------
!! read namelist for physics
!! (chose the physical package and subsequent parameters)
!!-----------------------------------------------------------------------
!!
  INTEGER :: i_status

  SELECT CASE (iforcing)
    !
  CASE (inoforcing)
    !
    ! nothing to be done
    !
  CASE (iheldsuarez, ildf_dry)
    !
    ldry_dycore       = .TRUE.
    !
  CASE (iecham,ildf_echam)
    !
! Temporarily removed by Hui
!   ldry_dycore       = .FALSE.
    !
  CASE (inwp)
    !
! Temporarily removed by Daniel
!    ldry_dycore     = .FALSE.

!    !> set default physics switches and values
!    CALL set_inwp_nml

!    !
!    !> final settings via namelist
!    CALL read_inwp_nml

    !
  CASE (impiom)
    !
    lmpiom_radiation  = .FALSE.
    lmpiom_convection = .FALSE.
    lmpiom_gentmcwill = .FALSE.
    !
    CALL position_nml ('mpiom_phy_ctl', status=i_status)
    !
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, mpiom_phy_ctl)
    END SELECT
    !  write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=mpiom_phy_ctl)
    !
   CASE DEFAULT
    !
  END SELECT
!
END SUBROUTINE setup_physics

!-------------------------------------------------------------------------
!
END MODULE mo_global_variables
