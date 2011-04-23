!>
!!  Module contains some constants relevant for implementational issues.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn (2005)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Peter Korn (2006-08):
!!  - added identifier for edges, cells, vertices (used for memory allocation)
!! @par
!!  Modification by Peter Korn (2006-12):
!!   - added identifier for variables (used for land-sea mask application)
!!  Modification by Hui Wan (2007-02):
!!  - added "nproma" and some comments
!!  Modification by Guenther Zaengl (2008-10-23):
!!  - added parameters defining the range of the refin_ctrl flags
!!    (these are also set in the patch generator and must match each other)
!! @par
!!  Modification by Hui Wan (2009-08-07)
!!  - added identifiers for time stepping methods.
!!  Modifications by Daniel Reinert (2010-10-06)
!!  - added identifiers for NWP physics time control
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
MODULE mo_impl_constants
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  INTEGER, PARAMETER :: MAX_CHAR_LENGTH     = 1024

  INTEGER, PARAMETER :: SUCCESS             = 0
  INTEGER,PARAMETER  :: CELLS               = 123
  INTEGER,PARAMETER  :: EDGES               = 345
  INTEGER,PARAMETER  :: VERTS               = 678

!-------------------------------------------------------------------------------
! Comments by Hui:
! According to Luis' explanation, the declarations above related to the blocking
! are not correct. There should be a single NPROMA, and different NBLKS values
! for edges, cells and vertices.
!
! Take the triangular cells for example. Considering that the number of cells is
! different from patch to patch, it may be a good idea NOT to declare "nblks_c" as
! a 2D array HERE, but as a component of the "patch" type. Since we want to do
! calculations ONLY for the internal cells, we also need to know how many blocks
! (e.g. "nblks_i_c") we have for the internal cells. Furthermore, it is also
! necessary to define a a parameter "npromz_i_c" for each patch, which stores the
! number of valid items in the last block of the internal cells.
!   At the very beginning of the execution of the program, "nproma" can be read
! from the namelist file. Then, after reading in the number of cells and external halo
! cells of a certain patch, we calculate the value of "nblks_c", "nblks_i_c"
! and "npromz_i_c" for that patch by:
!
!   n_patch_cell_all     = ptr_patch%ncells + ptr_patch%n_e_halo_cells
!   ptr_patch%nblks_c    = ( n_patch_cell_all - 1 )/nproma + 1
!
!   ptr_patch%nblks_i_c  = ( ptr_patch%ncell - 1 )/nproma + 1
!   ptr_patch%npromz_i_c = ptr_patch%ncell - (ptr_patch%nblks_i_c - 1)*nproma
!
!   (cf. echam5, mo_decomposition.f90,
!        the last few lines of SUBROUTINE decompose_grid)
!
! These calculations can be done in {\it mo\_model\_domain\_import} .
!-------------------------------------------------------------------------------

  !identifiers for prognostic and diagnostic variables
  ! required for mask application
  !
  !prognostic
  INTEGER, PARAMETER :: VELOCITY_NORMAL     = 1
  INTEGER, PARAMETER :: VELOCITY_TANGENTIAL = 2
  INTEGER, PARAMETER :: HEIGHT              = 3
  INTEGER, PARAMETER :: TRACER              = 4
  !
  !diagnostic
  INTEGER, PARAMETER :: VORTICITY           = 5
  INTEGER, PARAMETER :: DIVERGENCE          = 6
  INTEGER, PARAMETER :: KINETIC_ENERGY      = 7
  INTEGER, PARAMETER :: THICK_ED            = 8

  ! external
  ! #slo# changed 2010-08-03
  INTEGER,PARAMETER :: LAND          =  2  !  inner land
  INTEGER,PARAMETER :: LAND_BOUNDARY =  1  !  e.g. land cell with neighbouring wet cell
  INTEGER,PARAMETER :: BOUNDARY      =  0  !  edge where cells are differing
  INTEGER,PARAMETER :: SEA_BOUNDARY  = -1  !  e.g. wet cell with neighbouring land cell
  INTEGER,PARAMETER :: SEA           = -2  !  inner sea
  !INTEGER,PARAMETER :: LAND     = 0
  !INTEGER,PARAMETER :: SEA      = 1
  !INTEGER,PARAMETER :: BOUNDARY = 2

  ! dimensions of index list fields and index ranges for which grid points are reordered
  ! according to the refin_ctrl values
  ! Specifically:
  ! max_rl* indicates the number of cell/edge/vertex rows along the lateral boundary of nested
  !   domains for which grid points are reordered, i.e. moved to the beginning of the index lists;
  !   the number of cell rows for which the refin_ctrl flag is set is determined by the variable
  !   bdy_indexing_depth in prepare_gridref; it is in general larger than max_rlcell
  !  (the refin_ctrl flag here counts the distance from the lateral boundary in units of cell rows)
  ! ABS(min_rl*_int)-1 indicates the number of cell/edge/vertex rows overlapping with the lateral boundary
  !   of a nested domain for which grid points are reordered, i.e. moved to the end of the index lists;
  !   min_rl*_int refers to grid points overlapping with interior points of nested domains
  ! Finally, the indices between min_rl*_int-1 and min_rl*_int are reserved for halo points emerging
  !   from the MPI domain decomposition; these parts of the index lists are empty on exit of the
  !   grid generator and are filled only on mo_subdivision. However, the index list fields are always
  !   dimensioned with (min_rl*:max_rl*). The values set below are sufficient for a halo 
  !   width of two full cell rows; normally we use one, but stencils for high-order schemes may
  !   sometime require a halo width of two full rows
  !
  INTEGER, PARAMETER :: max_hw         = 2                         ! maximum halo width (n_ghost_rows)
  !
  INTEGER, PARAMETER :: min_rlcell_int = -4                        ! previously -6
  INTEGER, PARAMETER :: min_rlcell     = min_rlcell_int - 2*max_hw 
  INTEGER, PARAMETER :: max_rlcell     = 5                         ! previously 8
  INTEGER, PARAMETER :: min_rlvert_int = min_rlcell_int
  INTEGER, PARAMETER :: min_rlvert     = min_rlvert_int - (max_hw+1)
  INTEGER, PARAMETER :: max_rlvert     = max_rlcell
  INTEGER, PARAMETER :: min_rledge_int = 2*min_rlcell_int
  INTEGER, PARAMETER :: min_rledge     = min_rledge_int - (2*max_hw+1)
  INTEGER, PARAMETER :: max_rledge     = 2*max_rlcell

  ! maximum allowed number of model domains (10 should be enough for the time being)
  INTEGER, PARAMETER :: max_dom = 10

  ! Maximum allowed number of physical model domains
  INTEGER, PARAMETER :: max_phys_dom = 30

  ! maximum allowed number of tracers (20 should be enough for the time being)
  INTEGER, PARAMETER :: max_ntracer = 20

  ! identifiers for time stepping schemes

  INTEGER,PARAMETER :: TRACER_ONLY   = 1 ! pure tracer advection
  INTEGER,PARAMETER :: TWO_TL_SI     = 2 ! semi-implicit two time level
  INTEGER,PARAMETER :: LEAPFROG_EXPL = 3 ! explicit leapfrog
  INTEGER,PARAMETER :: LEAPFROG_SI   = 4 ! semi-implicit leapfrog
  INTEGER,PARAMETER :: RK4           = 5 ! standard 4th-order Runge-Kutta method
  INTEGER,PARAMETER :: SSPRK54       = 6 ! SSP RK(5,4)
  INTEGER,PARAMETER :: UNKNOWN       = 7

  ! Scheme for the "slow" component in the TWO_TL_SI time stepping

  INTEGER,PARAMETER :: EULER_FORWARD = 1
  INTEGER,PARAMETER :: AB2           = 2

  ! identifiers for NWP time control variables lcall_phy, t_elapsed_phy
  INTEGER, PARAMETER :: itupdate =  1
  INTEGER, PARAMETER :: itsatad  =  2
  INTEGER, PARAMETER :: itconv   =  3
  INTEGER, PARAMETER :: itccov   =  4
  INTEGER, PARAMETER :: itrad    =  5
  INTEGER, PARAMETER :: itsso    =  6
  INTEGER, PARAMETER :: itgscp   =  7
  INTEGER, PARAMETER :: itturb   =  8
  INTEGER, PARAMETER :: itradheat=  9 !calculation of radiative heating rates from radiative
                                      !fluxes with updated solar zenith angle
  INTEGER, PARAMETER :: itsfc    =  10

  INTEGER, PARAMETER :: iphysproc = 10! for NWP:
                                      ! number of slow physical processes:
                                      ! convection, sscloud cover, radiation, radheat, sso,
                                      ! microphysics, saturation adjustment, tracerupdate

!--------------------------------------------------------------------
END MODULE mo_impl_constants
