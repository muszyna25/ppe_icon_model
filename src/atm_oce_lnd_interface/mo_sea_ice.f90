!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters used for surface forcing
!! of the hydrostatic ocean model.
!!
!! @author Peter Korn, MPI
!! @author Dirk Notz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
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
MODULE mo_sea_ice
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: dtime
USE mo_dynamics_config,     ONLY: nold
USE mo_model_domain,        ONLY: t_patch
USE mo_exception,           ONLY: finish, message
USE mo_impl_constants,      ONLY: success, max_char_length
USE mo_loopindices,         ONLY: get_indices_c

IMPLICIT NONE

PRIVATE


CHARACTER(len=*), PARAMETER :: version = '$Id$'
! Public interface

! public subroutines
PUBLIC :: construct_sea_ice 
PUBLIC :: destruct_sea_ice

! public types
PUBLIC :: t_sea_ice


! Definition of forcing types


TYPE t_sea_ice

! The description of the sea-ice state, defined on cell-centers
! dimension: (nproma, nblks_c)

  LOGICAL, ALLOCATABLE :: &
    isice(:,:,:)    ! Logical field that marks ice-covered grid cells
  
  REAL(wp), ALLOCATABLE :: &
    alb(:,:,:)         ,  & ! Albedo of snow-ice system
    Tsurf(:,:,:)       ,  & ! Surface temperature                      [C]
    T1 (:,:,:)         ,  & ! Temperature upper layer                  [C]
    T2 (:,:,:)         ,  & ! Temperature lower layer                  [C]
    E1(:,:,:)          ,  & ! Energy content upper layer               [Jm/kg]
    E2(:,:,:)          ,  & ! Energy content lower layer               [Jm/kg]
    hi(:,:,:)          ,  & ! Ice thickness                            [m]
    hs(:,:,:)          ,  & ! Snow thickness                           [m]
    hiold(:,:,:)       ,  & ! Ice thickness at previous time step      [m]
    hsold(:,:,:)       ,  & ! Snow thickness at previous time step     [m]
    Qtop(:,:,:)        ,  & ! Energy flux available for surface melting[W/m2]
    Qbot(:,:,:)        ,  & ! Energy flux at ice-ocean interface       [W/m2]
    heatocei(:,:,:)    ,  & ! Energy to ocean when all ice is melted   [J]
    snow_to_ice(:,:,:) ,  & ! amount of snow that is transformed to ice[m]
    surfmelt(:,:,:)    ,  & ! surface melt water running into ocean    [m]
    surfmeltT(:,:,:)   ,  & ! Mean temperature of surface melt water   [C]
    evapwi(:,:,:)      ,  & ! amount of evaporated water if no ice left[kg/m2]
    conc(:,:,:)             ! ice concentration in each ice class
    
  REAL(wp), ALLOCATABLE :: &
    u(:,:)          ,  & ! Zonal velocity                            [m/s]
    v(:,:)          ,  & ! Meridional velocity                       [m/s]
    concSum(:,:)    ,  & ! Total ice concentration within a grid cell    
    newice(:,:)     ,  & ! New ice growth in open water              [m]
    zUnderIce(:,:)       ! water in upper ocean grid cell below ice  [m]

   INTEGER ::  kice = 1   ! Number of ice-thickness classes

  REAL(wp), ALLOCATABLE ::  hi_lim(:)   ! Thickness limits 

END TYPE t_sea_ice


CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of sea-ice model, allocates all components and assigns zero. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE construct_sea_ice(ppatch, p_ice, i_no_ice_thick_class)
  TYPE(t_patch),     INTENT(in)     :: ppatch
  TYPE (t_sea_ice),  INTENT (INOUT) :: p_ice
  INTEGER, INTENT(IN)               :: i_no_ice_thick_class

  !Local variables
  !INTEGER i

  INTEGER :: nblks_c, ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice: constructor'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )

  nblks_c = ppatch%nblks_c

  p_ice%kice = i_no_ice_thick_class

  ALLOCATE(p_ice%isice(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for isice failed')
  END IF

  ALLOCATE(p_ice%alb(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for alb failed')
  END IF

 ALLOCATE(p_ice%Tsurf(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Tsurf failed')
  END IF

 ALLOCATE(p_ice%T1(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for T1 failed')
  END IF

  ALLOCATE(p_ice%T2(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for T2 failed')
  END IF

 ALLOCATE(p_ice%E1(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for E1 failed')
  END IF

  ALLOCATE(p_ice%E2(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for E2 failed')
  END IF

  ALLOCATE(p_ice%hi(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hi failed')
  END IF

  ALLOCATE(p_ice%hs(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hs failed')
  END IF

  ALLOCATE(p_ice%hiold(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hiold failed')
  END IF

  ALLOCATE(p_ice%hsold(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hsold failed')
  END IF

  ALLOCATE(p_ice%Qtop(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Qtop failed')
  END IF

  ALLOCATE(p_ice%Qbot(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for Qbot failed')
  END IF

  ALLOCATE(p_ice%heatocei(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for heatocei failed')
  END IF

  ALLOCATE(p_ice%snow_to_ice(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for snow_to_ice failed')
  END IF

  ALLOCATE(p_ice%surfmelt(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surfmelt failed')
  END IF

  ALLOCATE(p_ice%surfmeltT(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for surfmeltT failed')
  END IF

  ALLOCATE(p_ice%evapwi(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for evapwi failed')
  END IF

  ALLOCATE(p_ice%conc(nproma,nblks_c,i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for conc failed')
  END IF

  ALLOCATE(p_ice%u(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for u failed')
  END IF
    
  ALLOCATE(p_ice%v(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for v failed')
  END IF

  ALLOCATE(p_ice%concSum(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for concSum failed')
  END IF

  ALLOCATE(p_ice%newice(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for newice failed')
  END IF

  ALLOCATE(p_ice%zUnderIce(nproma,nblks_c), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for zUnderIce failed')
  END IF

  ALLOCATE(p_ice%hi_lim(i_no_ice_thick_class), STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'allocation for hi_lim failed')
  END IF


  p_ice%alb        = 0.0_wp       
  p_ice%Tsurf      = 0.0_wp       
  p_ice%T1         = 0.0_wp        
  p_ice%T2         = 0.0_wp      
  p_ice%E1         = 0.0_wp        
  p_ice%E2         = 0.0_wp          
  p_ice%hi         = 0.0_wp          
  p_ice%hs         = 0.0_wp        
  p_ice%hiold      = 0.0_wp     
  p_ice%hsold      = 0.0_wp       
  p_ice%Qtop       = 0.0_wp       
  p_ice%Qbot       = 0.0_wp       
  p_ice%heatocei   = 0.0_wp    
  p_ice%snow_to_ice= 0.0_wp  
  p_ice%surfmelt   = 0.0_wp    
  p_ice%surfmeltT  = 0.0_wp   
  p_ice%evapwi     = 0.0_wp      
  p_ice%conc       = 0.0_wp         
    
  
 
  p_ice%snow_to_ice= 0.0_wp 
  p_ice%surfmelt   = 0.0_wp 
  p_ice%surfmeltT  = 0.0_wp 
  p_ice%evapwi     = 0.0_wp 
  p_ice%conc       = 0.0_wp 
  p_ice%u          = 0.0_wp 
  p_ice%v          = 0.0_wp 
  p_ice%concSum    = 0.0_wp 
  p_ice%newice     = 0.0_wp 
  p_ice%zUnderIce  = 0.0_wp 

  IF(p_ice%kice==1)THEN
    p_ice%hi_lim = 0.0_wp
  ELSEIF(p_ice%kice==8)THEN
    p_ice%hi_lim(:)=(/ 0.0_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp /)
  ENDIF


END SUBROUTINE construct_sea_ice
!-------------------------------------------------------------------------
  !
  !> Destructor of sea-ice model, deallocates all components. 
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !
  SUBROUTINE destruct_sea_ice(p_ice)
  TYPE (t_sea_ice),  INTENT (INOUT)    :: p_ice


  !Local variables
  INTEGER :: ist
  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice: destructor'
  !-------------------------------------------------------------------------
  CALL message(TRIM(routine), 'start' )
  DEALLOCATE(p_ice%alb, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for alb failed')
  END IF

 DEALLOCATE(p_ice%Tsurf, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Tsurf failed')
  END IF

 DEALLOCATE(p_ice%T1, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for T1 failed')
  END IF

  DEALLOCATE(p_ice%T2, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for T2 failed')
  END IF

 DEALLOCATE(p_ice%E1, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for E1 failed')
  END IF

  DEALLOCATE(p_ice%E2, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for E2 failed')
  END IF

  DEALLOCATE(p_ice%hi, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hi failed')
  END IF

  DEALLOCATE(p_ice%hs, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hs failed')
  END IF

  DEALLOCATE(p_ice%hiold, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hiold failed')
  END IF

  DEALLOCATE(p_ice%hsold, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for hsold failed')
  END IF

  DEALLOCATE(p_ice%Qtop, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Qtop failed')
  END IF

  DEALLOCATE(p_ice%Qbot, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for Qbot failed')
  END IF

  DEALLOCATE(p_ice%heatocei, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for heatocei failed')
  END IF

  DEALLOCATE(p_ice%snow_to_ice, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for snow_to_ice failed')
  END IF

  DEALLOCATE(p_ice%surfmelt, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for surfmelt failed')
  END IF

  DEALLOCATE(p_ice%surfmeltT, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for surfmeltT failed')
  END IF

  DEALLOCATE(p_ice%evapwi, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for evapwi failed')
  END IF

  DEALLOCATE(p_ice%conc, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for conc failed')
  END IF

  DEALLOCATE(p_ice%u, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for u failed')
  END IF
    
  DEALLOCATE(p_ice%v, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'adellocation for v failed')
  END IF

  DEALLOCATE(p_ice%concSum, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for concSum failed')
  END IF

  DEALLOCATE(p_ice%newice, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for newice failed')
  END IF

  DEALLOCATE(p_ice%zUnderIce, STAT=ist)
  IF (ist/=SUCCESS) THEN
    CALL finish(TRIM(routine),'deallocation for zUnderIce failed')
  END IF
 
END SUBROUTINE destruct_sea_ice
!-------------------------------------------------------------------------


END MODULE mo_sea_ice
