!>
!! Configuration of the ECHAM cumulus convection scheme. 
!! Includes switches and tuning parameters, as well as 
!! other control variables.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2011-07)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_echam_conv_config

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_conv_config, echam_conv_config

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing main swithes for configuring 
  !! the cumulus convection scheme of ECHAM
  !!
  TYPE t_echam_conv_config

    ! Namelist variables

    INTEGER :: iconv     !< 1,2,3 for different convection schemes
    INTEGER :: ncvmicro  !< 0 or 1. Scheme for convective microphysics
  
    LOGICAL :: lmfpen    !< true when penetrative convection is switched on
    LOGICAL :: lmfmid    !< true when midlevel    convection is switched on
    LOGICAL :: lmfdd     !< true when cumulus downdraft      is switched on
    LOGICAL :: lmfdudv   !< true when cumulus friction       is switched on

    REAL(wp) :: dlev     !< "zdlev" in subroutine "cuasc". 
                         !< Critical thickness (unit: Pa) necessary for 
                         !< the onset of convective precipitation
  
    REAL(wp) :: cmftau   !< characteristic adjustment time scale
                         !< (replaces "ztau" in "cumastr"
    REAL(wp) :: cmfctop  !< fractional convective mass flux across the 
                         !< top of cloud
    REAL(wp) :: cprcon   !< coefficient for determining conversion
                         !< from cloud water to rain
  
    REAL(wp) :: cminbuoy !< minimum excess buoyancy
    REAL(wp) :: entrpen  !< entrainment rate for penetrative convection
 
    ! Currently unused namelist variables 
    !INTEGER :: nauto        !< autoconversion scheme. 1 or 2.
    !LOGICAL :: lconvmassfix !< aerosol mass fixer in convection
    !LOGICAL :: lmfscv       !< true when shallow convection is switched on

    ! Derived variables

    INTEGER :: nmctop    !< max. level for cloud base of mid level conv.

    REAL(wp),ALLOCATABLE :: cevapcu(:)  !< evaporation coefficient for kuo0
                                        !< In ECHAM6 it is declared in 
                                        !< mo_physc2, allocated in subroutine 
                                        !< alloc_mods, and initialized in 
                                        !< subroutine iniphy.
  
  END TYPE t_echam_conv_config

  !>
  !! The configuration state (variable).
  !! So far we have not yet tried to use different configurations for different
  !! domains (grid levels) in experiments with local grid refinement (nesting),
  !! thus the variable is declared as a scalar. Later it might be changed into
  !! an array of shape (/n_dom/) or (/MAX_DOM/).
  !!
  TYPE(t_echam_conv_config) :: echam_conv_config

END MODULE mo_echam_conv_config
