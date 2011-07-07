!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_lnd_nwp_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for NWP land
  !!--------------------------------------------------------------------------


  TYPE t_nwp_lnd_config

    INTEGER ::  nlev_soil, nztlev  !! number of soil layers, time integration scheme
    INTEGER ::  nlev_snow        !! number of snow layers
    INTEGER ::  nsfc_subs        !! number of TILES

    LOGICAL ::       &
         lseaice,    & !> forecast with sea ice model
         llake,      & !! forecst with lake model FLake
         lmelt     , & !! soil model with melting process
         lmelt_var , & !! freezing temperature dependent on water content
         lmulti_snow   !! run the multi-layer snow model

  END TYPE t_nwp_lnd_config

  !>
  !!
  TYPE(t_nwp_lnd_config) :: nwp_lnd_config(max_dom)

END MODULE mo_lnd_nwp_config
