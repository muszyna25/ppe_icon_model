!>
!! Namelist for ART-package
!!
!! Subroutine is called by read_atmo_namelists for setting up the ART-package 
!!
!! @author Daniel Reinert, DWD
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-12-08)
!! Modification of namelist parameters, Kristina Lundgren (2012-03-21)
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
MODULE mo_art_nml

  USE mo_kind                ,ONLY: wp
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_exception           ,ONLY: finish
  USE mo_io_units            ,ONLY: nnml, nnml_output
  USE mo_master_control      ,ONLY: is_restart_run
  USE mo_impl_constants      ,ONLY: max_dom
  USE mo_namelist            ,ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi                 ,ONLY: my_process_is_stdio
  USE mo_io_restart_namelist ,ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_art_config          ,ONLY: art_config ,t_volc_list,MAX_NUM_VOLC 

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!

 TYPE t_list_volcanoes
    REAL(wp)           :: lat, lon
    CHARACTER (LEN=20) :: zname
!    INTEGER            :: nstart
    
 END TYPE t_list_volcanoes

 TYPE (t_list_volcanoes) :: art_volclist_tot(MAX_NUM_VOLC) !>list of volcanoes
 
    LOGICAL :: lart                        !< main switch for using the ART-package
                                           !< .TRUE.: switch ON
                                           !< .FALSE.: switch OFF

    LOGICAL :: lart_volcano                !< Treatment of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_emis_volcano           !< Emission of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_conv_volcano           !< Convection of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_wash_volcano           !< Washout of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_rad_volcano            !< Radiative impact of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_cloud_volcano          !< Cloud volcanic ash interaction (TRUE/FALSE)

    INTEGER :: nart_emis_volcano_update    !< Time interval for reading volcano emission file

  LOGICAL :: lart_volclist  !< Input of Volcano coordinates. TRUE:use nml, FALSE:external data file is used.
  NAMELIST/art_nml/ lart, lart_volcano, lart_emis_volcano,lart_conv_volcano, lart_wash_volcano, &
    &               lart_rad_volcano, lart_cloud_volcano, nart_emis_volcano_update,art_volclist_tot, lart_volclist


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for ART-package. 
  !!
  !! This subroutine 
  !! - reads the Namelist for the ART-package
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-12-08)
  !!
  SUBROUTINE read_art_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: jb, jc, nblks, npromz, nvolc, ivolc   
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_art_nml: read_art_nml'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    lart                = .FALSE.        ! ART-package switched off
    lart_volcano        = .FALSE.        ! Treatment of volcanic ash
    lart_emis_volcano   = .FALSE.        ! Emission of volcanic ash
    lart_conv_volcano   = .FALSE.        ! Convection of volcanic ash
    lart_wash_volcano   = .FALSE.        ! Washout of volcanic ash 
    lart_rad_volcano    = .FALSE.        ! Radiative impact of volcanic ash
    lart_cloud_volcano  = .FALSE.        ! Impact on clouds

    nart_emis_volcano_update= 0          ! Time interval for reading emission file
    lart_volclist=.FALSE.
    art_volclist_tot(:)%lon   = -1._wp     ! Longitude coordinate of each volcano. 
                                           !-1 is used for creating the list of volcanoes.  
    art_volclist_tot(:)%lat   = 0._wp     ! Latitude coordinate of each volcano
    art_volclist_tot(:)%zname = ""        ! Name of volcanoes
!   art_volclist_tot(:)%nstart= 1         ! Start time of volcanic eruption
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('art_nml')
      READ(funit,NML=art_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('art_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, art_nml)
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    
    ! Determine lenght of volcano list, i.e. number ov volcanoes
      nvolc=0
      DO
        nvolc = nvolc + 1
        IF (nvolc > MAX_NUM_VOLC) EXIT           ! maximum number reached.
        IF (art_volclist_tot(nvolc)%lon == -1._wp) EXIT   ! default value --> this component is not filled. Max. comp. reached.
      ENDDO
      nvolc=nvolc-1
    DO jg= 0,max_dom
      art_config(jg)%lart          = lart
      art_config(jg)%lart_volcano     = lart_volcano
      art_config(jg)%lart_emis_volcano     = lart_emis_volcano
      art_config(jg)%lart_conv_volcano  = lart_conv_volcano
      art_config(jg)%lart_wash_volcano  = lart_wash_volcano
      art_config(jg)%lart_rad_volcano     = lart_rad_volcano
      art_config(jg)%lart_cloud_volcano   = lart_cloud_volcano
      art_config(jg)%nart_emis_volcano_update   = nart_emis_volcano_update 
      art_config(jg)%lart_volclist=lart_volclist
      art_config(jg)%nvolc         = nvolc
      
      nblks=nvolc/nproma+1
      npromz=nvolc-nproma*(nblks-1)
      art_config(jg)%nblks         = nblks   ! 
      art_config(jg)%npromz        = npromz  !  

      ALLOCATE(art_config(jg)%volclist(nproma,nblks))

      jc=0
      jb=1
       DO ivolc=1,nvolc
          jc=jc+1
          IF (jc>nproma) THEN
            jc = 1
            jb = jb + 1
          ENDIF
          art_config(jg)%volclist(jc,jb)%zname=art_volclist_tot(ivolc)%zname
          art_config(jg)%volclist(jc,jb)%location%lat=art_volclist_tot(ivolc)%lat
          art_config(jg)%volclist(jc,jb)%location%lon=art_volclist_tot(ivolc)%lon
       ENDDO !nvolc
    ENDDO !jg



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=art_nml)                    
      CALL store_and_close_namelist(funit, 'art_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=art_nml)


  END SUBROUTINE read_art_namelist

END MODULE mo_art_nml
