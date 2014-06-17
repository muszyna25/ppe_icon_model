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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_nml

  USE mo_exception,           ONLY: message, finish, message_text
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
  USE mo_art_config          ,ONLY: art_config ,t_volc_list,max_volc_input
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!

 TYPE t_list_volcanoes
    REAL(wp)           :: lat, lon
    CHARACTER (LEN=20) :: zname
!    INTEGER            :: nstart
    
 END TYPE t_list_volcanoes

    TYPE (t_list_volcanoes) :: art_volclist_tot(max_volc_input) !>list of volcanoes
                                   
    LOGICAL :: lart_seasalt        !< Treatment of sea salt aerosol (TRUE/FALSE)
    
    LOGICAL :: lart_dust           !< Treatment of mineral dust aerosol (TRUE/FALSE)

    LOGICAL :: lart_volcano        !< Treatment of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_emiss          !< Emission of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_conv           !< Convection of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_wash           !< Washout of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_rad            !< Radiative impact of volcanic ash (TRUE/FALSE)

    LOGICAL :: lart_cloud          !< Cloud volcanic ash interaction (TRUE/FALSE)

    INTEGER :: nart_emis_volcano_update    !< Time interval for reading volcano emission file

    LOGICAL :: lart_volclist  !< Input of Volcano coordinates. TRUE:use nml, FALSE:external data file is used.

    CHARACTER (LEN=120) :: volcanofile_path
    
    CHARACTER (LEN=120) :: art_folder

    LOGICAL :: lart_radioact                !< Treatment of volcanic ash (TRUE/FALSE)
   
    LOGICAL :: lart_decay_radioact          !< Treatment of volcanic ash (TRUE/FALSE)

    CHARACTER (LEN=120) :: radioactfile_path

    LOGICAL :: lart_chemtracer                !< Treatment of chemical tracer (TRUE/FALSE)

    LOGICAL :: lart_loss_chemtracer           !< Treatment of chemical loss (TRUE/FALSE)

    NAMELIST/art_nml/ lart_seasalt, lart_dust, lart_volcano, lart_emiss,        &
   &               lart_conv, lart_wash, lart_rad, lart_cloud,                        &
   &               nart_emis_volcano_update,art_volclist_tot, lart_volclist,          &
   &               volcanofile_path,lart_radioact,lart_decay_radioact,                &
   &               radioactfile_path, lart_chemtracer,lart_loss_chemtracer, art_folder

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
    INTEGER :: iunit

    !-----------------------
    ! 1. default settings   
    !-----------------------
    lart_seasalt        = .FALSE.        ! Treatment of sea salt aerosol
    lart_dust           = .FALSE.        ! Treatment of mineral dust aerosol
    lart_volcano        = .FALSE.        ! Treatment of volcanic ash
    lart_emiss   = .FALSE.        ! Emission of volcanic ash
    lart_conv   = .FALSE.        ! Convection of volcanic ash
    lart_wash   = .FALSE.        ! Washout of volcanic ash 
    lart_rad    = .FALSE.        ! Radiative impact of volcanic ash
    lart_cloud  = .FALSE.        ! Impact on clouds

    nart_emis_volcano_update= 0          ! Time interval for reading emission file
    lart_volclist=.FALSE.
    art_volclist_tot(:)%lon   = -1000._wp     ! Longitude coordinate of each volcano. 
                                           !-1000 is used for creating the list of volcanoes.  
    art_volclist_tot(:)%lat   = 0._wp     ! Latitude coordinate of each volcano
    art_volclist_tot(:)%zname = ""        ! Name of volcanoes
!   art_volclist_tot(:)%nstart= 1         ! Start time of volcanic eruption

    volcanofile_path          =  "./volcanofile"

    lart_radioact             = .FALSE.    !< Treatment of radioactive nuclides
    lart_decay_radioact       = .FALSE.    !< Treatment of radioactive decay
    radioactfile_path         =  "./radioactfile"

    lart_chemtracer           = .FALSE.     !< Treatment of chemical tracer (TRUE/FALSE)
    lart_loss_chemtracer      = .FALSE.     !< Treatment of chemical loss (TRUE/FALSE)

    art_folder                =  "./art/"

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
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, art_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, art_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, art_nml)    ! write settings to temporary text file
      END IF
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
        IF (nvolc > max_volc_input) THEN
          CALL message(TRIM(routine), 'number of volcanos in input file &
        & exeedes maximum allowed number of 20, please use file input')
          EXIT           ! maximum number reached.
        ENDIF
        IF (art_volclist_tot(nvolc)%lon == -1000._wp) EXIT   ! default value --> this component is not filled. Max. comp. reached.
      ENDDO
      nvolc=nvolc-1
    DO jg= 0,max_dom
      art_config(jg)%lart_seasalt             = lart_seasalt
      art_config(jg)%lart_dust                = lart_dust
      art_config(jg)%lart_volcano             = lart_volcano
      art_config(jg)%lart_emiss               = lart_emiss
      art_config(jg)%lart_conv                = lart_conv
      art_config(jg)%lart_wash                = lart_wash
      art_config(jg)%lart_rad                 = lart_rad
      art_config(jg)%lart_cloud               = lart_cloud
      art_config(jg)%nart_emis_volcano_update = nart_emis_volcano_update 
      art_config(jg)%lart_volclist            = lart_volclist
      art_config(jg)%nvolc                    = nvolc
      art_config(jg)%volcanofile_path         = volcanofile_path
      art_config(jg)%lart_radioact            = lart_radioact
      art_config(jg)%lart_decay_radioact      = lart_radioact
      art_config(jg)%radioactfile_path        = radioactfile_path
      art_config(jg)%lart_chemtracer          = lart_chemtracer
      art_config(jg)%lart_loss_chemtracer     = lart_loss_chemtracer
      art_config(jg)%art_folder               = art_folder
      
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
