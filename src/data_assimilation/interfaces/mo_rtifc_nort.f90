!
!+ RTTOV interface for the case, that RTTOV is not available
!
MODULE mo_rtifc_nort
!
! Description:
!   This module contains version specific stuff for the mo_rtifc module for
!   the case that RTTOV is not available. All subroutines just set
!   the exit status to "ERR_NO_RTTOV_LIB".
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     Robin Faulwetter
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!=================================================

!---------------
! MACRO SETTINGS
!---------------

#include "mo_rtifc_macros.incf"

!========================
#if (_RTTOV_VERSION <= 0)
!========================
  
!-------------
! Modules used
!-------------
  use mo_rtifc_base

#if defined(_DACE_)
  use mo_rad,             only: t_radv             ! derived type to store radiance obs.
#endif


  implicit none

  private


  ! subroutines
  public :: rtifc_version          ! Version string
  public :: rtifc_init             ! Initialise RTTOV modules, read coeffs
  public :: rtifc_cleanup          ! frees memory allocated by rtifc_init
  public :: rtifc_get_preslev      ! Get pressure levels
  public :: rtifc_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rtifc_direct           ! calls RTTOV direct routine
  public :: rtifc_k                ! calls RTTOV K routine
  public :: rtifc_dealloc_profiles ! deallocate the profiles
  public :: rtifc_print_profiles

  ! Fake RTTOV options
  public :: rttov_options
  public :: rttov_opts_def

  ! parameters/variables
  public :: rtifc_vers

  interface rtifc_fill_input
    module procedure rtifc_fill_input_var
#if defined(_DACE_)
    module procedure rtifc_fill_input_rad
#endif
  end interface


  type rttov_options
    ! Fake RTTOV options
  end type rttov_options


  ! ifc version
  integer, parameter :: rtifc_vers = 0

  ! RTTOV default options
  type(rttov_options) :: rttov_opts_def

contains


  function rtifc_version() result(vers)
    character(len=17) :: vers
    
    vers = rttov_version()
    write(vers(13:),'("IFC",I2.2)') rtifc_vers
  end function rtifc_version


  subroutine rtifc_init(instruments,channels,nchans_inst,ropts,my_proc_id,n_proc,  &
                        io_proc_id,mpi_comm_type,status,&
                        path_coefs)
    integer,             intent(in)          :: instruments(:,:)
    integer,             intent(in)          :: channels(:,:)
    integer,             intent(in)          :: nchans_inst(:)
    type(rttov_options), intent(in)          :: ropts(:)
    integer,             intent(in)          :: my_proc_id             
    integer,             intent(in)          :: n_proc                 
    integer,             intent(in)          :: io_proc_id             
    integer,             intent(in)          :: mpi_comm_type          
    integer,             intent(out)         :: status                ! exit status
    character(*),        intent(in), optional:: path_coefs 

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_init



  subroutine rtifc_cleanup(status)
    integer, intent(out) :: status  ! exit status

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_cleanup



  subroutine rtifc_get_preslev(instrIdx, preslev, status)
    integer,  intent(in)  :: instrIdx
    real(wp), intent(out) :: preslev(:)
    integer, intent(out) :: status  ! exit status

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_get_preslev

#if defined(_DACE_)
  subroutine rtifc_fill_input_rad (status,rad,ropts,ivect,istart,iend, pe, &
                                   wr_profs, wr_profs_fmt)
    integer,             intent(out)          :: status   
    type(t_radv),        intent(in)           :: rad      
    type(rttov_options), intent(in)           :: ropts
    integer,             intent(in), optional :: ivect           
    integer,             intent(in), optional :: istart   
    integer,             intent(in), optional :: iend     
    integer,             intent(in), optional :: pe
    integer,             intent(in), optional :: wr_profs(:)
    character(len=*),    intent(in), optional :: wr_profs_fmt
  
    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_fill_input_rad
#endif


  subroutine rtifc_fill_input_var (status,ropts,press,temp,humi,t2m,q2m,psurf,hsurf,&
                                   u10m,v10m,stemp,stype,lat,lon,sat_zen,sun_zen,   &
                                   sat_azi,sun_azi,cloud,cfrac,id,ctp,cfraction,    &
                                   ivect,istart,iend,pe,wr_profs, wr_profs_fmt)
    integer,             intent(out)          :: status   
    type(rttov_options), intent(in)           :: ropts
    real(wp),            intent(in)           :: press     (:,:) 
    real(wp),            intent(in)           :: temp      (:,:) 
    real(wp),            intent(in)           :: humi      (:,:) 
    real(wp),            intent(in)           :: t2m         (:) 
    real(wp),            intent(in)           :: q2m         (:) 
    real(wp),            intent(in)           :: psurf       (:) 
    real(wp),            intent(in)           :: hsurf       (:) 
    real(wp),            intent(in)           :: u10m        (:) 
    real(wp),            intent(in)           :: v10m        (:) 
    real(wp),            intent(in)           :: stemp       (:) 
    integer,             intent(in)           :: stype       (:) 
    real(wp),            intent(in)           :: lat         (:) 
    real(wp),            intent(in)           :: lon         (:) 
    real(wp),            intent(in)           :: sat_zen     (:) 
    real(wp),            intent(in), optional :: sun_zen     (:) 
    real(wp),            intent(in), optional :: ctp         (:) 
    real(wp),            intent(in), optional :: cfraction   (:) 
    real(wp),            intent(in), optional :: sat_azi     (:) 
    real(wp),            intent(in), optional :: sun_azi     (:) 
    real(wp),            intent(in), optional :: cloud   (:,:,:) 
    real(wp),            intent(in), optional :: cfrac     (:,:) 
    integer,             intent(in), optional :: id          (:) 
    integer,             intent(in), optional :: ivect           
    integer,             intent(in), optional :: istart   
    integer,             intent(in), optional :: iend     
    integer,             intent(in), optional :: pe
    integer,             intent(in), optional :: wr_profs(:)
    character(len=*),    intent(in), optional :: wr_profs_fmt
  
    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_fill_input_var

  

  subroutine rtifc_direct (instrIdx,lprofs,chans,ropts,emissiv,t_b,status,t_b_clear,rad, &
                             radclear,radupclear,raddnclear,refdnclear,radovercast,radtotal,    &
                             transm,transmtotal,opdep,height,istore,errorstatus,reg_lim,rflag,  &
                             dealloc,iprint,rad_out_flg,pe,l_pio, l_opdep)
    integer,             intent(in)          :: instrIdx           
    integer,             intent(in)          :: lprofs         (:) 
    integer,             intent(in)          :: chans          (:) 
    type(rttov_options), intent(in)          :: ropts
    real(wp),            intent(inout)       :: emissiv      (:,:) 
    real(wp),            intent(out)         :: t_b          (:,:) 
    integer,             intent(out)         :: status             
    real(wp),            intent(out),optional:: t_b_clear    (:,:) 
    real(wp),            intent(out),optional:: rad          (:,:) 
    real(wp),            intent(out),optional:: radclear     (:,:) 
    real(wp),            intent(out),optional:: radupclear   (:,:) 
    real(wp),            intent(out),optional:: raddnclear   (:,:) 
    real(wp),            intent(out),optional:: refdnclear   (:,:) 
    real(wp),            intent(out),optional:: radovercast(:,:,:) 
    real(wp),            intent(out),optional:: radtotal     (:,:) 
    real(wp),            intent(out),optional:: transm     (:,:,:) 
    real(wp),            intent(out),optional:: transmtotal  (:,:) 
    real(wp),            intent(out),optional:: opdep      (:,:,:) 
    real(wp),            intent(out),optional:: height       (:,:) 
    integer,             intent(in) ,optional:: istore       (:,:) 
    integer,             intent(out),optional:: errorstatus    (:) 
    integer,             intent(out),optional:: reg_lim    (:,:,:) 
    integer,             intent(out),optional:: rflag        (:,:) 
    integer,             intent(in) ,optional:: iprint             
    logical,             intent(in) ,optional:: dealloc
    integer,             intent(in) ,optional:: rad_out_flg        
    integer,             intent(in) ,optional:: pe
    logical,             intent(in) ,optional:: l_pio
    logical,             intent(in) ,optional:: l_opdep
  
    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_direct


  subroutine rtifc_k (instrIdx,lprofs,chans,ropts,emissiv,emissiv_k,temp_k,     &
                      humi_k,t2m_k,q2m_k,stemp_k,t_b,status,t_b_clear,rad,       &
                      radclear,radtotal,radovercast,transm,opdep,psurf_k,u10m_k,v10m_k,&
                      o3_surf_k,wfetc_k,ctp_k,cfraction_k,clw_k,o3_k,     &
                      co2_k,n2o_k,co_k,ch4_k,istore,reg_lim,rflag,dealloc, &
                      iprint,rad_out_flg,pe,l_pio,l_opdep)
    integer ,            intent(in)          :: instrIdx           
    integer ,            intent(in)          :: lprofs         (:) 
    integer ,            intent(in)          :: chans          (:) 
    type(rttov_options), intent(in)          :: ropts
    real(wp),            intent(inout)       :: emissiv      (:,:) 
    real(wp),            intent(inout)       :: emissiv_k    (:,:) 
    real(wp),            intent(out)         :: temp_k     (:,:,:) 
    real(wp),            intent(out)         :: humi_k     (:,:,:) 
    real(wp),            intent(out)         :: t2m_k        (:,:) 
    real(wp),            intent(out)         :: q2m_k        (:,:) 
    real(wp),            intent(out)         :: stemp_k      (:,:) 
    real(wp),            intent(out)         :: t_b          (:,:) 
    integer,             intent(out)         :: status             
    real(wp),            intent(out),optional:: t_b_clear    (:,:) 
    real(wp),            intent(out),optional:: rad          (:,:) 
    real(wp),            intent(out),optional:: radclear     (:,:) 
    real(wp),            intent(out),optional:: radtotal     (:,:) 
    real(wp),            intent(out),optional:: radovercast(:,:,:) 
    real(wp),            intent(out),optional:: transm     (:,:,:) 
    real(wp),            intent(out),optional:: opdep      (:,:,:) 
    real(wp),            intent(out),optional:: psurf_k      (:,:) 
    real(wp),            intent(out),optional:: u10m_k       (:,:) 
    real(wp),            intent(out),optional:: v10m_k       (:,:) 
    real(wp),            intent(out),optional:: o3_surf_k    (:,:) 
    real(wp),            intent(out),optional:: wfetc_k      (:,:) 
    real(wp),            intent(out),optional:: ctp_k        (:,:) 
    real(wp),            intent(out),optional:: cfraction_k  (:,:) 
    real(wp),            intent(out),optional:: clw_k      (:,:,:) 
    real(wp),            intent(out),optional:: o3_k       (:,:,:) 
    real(wp),            intent(out),optional:: co2_k      (:,:,:) 
    real(wp),            intent(out),optional:: n2o_k      (:,:,:) 
    real(wp),            intent(out),optional:: co_k       (:,:,:) 
    real(wp),            intent(out),optional:: ch4_k      (:,:,:) 
    integer,             intent(in) ,optional:: istore       (:,:) 
    integer,             intent(out),optional:: reg_lim    (:,:,:) 
    integer,             intent(out),optional:: rflag        (:,:) 
    logical,             intent(in) ,optional:: dealloc
    integer,             intent(in) ,optional:: rad_out_flg        
    integer,             intent(in) ,optional:: iprint
    integer,             intent(in) ,optional:: pe
    logical,             intent(in) ,optional:: l_pio
    logical,             intent(in) ,optional:: l_opdep
  
    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_k


  subroutine rtifc_dealloc_profiles(status)
    integer, intent(out)         :: status             ! exit status

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_dealloc_profiles


  subroutine rtifc_print_profiles(unit)
    integer, intent(in) :: unit
  end subroutine rtifc_print_profiles


!===============================
#endif /* _RTTOV_VERSION <= 0 */
!===============================

end module mo_rtifc_nort
