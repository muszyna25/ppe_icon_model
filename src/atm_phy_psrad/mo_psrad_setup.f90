module mo_psrad_setup

  use mo_psrad_general, only : wp 
  public :: setup_cloud_inhomogeneity, psrad_basic_setup

contains

  subroutine setup_cloud_inhomogeneity(zinhoml1_, zinhoml2_, &
    zinhoml3_, zinhomi_)
    use mo_psrad_cloud_optics, only : zinhoml1, zinhoml2, zinhoml3, zinhomi
    real(wp), intent(in) :: zinhoml1_, zinhoml2_, zinhoml3_, zinhomi_
    zinhoml1 = zinhoml1_
    zinhoml2 = zinhoml2_
    zinhoml3 = zinhoml3_
    zinhomi  = zinhomi_
  end subroutine setup_cloud_inhomogeneity

  subroutine psrad_basic_setup

    use mo_psrad_cloud_optics, only : setup_cloud_optics
    use mo_psrad_lrtm_setup, only : setup_lrtm
    use mo_psrad_srtm_setup, only : setup_srtm
    use mo_psrad_general, only : finish_cb, default_finish

    implicit none

    finish_cb => default_finish

    call setup_cloud_optics
    call setup_lrtm
    call setup_srtm
  end subroutine
  
end module mo_psrad_setup
