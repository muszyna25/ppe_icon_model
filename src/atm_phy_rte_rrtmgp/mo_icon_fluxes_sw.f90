
module mo_icon_fluxes_sw

  use mo_kind,          only: wp
  use mo_optical_props, only: ty_optical_props
  use mo_fluxes, only: ty_fluxes_broadband, &
    are_desired_broadband, reduce_broadband
  implicit none

  type, extends(ty_fluxes_broadband), public :: ty_icon_fluxes_sw
    real(wp), dimension(:), pointer :: &
      vis_dn_dir_sfc => NULL(), &
      par_dn_dir_sfc => NULL(), &
      nir_dn_dir_sfc => NULL(), &
      vis_dn_dff_sfc => NULL(), &
      par_dn_dff_sfc => NULL(), &
      nir_dn_dff_sfc => NULL(), &
      vis_up_sfc => NULL(), &
      par_up_sfc => NULL(), &
      nir_up_sfc => NULL()
    REAL(wp), ALLOCATABLE :: &
      band_weight(:), & ! adjustment for current Earth/Sun distance
      frc_par(:), &
      frc_vis(:)
    
  contains
    procedure, public :: reduce      => reduce_icon
    procedure, public :: are_desired => are_desired_icon
    final             :: del
  end type ty_icon_fluxes_sw

  public :: set_fractions

contains

  subroutine set_fractions(this, optical_props, solar_constant, ssi_fraction)
    TYPE(ty_icon_fluxes_sw), INTENT(INOUT) :: this
    CLASS(ty_optical_props), INTENT(IN) :: optical_props
    REAL(wp), INTENT(IN) :: solar_constant, ssi_fraction(:)

    INTEGER :: i, nbndsw
    REAL(wp), PARAMETER :: nir_vis_boundary   = 14500._wp
    REAL(wp), ALLOCATABLE :: wavenum(:,:), delwave(:)

    nbndsw = optical_props%get_nband()
    
    IF (.not. ALLOCATED(this%frc_par)) THEN
      ALLOCATE(this%frc_par(nbndsw))
      ALLOCATE(this%frc_vis(nbndsw))
      ALLOCATE(this%band_weight(nbndsw))
      ALLOCATE(wavenum(2,nbndsw))
      ALLOCATE(delwave(nbndsw))
      
      wavenum = optical_props%get_band_lims_wavenumber()
      delwave = wavenum(2,:) - wavenum(1,:)

      this%frc_par(1:nbndsw) = 0.0
      this%frc_par(9) = 0.533725_wp
      this%frc_par(10) = 1.0_wp
      this%frc_par(11) = 0.550164_wp

      DO i = 1, nbndsw 
        this%frc_vis(i) = MAX(0.0_wp, MIN(1.0_wp, &
          (wavenum(2,i) - nir_vis_boundary) / delwave(i) ))
      ENDDO

      ! DA TODO: move to the GPU
      !$ACC enter data copyin(this, this%frc_par, this%frc_vis, this%band_weight)
    ENDIF

    ! --- weight radiation within a band for the solar cycle ---
    ! solar_constant contains TSI (the "solar constant") scaled with the
    ! Sun-Earth distance. ssi_fraction contains the relative contribution
    ! of each band to TSI. ssi_default is the (originally only
    ! implicitly defined) solar flux in the 14 bands.

    ! This routine is called from within RRTMGP, so it shouldn't be async
    !$ACC parallel loop default(none) present(this) copyin(ssi_fraction) gang vector
    DO i = 1, nbndsw
      this%band_weight(i) = solar_constant*ssi_fraction( MOD(i, nbndsw)+1 ) ! / ssi_default(:)
    END DO
  END SUBROUTINE set_fractions

  function reduce_icon(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_icon_fluxes_sw),        intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: nlev, ncol, ngpt, nbndsw, isfc, band, gpt, limits(2), jl
    integer :: band2gpt(2, spectral_disc%get_nband())

    error_msg = reduce_broadband(this, gpt_flux_up, gpt_flux_dn, &
      spectral_disc, top_at_1, gpt_flux_dn_dir)
    if (TRIM(error_msg) /= '') return

    ncol = size(gpt_flux_up,1)
    nlev = size(gpt_flux_up,2)
    ngpt = size(gpt_flux_up,3)
    nbndsw = spectral_disc%get_nband()
    if (top_at_1) then
      isfc = nlev
    else
      isfc = 1
    endif

    !DA TODO: this has to run on GPU
    band2gpt(:,:) = spectral_disc%band2gpt(:,:)

    ! This routine is called from within RRTMGP, so it shouldn't be async
    !$ACC parallel default(present) copyin(band2gpt) vector_length(64)
    !$ACC loop gang vector private(limits)
    DO jl = 1, ncol
      this%vis_dn_dir_sfc(jl) = 0.0_wp
      this%par_dn_dir_sfc(jl) = 0.0_wp
      this%nir_dn_dir_sfc(jl) = 0.0_wp
      this%vis_dn_dff_sfc(jl) = 0.0_wp
      this%par_dn_dff_sfc(jl) = 0.0_wp
      this%nir_dn_dff_sfc(jl) = 0.0_wp
      this%vis_up_sfc(jl) = 0.0_wp
      this%par_up_sfc(jl) = 0.0_wp
      this%nir_up_sfc(jl) = 0.0_wp

      !$ACC loop seq
      DO band = 1, nbndsw
        limits(:) = band2gpt(:, band)
        !$ACC loop seq
        DO gpt = limits(1), limits(2)
          this%vis_dn_dir_sfc(jl) = this%vis_dn_dir_sfc(jl) + &
            this%frc_vis(band) * gpt_flux_dn_dir(jl,isfc,gpt)
          this%par_dn_dir_sfc(jl) = this%par_dn_dir_sfc(jl) + &
            this%frc_par(band) * gpt_flux_dn_dir(jl,isfc,gpt)
          this%nir_dn_dir_sfc(jl) = this%nir_dn_dir_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * gpt_flux_dn_dir(jl,isfc,gpt)

          this%vis_dn_dff_sfc(jl) = this%vis_dn_dff_sfc(jl) + &
            this%frc_vis(band) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))
          this%par_dn_dff_sfc(jl) = this%par_dn_dff_sfc(jl) + &
            this%frc_par(band) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))
          this%nir_dn_dff_sfc(jl) = this%nir_dn_dff_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * ( &
            gpt_flux_dn(jl,isfc,gpt) - gpt_flux_dn_dir(jl,isfc,gpt))

          this%vis_up_sfc(jl) = this%vis_up_sfc(jl) + &
            this%frc_vis(band) * gpt_flux_up(jl,isfc,gpt)
          this%par_up_sfc(jl) = this%par_up_sfc(jl) + &
            this%frc_par(band) * gpt_flux_up(jl,isfc,gpt)
          this%nir_up_sfc(jl) = this%nir_up_sfc(jl) + &
            (1.0_wp - this%frc_vis(band)) * gpt_flux_up(jl,isfc,gpt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC end parallel
  
  end function reduce_icon

  function are_desired_icon(this)
    class(ty_icon_fluxes_sw), intent(in) :: this
    logical                              :: are_desired_icon

    are_desired_icon = are_desired_broadband(this)
  end function are_desired_icon

  subroutine del(this)
    type(ty_icon_fluxes_sw), intent(inout) :: this
    !$ACC exit data delete(this%frc_par, this%frc_vis, this%band_weight, this)
  end subroutine del
  ! --------------------------------------------------------------------------------------
end module mo_icon_fluxes_sw
