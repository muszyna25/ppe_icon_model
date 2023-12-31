 module cvmix_tidal

!BOP
!\newpage
! !MODULE: cvmix_tidal
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to initialize the derived types needed for
!  tidal mixing (currently just the Simmons scheme) and to set the viscosity
!  and diffusivity coefficients accordingly.
!\\
!\\
!  References:\\
!  * HL Simmons, SR Jayne, LC St. Laurent, and AJ Weaver.
!  Tidally Driven Mixing in a Numerical Model of the Ocean General Circulation.
!  Ocean Modelling, 2004.
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_zero,                               &
                                    cvmix_data_type,                          &
                                    cvmix_strlen,                             &
                                    cvmix_global_params_type,                 &
                                    CVMIX_OVERWRITE_OLD_VAL,                  &
                                    CVMIX_SUM_OLD_AND_NEW_VALS,               &
                                    CVMIX_MAX_OLD_AND_NEW_VALS
  use cvmix_utils,           only : cvmix_update_wrap
  use cvmix_put_get,         only : cvmix_put

  USE mo_exception,          ONLY: finish

!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_tidal
  public :: cvmix_compute_vert_dep
  public :: cvmix_coeffs_tidal
  public :: cvmix_compute_Simmons_invariant
  public :: cvmix_put_tidal
  public :: cvmix_get_tidal_real
  public :: cvmix_get_tidal_str

  interface cvmix_coeffs_tidal
    module procedure cvmix_coeffs_tidal_low
    module procedure cvmix_coeffs_tidal_wrap
  end interface cvmix_coeffs_tidal

  interface cvmix_compute_Simmons_invariant
    module procedure cvmix_compute_Simmons_invariant_low
    module procedure cvmix_compute_Simmons_invariant_wrap
  end interface cvmix_compute_Simmons_invariant

  interface cvmix_put_tidal
    module procedure cvmix_put_tidal_int
    module procedure cvmix_put_tidal_real
    module procedure cvmix_put_tidal_str
  end interface cvmix_put_tidal

! !PUBLIC TYPES:

  ! cvmix_tidal_params_type contains the necessary parameters for tidal mixing
  ! (currently just Simmons)
  type, public :: cvmix_tidal_params_type
    private
      ! Tidal mixing scheme being used (currently only support Simmons et al)
      character(len=cvmix_strlen) :: mix_scheme

      ! efficiency is the mixing efficiency (Gamma in Simmons)
      real(cvmix_r8) :: efficiency           ! units: unitless (fraction)

      ! local_mixing_frac is the tidal dissipation efficiency (q in Simmons)
      real(cvmix_r8) :: local_mixing_frac    ! units: unitless (fraction)

      ! vertical_decay_scale is zeta in the Simmons paper (used to compute the
      ! vertical deposition function)
      real(cvmix_r8) :: vertical_decay_scale ! units: m

      ! depth_cutoff is depth of the shallowest column where tidal mixing is
      ! computed (like all depths, positive => below the surface)
      real(cvmix_r8) :: depth_cutoff         ! units: m

      ! max_coefficient is the largest acceptable value for diffusivity
      real(cvmix_r8) :: max_coefficient      ! units: m^2/s

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals

      ! Note: need to include some logic to avoid excessive memory use
  end type cvmix_tidal_params_type
!EOP

  type(cvmix_tidal_params_type), target :: CVmix_tidal_params_saved

  CHARACTER(LEN=*), PARAMETER :: module_name = 'cvmix_tidal'

contains

!BOP

! !IROUTINE: cvmix_init_tidal
! !INTERFACE:

  subroutine cvmix_init_tidal(CVmix_tidal_params_user, mix_scheme, efficiency,&
                              vertical_decay_scale, max_coefficient,          &
                              local_mixing_frac, depth_cutoff, old_vals)

! !DESCRIPTION:
!  Initialization routine for tidal mixing. There is currently just one
!  supported schemes - set \verb|mix_scheme = 'simmons'| to use the Simmons
!  mixing scheme.
!
! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), optional, intent(in) :: mix_scheme, old_vals
    real(cvmix_r8),   optional, intent(in) :: efficiency
    real(cvmix_r8),   optional, intent(in) :: vertical_decay_scale
    real(cvmix_r8),   optional, intent(in) :: max_coefficient
    real(cvmix_r8),   optional, intent(in) :: local_mixing_frac
    real(cvmix_r8),   optional, intent(in) :: depth_cutoff

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), optional, target, intent(inout) ::         &
                                              CVmix_tidal_params_user

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_init_tidal'

!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_out

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_out => CVmix_tidal_params_user
    else
      CVmix_tidal_params_out => CVmix_tidal_params_saved
    end if

    if (present(mix_scheme)) then
      call cvmix_put_tidal("mix_scheme", trim(mix_scheme),                    &
                           CVmix_tidal_params_user)
    else
      call cvmix_put_tidal("mix_scheme", "Simmons", CVmix_tidal_params_user)
    end if

    select case (trim(CVmix_tidal_params_out%mix_scheme))
      case ('simmons','Simmons')
        ! Unitless parameters
        if (present(efficiency)) then
          call cvmix_put_tidal("efficiency", efficiency,                      &
                               CVmix_tidal_params_user)
        else
          call cvmix_put_tidal("efficiency", 0.2_cvmix_r8,                    &
                               CVmix_tidal_params_user)
        end if

        if (present(local_mixing_frac)) then
          call cvmix_put_tidal("local_mixing_frac", local_mixing_frac,        &
                               CVmix_tidal_params_user)
        else
          call cvmix_put_tidal("local_mixing_frac", 3, CVmix_tidal_params_user)
        end if

        ! Parameters with units
        if (present(vertical_decay_scale)) then
          call cvmix_put_tidal("vertical_decay_scale", vertical_decay_scale,  &
                               CVmix_tidal_params_user)
        else
          call cvmix_put_tidal("vertical_decay_scale", 500,                   &
                               CVmix_tidal_params_user)
        end if

        if (present(depth_cutoff)) then
          call cvmix_put_tidal("depth_cutoff", depth_cutoff,                  &
                               CVmix_tidal_params_user)
        else
          ! Default: no cutoff depth => 0 m
          call cvmix_put_tidal("depth_cutoff", 0, CVmix_tidal_params_user)
        end if

        if (present(max_coefficient)) then
          call cvmix_put_tidal("max_coefficient", max_coefficient,            &
                               CVmix_tidal_params_user)
        else
          call cvmix_put_tidal("max_coefficient", 50e-4_cvmix_r8,             &
                               CVmix_tidal_params_user)
        end if

      case DEFAULT
!        print*, "ERROR: ", trim(mix_scheme), " is not a valid choice for ", &
!                "tidal mixing."
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(mix_scheme)//' is not a valid choice for &
             tidal mixing.')


    end select

    if (present(old_vals)) then
      select case (trim(old_vals))
        case ("overwrite")
          call cvmix_put_tidal('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,    &
                               cvmix_tidal_params_user)
        case ("sum")
          call cvmix_put_tidal('handle_old_vals', CVMIX_SUM_OLD_AND_NEW_VALS, &
                               cvmix_tidal_params_user)
        case ("max")
          call cvmix_put_tidal('handle_old_vals', CVMIX_MAX_OLD_AND_NEW_VALS, &
                               cvmix_tidal_params_user)
        case DEFAULT
!          print*, "ERROR: ", trim(old_vals), " is not a valid option for ",   &
!                  "handling old values of diff and visc."
!          stop 1
          CALL finish(method_name,'ERROR: '//TRIM(old_vals)//' is not a valid option for &
               handling old values of diff and visc.')
      end select
    else
      call cvmix_put_tidal('handle_old_vals', CVMIX_OVERWRITE_OLD_VAL,        &
                               cvmix_tidal_params_user)
    end if

!EOC

  end subroutine cvmix_init_tidal

!BOP

! !IROUTINE: cvmix_coeffs_tidal_wrap
! !INTERFACE:

  subroutine cvmix_coeffs_tidal_wrap(CVmix_vars, CVmix_params,                &
                                     CVmix_tidal_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for tidal mixing
!  parameterizations.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type),  target, optional, intent(in) ::           &
                                            CVmix_tidal_params_user
    type(cvmix_global_params_type), intent(in) :: CVmix_params

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP
!BOC

    ! Local variables
    real(cvmix_r8), dimension(CVmix_vars%max_nlev+1) :: new_Mdiff, new_Tdiff
    type(cvmix_tidal_params_type),  pointer :: CVmix_tidal_params_in
    integer :: nlev, max_nlev

    CVmix_tidal_params_in => CVmix_tidal_params_saved
    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_in => CVmix_tidal_params_user
    end if
    nlev = CVmix_vars%nlev
    max_nlev = CVmix_vars%max_nlev

    call cvmix_coeffs_tidal(new_Mdiff, new_Tdiff,                             &
                            CVmix_vars%SqrBuoyancyFreq_iface,                 &
                            CVmix_vars%OceanDepth, CVmix_vars%SimmonsCoeff,   &
                            CVmix_vars%VertDep_iface, nlev, max_nlev,         &
                            CVMix_params, CVmix_tidal_params_user)
    call cvmix_update_wrap(CVmix_tidal_params_in%handle_old_vals, max_nlev,   &
                           Mdiff_out = CVmix_vars%Mdiff_iface,                &
                           Tdiff_out = CVmix_vars%Tdiff_iface,                &
                           new_Mdiff = new_Mdiff,                             &
                           new_Tdiff = new_Tdiff)

!EOC

  end subroutine cvmix_coeffs_tidal_wrap

!BOP

! !IROUTINE: cvmix_coeffs_tidal_low
! !INTERFACE:

  subroutine cvmix_coeffs_tidal_low(Mdiff_out, Tdiff_out, Nsqr, OceanDepth,   &
                                    SimmonsCoeff, vert_dep, nlev, max_nlev,   &
                                    CVmix_params, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Computes vertical diffusion coefficients for tidal mixing
!  parameterizations.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type),  target, optional, intent(in) ::           &
                                            CVmix_tidal_params_user
    integer,                               intent(in) :: nlev, max_nlev
    type(cvmix_global_params_type),        intent(in) :: CVmix_params
    real(cvmix_r8),                        intent(in) :: OceanDepth
    real(cvmix_r8),                        intent(in) :: SimmonsCoeff
    real(cvmix_r8), dimension(max_nlev+1), intent(in) :: Nsqr, vert_dep

! !INPUT/OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Mdiff_out
    real(cvmix_r8), dimension(max_nlev+1), intent(inout) :: Tdiff_out

!EOP
!BOC

    ! Local variables
    integer        :: k

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_coeffs_tidal_low'

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params => CVmix_tidal_params_user
    else
      CVmix_tidal_params => CVmix_tidal_params_saved
    end if

    select case (trim(CVmix_tidal_params%mix_scheme))
      case ('simmons','Simmons')
        Tdiff_out = cvmix_zero
        if (OceanDepth.ge.CVmix_tidal_params%depth_cutoff) then
          do k=1, nlev+1
            if (Nsqr(k).gt.cvmix_zero) &
              Tdiff_out(k) = SimmonsCoeff*vert_dep(k)/Nsqr(k)
            if (Tdiff_out(k).gt.CVmix_tidal_params%max_coefficient) &
              Tdiff_out(k) = CVmix_tidal_params%max_coefficient
          end do
        end if

      case DEFAULT
        ! Note: this error should be caught in cvmix_init_tidal
!        print*, "ERROR: invalid choice for type of tidal mixing."
!        stop 1
        CALL finish(method_name,'ERROR: invalid choice for type of tidal mixing.')

    end select
    Mdiff_out = CVmix_params%Prandtl*Tdiff_out

!EOC

  end subroutine cvmix_coeffs_tidal_low

!BOP

! !IROUTINE: cvmix_compute_vert_dep
! !INTERFACE:

  function cvmix_compute_vert_dep(zw, zt, nlev, CVmix_tidal_params)

! !DESCRIPTION:
!  Computes the vertical deposition function needed for Simmons et al tidal
!  mixing.
!\\
!\\
!
! !USES:
!  only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_tidal_params_type),     intent(in) :: CVmix_tidal_params
    integer,                           intent(in) :: nlev
    real(cvmix_r8), dimension(nlev+1), intent(in) :: zw
    real(cvmix_r8), dimension(nlev),   intent(in) :: zt

! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(nlev+1) :: cvmix_compute_vert_dep

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: tot_area, num, thick
    integer        :: k

    ! Compute vertical deposition
    tot_area = cvmix_zero
    cvmix_compute_vert_dep(1) = cvmix_zero
    cvmix_compute_vert_dep(nlev+1) = cvmix_zero
    do k=2,nlev
      num = -zw(k)/CVmix_tidal_params%vertical_decay_scale
      ! Simmons vertical deposition
      ! Note that it is getting normalized (divide through by tot_area)
      ! So multiplicative constants that are independent of z are omitted
      cvmix_compute_vert_dep(k) = exp(num)

      ! Compute integral of vert_dep via trapezoid rule
      ! (looks like midpoint rule, but vert_dep = 0 at z=0 and z=-ocn_depth)
      thick = zt(k-1) - zt(k)
      tot_area = tot_area + cvmix_compute_vert_dep(k)*thick
    end do
    ! Normalize vert_dep (need integral = 1.0D0)
    cvmix_compute_vert_dep = cvmix_compute_vert_dep/tot_area

!EOC

  end function cvmix_compute_vert_dep

!BOP

! !IROUTINE: cvmix_compute_Simmons_invariant_wrap
! !INTERFACE:

  subroutine cvmix_compute_Simmons_invariant_wrap(CVmix_vars, CVmix_params,   &
                                                  energy_flux,                &
                                                  CVmix_tidal_params_user)

! !DESCRIPTION:
!  Compute the time-invariant portion of the tidal mixing coefficient using
!  the Simmons, et al., scheme.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    type(cvmix_global_params_type), intent(in) :: CVmix_params
    real(cvmix_r8), intent(in) :: energy_flux
    type(cvmix_tidal_params_type),  target, optional, intent(in) ::           &
                                            CVmix_tidal_params_user

! !INPUT/OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

!EOP

    ! local variables
    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params => CVmix_tidal_params_user
    else
      CVmix_tidal_params => CVmix_tidal_params_saved
    end if

    call cvmix_put(CVmix_vars, 'SimmonsCoeff', cvmix_zero)
    call cvmix_put(CVmix_vars, 'VertDep', cvmix_zero)
    call cvmix_compute_Simmons_invariant_low(CVmix_vars%nlev,                 &
                                             energy_flux,                     &
                                             CVmix_params%FreshWaterDensity,  &
                                             CVmix_vars%SimmonsCoeff,         &
                                             CVmix_vars%VertDep_iface,        &
                                             CVmix_vars%zw_iface,             &
                                             CVmix_vars%zt_cntr,              &
                                             CVMix_tidal_params_user)

!EOC

  end subroutine cvmix_compute_Simmons_invariant_wrap

!BOP

! !IROUTINE: cvmix_compute_Simmons_invariant_low
! !INTERFACE:

  subroutine cvmix_compute_Simmons_invariant_low(nlev, energy_flux, rho,      &
                                                 SimmonsCoeff, VertDep, zw,   &
                                                 zt, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Compute the time-invariant portion of the tidal mixing coefficient using
!  the Simmons, et al., scheme.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    integer,        intent(in) :: nlev
    real(cvmix_r8), intent(in) :: energy_flux, rho
    real(cvmix_r8), dimension(:), intent(in) :: zw, zt
    type(cvmix_tidal_params_type),  target, optional, intent(in) ::           &
                                            CVmix_tidal_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8), intent(out) :: SimmonsCoeff
    real(cvmix_r8), dimension(nlev+1), intent(inout) :: VertDep

!EOP

    ! local variables
    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params => CVmix_tidal_params_user
    else
      CVmix_tidal_params => CVmix_tidal_params_saved
    end if

    SimmonsCoeff = CVmix_tidal_params%local_mixing_frac *                     &
                   CVmix_tidal_params%efficiency *                            &
                   energy_flux/rho
    VertDep = cvmix_compute_vert_dep(zw(1:nlev+1), zt(1:nlev), nlev,          &
                                     CVmix_tidal_params)
!BOC

!EOC

  end subroutine cvmix_compute_Simmons_invariant_low

!BOP

! !IROUTINE: cvmix_put_tidal_int
! !INTERFACE:

  subroutine cvmix_put_tidal_int(varname, val, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), optional, target, intent(inout) ::         &
                                              CVmix_tidal_params_user

!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_out

    CVmix_tidal_params_out => CVmix_tidal_params_saved
    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_out => CVmix_tidal_params_user
    end if

    select case (trim(varname))
      case ('old_vals', 'handle_old_vals')
        CVmix_tidal_params_out%handle_old_vals = val
      case DEFAULT
        call cvmix_put_tidal(varname, real(val,cvmix_r8),                     &
                             CVmix_tidal_params_user)
    end select

!EOC

  end subroutine cvmix_put_tidal_int

!BOP

! !IROUTINE: cvmix_put_tidal_real
! !INTERFACE:

  subroutine cvmix_put_tidal_real(varname, val, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Write a real value into a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), optional, target, intent(inout) ::         &
                                              CVmix_tidal_params_user

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_tidal_real'
!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_out

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_out => CVmix_tidal_params_user
    else
      CVmix_tidal_params_out => CVmix_tidal_params_saved
    end if

    select case (trim(varname))
      case ('efficiency')
        CVmix_tidal_params_out%efficiency = val
      case ('vertical_decay_scale')
        CVmix_tidal_params_out%vertical_decay_scale = val
      case ('max_coefficient')
        CVmix_tidal_params_out%max_coefficient = val
      case ('local_mixing_frac')
        CVmix_tidal_params_out%local_mixing_frac = val
      case ('depth_cutoff')
        CVmix_tidal_params_out%depth_cutoff = val
      case DEFAULT
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//trim(varname)//' not a valid choice!')

    end select

!EOC

  end subroutine cvmix_put_tidal_real

!BOP

! !IROUTINE: cvmix_put_tidal_str
! !INTERFACE:

  subroutine cvmix_put_tidal_str(varname, val, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Write a string into a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_tidal_params_type), optional, target, intent(inout) ::         &
                                              CVmix_tidal_params_user

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_tidal_str'
!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_out

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_out => CVmix_tidal_params_user
    else
      CVmix_tidal_params_out => CVmix_tidal_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        CVmix_tidal_params_out%mix_scheme = val
      case DEFAULT
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice!')
    end select

!EOC

  end subroutine cvmix_put_tidal_str

!BOP

! !IROUTINE: cvmix_get_tidal_real
! !INTERFACE:

  function cvmix_get_tidal_real(varname, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Returns the real value of a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_tidal_params_type), optional, target, intent(in) ::            &
                                           CVmix_tidal_params_user

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_get_tidal_real

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_get_tidal_real'

!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_in

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_in => CVmix_tidal_params_user
    else
      CVmix_tidal_params_in => CVmix_tidal_params_saved
    end if

    cvmix_get_tidal_real = cvmix_zero
    select case (trim(varname))
      case ('efficiency')
        cvmix_get_tidal_real = CVmix_tidal_params_in%efficiency
      case ('vertical_decay_scale')
        cvmix_get_tidal_real = CVmix_tidal_params_in%vertical_decay_scale
      case ('max_coefficient')
        cvmix_get_tidal_real = CVmix_tidal_params_in%max_coefficient
      case ('local_mixing_frac')
        cvmix_get_tidal_real = CVmix_tidal_params_in%local_mixing_frac
      case ('depth_cutoff')
        cvmix_get_tidal_real = CVmix_tidal_params_in%depth_cutoff
      case DEFAULT
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice!')

    end select

!EOC

  end function cvmix_get_tidal_real

!BOP

! !IROUTINE: cvmix_get_tidal_str
! !INTERFACE:

  function cvmix_get_tidal_str(varname, CVmix_tidal_params_user)

! !DESCRIPTION:
!  Returns the string value of a cvmix\_tidal\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),                                intent(in) :: varname
    type(cvmix_tidal_params_type), optional, target, intent(in) ::            &
                                           CVmix_tidal_params_user

! !OUTPUT PARAMETERS:
    character(len=cvmix_strlen) :: cvmix_get_tidal_str

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_get_tidal_str'

!EOP
!BOC

    type(cvmix_tidal_params_type), pointer :: CVmix_tidal_params_in

    if (present(CVmix_tidal_params_user)) then
      CVmix_tidal_params_in => CVmix_tidal_params_user
    else
      CVmix_tidal_params_in => CVmix_tidal_params_saved
    end if

    select case (trim(varname))
      case ('mix_scheme')
        cvmix_get_tidal_str = trim(CVmix_tidal_params_in%mix_scheme)
      case DEFAULT
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice!')
    end select

!EOC

  end function cvmix_get_tidal_str

end module cvmix_tidal
