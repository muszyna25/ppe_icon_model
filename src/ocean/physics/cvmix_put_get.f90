module cvmix_put_get

!BOP
!\newpage
! !MODULE: cvmix_put_get
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to pack data into the cvmix datatypes
!  (allocating memory as necessary) and then unpack the data out. If we switch
!  to pointers, the pack will just point at the right target and the unpack
!  will be un-necessary.
!\\
!\\

! !USES:

   use cvmix_kinds_and_types, only : cvmix_r8,                  &
                                     cvmix_data_type,           &
                                     cvmix_global_params_type
   use cvmix_utils,           only : cvmix_att_name


   USE mo_exception,          ONLY: finish, message

!EOP

  implicit none
  private
  save

!BOP

! !PUBLIC MEMBER FUNCTIONS:
  public :: cvmix_put

  interface cvmix_put
    module procedure cvmix_put_int
    module procedure cvmix_put_real
    module procedure cvmix_put_real_1D
    module procedure cvmix_put_global_params_int
    module procedure cvmix_put_global_params_real
  end interface cvmix_put

  CHARACTER(LEN=*), PARAMETER :: module_name = 'cvmix_put_get'


!EOP

contains

!BOP

! !IROUTINE: cvmix_put_int
! !INTERFACE:

  subroutine cvmix_put_int(CVmix_vars, varname, val, nlev_in)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
    integer,          optional, intent(in) :: nlev_in

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_int'

!EOP
!BOC

    ! Local variables
    integer :: nlev

    if (present(nlev_in)) then
      nlev = nlev_in
    else
      nlev = CVmix_vars%max_nlev
    end if

    if ((trim(varname).ne.'nlev').and.(nlev.eq.-1)) then
!      print*, "ERROR: you must specify the number of levels before ", &
!              "you can pack data into a cvmix_data_type!"
!      print*, "You tried to set ", trim(varname)
!      stop 1
      CALL message(method_name,'ERROR: you must specify the number of levels before &
                                  you can pack data into a cvmix_data_type!')
      CALL finish(method_name,'You tried to set '//trim(varname))
    end if

    select case (trim(cvmix_att_name(varname)))
      case ('nlev')
        CVmix_vars%nlev = val
        if (CVmix_vars%max_nlev.eq.-1) then
          CVmix_vars%max_nlev= val
        end if
      case ('max_nlev')
        CVmix_vars%max_nlev = val
        if (CVmix_vars%nlev.eq.-1) then
          CVmix_vars%nlev= val
        end if
      case default
        ! All other scalars are real(cvmix_r8)
        call cvmix_put_real(CVmix_vars, varname, real(val,cvmix_r8), nlev_in)
    end select
!EOC

  end subroutine cvmix_put_int

!BOP

! !IROUTINE: cvmix_put_real
! !INTERFACE:

  subroutine cvmix_put_real(CVmix_vars, varname, val, nlev_in)

! !DESCRIPTION:
!  Write a real value into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),           intent(in) :: varname
    real(cvmix_r8),             intent(in) :: val
    integer,          optional, intent(in) :: nlev_in

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_real'

!EOP
!BOC

    ! Local variables
    integer :: nlev

    if (present(nlev_in)) then
      nlev = nlev_in
    else
      nlev = CVmix_vars%max_nlev
    end if

    if (nlev.eq.-1) then
!      print*, "ERROR: you must specify the number of levels before ", &
!              "you can pack data into a cvmix_data_type!"
!      print*, "You tried to set ", trim(varname)
!      stop 1
      CALL message(method_name,'ERROR: you must specify the number of levels before &
                                  you can pack data into a cvmix_data_type!')
      CALL finish(method_name,'You tried to set '//trim(varname))
    end if

    select case (trim(cvmix_att_name(varname)))
      case ('OceanDepth')
        CVmix_vars%OceanDepth = val
      case ('BoundaryLayerDepth')
        CVmix_vars%BoundaryLayerDepth = val
      case ('SeaSurfaceHeight')
        CVmix_vars%SeaSurfaceHeight = val
      case ('SurfaceFriction')
        CVmix_vars%SurfaceFriction = val
      case ("SurfaceBuoyancyForcing")
        CVmix_vars%SurfaceBuoyancyForcing = val
      case ("lat")
        CVmix_vars%lat = val
      case ("lon")
        CVmix_vars%lon = val
      case ("Coriolis")
        CVmix_vars%Coriolis = val
      case ("kOBL_depth")
        CVmix_vars%kOBL_depth = val
      ! QL, 150610
      case ("LangmuirEnhancementFactor")
        CVmix_vars%LangmuirEnhancementFactor = val
      case ("SurfaceStokesDrift")
        CVmix_vars%SurfaceStokesDrift = val
      case ('SimmonsCoeff')
        CVmix_vars%SimmonsCoeff = val
! IDEMIX / TKE
      case ('forc_iw_bottom')
        CVmix_vars%forc_iw_bottom = val
      case ('forc_iw_surface')
        CVmix_vars%forc_iw_surface = val
      case ('forc_tke_surf')
        CVmix_vars%forc_tke_surf = val
      case ("forc_rho_surf")
        CVmix_vars%forc_rho_surf = val
      case ('dtime')
        CVmix_vars%dtime = val

      case ("dzw")
!        print*, "WARNING: you are setting the cell midpoint to midpoint ",    &
!                "distance in all levels to a constant value"
        if (.not.associated(CVmix_vars%dzw)) then
          allocate(CVmix_vars%dzw(nlev+1))
        end if
        CVmix_vars%dzw(:) = val
      case ("Mdiff_iface")
        if (.not.associated(CVmix_vars%Mdiff_iface)) then
          allocate(CVmix_vars%Mdiff_iface(nlev+1))
        end if
        CVmix_vars%Mdiff_iface(:) = val
      case ("Tdiff_iface")
        if (.not.associated(CVmix_vars%Tdiff_iface)) then
          allocate(CVmix_vars%Tdiff_iface(nlev+1))
        end if
        CVmix_vars%Tdiff_iface(:) = val
      case ("Sdiff_iface")
        if (.not.associated(CVmix_vars%Sdiff_iface)) then
          allocate(CVmix_vars%Sdiff_iface(nlev+1))
        end if
        CVmix_vars%Sdiff_iface(:) = val
      case ("ShearRichardson_iface")
!        print*, "WARNING: you are setting the Richardson number in all ",     &
!                "levels to a constant value"
        if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
          allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
        end if
        CVmix_vars%ShearRichardson_iface(:) = val
      case ("SqrBuoyancyFreq_iface")
!        print*, "WARNING: you are setting the buoyancy in all levels to a ",  &
!                "constant value"
        if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
          allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
        end if
        CVmix_vars%SqrBuoyancyFreq_iface(:) = val
      case ("kpp_nonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
        CVmix_vars%kpp_Snonlocal_iface(:) = val
      case ("kpp_Tnonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
      case ("kpp_Snonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Snonlocal_iface(:) = val

      case ("dzt")
!        print*, "WARNING: you are setting the cell thickness in all levels ", &
!                "to a constant value"
        if (.not.associated(CVmix_vars%dzt)) then
          allocate(CVmix_vars%dzt(nlev))
        end if
        CVmix_vars%dzt(:) = val
      case ("WaterDensity_cntr")
!        print*, "WARNING: you are setting the density in all levels to a ",   &
!                "constant value"
        if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
          allocate(CVmix_vars%WaterDensity_cntr(nlev))
        end if
        CVmix_vars%WaterDensity_cntr(:) = val
      case ("AdiabWaterDensity_cntr")
!        print*, "WARNING: you are setting the adiabatic density in all ",     &
!                "levels to a constant value"
        if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
          allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
        end if
        CVmix_vars%AdiabWaterDensity_cntr(:) = val
      case ("BulkRichardson_cntr")
!        print*, "WARNING: you are setting the bulk Richardson number in all", &
!                " levels to a constant value"
        if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
          allocate(CVmix_vars%BulkRichardson_cntr(nlev))
        end if
        CVmix_vars%BulkRichardson_cntr(:) = val
      case ('strat_param_num')
!        print*, "WARNING: you are setting the numerator of the ",             &
!                "stratification parameter in all levels to a constant value"
        if (.not.associated(CVmix_vars%strat_param_num)) then
          allocate(CVmix_vars%strat_param_num(nlev))
        end if
        CVmix_vars%strat_param_num(:) = val
      case ('strat_param_denom')
!        print*, "WARNING: you are setting the denominator of the ",           &
!                "stratification parameter in all levels to a constant value"
        if (.not.associated(CVmix_vars%strat_param_denom)) then
          allocate(CVmix_vars%strat_param_denom(nlev))
        end if
        CVmix_vars%strat_param_denom(:) = val
      case ("VertDep_iface")
        if (.not.associated(CVmix_vars%VertDep_iface)) then
          allocate(CVmix_vars%VertDep_iface(nlev+1))
        end if
        CVmix_vars%VertDep_iface(:) = val

      case default
!       print*, "ERROR: ", trim(varname), " not a valid choice for cvmix_put!"
!       stop 1
        CALL finish(method_name,'ERROR: '//trim(varname)//' not a valid choice for cvmix_put!')
    end select
!EOC

  end subroutine cvmix_put_real

!BOP

! !IROUTINE: cvmix_put_real_1D
! !INTERFACE:

  subroutine cvmix_put_real_1D(CVmix_vars, varname, val, nlev_in)

! !DESCRIPTION:
!  Write an array of real values into a cvmix\_data\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*),             intent(in) :: varname
    real(cvmix_r8), dimension(:), intent(in) :: val
    integer,        optional,     intent(in) :: nlev_in

! !OUTPUT PARAMETERS:
    type(cvmix_data_type), intent(inout) :: CVmix_vars
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_real_1D'
!EOP
!BOC

    ! Local variables
    integer :: nlev

    if (present(nlev_in)) then
      nlev = nlev_in
    else
      nlev = CVmix_vars%max_nlev
    end if

    if (nlev.eq.-1) then
!      print*, "ERROR: you must specify the number of levels before ", &
!              "you can pack data into a cvmix_data_type!"
!      print*, "You tried to set ", trim(varname)
!      stop 1
      CALL message(method_name,'ERROR: you must specify the number of levels before &
                                  you can pack data into a cvmix_data_type!')
      CALL finish(method_name,'You tried to set '//trim(varname))
    end if

    select case (trim(cvmix_att_name(varname)))
      case ("zw_iface")
        if (.not.associated(CVmix_vars%zw_iface)) then
          allocate(CVmix_vars%zw_iface(nlev+1))
        end if
        CVmix_vars%zw_iface(:) = val
      case ("dzw")
        if (.not.associated(CVmix_vars%dzw)) then
          allocate(CVmix_vars%dzw(nlev+1))
        end if
        CVmix_vars%dzw(:) = val
      case ("Mdiff_iface")
        if (.not.associated(CVmix_vars%Mdiff_iface)) then
          allocate(CVmix_vars%Mdiff_iface(nlev+1))
        end if
        CVmix_vars%Mdiff_iface(:) = val
      case ("Tdiff_iface")
        if (.not.associated(CVmix_vars%Tdiff_iface)) then
          allocate(CVmix_vars%Tdiff_iface(nlev+1))
        end if
        CVmix_vars%Tdiff_iface(:) = val
      case ("Sdiff_iface")
        if (.not.associated(CVmix_vars%Sdiff_iface)) then
          allocate(CVmix_vars%Sdiff_iface(nlev+1))
        end if
        CVmix_vars%Sdiff_iface(:) = val
      case ("ShearRichardson_iface")
        if (.not.associated(CVmix_vars%ShearRichardson_iface)) then
          allocate(CVmix_vars%ShearRichardson_iface(nlev+1))
        end if
        CVmix_vars%ShearRichardson_iface(:) = val
      case ("SqrBuoyancyFreq_iface")
        if (.not.associated(CVmix_vars%SqrBuoyancyFreq_iface)) then
          allocate(CVmix_vars%SqrBuoyancyFreq_iface(nlev+1))
        end if
        CVmix_vars%SqrBuoyancyFreq_iface(:) = val
      case ("kpp_nonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
        CVmix_vars%kpp_Snonlocal_iface(:) = val
      case ("kpp_Tnonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Tnonlocal_iface)) then
          allocate(CVmix_vars%kpp_Tnonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Tnonlocal_iface(:) = val
      case ("kpp_Snonlocal_iface")
        if (.not.associated(CVmix_vars%kpp_Snonlocal_iface)) then
          allocate(CVmix_vars%kpp_Snonlocal_iface(nlev+1))
        end if
        CVmix_vars%kpp_Snonlocal_iface(:) = val
      case ("VertDep_iface")
        if (.not.associated(CVmix_vars%VertDep_iface)) then
          allocate(CVmix_vars%VertDep_iface(nlev+1))
        end if
        CVmix_vars%VertDep_iface(:) = val
      case ("zt_cntr")
        if (.not.associated(CVmix_vars%zt_cntr)) then
          allocate(CVmix_vars%zt_cntr(nlev))
        end if
        CVmix_vars%zt_cntr(:) = val
      case ("dzt")
        if (.not.associated(CVmix_vars%dzt)) then
          allocate(CVmix_vars%dzt(nlev))
        end if
        CVmix_vars%dzt(:) = val
      case ("WaterDensity_cntr")
        if (.not.associated(CVmix_vars%WaterDensity_cntr)) then
          allocate(CVmix_vars%WaterDensity_cntr(nlev))
        end if
        CVmix_vars%WaterDensity_cntr(:) = val
      case ("AdiabWaterDensity_cntr")
        if (.not.associated(CVmix_vars%AdiabWaterDensity_cntr)) then
          allocate(CVmix_vars%AdiabWaterDensity_cntr(nlev))
        end if
        CVmix_vars%AdiabWaterDensity_cntr(:) = val
      case ("BulkRichardson_cntr")
        if (.not.associated(CVmix_vars%BulkRichardson_cntr)) then
          allocate(CVmix_vars%BulkRichardson_cntr(nlev))
        end if
        CVmix_vars%BulkRichardson_cntr(:) = val
      case ('strat_param_num')
        if (.not.associated(CVmix_vars%strat_param_num)) then
          allocate(CVmix_vars%strat_param_num(nlev))
        end if
        CVmix_vars%strat_param_num(:) = val
      case ('strat_param_denom')
        if (.not.associated(CVmix_vars%strat_param_denom)) then
          allocate(CVmix_vars%strat_param_denom(nlev))
        end if
        CVmix_vars%strat_param_denom(:) = val
      case ("Vx_cntr")
        if (.not.associated(CVmix_vars%Vx_cntr)) then
          allocate(CVmix_vars%Vx_cntr(nlev))
        end if
        CVmix_vars%Vx_cntr(:) = val
      case ("Vy_cntr")
        if (.not.associated(CVmix_vars%Vy_cntr)) then
          allocate(CVmix_vars%Vy_cntr(nlev))
        end if
        CVmix_vars%Vy_cntr(:) = val

! IDEMIX / TKE
      case ("KappaM_iface")
        if (.not.associated(CVmix_vars%KappaM_iface)) then
          allocate(CVmix_vars%KappaM_iface(nlev+1))
        end if
        CVmix_vars%KappaM_iface(:) = val

      case ("KappaH_iface")
        if (.not.associated(CVmix_vars%KappaH_iface)) then
          allocate(CVmix_vars%KappaH_iface(nlev+1))
        end if
        CVmix_vars%KappaH_iface(:) = val

      case ("Ssqr_iface")
        if (.not.associated(CVmix_vars%Ssqr_iface)) then
          allocate(CVmix_vars%Ssqr_iface(nlev+1))
        end if
        CVmix_vars%Ssqr_iface(:) = val

      case ("Nsqr_iface")
        if (.not.associated(CVmix_vars%Nsqr_iface)) then
          allocate(CVmix_vars%Nsqr_iface(nlev+1))
        end if
        CVmix_vars%Nsqr_iface(:) = val

      case ("TKE")
        if (.not.associated(CVmix_vars%TKE)) then
          allocate(CVmix_vars%TKE(nlev+1))
        end if
        CVmix_vars%TKE(:) = val

      case ("E_iw")
        if (.not.associated(CVmix_vars%E_iw)) then
          allocate(CVmix_vars%E_iw(nlev+1))
        end if
        CVmix_vars%E_iw(:) = val

      case ("iw_diss")
        if (.not.associated(CVmix_vars%iw_diss)) then
          allocate(CVmix_vars%iw_diss(nlev+1))
        end if
        CVmix_vars%iw_diss(:) = val

      case ("alpha_c")
        if (.not.associated(CVmix_vars%alpha_c)) then
          allocate(CVmix_vars%alpha_c(nlev+1))
        end if
        CVmix_vars%alpha_c(:) = val

      case default
!        print*, "ERROR: ", trim(varname), " not a valid choice for cvmix_put!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice for cvmix_put!')

    end select

!EOC

  end subroutine cvmix_put_real_1D

!BOP

! !IROUTINE: cvmix_put_global_params_int
! !INTERFACE:

  subroutine cvmix_put_global_params_int(CVmix_params, varname, val)

! !DESCRIPTION:
!  Write an integer value into a cvmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: val

! !OUTPUT PARAMETERS:
    type (cvmix_global_params_type), intent(inout) :: CVmix_params

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_global_params_int'

!EOP
!BOC

    select case (trim(varname))
      case ('max_nlev')
        CVmix_params%max_nlev = val

      case default
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice!')

    end select
!EOC

  end subroutine cvmix_put_global_params_int

!BOP

! !IROUTINE: cvmix_put_global_params_real
! !INTERFACE:

  subroutine cvmix_put_global_params_real(CVmix_params, varname, val)

! !DESCRIPTION:
!  Write a real value into a cvmix\_global\_params\_type variable.
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    character(len=*), intent(in) :: varname
    real(cvmix_r8),   intent(in) :: val

! !OUTPUT PARAMETERS:
    type(cvmix_global_params_type), intent(inout) :: CVmix_params
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':cvmix_put_global_params_type'
!EOP
!BOC

    select case (trim(varname))
      case ('prandtl','Prandtl')
        CVmix_params%prandtl = val
      case ('fw_rho','FreshWaterDensity')
        CVmix_params%FreshWaterDensity = val
      case ('sw_rho','SaltWaterDensity')
        CVmix_params%SaltWaterDensity = val
      ! QL, 160708, add gravity
      case ('g','Gravity')
        CVmix_params%Gravity = val
      case default
!        print*, "ERROR: ", trim(varname), " not a valid choice!"
!        stop 1
        CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice!')

    end select
!EOC

  end subroutine cvmix_put_global_params_real

end module cvmix_put_get

