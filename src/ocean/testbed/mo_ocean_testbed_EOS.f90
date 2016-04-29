!>
!! Contains the main stepping method_name the 3-dim hydrostatic ocean model.
!!
!! @author Leonidas Linardakis, MPI
!!
!! @par Revision History
!
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_testbed_EOS
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, sp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ocean_thermodyn,   ONLY: calculate_density_jmdwfg06_onColumn, calculate_density_mpiom_onColumn, convert_insitu2pot_temp_func
  USE mo_physical_constants,  ONLY: grav, sal_ref, b_s, &
    & sitodbar, sfc_press_bar
  USE mo_ocean_nml,           ONLY: OceanReferenceDensity, ReferencePressureIndbars

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_EOS
  
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE ocean_test_EOS()

    INTEGER, PARAMETER :: temperature_size = 14
!     REAL(wp) :: temperature(temperature_size) = (/-1.7_wp, -1.5_wp, -1.2_wp, -1.0_wp, -0.8_wp, -0.6_wp, -0.4_wp, -0.2_wp, 0.0_wp, 0.2_wp, 0.4_wp, 1.0_wp, 2.0_wp/)
    REAL(wp) :: temperature(temperature_size) = (/-2.0_wp, 0.0_wp, 2.0_wp, 4.0_wp, 7.0_wp, 10.0_wp, 13.0_wp, 16.0_wp, 19.0_wp, 22.0_wp, 25.4_wp, 28.0_wp, 31.0_wp, &
      & 35.0_wp/)
    INTEGER, PARAMETER :: salinity_size = 1
!     REAL(wp) :: salinity(salinity_size) = (/33.70_wp, 33.75_wp, 33.8_wp, 33.85_wp, 33.9_wp, 33.95_wp, 34.0_wp, 34.5_wp /)
    REAL(wp) :: salinity(salinity_size) = (/35._wp /)
    INTEGER, PARAMETER :: depth_size = 20
    REAL(wp) :: depth(depth_size) = (/0.0_wp, 100.0_wp,  200.0_wp, 300.0_wp, 500.0_wp, 750.0_wp, 1000.0_wp, 1500.0_wp, 2000.0_wp,  &
     & 2500.0_wp, 3000.0_wp, 3500.0_wp, 4000.0_wp, 4500.0_wp, 5000.0_wp, 6000.0_wp, 7000.0_wp, 8000._wp, 9000._wp, 10000._wp /)

    INTEGER, PARAMETER :: columnn_size = temperature_size ! depth_size
    REAL(wp) :: temperature_column(columnn_size)
    REAL(wp) :: salinity_column(columnn_size)
    REAL(wp) :: pressure_column(columnn_size)
    REAL(wp) :: rho_MPIOM(columnn_size)
    REAL(wp) :: rho_EOS3(columnn_size)

    INTEGER :: t,s,p,l

!     DO t=1,temperature_size
!       write(0,*) "temperature=", temperature(t)
!       temperature_column = temperature(t)

      DO p=1, depth_size
        write(0,*) "depth=", depth(p)

      DO s=1,salinity_size
        salinity_column = salinity(s)
        write(0,*) "salinity=", salinity(s)

!         pressure_column = depth * OceanReferenceDensity * grav * sitodbar

        DO l=1,columnn_size
!           temperature_column(l) = convert_insitu2pot_temp_func(temperature(t), salinity(s), depth(l) / 10.0_wp)
          temperature_column(l) = convert_insitu2pot_temp_func(temperature(l), salinity(s), depth(p) / 10.0_wp)
        ENDDO

        pressure_column = depth(p) !* ReferencePressureIndbars
        rho_EOS3 = calculate_density_jmdwfg06_onColumn( &
          & temperature_column,  salinity_column, pressure_column, columnn_size)

        pressure_column = depth(p) / 10.0_wp !* OceanReferenceDensity * sitodbar
        rho_MPIOM = calculate_density_mpiom_onColumn( &
          & temperature_column,  salinity_column, pressure_column, columnn_size)

        

        DO l=1,columnn_size
          write(0,*) temperature(l), rho_MPIOM(l), rho_EOS3(l), (rho_EOS3(l)-rho_MPIOM(l)) ! / rho_MPIOM(p)

        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE ocean_test_EOS
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_EOS
