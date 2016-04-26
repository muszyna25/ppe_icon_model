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
  USE mo_ocean_thermodyn,   ONLY: calculate_density_jmdwfg06_onColumn, calculate_density_mpiom_onColumn

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_EOS
  
  
  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE ocean_test_EOS()

    INTEGER, PARAMETER :: temperature_size = 13
    REAL(wp) :: temperature(temperature_size) = (/-1.7_wp, -1.5_wp, -1.2_wp, -1.0_wp, -0.8_wp, -0.6_wp, -0.4_wp, -0.2_wp, 0.0_wp, 0.2_wp, 0.4_wp, 1.0_wp, 2.0_wp/)
    INTEGER, PARAMETER :: salinity_size = 8
    REAL(wp) :: salinity(salinity_size) = (/33.70_wp, 33.75_wp, 33.8_wp, 33.85_wp, 33.9_wp, 33.95_wp, 34.0_wp, 34.5_wp /)
    INTEGER, PARAMETER :: pressure_size = 15
    REAL(wp) :: pressure(pressure_size) = (/0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp,  50.0_wp, 100.0_wp, 200.0_wp, 300.0_wp, 500.0_wp, 750.0_wp, 1000.0_wp,  &
     & 1250.0_wp, 1500.0_wp, 1750.0_wp, 2000.0_wp /)

    INTEGER, PARAMETER :: columnn_size = pressure_size
    REAL(wp) :: temperature_column(columnn_size)
    REAL(wp) :: salinity_column(columnn_size)
    REAL(wp) :: pressure_column(columnn_size)
    REAL(wp) :: rho_MPIOM(columnn_size)
    REAL(wp) :: rho_EOS3(columnn_size)


    INTEGER :: t,s, p

    DO t=1,temperature_size
      write(0,*) "temperature=", temperature(t)
      temperature_column = temperature(t)

      DO s=1,salinity_size
        salinity_column = salinity(s)
        write(0,*) "salinity=", salinity(s)

        pressure_column = pressure

        rho_EOS3 = calculate_density_jmdwfg06_onColumn( &
          & temperature_column,  salinity_column, pressure_column, columnn_size)

        rho_MPIOM = calculate_density_mpiom_onColumn( &
          & temperature_column,  salinity_column, pressure_column, columnn_size)

        DO p=1,pressure_size
          write(0,*) pressure(p), rho_MPIOM(p), rho_EOS3(p)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ocean_test_EOS
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_EOS
