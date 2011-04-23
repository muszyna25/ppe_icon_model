PROGRAM atm_hyb_sp_gw
  ! This program computes NLEV Level with total height =10km  for a given
  ! Brunt Vaisala frequency as needed for the gravity wave test case.
  IMPLICIT NONE
  INTEGER,  PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  CHARACTER (LEN=50)  :: hyb='atm_hyb_sp_gw'

  INTEGER,  PARAMETER :: nlev   = 20
  INTEGER,  PARAMETER :: nlevp1 = nlev+1

  REAL(wp), PARAMETER :: bruntvais = 0.01_wp ! change according to namelist parameter
                                             ! gw_brunt_vais if applicable

  REAL(wp), PARAMETER :: rd    = 287.04_wp
  REAL(wp), PARAMETER :: cpd   = 1004.64_wp
  REAL(wp), PARAMETER :: kappa = rd/cpd
  REAL(wp), PARAMETER :: grav  = 9.80616_wp

  REAL(wp), PARAMETER :: t0    = 300.0_wp
  REAL(wp), PARAMETER :: ps0   = 100000.0_wp
  REAL(wp), PARAMETER :: ztop  = 10000.0_wp  ! total height =10km
  REAL(wp), PARAMETER :: c     = 1.0_wp

  INTEGER  :: i
  REAL(wp) :: sot,zdel,z,pz
  REAL(wp) :: a(nlevp1),b(nlevp1),eta(nlevp1)

  sot   = (grav/bruntvais)**2/cpd/t0
  zdel  = ztop/(nlevp1-1)
  DO i= 1,nlevp1
    z   = zdel * (nlevp1-i)
    pz  = ps0*(1.0_wp-sot+sot*EXP(-bruntvais*bruntvais*z/grav))**(1.0_wp/kappa)
    eta(i)=pz/ps0
    b(i) = ( (eta(i)-eta(1))/(1.0_wp-eta(1)))**c
    a(i) = ps0*(eta(i)-b(i))
  ENDDO

  OPEN(UNIT=10,FILE=hyb)
  WRITE(10,'(a)') '    k     vct_a(k) [Pa]   vct_b(k) []'
  DO i=1,nlevp1
    WRITE(10,'(i5,f18.10,f14.10)') i,a(i),b(i)
  ENDDO
  WRITE(10,'(a)') '====================================='
  WRITE(10,'(a)') 'Source:'
  WRITE(10,'(a)') '- atm_hyb_sp_gw.f90'
  WRITE(10,'(a)') 'Comments:'
  WRITE(10,'(a)') '- for gravity-wave test'
  WRITE(10,'(a)') 'References:'
  WRITE(10,'(a)') '- n.a.'

  CLOSE(10)
END PROGRAM atm_hyb_sp_gw
