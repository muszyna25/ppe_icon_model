!>
!! @brief Read and apply monthly aerosol optical properties of S. Kinne
!! from yearly files.
!!
!! @author B. Stevens, K. Peters, J.S. Rast (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_bc_aeropt_splumes

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_read_interface,       ONLY: openInputFile, read_1D, &
                                   & read_bcast_real_2D, read_bcast_real_3D, &
                                   & closeFile
  USE mo_datetime,             ONLY: t_datetime
  USE mo_physical_constants,   ONLY: rgrav
  USE mo_model_domain,         ONLY: p_patch
  USE mo_psrad_srtm_setup,     ONLY: &
      &  sw_wv1 => wavenum1     ,&     !< smallest wave number in each of the sw bands
      &  sw_wv2 => wavenum2            !< largest wave number in each of the sw bands
  USE mo_rrtm_params,          ONLY: &
      & jpb1                    ,&     !< index for lower sw band
      & jpb2                           !< index for upper sw band
  USE mo_math_constants,       ONLY: rad2deg

!!$, on_cells, &
!!$    &                                t_stream_id, read_0D_real, read_3D_time

  IMPLICIT NONE

  PRIVATE
  PUBLIC                  :: setup_bc_aeropt_splumes, add_bc_aeropt_splumes

  INTEGER, PARAMETER      ::     &
       nplumes   = 9            ,& !< Number of plumes
       nfeatures = 2            ,& !< Number of features per plume
       ntimes    = 52           ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251             !< Number of years of available forcing
  REAL(wp), POINTER ::                    &
       plume_lat   (:)  ,& !< (nplumes) latitude where plume maximizes
       plume_lon   (:)  ,& !< (nplumes) longitude where plume maximizes
       beta_a      (:)  ,& !< (nplumes) parameter a for beta function 
                           !< vertical profile
       beta_b      (:)  ,& !< (nplumes) parameter b for beta function 
                           !< vertical profile
       aod_spmx    (:)  ,& !< (nplumes) aod at 550 for simple plume (maximum)
       aod_fmbg    (:)  ,& !< (nplumes) aod at 550 for fine mode 
                           !< background (for twomey effect)
       asy550      (:)  ,& !< (nplumes) asymmetry parameter for plume at 550nm
       ssa550      (:)  ,& !< (nplumes) single scattering albedo for 
                           !< plume at 550nm
       angstrom    (:)  ,& !< (nplumes) angstrom parameter for plume 
       sig_lon_E   (:,:),& !< (nfeatures,nplumes) Eastward extent of 
                           !< plume feature
       sig_lon_W   (:,:),& !< (nfeatures,nplumes) Westward extent of 
                           !< plume feature
       sig_lat_E   (:,:),& !< (nfeatures,nplumes) Southward extent of 
                           !< plume feature
       sig_lat_W   (:,:),& !< (nfeatures,nplumes) Northward extent of 
                           !< plume feature
       theta       (:,:),& !< (nfeatures,nplumes) Rotation angle of feature
       ftr_weight  (:,:),& !< (nfeatures,nplumes) Feature weights = 
                           !< (nfeatures + 1) to account for BB background
       year_weight (:,:)    ,& !< (nyear,nplumes) Yearly weight for plume
       ann_cycle   (:,:,:)     !< (nfeatures,ntimes,nplumes) annual cycle for feature
  REAL(wp)                 :: &
    time_weight (nfeatures,nplumes)       !< Time-weights to account for BB background

  CHARACTER(LEN=256)       :: cfname
  LOGICAL                  :: sp_initialized

  CONTAINS

  ! -----------------------------------------------------------------
  ! SETUP_BC_AEROPT_SPLUMES:  This subroutine should be called at initialization to 
  !            read the netcdf data that describes the simple plume
  !            climatology.  The information needs to be either read 
  !            by each processor or distributed to processors.
  !
  SUBROUTINE setup_bc_aeropt_splumes
    !
    ! ---------- 
    !
    INTEGER           :: ifile_id

    cfname='MACv2.0-SP_v1.nc'
    ifile_id=openInputFile(cfname)

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lat',&
                       & return_pointer=plume_lat, file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lon',&
                       & return_pointer=plume_lon, file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_a',   &
                       & return_pointer=beta_a, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_b',   &
                       & return_pointer=beta_b, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_spmx', &
                       & return_pointer=aod_spmx, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_fmbg', &
                       & return_pointer=aod_fmbg, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='ssa550',   &
                       & return_pointer=ssa550, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='asy550',   &
                       & return_pointer=asy550, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='angstrom', &
                       & return_pointer=angstrom, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_W',&
                       & return_pointer=sig_lat_W, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_E',&
                       & return_pointer=sig_lat_E, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_W',&
                       & return_pointer=sig_lon_W, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_E',&
                       & return_pointer=sig_lon_E, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='theta',    &
                       & return_pointer=theta, file_name=cfname,             &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='ftr_weight',&
                       & return_pointer=ftr_weight, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),                &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='setup_bc_aeropt_splumes'              )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='year_weight',&
                       & return_pointer=year_weight, file_name=cfname,         &
                       & variable_dimls=(/nyears,nplumes/),                    &
                       & module_name='mo_bc_aeropt_splumes',                   &
                       & sub_prog_name='setup_bc_aeropt_splumes'               )
    CALL read_3d_wrapper(ifile_id=ifile_id,        variable_name='ann_cycle', &
                       & return_pointer=ann_cycle, file_name=cfname,          &
                       & variable_dimls=(/nfeatures,ntimes,nplumes/),         &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='setup_bc_aeropt_splumes'               )
    CALL closeFile(ifile_id)
    sp_initialized = .TRUE.
    RETURN
  END SUBROUTINE SETUP_BC_AEROPT_SPLUMES
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight. 
  !
  SUBROUTINE set_time_weight(year_fr)
    !
    ! ---------- 
    !
    REAL(wp), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    INTEGER          ::  &
         idate(8)       ,& !< integer array of system clock date yyyy, mm, dd
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100) 
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    CHARACTER(LEN=1024) :: message
    !
    ! ---------- 
    !
    iyear = FLOOR(year_fr) - 1850 + 1
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
    END DO    
    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properites based on the
  ! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
  ! optical properties on a host models vertical grid. 
  !
  SUBROUTINE sp_aop_profile( &
     & nlevels        ,ncol           ,ncol_max       ,lambda         ,oro            ,lon            , &
     & lat            ,year_fr        ,z              ,dz             ,dNovrN         ,aod_prof       , &
     & ssa_prof       ,asy_prof       )
    !
    ! ---------- 
    !
    INTEGER, INTENT(IN)        :: &
       & nlevels,                 & !< number of levels
       & ncol,                    & !< number of columns
       & ncol_max                   !< first dimension of 2d-vars as declared in calling (sub)program

    REAL(wp), INTENT(IN)       :: &
       & lambda,                  & !< wavelength
       & year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
       & oro(ncol),               & !< orographic height (m)
       & lon(ncol),               & !< longitude in degrees E
       & lat(ncol),               & !< latitude in degrees N
       & z (ncol_max,nlevels),    & !< height above sea-level (m)
       & dz(ncol_max,nlevels)       !< level thickness (difference between half levels)

    REAL(wp), INTENT(OUT)      ::     &
       & dNovrN(ncol)               , & !< anthropogenic incroment to cloud drop number concentration 
       & aod_prof(ncol_max,nlevels) , & !< profile of aerosol optical depth
       & ssa_prof(ncol_max,nlevels) , & !< profile of single scattering albedo
       & asy_prof(ncol_max,nlevels)     !< profile of asymmetry parameter

    INTEGER                    :: iplume, icol, k

    REAL(wp)                   ::  &
       & eta(ncol_max,nlevels),    & !< normalized height (by 15 km)
       & z_beta(ncol_max,nlevels), & !< profile for scaling column optical depth
       & prof(ncol_max,nlevels),   & !< scaled profile (by beta function)
       & beta_sum(ncol),           & !< vertical sum of beta function
       & ssa(ncol),                & !< aerosol optical depth 
       & asy(ncol),                & !< aerosol optical depth 
       & cw_an(ncol),              & !< column weight for simple plume (anthropogenic) aod at 550 nm
       & cw_bg(ncol),              & !< column weight for fine-mode background aod at 550 nm
       & caod_sp(ncol),            & !< column simple plume (anthropogenic) aod at 550 nm
       & caod_bg(ncol),            & !< column fine-mode background aod at 550 nm
       & a_plume1,                 & !< gaussian longitude factor for feature 1
       & a_plume2,                 & !< gaussian longitude factor for feature 2
       & b_plume1,                 & !< gaussian latitude factor for feature 1
       & b_plume2,                 & !< gaussian latitude factor for feature 2
       & delta_lat,                & !< latitude offset
       & delta_lon,                & !< longitude offset
       & lon1,                     & !< rotated longitude for feature 1
       & lat1,                     & !< rotated latitude for feature 2
       & lon2,                     & !< rotated longitude for feature 1
       & lat2,                     & !< rotated latitude for feature 2
       & f1,                       & !< contribution from feature 1
       & f2,                       & !< contribution from feature 2
       & aod_550,                  & !< aerosol optical depth at 550nm
       & aod_lmd,                  & !< aerosol optical depth at input wavelength
       & lfactor                     !< factor to compute wavelength dependence of optical properties
    !
    ! ---------- 
    !
    ! initialize input data (by calling setup at first instance) 
    !
    IF (.NOT.sp_initialized) CALL setup_bc_aeropt_splumes
    !
    ! get time weights
    !
    CALL set_time_weight(year_fr)
    !
    ! initialize variables, including output
    !
    DO k=1,nlevels
      DO icol=1,ncol
        aod_prof(icol,k) = 0.0
        ssa_prof(icol,k) = 0.0
        asy_prof(icol,k) = 0.0
        z_beta(icol,k)   = MERGE(1.0_wp, 0.0_wp, z(icol,k) >= oro(icol))
        eta(icol,k)      = MAX(0.0_wp,MIN(1.0_wp,z(icol,k)/15000.))
      END DO
    END DO
    DO icol=1,ncol
      dNovrN(icol)   = 1.0_wp
      caod_sp(icol)  = 0.00_wp
      caod_bg(icol)  = 0.02_wp
    END DO
    !
    ! sum contribution from plumes to construct composite profiles of aerosol otpical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO icol=1,ncol
        beta_sum(icol) = 0._wp
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1._wp) * (1._wp-eta(icol,k))**(beta_b(iplume)-1._wp))*dz(icol,k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k)
        END DO
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = prof(icol,k) / beta_sum(icol) * z_beta(icol,k)
        END DO
      END DO
      !
      ! calculate plume weights
      !
      DO icol=1,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat = lat(icol) - plume_lat(iplume)
        delta_lon = lon(icol) - plume_lon(iplume)
        delta_lon = MERGE ( delta_lon-SIGN(360._wp,delta_lon) , delta_lon , ABS(delta_lon) > 180._wp )
        a_plume1  = 0.5_wp / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0)**2)
        b_plume1  = 0.5_wp / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0)**2)
        a_plume2  = 0.5_wp / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0)**2)
        b_plume2  = 0.5_wp / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-1._wp* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-1._wp* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))) 

        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)  
        cw_bg(icol) = f1 * aod_fmbg(iplume) + f2 * aod_fmbg(iplume) 
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_wp,700.0_wp/lambda)
        ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
        asy(icol) =  asy550(iplume) * SQRT(lfactor)
      END DO
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the anstrom parameter. 
      !      
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_wp))
      DO k=1,nlevels
        DO icol = 1,ncol
          aod_550          = prof(icol,k)     * cw_an(icol)
          aod_lmd          = aod_550          * lfactor
          caod_sp(icol)    = caod_sp(icol)    + prof(icol,k) * cw_an(icol)
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)
          asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
          ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
          aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
        END DO
      END DO
    END DO
    !
    ! complete optical depth weighting
    !
    DO k=1,nlevels
      DO icol = 1,ncol
        asy_prof(icol,k) = MERGE(asy_prof(icol,k)/ssa_prof(icol,k), 0.0_wp, ssa_prof(icol,k) > TINY(1._wp))
        ssa_prof(icol,k) = MERGE(ssa_prof(icol,k)/aod_prof(icol,k), 1.0_wp, aod_prof(icol,k) > TINY(1._wp))
      END DO
    END DO
    !
    ! calcuate effective radius normalization (divisor) factor
    !
    DO icol=1,ncol
      dNovrN(icol) = LOG((1000.0_wp * (caod_sp(icol) + caod_bg(icol))) + 3.0_wp)/LOG((1000.0_wp * caod_bg(icol)) + 3.0_wp)
    END DO

    RETURN
  END SUBROUTINE sp_aop_profile
  ! ------------------------------------------------------------------------------------------------------------------------
  ! ADD_BC_AEROPT_SPLUMES:  This subroutine provides the interface to simple plume (sp) fit to the MPI Aerosol Climatology (Version 2).
  ! It does so by collecting or deriving spatio-temporal information and calling the simple plume aerosol subroutine and
  ! incrementing the background aerosol properties (and effective radisu) with the anthropogenic plumes.
  !
  SUBROUTINE add_bc_aeropt_splumes                                                ( &
     & jg             ,kproma         ,kbdim          ,klev           ,krow        ,&
     & nb_sw          ,datetime       ,geoi           ,geom           ,oromea      ,&
     & aod_sw_vr      ,ssa_sw_vr      ,asy_sw_vr      ,x_cdnc                      )
    !
    ! --- 0.1 Variables passed through argument list
    INTEGER, INTENT(IN) ::            &
         jg                          ,& !< domain index
         kproma                      ,& !< number of elements in current block
         kbdim                       ,& !< block dimension (greater than or equal to kproma)
         klev                        ,& !< number of full levels
         krow                        ,& !< index for current block
         nb_sw                          !< number of bands in short wave

    TYPE(t_datetime), INTENT(IN) :: datetime

    REAL(wp), INTENT (IN)        :: &
         geoi(kbdim,klev+1),        & !< geopotential wrt surface at layer interfaces
         geom(kbdim,klev),          & !< geopotential wrt surface at layer centres
         oromea(kbdim)                !< orography in metres

    REAL(wp), INTENT (INOUT) ::       &
         aod_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol shortwave optical depth
         ssa_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol single scattering albedo
         asy_sw_vr(kbdim,klev,nb_sw) ,& !< Aerosol asymmetry parameter
         x_cdnc(kbdim)                  !< Scale factor for Cloud Droplet Number Concentration
  
    !
    ! --- 0.2 Dummy variables
    !
    INTEGER ::                        &
         jk                          ,& !< index for looping over vertical dimension
         jki                         ,& !< index for looping over vertical dimension for reversing
         jl                          ,& !< index for looping over block
         jwl                         ,& !< index for looping over wavelengths
         j_sw                           !< index for looping over wavelengths
    
    REAL(wp) ::                       &
         year_fr                     ,& !< time in year fraction (1989.0 is 0Z on Jan 1 1989)
         lambda                      ,& !< wavelength at central band wavenumber [nm]
         z_sfc(kproma)               ,& !< surface height [m]
         lon_sp(kproma)              ,& !< longitude passed to sp
         lat_sp(kproma)              ,& !< latitude passed to sp
         z_fl_vr(kbdim,klev)         ,& !< level height [m], vertically reversed indexing (1=lowest level)
         dz_vr(kbdim,klev)           ,& !< level thickness [m], vertically reversed 
         sp_aod_vr(kbdim,klev)       ,& !< simple plume aerosol optical depth, vertically reversed 
         sp_ssa_vr(kbdim,klev)       ,& !< simple plume single scattering albedo, vertically reversed
         sp_asy_vr(kbdim,klev)       ,& !< simple plume asymmetry factor, vertically reversed indexing
         sp_xcdnc(kproma)               !< drop number scale factor

    year_fr=datetime%yeafrc+REAL(datetime%year)
    IF (datetime%year > 1850) THEN
      ! 
      ! --- 1.1 geographic information
      !
      DO jk=1,klev
        jki=klev-jk+1
        DO jl=1,kproma
          dz_vr  (jl,jk) = rgrav*(geoi(jl,jki+1)-geoi(jl,jki))
          z_fl_vr(jl,jk) = rgrav*geom(jl,jki)+oromea(jl)
        END DO
      END DO
      z_sfc(1:kproma)  = oromea(1:kproma)
      lon_sp(1:kproma) = p_patch(jg)%cells%center(1:kproma,krow)%lon*rad2deg
      lat_sp(1:kproma) = p_patch(jg)%cells%center(1:kproma,krow)%lat*rad2deg
      ! 
      ! --- 1.2 Aerosol Shortwave properties
      !
      ! get aerosol optical properties in each band, and adjust effective radius
      !
      DO jwl = 1,nb_sw
        j_sw   = jpb1 + jwl - 1
        lambda = 1.e7_wp/ (0.5_wp * (sw_wv1(j_sw) + sw_wv2(j_sw)))
        CALL sp_aop_profile                                                                   ( &
           & klev               ,kproma             ,kbdim               ,lambda              , &
           & z_sfc(:)           ,lon_sp(:)          ,lat_sp(:)           ,year_fr             , &
           & z_fl_vr(:,:)       ,dz_vr(:,:)         ,sp_xcdnc(:)         ,sp_aod_vr(:,:)      , &
           & sp_ssa_vr(:,:)     ,sp_asy_vr(:,:)                                               )

        sp_xcdnc(:)=0.01_wp
        sp_aod_vr(:,:)=0.01_wp
        sp_ssa_vr(:,:)=0.9_wp
        sp_asy_vr(:,:)=0.9_wp

        DO jk=1,klev
          DO jl=1,kproma
            asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl) * ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)    &
                 + sp_asy_vr(jl,jk)   * sp_ssa_vr(jl,jk)    * sp_aod_vr(jl,jk)
            ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)                           &
                 + sp_ssa_vr(jl,jk)   * sp_aod_vr(jl,jk)
            aod_sw_vr(jl,jk,jwl) = aod_sw_vr(jl,jk,jwl) + sp_aod_vr(jl,jk)
            asy_sw_vr(jl,jk,jwl) = MERGE(asy_sw_vr(jl,jk,jwl)/ssa_sw_vr(jl,jk,jwl),asy_sw_vr(jl,jk,jwl), &
                 ssa_sw_vr(jl,jk,jwl) > TINY(1.0_wp))
            ssa_sw_vr(jl,jk,jwl) = MERGE(ssa_sw_vr(jl,jk,jwl)/aod_sw_vr(jl,jk,jwl),ssa_sw_vr(jl,jk,jwl), &
                 aod_sw_vr(jl,jk,jwl) > TINY(1.0_wp))
          END DO
        END DO
      END DO

      DO jl=1,kproma
        x_cdnc(jl) = sp_xcdnc(jl)
        !*sp_xre(jl)
      END DO
      RETURN
    END IF
 
  END SUBROUTINE add_bc_aeropt_splumes

  SUBROUTINE read_1d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(1)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length, cj_length
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_1D(file_id=ifile_id,         variable_name=variable_name,      &
                 return_pointer=return_pointer                               )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1)) THEN
         WRITE(ci_length,*) SIZE(return_pointer,1)
         WRITE(cj_length,*) variable_dimls(1)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('// &
                & TRIM(ADJUSTL(cj_length))//') has wrong dimension length '// &
                & TRIM(ADJUSTL(ci_length))//' in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_1d_wrapper
  SUBROUTINE read_2d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:,:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(2)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(2), cj_length(2)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_2D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  return_pointer=return_pointer                    )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1) .OR. &
         & SIZE(return_pointer,2)/=variable_dimls(2)) THEN
         WRITE(ci_length(1),*) SIZE(return_pointer,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(return_pointer,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_2d_wrapper
  SUBROUTINE read_3d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:,:,:) !< values of 
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(3)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(3), cj_length(3)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_3D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  return_pointer=return_pointer                    )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1) .OR. &
         & SIZE(return_pointer,2)/=variable_dimls(2) .OR. &
         & SIZE(return_pointer,3)/=variable_dimls(3)) THEN
         WRITE(ci_length(1),*) SIZE(return_pointer,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(return_pointer,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         WRITE(ci_length(3),*) SIZE(return_pointer,3)
         WRITE(cj_length(3),*) variable_dimls(3)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//','//&
                & TRIM(ADJUSTL(cj_length(3)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//','//&
                & TRIM(ADJUSTL(ci_length(3)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
!!$         WRITE(0,*) TRIM(ADJUSTL(message1))
!!$         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_3d_wrapper
  
END MODULE mo_bc_aeropt_splumes
