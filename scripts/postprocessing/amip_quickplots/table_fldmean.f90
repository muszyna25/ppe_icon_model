! ECHAM post processing
! Generates tables of averages
!
MODULE mo_tables

  IMPLICIT NONE

  INTEGER, PARAMETER :: kmax = 276

  CHARACTER(len=60) :: string(kmax)

  REAL :: offset(kmax), factor(kmax)

CONTAINS

  SUBROUTINE init_tables 

write(6,*) "kmax",kmax

    string( 34) = 'LOW_CLD       %      low cloud                              ';  factor( 34) =     100.00; offset( 34) =    0.00;
    string( 35) = 'MID_CLD       %      mid cloud                              ';  factor( 35) =     100.00; offset( 35) =    0.00;
    string( 36) = 'HIH_CLD       %      high cloud                             ';  factor( 36) =     100.00; offset( 36) =    0.00;
    string( 61) = 'prlr          mm/d   precip large scale rain                ';  factor( 61) =   86400.00; offset( 61) =    0.00;
    string( 62) = 'prls          mm/d   precip large scale snow                ';  factor( 62) =   86400.00; offset( 62) =    0.00;
    string( 63) = 'prcr          mm/d   precip convective rain                 ';  factor( 63) =   86400.00; offset( 63) =    0.00;
    string( 64) = 'prcs          mm/d   precip convective snow                 ';  factor( 64) =   86400.00; offset( 64) =    0.00;
    string( 68) = 'FAGE          %      aging factor of snow on ice          I ';  factor( 68) =     100.00; offset( 68) =    0.00;
    string( 69) = 'SNIFRAC       %      fraction of ice covered with snow    I ';  factor( 69) =     100.00; offset( 69) =    0.00;
    string( 70) = 'BAREFRAC      %      bare ice fraction                    I ';  factor( 70) =     100.00; offset( 70) =    0.00;
    string( 71) = 'ALSOM         %      albedo of melt ponds                 I ';  factor( 71) =     100.00; offset( 71) =    0.00;
    string( 72) = 'ALSOBS        %      albedo bare ice and snow             I ';  factor( 72) =     100.00; offset( 72) =    0.00;
    string( 73) = 'SICEPDW       cm     melt pond depth on sea-ice           I ';  factor( 73) =     100.00; offset( 73) =    0.00;
    string( 74) = 'SICEPDI       cm     ice thickness on melt pond           I ';  factor( 74) =     100.00; offset( 74) =    0.00;
    string( 75) = 'TSICEPDI      C      ice temperature on frozen melt pond  I ';  factor( 75) =       1.00; offset( 75) = -273.15;
    string( 76) = 'SICEPRES      W/m2   residual heat flux                   I ';  factor( 76) =       1.00; offset( 76) =    0.00;
    string( 77) = 'AMELTDEPTH    cm     total melt pond depth                I ';  factor( 77) =     100.00; offset( 77) =    0.00;
    string( 78) = 'AMELTFRAC     %      frac. area of melt ponds on sea-ice  I ';  factor( 78) =     100.00; offset( 78) =    0.00;
    
    string( 79) = 'ALBEDO_VIS_DIR %     surface albedo visible range direct    ';  factor( 79) =     100.00; offset( 79) =    0.00;
    string( 80) = 'ALBEDO_NIR_DIR %     surface albedo NIR range direct        ';  factor( 80) =     100.00; offset( 80) =    0.00;
    string( 81) = 'ALBEDO_VIS_DIF %     surface albedo visible range diffuse   ';  factor( 81) =     100.00; offset( 81) =    0.00;
    string( 82) = 'ALBEDO_NIR_DIF %     surface albedo NIR range diffuse       ';  factor( 82) =     100.00; offset( 82) =    0.00;
    string( 83) = 'OCU           m/s    ocean eastw. velocity  -coupled mode W ';  factor( 83) =       1.00; offset( 83) =    0.00;
    string( 84) = 'OCV           m/s    ocean northw. velocity -coupled mode W ';  factor( 84) =       1.00; offset( 84) =    0.00;
    string( 85) = 'TRADL         W/m2   net LW radiation 200mb                 ';  factor( 85) =       1.00; offset( 85) =    0.00;
    string( 86) = 'SRADL         W/m2   net SW radiation 200mb                 ';  factor( 86) =       1.00; offset( 86) =    0.00;
    string( 87) = 'TRAFL         W/m2   net LW radiation 200mb (clear sky)     ';  factor( 87) =       1.00; offset( 87) =    0.00;
    string( 88) = 'SRAFL         W/m2   net LW radiation 200mb (clear sky)     ';  factor( 88) =       1.00; offset( 88) =    0.00;
    string( 89) = 'AMLCORAC      W/m2   mixed layer correction               WA';  factor( 89) =       1.00; offset( 89) =    0.00;
    string( 90) = 'AMLHEATAC     J/m2   mixed layer heat content             WA';  factor( 90) = 1.7889E-11; offset( 90) =    0.00;
    string( 91) = 'rlns_ice      W/m2   net LW radiation over ice              ';  factor( 91) =       1.00; offset( 91) =    0.00;
    string( 92) = 'rlns_wrt      W/m2   net LW radiation over water            ';  factor( 92) =       1.00; offset( 92) =    0.00;
    string( 93) = 'rlns_lnd      W/m2   net LW radiation over land             ';  factor( 93) =       1.00; offset( 93) =    0.00;
    string( 94) = 'rsns_ice      W/m2   net SW radiation over ice              ';  factor( 94) =       1.00; offset( 94) =    0.00;
    string( 95) = 'rsns_wrt      W/m2   net SW radiation over water            ';  factor( 95) =       1.00; offset( 95) =    0.00;
    string( 96) = 'rsns_lnd      W/m2   net SW radiation over land             ';  factor( 96) =       1.00; offset( 96) =    0.00;
    string( 97) = 'sic           %      ice cover (fraction of grid box)       ';  factor( 97) =     100.00; offset( 97) =    0.00;
    string( 98) = '              -----                              -----      ';  factor( 98) =       1.00; offset( 98) =    0.00;
    string( 99) = '              -----                              -----      ';  factor( 99) =       1.00; offset( 99) =    0.00;
    string(100) = 'ALBEDO_VIS    %      surface albedo visible range           ';  factor(100) =       1.00; offset(100) =    0.00;
    string(101) = 'ALBEDO_NIR    %      surface albedo NIR range               ';  factor(101) =       1.00; offset(101) =    0.00;
    string(102) = 'ts_ice        C      surface temperature of ice             ';  factor(102) =       1.00; offset(102) = -273.15;
    string(103) = 'ts_wtr        C      surface temperature of water           ';  factor(103) =       1.00; offset(103) = -273.15;
    string(104) = 'tauu_ice      mN/m2  zonal wind stress over ice             ';  factor(104) =    1000.00; offset(104) =    0.00;
    string(105) = 'tauv_ice      mN/m2  meridional wind stress over ice        ';  factor(105) =    1000.00; offset(105) =    0.00;
    string(106) = 'tauu_wrt      mN/m2  zonal wind stress over water           ';  factor(106) =    1000.00; offset(106) =    0.00;
    string(107) = 'tauv_wrt      mN/m2  meridional wind stress over water      ';  factor(107) =    1000.00; offset(107) =    0.00;
    string(108) = 'tauu_lnd      mN/m2  zonal wind stress over land            ';  factor(108) =    1000.00; offset(108) =    0.00;
    string(109) = 'tauv_lnd      mN/m2  meridional wind stress over land       ';  factor(109) =    1000.00; offset(109) =    0.00;
    string(110) = 'hfls_ice      W/m2   latent heat flux over ice              ';  factor(110) =       1.00; offset(110) =    0.00;
    string(111) = 'hfls_wrt      W/m2   latent heat flux over water            ';  factor(111) =       1.00; offset(111) =    0.00;
    string(112) = 'hfls_ice      W/m2   latent heat flux over land             ';  factor(112) =       1.00; offset(112) =    0.00;
    string(113) = 'EVAPIAC       mm/d   evaporation over ice                 IA';  factor(113) =   86400.00; offset(113) =    0.00;
    string(114) = 'EVAPWAC       mm/d   evaporation over water               WA';  factor(114) =   86400.00; offset(114) =    0.00;
    string(115) = 'EVAPLAC       mm/d   evaporation over land                LA';  factor(115) =   86400.00; offset(115) =    0.00;
    string(116) = 'AZ0I          cm     roughness length over ice            I ';  factor(116) =     100.00; offset(116) =    0.00;
    string(117) = 'AZ0W          cm     roughness length over water          W ';  factor(117) =     100.00; offset(117) =    0.00;
    string(118) = 'AZ0L          cm     roughness length over land           L ';  factor(118) =     100.00; offset(118) =    0.00;
    string(119) = 'hfss_ice      W/m2   sensible heat flux over ice            ';  factor(119) =       1.00; offset(119) =    0.00;
    string(120) = 'hfss_wrt      W/m2   sensible heat flux over water          ';  factor(120) =       1.00; offset(120) =    0.00;
    string(121) = 'hfss_lnd      W/m2   sensible heat flux over land           ';  factor(121) =       1.00; offset(121) =    0.00;
    string(122) = 'ALSOI         %      albedo of ice                        I ';  factor(122) =     100.00; offset(122) =    0.00;
    string(123) = 'ALSOW         %      albedo of water                      W ';  factor(123) =     100.00; offset(123) =    0.00;
    string(124) = 'ALSOL         %      albedo of land                       L ';  factor(124) =     100.00; offset(124) =    0.00;
    string(125) = 'AHFICE        W/m2   conductive heat flux through ice     I ';  factor(125) =       1.00; offset(125) =    0.00;
    string(126) = 'QRES          W/m2   residual heat flux for melt. sea ice I ';  factor(126) =       1.00; offset(126) =    0.00;
    string(127) = 'ALAKE         %      lake fraction                          ';  factor(127) =     100.00; offset(127) =    0.00;
    string(128) = 'RINTOP        %      frequency of low level inversion       ';  factor(128) =       1.00; offset(128) =    0.00;
    string(129) = 'GEOSP (Z0)    m      surface geopotential (orography)       ';  factor(129) =    .101937; offset(129) =    0.00;
    string(130) = 'STP           C      temperature                            ';  factor(130) =       1.00; offset(130) = -273.15;
    string(131) = 'U             m/s    u-velocity                             ';  factor(131) =       1.00; offset(131) =    0.00;
    string(132) = 'V             m/s    v-velocity                             ';  factor(132) =       1.00; offset(132) =    0.00;
    string(133) = 'Q             kg/kg  specific humidity                      ';  factor(133) =       1.00; offset(133) =    0.00;
    string(134) = 'APS           hPa    surface pressure                       ';  factor(134) =       0.01; offset(134) =    0.00;
    string(135) = 'OMEGA         Pa/s   vertical velocity                      ';  factor(135) =       1.00; offset(135) =    0.00;
    string(136) = 'ACDNC         1/m3   cloud droplet number concentration     ';  factor(136) =       1.00; offset(136) =    0.00;
    string(137) = 'APMEB         mm/d   (P-E) error                            ';  factor(137) =   86400.00; offset(137) =    0.00;
    string(138) = 'SVO           1/s    vorticity                              ';  factor(138) =       1.00; offset(138) =    0.00;
    string(139) = 'ts_lnd        C      surface temperature of land            ';  factor(139) =       1.00; offset(139) = -273.15;
    string(140) = 'WS            cm     soil wetness                         L ';  factor(140) =     100.00; offset(140) =    0.00;
    string(141) = 'SN            cm     water equivalent of snow depth       L ';  factor(141) =     100.00; offset(141) =    0.00;
    string(142) = 'APRL          mm/d   large scale precipitation              ';  factor(142) =   86400.00; offset(142) =    0.00;
    string(143) = 'APRC          mm/d   convective precipitation               ';  factor(143) =   86400.00; offset(143) =    0.00;
    string(144) = 'APRS          mm/d   snow fall                              ';  factor(144) =   86400.00; offset(144) =    0.00;
    string(145) = 'VDIS          W/m2   boundary layer dissipation             ';  factor(145) =       1.00; offset(145) =    0.00;
    string(146) = 'hfss          W/m2   sensible heat flux                     ';  factor(146) =       1.00; offset(146) =    0.00;
    string(147) = 'hfls          W/m2   latent heat flux                       ';  factor(147) =       1.00; offset(147) =    0.00;
    string(148) = 'STREAM        m2/s   streamfunction                         ';  factor(148) =       1.00; offset(148) =    0.00;
    string(149) = 'VELOPOT       m2/s   velocity potential                     ';  factor(149) =       1.00; offset(149) =    0.00;
    string(150) = 'clivi         g/m2   vertically integrated cloud ice        ';  factor(150) =    1000.00; offset(150) =    0.00;
    string(151) = 'SLP           hPa    mean sea level pressure -1000 hPa      ';  factor(151) =       0.01; offset(151) = -1000.0;
    string(152) = 'STP                  log surface pressure                   ';  factor(152) =       1.00; offset(152) =    0.00;
    string(153) = 'XL            kg/kg  cloud water                            ';  factor(153) =       1.00; offset(153) =    0.00;
    string(154) = 'XI            kg/kg  cloud ice                              ';  factor(154) =       1.00; offset(154) =    0.00;
    string(155) = 'SD             1/s   divergence                             ';  factor(155) =       1.00; offset(155) =    0.00;
    string(156) = 'GEOPOTH       m      geopotential height                    ';  factor(156) =       1.00; offset(156) =    0.00;
    string(157) = 'RHUMIDITY     %      relative humidity                      ';  factor(157) =     100.00; offset(157) =    0.00;
    string(158) = '              -----                                 -----   ';  factor(158) =       1.00; offset(158) =    0.00;
    string(159) = 'WIND10W       m/s    windspeed over water                   ';  factor(159) =       1.00; offset(159) =    0.00;
    string(160) = 'RUNOFF        mm/d   surface runoff and drainage          LA';  factor(160) =   86400.00; offset(160) =    0.00;
    string(161) = 'DRAIN         mm/d   drainage                             LA';  factor(161) =   86400.00; offset(161) =    0.00;
    string(163) = 'ACLCV         %      total cloud cover                      ';  factor(163) =     100.00; offset(163) =    0.00;
    string(164) = 'clt           %      total cloud cover                      ';  factor(164) =     100.00; offset(164) =    0.00;
    string(165) = 'U10           m/s    10 m u-velocity                        ';  factor(165) =       1.00; offset(165) =    0.00;
    string(166) = 'V10           m/s    10 m v-velocity                        ';  factor(166) =       1.00; offset(166) =    0.00;
    string(167) = 'TEMP2         C      2 m temperature                        ';  factor(167) =       1.00; offset(167) = -273.15;
    string(168) = 'DEW2          C      2 m dew-point temperature              ';  factor(168) =       1.00; offset(168) = -273.15;
    string(169) = 'ts            C      surface temperature                    ';  factor(169) =       1.00; offset(169) = -273.15;
    string(170) = 'XVAR          kg/kg  variance of total water amount         ';  factor(170) =       1.00; offset(170) =    0.00;
    string(171) = 'WIND10        m/s    10m windspeed                          ';  factor(171) =       1.00; offset(171) =    0.00;
    string(172) = 'SLM           %      land sea mask                          ';  factor(172) =     100.00; offset(172) =    0.00;
    string(173) = 'AZ0           cm     roughness length                       ';  factor(173) =     100.00; offset(173) =    0.00;
    string(174) = 'ALB           %      surface background albedo              ';  factor(174) =     100.00; offset(174) =    0.00;
    string(175) = 'alb           %      surface albedo                         ';  factor(175) =     100.00; offset(175) =    0.00;
    string(176) = 'rsns          W/m2   net surface SW radiation               ';  factor(176) =       1.00; offset(176) =    0.00;
    string(177) = 'rlns          W/m2   net surface LW radiation               ';  factor(177) =       1.00; offset(177) =    0.00;
    string(178) = 'rsnt          W/m2   net top SW radiation                   ';  factor(178) =       1.00; offset(178) =    0.00;
    string(179) = 'rlnt          W/m2   net top LW radiation  (-OLR)           ';  factor(179) =       1.00; offset(179) =    0.00;
    string(180) = 'tauu          mN/m2  u-stress                               ';  factor(180) =    1000.00; offset(180) =    0.00;
    string(181) = 'tauv          mN/m2  v-stress                               ';  factor(181) =    1000.00; offset(181) =    0.00;
    string(182) = 'evspsbl       mm/d   evaporation                            ';  factor(182) =   86400.00; offset(182) =    0.00;
    string(183) = 'XSKEW                skewness of total water amount         ';  factor(183) =       1.00; offset(183) =    0.00;
    string(184) = 'rsdt          W/m2   top incoming SW radiation              ';  factor(184) =       1.00; offset(184) =    0.00;
    string(185) = 'SRAFS         W/m2   net surface SW radiation (clear sky)   ';  factor(185) =       1.00; offset(185) =    0.00;
    string(186) = 'TRAFS         W/m2   net surface LW radiation (clear sky)   ';  factor(186) =       1.00; offset(186) =    0.00;
    string(187) = 'SRAF0         W/m2   net top SW radiation (clear sky)       ';  factor(187) =       1.00; offset(187) =    0.00;
    string(188) = 'TRAF0         W/m2   net top LW radiation (clear sky)       ';  factor(188) =       1.00; offset(188) =    0.00;
    string(189) = 'SCLFS         W/m2   net surface SW cloud forcing(176-185)  ';  factor(189) =       1.00; offset(189) =    0.00;
    string(190) = 'TCLFS         W/m2   net surface LW cloud forcing(177-186)  ';  factor(190) =       1.00; offset(190) =    0.00;
    string(191) = 'SCLF0         W/m2   net top SW cloud forcing (178-187)     ';  factor(191) =       1.00; offset(191) =    0.00;
    string(192) = 'TCLF0         W/m2   net top LW cloud forcing (179-188)     ';  factor(192) =       1.00; offset(192) =    0.00;
    string(193) = 'WL            mm     skin reservoir                       L ';  factor(193) =    1000.00; offset(193) =    0.00;
    string(194) = 'SLF           %      fractional land cover                  ';  factor(194) =     100.00; offset(194) =    0.00;
    string(195) = 'tauu_sso      mN/m2  u-gravity wave stress                  ';  factor(195) =    1000.00; offset(195) =    0.00;
    string(196) = 'tauv_sso      mN/m2  v-gravity wave stress                  ';  factor(196) =    1000.00; offset(196) =    0.00;
    string(197) = 'diss_sso      W/m2   gravity wave dissipation               ';  factor(197) =       1.00; offset(197) =    0.00;
    string(198) = 'VGRAT                vegetation ratio                     L ';  factor(198) =       1.00; offset(198) =    0.00;
    string(199) = 'OROSTD        m      orographic standard deviation        L ';  factor(199) =       1.00; offset(199) =    0.00;
    string(200) = 'VLT                  leaf area index                      L ';  factor(200) =       1.00; offset(200) =    0.00;
    string(201) = 'T2MAX         C      maximum 2-meter temperature            ';  factor(201) =       1.00; offset(201) = -273.15;
    string(202) = 'T2MIN         C      minimum 2-meter temperature            ';  factor(202) =       1.00; offset(202) = -273.15;
    string(203) = 'SRADOU        W/m2   top SW radiation upward                ';  factor(203) =       1.00; offset(203) =    0.00;
    string(204) = 'SRADSU        W/m2   surface SW radiation upward            ';  factor(204) =       1.00; offset(204) =    0.00;
    string(205) = 'TRADSU        W/m2   surface LW radiation upward            ';  factor(205) =       1.00; offset(205) =    0.00;
    string(206) = 'GRNDFLUX      W/m2   surface ground heat flux             LA';  factor(206) =       1.00; offset(206) =    0.00;
    string(207) = 'TSOIL         C      deep soil temperature (layer 1)      L ';  factor(207) =       1.00; offset(207) = -273.15;
    string(208) = 'AHFCON        W/m2   conductive heat flux through ice     IA';  factor(208) =       1.00; offset(208) =    0.00;
    string(209) = 'AHFRES        W/m2   res. heat flux for melting ice       IA';  factor(209) =       1.00; offset(209) =    0.00;
    string(210) = 'SEAICE        %      ice cover (fraction of ice+water)    S ';  factor(210) =     100.00; offset(210) =    0.00;
    string(211) = 'sit           m      ice thickness                          ';  factor(211) =       1.00; offset(211) =    0.00;
    string(212) = 'FOREST        %      forest fraction                      L ';  factor(212) =     100.00; offset(212) =    0.00;
    string(213) = 'GLD           m      glacier thickness                    S ';  factor(213) =       1.00; offset(213) =    0.00;
    string(214) = 'SNI           cm     water equivalent of snow on ice      S ';  factor(214) =     100.00; offset(214) =    0.00;
    string(215) = 'ROGL          mm/d   glacier runoff                       LA';  factor(215) =   86400.00; offset(215) =    0.00;
    string(216) = 'WIMAX         m/s    maximum 10-m wind                      ';  factor(216) =       1.00; offset(216) =    0.00;
    string(217) = 'TOPMAX        hPa    maximum height convective cloud tops   ';  factor(217) =       0.01; offset(217) =    0.00;
    string(218) = 'SNMEL         mm/d   snow melt                            LA';  factor(218) =   86400.00; offset(218) =    0.00;
    string(219) = 'RUNTOC        mm/d   surface runoff into ocean            SA';  factor(219) =   86400.00; offset(219) =    0.00;
    string(220) = 'RUNLND        W/m2   surface runoff not runn. into ocean    ';  factor(220) =   86400.00; offset(220) =    0.00;
    string(221) = 'APMEGL        mm/d   P-E over land ice                    LA';  factor(221) =   86400.00; offset(221) =    0.00;
    string(222) = 'SNACL         mm/d   snow accumulation over land          LA';  factor(222) =   86400.00; offset(222) =    0.00;
    string(223) = 'ACLCAC        %      cloud cover                            ';  factor(223) =     100.00; offset(223) =    0.00;
    string(224) = 'TKE           m2/s2  turbulent kinetic energy               ';  factor(224) =       1.00; offset(224) =    0.00;
    string(225) = 'TKEM1         m2/s2  turbulent kinetic energy (t-1)         ';  factor(225) =       1.00; offset(225) =    0.00;
    string(226) = 'FAO           0...5  "FAO" data set (soil data flags)     L ';  factor(226) =       1.00; offset(226) =    0.00;
    string(227) = 'RGCGN                heat capacity of soil/100000.        L ';  factor(227) =     1.E-06; offset(227) =    0.00;
    string(228) = 'SODIF                soil diffusivity*1000000.            L ';  factor(228) = 1000000.00; offset(228) =    0.00;
    string(229) = 'WSMX          cm     field capacity of soil               L ';  factor(229) =     100.00; offset(229) =    0.00;
    string(230) = 'prw           kg/m2  vertically integrated water vapor      ';  factor(230) =       1.00; offset(230) =    0.00;
    string(231) = 'cllvi         g/m2   vertically integrated cloud water      ';  factor(231) =    1000.00; offset(231) =    0.00;
    string(232) = 'GLAC          %      fraction of land covered by glaciers L ';  factor(232) =     100.00; offset(232) =    0.00;
    string(233) = 'SNC           mm     snow depth at the canopy             L ';  factor(233) =    1000.00; offset(233) =    0.00;
    string(234) = 'RTYPE         0...3  type of covection                      ';  factor(234) =       1.00; offset(234) =    0.00;
    string(235) = 'ABSO4         mg/m2  antropogenic sulfur burden             ';  factor(235) = 1000000.00; offset(235) =    0.00;
    string(236) = 'AO3           kg/m2  ipcc ozone                             ';  factor(236) =       1.00; offset(236) =    0.00;
    string(237) = 'TROPO         Pa     WMO defined tropopause height          ';  factor(237) =       1.00; offset(237) =    0.00;
    string(238) = 'THVSIG        C      stddev virt. pot temp                  ';  factor(238) =       1.00; offset(238) = -273.15;
    string(239) = '              -----                                 -----   ';  factor(239) =       1.00; offset(239) =    0.00;
    string(240) = '              -----                                 -----   ';  factor(240) =       1.00; offset(240) =    0.00;
    string(241) = 'XT(1)                tracer gas                     -----   ';  factor(241) =       1.00; offset(241) =    0.00;
    string(242) = 'XT(2)                tracer gas                     -----   ';  factor(242) =       1.00; offset(242) =    0.00;
    string(243) = '              -----                                 -----   ';  factor(243) =       1.00; offset(243) =    0.00;
    string(244) = '              -----                                 -----   ';  factor(244) =       1.00; offset(244) =    0.00;
    string(245) = '              -----                                 -----   ';  factor(245) =       1.00; offset(245) =    0.00;
    string(246) = '              -----                                 -----   ';  factor(246) =       1.00; offset(246) =    0.00;
    string(247) = '              -----                                 -----   ';  factor(247) =       1.00; offset(247) =    0.00;
    string(248) = '              -----                                 -----   ';  factor(248) =       1.00; offset(248) =    0.00;
    string(249) = '              -----                                 -----   ';  factor(249) =       1.00; offset(249) =    0.00;
    string(250) = '              -----                                 -----   ';  factor(250) =       1.00; offset(250) =    0.00;
    string(251) = '              -----                                 -----   ';  factor(251) =       1.00; offset(251) =    0.00;
    string(252) = '              -----                                 -----   ';  factor(252) =       1.00; offset(252) =    0.00;
    string(253) = '              -----                                 -----   ';  factor(253) =       1.00; offset(253) =    0.00;
    string(254) = '              -----                                 -----   ';  factor(254) =       1.00; offset(254) =    0.00;
    string(255) = '              -----                                 -----   ';  factor(255) =       1.00; offset(255) =    0.00;
    string(256) = '              -----                                 -----   ';  factor(256) =       1.00; offset(256) =    0.00;
    string(257) = '              -----                                 -----   ';  factor(257) =       1.00; offset(257) =    0.00;
    string(258) = '              -----                                 -----   ';  factor(258) =       1.00; offset(258) =    0.00;
    string(259) = '              -----                                 -----   ';  factor(259) =       1.00; offset(259) =    0.00;
    string(260) = 'PRECIP        mm/d   total precipitation (142+143)          ';  factor(260) =   86400.00; offset(260) =    0.00;
    string(261) = 'net_top       W/m2   total top radiation (178+179)          ';  factor(261) =       1.00; offset(261) =    0.00;
    string(262) = '              -----                                 -----   ';  factor(262) =       1.00; offset(262) =    0.00;
    string(263) = '              W/m2   total surface heat flux                ';  factor(263) =       1.00; offset(263) =    0.00;
    string(264) = '              -----                                 -----   ';  factor(264) =       1.00; offset(264) =    0.00;
    string(265) = '              %      WS/WSMX                              L ';  factor(265) =     100.00; offset(265) =    0.00;
    string(266) = '              mm/d   total freshwater flux (260+182)        ';  factor(266) =   86400.00; offset(266) =    0.00;
    string(267) = '              %      water (fraction of grid box)           ';  factor(267) =     100.00; offset(267) =    0.00;
    string(268) = '              mm/d   freshwater flux - ICE                IA';  factor(268) =   86400.00; offset(268) =    0.00;
    string(269) = '              mm/d   freshwater flux - WATER              WA';  factor(269) =   86400.00; offset(269) =    0.00;
    string(270) = '              mm/d   water budget - LAND                  LA';  factor(270) =   86400.00; offset(270) =    0.00;
    string(271) = '              W/m2   heat budget         - LAND           LA';  factor(271) =       1.00; offset(271) =    0.00;
    string(272) = '              W/m2   heat budget         - WATER          WA';  factor(272) =       1.00; offset(272) =    0.00;
    string(273) = '              W/m2   heat budget         - ICE            IA';  factor(273) =       1.00; offset(273) =    0.00;
    string(274) = '              mm/d   total precipitation - ICE            I ';  factor(274) =   86400.00; offset(274) =    0.00;
    string(275) = '              mm/d   total precipitation - WATER          W ';  factor(275) =   86400.00; offset(275) =    0.00;
    string(276) = '              mm/d   total precipitation - LAND           L ';  factor(276) =   86400.00; offset(276) =    0.00;

  END SUBROUTINE init_tables

END MODULE mo_tables

MODULE mo_util_string

CONTAINS

  FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(len=*)              ,INTENT(in) :: upper
    CHARACTER(len=LEN_TRIM(upper))            :: tolower
    
    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')
    
    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
           ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO
    
  END FUNCTION tolower

END MODULE mo_util_string

PROGRAM momitt

  USE mo_tables
  USE mo_util_string

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  REAL(dp) , PARAMETER :: m_pi = 3.14159265358979323846D0
  REAL(dp), ALLOCATABLE :: phi(:), pw(:)

  INTEGER, PARAMETER :: npp = 15



  REAL, ALLOCATABLE :: field(:,:), dat(:,:,:), zm(:,:), xzm(:)

  REAL :: nh(kmax), sh(kmax), gm(kmax), faktor(kmax), addi(kmax)

  REAL :: xnh, xsh, xgm, sumpnh, sumpsh, sumpgm 

  INTEGER :: year1, year2, iLONG, landcode
  INTEGER :: ihead(8), icode(kmax), imean

  INTEGER :: nlon, nlat, nlath

  INTEGER :: i, j, k, kk, n, ic, ida, kend, ja, jh, jo, ilayer207


  
  LOGICAL :: linit = .TRUE.

  CHARACTER(len=82) :: titel

  CHARACTER(len=80) :: average

  CHARACTER(len=49) :: trange
  CHARACTER(len=60) :: stringo(kmax)

  CHARACTER(len=40) :: half
  CHARACTER(len=20) :: expnam

  CHARACTER(len= 3) :: cmean 

  CHARACTER(len=14) :: month(17) = (/ &
       '  January     ','  February    ','  March       ', &
       '  April       ','  May         ','  June        ', &
       '  July        ','  August      ','  September   ', &
       '  October     ','  November    ','  December    ', &
       's DEC/JAN/FEB ','s MAR/APR/MAY ','s JUN/JUL/AUG ', &
       's SEP/OCT/NOV ','s JAN...DEC   ' /)

  NAMELIST / exper / expnam, average, year1, year2, iLONG, landcode

  !     -----------------------------------------------------
  !     Read the experimentname from namelist experim


  write (6,*) 'namelist:'
  READ (*,exper)
  write (6,*) 'namelist',expnam, average, year1, year2, iLONG, landcode


  ! Initilaize tables

  CALL init_tables 

  cmean(1:3) = tolower(average)

average_periode: SELECT CASE (cmean)
  CASE ('jan') ! Januray
    imean =  1
  CASE ('feb') ! February
    imean =  2
  CASE ('mar') ! March
    imean =  3
  CASE ('apr') ! April
    imean =  4
  CASE ('may') ! May
    imean =  5
  CASE ('jun') ! June
    imean =  6
  CASE ('jul') ! July
    imean =  7
  CASE ('aug') ! August
    imean =  8
  CASE ('sep') ! September
    imean =  9
  CASE ('oct') ! October
    imean = 10
  CASE ('nov') ! November
    imean = 11
  CASE ('dec') ! December
    imean = 12
  CASE ('djf') ! dec/jan/feb
    imean = 13
  CASE ('mam') ! mar/apr/may
    imean = 14
  CASE ('jja') ! jun/jul/au
    imean = 15
  CASE ('son') ! sep/oct/nov
    imean = 16
  CASE ('ann') ! a whole year 
    imean = 17
  CASE DEFAULT
    WRITE (0,*) 'Unsupported average peride selected ...'
    STOP 'abort ...'
  END SELECT average_periode


  
  CALL NametoCode (kmax,string)
   

END PROGRAM momitt


SUBROUTINE NametoCode(kmax,string)

  ! 

  IMPLICIT NONE

  INTEGER, INTENT(in) :: kmax
  INTEGER :: echamCode(kmax),i
  REAL    :: globalVar(kmax)
  CHARACTER(len=10) :: iconVar(kmax)
  CHARACTER(len=60), INTENT(in) :: string(kmax)

  OPEN (unit=20,file='global.txt',form='FORMATTED')
  OPEN (unit=21,file='codetext.txt',form='FORMATTED')
    
  WRITE(21,*)' code   name      global    unit   description'
  DO i = 1, kmax
    READ(20,'(f9.4,1x,i4,1x,a10)',END=100)globalVar(i),echamCode(i),iconVar(i)
    write(6,*) "i: ",i,globalVar(i),echamCode(i),iconVar(i)
    WRITE(21,'(1x,1i5,2x,a10,1x,f9.4,1x,1a44)')echamCode(i),iconVar(i),globalVar(i),string(echamCode(i))(15:58)
  ENDDO
!    WRITE(9,'(1x,1i5,3x,1a46)') icode(kk), string(kk)(15:60)
  
100 Continue  

    
END SUBROUTINE
