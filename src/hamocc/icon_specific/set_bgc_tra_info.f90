      SUBROUTINE set_bgc_tracer_info(max_oce_tracer,bgc_tracer_names,&
      & bgc_tracer_longnames,&
      & bgc_tracer_codes,&
      & bgc_tracer_units, suffix)
      USE mo_impl_constants,         ONLY: max_char_length

      INTEGER, INTENT(in)            :: max_oce_tracer
      CHARACTER(LEN=max_char_length) :: bgc_tracer_names(max_oce_tracer),&
      & bgc_tracer_units(max_oce_tracer),&
      & bgc_tracer_longnames(max_oce_tracer)
      INTEGER :: bgc_tracer_codes(max_oce_tracer)
      CHARACTER(LEN=max_char_length), OPTIONAL :: suffix
    
       bgc_tracer_names(3)     = 'dic'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(3) = 'dic'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(3) = 'dissolved inorganic carbon'
       bgc_tracer_units(3)     = 'kmol C m-3'
       bgc_tracer_codes(3)     = 203
    
       bgc_tracer_names(4)     = 'alk'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(4) = 'alk_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(4) = 'alkalinity'
       bgc_tracer_units(4)     = 'kmol m-3'
       bgc_tracer_codes(4)     = 204
    
       bgc_tracer_names(5)     = 'phosph'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(5) = 'phosph_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(5) = 'phosphate'
       bgc_tracer_units(5)     = 'kmol P m-3'
       bgc_tracer_codes(5)     = 205

       bgc_tracer_names(6)     = 'nitrate'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(6) = 'nitrate_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(6) = 'nitrate'
       bgc_tracer_units(6)     = 'kmol P m-3'
       bgc_tracer_codes(6)     = 206

       bgc_tracer_names(7)     = 'gasnit'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(7) = 'gasnit_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(7) = 'gaseous nitrogen'
       bgc_tracer_units(7)     = 'kmol N m-3'
       bgc_tracer_codes(7)     = 207

       bgc_tracer_names(8)     = 'phy'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(8) = 'phy_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(8) = 'phytoplankton'
       bgc_tracer_units(8)     = 'kmol P m-3'
       bgc_tracer_codes(8)     = 208

       bgc_tracer_names(9)     = 'zoo'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(9) = 'zoo_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(9) = 'zooplankton'
       bgc_tracer_units(9)     = 'kmol P m-3'
       bgc_tracer_codes(9)     = 209

       bgc_tracer_names(10)     = 'cyano'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(10) = 'cyano_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(10) = 'cyanobacteria'
       bgc_tracer_units(10)     = 'kmol P m-3'
       bgc_tracer_codes(10)     = 210

       bgc_tracer_names(11)     = 'oxygen'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(11) = 'oxygen_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(11) = 'oxygen'
       bgc_tracer_units(11)     = 'kmol m-3'
       bgc_tracer_codes(11)     = 211

       bgc_tracer_names(12)     = 'silica'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(12) = 'silica_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(12) = 'silicate'
       bgc_tracer_units(12)     = 'kmol m-3'
       bgc_tracer_codes(12)     = 212

       bgc_tracer_names(13)     = 'doc'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(13) = 'doc_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(13) = 'dissolved_organic_carbon'
       bgc_tracer_units(13)     = 'kmol P m-3'
       bgc_tracer_codes(13)     = 213

       bgc_tracer_names(14)     = 'an2o'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(14) = 'an2o_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(14) = 'an2o'
       bgc_tracer_units(14)     = 'kmol  m-3'
       bgc_tracer_codes(14)     = 214

       bgc_tracer_names(15)     = 'det'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(15) = 'det_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(15) = 'particulate_organic_carbon'
       bgc_tracer_units(15)     = 'kmol P m-3'
       bgc_tracer_codes(15)     = 215
      
       bgc_tracer_names(16)     = 'doccya'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(16) = 'doccya_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(16) = 'dissolved_organic_carbon of cyanos'
       bgc_tracer_units(16)     = 'kmol P m-3'
       bgc_tracer_codes(16)     = 216

       bgc_tracer_names(17)     = 'iron'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(17) = 'iron_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(17) = 'iron'
       bgc_tracer_units(17)     = 'kmol Fe m-3'
       bgc_tracer_codes(17)     = 217

       bgc_tracer_names(18)     = 'dms'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(18) = 'dms_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(18) = 'dms'
       bgc_tracer_units(18)     = 'kmol DMS m-3'
       bgc_tracer_codes(18)     = 218

       bgc_tracer_names(19)     = 'calc'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(19) = 'calc_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(19) = 'calcium_carbonate'
       bgc_tracer_units(19)     = 'kmol C m-3'
       bgc_tracer_codes(19)     = 219

       bgc_tracer_names(20)     = 'opal'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(20) = 'opal_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(20) = 'opal shells'
       bgc_tracer_units(20)     = 'kmol m-3'
       bgc_tracer_codes(20)     = 220


       bgc_tracer_names(21)     = 'dust'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(21) = 'dust_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(21) = 'free dust in seawater'
       bgc_tracer_units(21)     = 'kmol m-3'
       bgc_tracer_codes(21)     = 221
       END SUBROUTINE
