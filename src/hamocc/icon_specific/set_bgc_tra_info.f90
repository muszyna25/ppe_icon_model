      SUBROUTINE set_bgc_tracer_info(no_tracer,max_oce_tracer,bgc_tracer_names,&
      & bgc_tracer_longnames,&
      & bgc_tracer_codes,&
      & bgc_tracer_units, suffix)
      USE mo_impl_constants,         ONLY: max_char_length

      INTEGER, INTENT(in)            :: max_oce_tracer, no_tracer
      CHARACTER(LEN=max_char_length) :: bgc_tracer_names(max_oce_tracer),&
      & bgc_tracer_units(max_oce_tracer),&
      & bgc_tracer_longnames(max_oce_tracer)
      INTEGER :: bgc_tracer_codes(max_oce_tracer),itrac, last_oce_code
      CHARACTER(LEN=max_char_length), OPTIONAL :: suffix
    
       last_oce_code=bgc_tracer_codes(no_tracer)

       ! Note: tracers are ordered according to HAMOCC ids in mo_parameter1_bgc

       itrac = no_tracer + 1
       bgc_tracer_names(itrac)     = 'dic'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'dic'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'dissolved inorganic carbon'
       bgc_tracer_units(itrac)     = 'kmol C m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac
    
       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'alk'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'alk_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'alkalinity'
       bgc_tracer_units(itrac)     = 'kmol m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac
    
       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'phosph'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'phosph_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'phosphate'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'nitrate'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'nitrate_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'nitrate'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'gasnit'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'gasnit_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'gaseous nitrogen'
       bgc_tracer_units(itrac)     = 'kmol N m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'phy'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'phy_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'phytoplankton'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'zoo'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'zoo_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'zooplankton'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'cyano'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'cyano_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'cyanobacteria'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'oxygen'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'oxygen_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'oxygen'
       bgc_tracer_units(itrac)     = 'kmol m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'silica'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'silica_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'silicate'
       bgc_tracer_units(itrac)     = 'kmol m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'doc'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'doc_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'dissolved_organic_carbon'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'an2o'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'an2o_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'an2o'
       bgc_tracer_units(itrac)     = 'kmol  m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'det'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'det_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'particulate_organic_carbon'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac
      
       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'doccya'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'doccya_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'dissolved_organic_carbon of cyanos'
       bgc_tracer_units(itrac)     = 'kmol P m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'iron'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'iron_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'iron'
       bgc_tracer_units(itrac)     = 'kmol Fe m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'dms'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'dms_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'dms'
       bgc_tracer_units(itrac)     = 'kmol DMS m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'calc'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'calc_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'calcium_carbonate'
       bgc_tracer_units(itrac)     = 'kmol C m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'opal'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'opal_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'opal shells'
       bgc_tracer_units(itrac)     = 'kmol m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac


       itrac=itrac+1
       bgc_tracer_names(itrac)     = 'dust'
       IF (PRESENT(suffix)) THEN
        bgc_tracer_names(itrac) = 'dust_'//TRIM(suffix)
       END IF
       bgc_tracer_longnames(itrac) = 'free dust in seawater'
       bgc_tracer_units(itrac)     = 'kmol m-3'
       bgc_tracer_codes(itrac)     = last_oce_code+itrac

       END SUBROUTINE
