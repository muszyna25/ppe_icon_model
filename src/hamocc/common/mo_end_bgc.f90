MODULE mo_end_bgc

  USE mo_memory_bgc
  USE mo_sedmnt
  USE mo_hamocc_nml, ONLY: l_cpl_co2 
  USE mo_bgc_memory_types, ONLY : deallocate_all_bgc_memory_types

  IMPLICIT NONE

  PUBLIC


CONTAINS

  SUBROUTINE CLEANUP_HAMOCC

    CALL deallocate_all_bgc_memory_types()
  
    DEALLOCATE (dzs)
    DEALLOCATE (seddzi)
    DEALLOCATE (z_sed)
    DEALLOCATE (seddw)
    DEALLOCATE (porsol)
    DEALLOCATE (porwah)
    DEALLOCATE (porwat)
    DEALLOCATE (pors2w)
    
  END SUBROUTINE CLEANUP_HAMOCC

END MODULE mo_end_bgc
