MODULE mo_end_bgc

  USE mo_memory_bgc
  USE mo_sedmnt
  USE mo_hamocc_nml, ONLY: l_cpl_co2 

  IMPLICIT NONE

  PUBLIC



CONTAINS

  SUBROUTINE CLEANUP_HAMOCC

    DEALLOCATE (wpoc)
    DEALLOCATE (bgctra)
    DEALLOCATE (bgctend)
    DEALLOCATE (bgcflux)
    DEALLOCATE (hi)
    DEALLOCATE (co3)
    DEALLOCATE (solco2)
    DEALLOCATE (satn2o)
    DEALLOCATE (satoxy)
    DEALLOCATE (satn2)
    DEALLOCATE (sedfluxo)
    DEALLOCATE (aksp)
    DEALLOCATE (aks3)
    DEALLOCATE (akf3)
    DEALLOCATE (ak1p3)
    DEALLOCATE (ak2p3)
    DEALLOCATE (ak3p3)
    DEALLOCATE (aksi3)
    DEALLOCATE (ak23)
    DEALLOCATE (ak13)
    DEALLOCATE (akb3)
    DEALLOCATE (akw3)
    DEALLOCATE (swr_frac)
    DEALLOCATE (meanswr)
    DEALLOCATE (aksurf)
    DEALLOCATE (atm)
    DEALLOCATE (atdifv)
    DEALLOCATE (suppco2)
    DEALLOCATE (co2flux)
    DEALLOCATE (o2flux)
    DEALLOCATE (n2flux)
    DEALLOCATE (n2oflux)

    IF (l_cpl_co2) THEN
      DEALLOCATE (co2trans)
      DEALLOCATE (co2conc)
      DEALLOCATE (co2flux_cpl)
    ENDIF


    DEALLOCATE (expoor)
    DEALLOCATE (expoca)
    DEALLOCATE (exposi)
    DEALLOCATE (kbo)
    DEALLOCATE (strahl)
    DEALLOCATE (bolay)
    DEALLOCATE (alar1max)
    DEALLOCATE (TSFmax)
    DEALLOCATE (TMFmax)
    DEALLOCATE (sedhpl)
    DEALLOCATE (burial)
    DEALLOCATE (powtra)
    DEALLOCATE (silpro)
    DEALLOCATE (prorca)
    DEALLOCATE (prcaca)
    DEALLOCATE (produs)
    DEALLOCATE (dzs)
    DEALLOCATE (seddzi)
    DEALLOCATE (z_sed)
    DEALLOCATE (seddw)
    DEALLOCATE (porsol)
    DEALLOCATE (porwah)
    DEALLOCATE (porwat)
    DEALLOCATE (pors2w)
    DEALLOCATE (pown2bud)
    DEALLOCATE (powh2obud)
    DEALLOCATE (sedtend)
    DEALLOCATE (seddenit)
    DEALLOCATE(sedlay)   
    
  END SUBROUTINE CLEANUP_HAMOCC

END MODULE mo_end_bgc
