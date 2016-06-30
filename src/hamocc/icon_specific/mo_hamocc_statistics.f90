!!
!! @author Irene Stemmler, MPI
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_hamocc_statistics
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_ocean_nml,              ONLY: n_zlev 
  USE mo_hamocc_types,           ONLY: t_hamocc_state, t_hamocc_acc, t_hamocc_tend
  USE mo_statistics
  USE mo_model_domain,           ONLY: t_subset_range
  USE mo_sedmnt,                ONLY: ks
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: update_hamocc_statistics
  PUBLIC :: compute_mean_hamocc_statistics
  PUBLIC :: reset_hamocc_statistics
  !-------------------------------------------------------------------------
  
CONTAINS
    
  !---------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_hamocc_statistics(hamocc_state,cells,edges,verts,max_zlev)
    TYPE(t_hamocc_state), INTENT(inout) :: hamocc_state
    TYPE(t_subset_range),      INTENT(IN)    :: cells,edges,verts
    INTEGER, INTENT(in)                      :: max_zlev
    
    INTEGER :: i
    
    ! update hamocc state accumulated values
    CALL add_fields(hamocc_state%p_acc%npp             , hamocc_state%p_tend%npp             , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%graz            , hamocc_state%p_tend%graz            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%graton          , hamocc_state%p_tend%graton          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%exud            , hamocc_state%p_tend%exud            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%exudz           , hamocc_state%p_tend%exudz           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%zoomor          , hamocc_state%p_tend%zoomor           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%phymor          , hamocc_state%p_tend%phymor           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%nfix            , hamocc_state%p_tend%nfix            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%phoc            , hamocc_state%p_tend%phoc            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%cyloss          , hamocc_state%p_tend%cyloss          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%remina          , hamocc_state%p_tend%remina          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%remins          , hamocc_state%p_tend%remins          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%reminn          , hamocc_state%p_tend%reminn          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%bacfra          , hamocc_state%p_tend%bacfra          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%delsil          , hamocc_state%p_tend%delsil          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%delcar          , hamocc_state%p_tend%delcar          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%dmsprod         , hamocc_state%p_tend%dmsprod         , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%dmsbac          , hamocc_state%p_tend%dmsbac          , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%dmsuv           , hamocc_state%p_tend%dmsuv           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%euexp           , hamocc_state%p_tend%euexp           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%flim            , hamocc_state%p_tend%flim            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%nlim            , hamocc_state%p_tend%nlim            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%plim            , hamocc_state%p_tend%plim            , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%aou             , hamocc_state%p_tend%aou             , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%cTlim           , hamocc_state%p_tend%cTlim           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%cLlim           , hamocc_state%p_tend%cLlim           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%cPlim           , hamocc_state%p_tend%cPlim           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%cFlim           , hamocc_state%p_tend%cFlim           , cells,levels=max_zlev)
    CALL add_fields(hamocc_state%p_acc%nfixd           , hamocc_state%p_tend%nfixd           , cells)
    CALL add_fields(hamocc_state%p_acc%cflux           , hamocc_state%p_tend%cflux           , cells)
    CALL add_fields(hamocc_state%p_acc%oflux           , hamocc_state%p_tend%oflux           , cells)
    CALL add_fields(hamocc_state%p_acc%dmsflux         , hamocc_state%p_tend%dmsflux         , cells)
    CALL add_fields(hamocc_state%p_acc%nflux           , hamocc_state%p_tend%nflux           , cells)
    CALL add_fields(hamocc_state%p_acc%n2oflux         , hamocc_state%p_tend%n2oflux         , cells)
    CALL add_fields(hamocc_state%p_acc%orginp          , hamocc_state%p_tend%orginp          , cells)
    CALL add_fields(hamocc_state%p_acc%silinp          , hamocc_state%p_tend%silinp          , cells)
    CALL add_fields(hamocc_state%p_acc%calinp          , hamocc_state%p_tend%calinp          , cells)
    CALL add_fields(hamocc_state%p_acc%prcaca          , hamocc_state%p_tend%prcaca          , cells)
    CALL add_fields(hamocc_state%p_acc%prorca          , hamocc_state%p_tend%prorca          , cells)
    CALL add_fields(hamocc_state%p_acc%silpro          , hamocc_state%p_tend%silpro          , cells)
    CALL add_fields(hamocc_state%p_acc%produs          , hamocc_state%p_tend%produs          , cells)
    CALL add_fields(hamocc_state%p_acc%coex90          , hamocc_state%p_tend%coex90          , cells)
    CALL add_fields(hamocc_state%p_acc%calex90         , hamocc_state%p_tend%calex90         , cells)
    CALL add_fields(hamocc_state%p_acc%opex90          , hamocc_state%p_tend%opex90          , cells)
    CALL add_fields(hamocc_state%p_acc%coex1000          , hamocc_state%p_tend%coex1000          , cells)
    CALL add_fields(hamocc_state%p_acc%calex1000         , hamocc_state%p_tend%calex1000         , cells)
    CALL add_fields(hamocc_state%p_acc%opex1000          , hamocc_state%p_tend%opex1000          , cells)
    CALL add_fields(hamocc_state%p_acc%coex2000          , hamocc_state%p_tend%coex2000          , cells)
    CALL add_fields(hamocc_state%p_acc%calex2000         , hamocc_state%p_tend%calex2000         , cells)
    CALL add_fields(hamocc_state%p_acc%opex2000          , hamocc_state%p_tend%opex2000          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflic         , hamocc_state%p_tend%sedflic          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflal         , hamocc_state%p_tend%sedflal          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflph         , hamocc_state%p_tend%sedflph          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflox         , hamocc_state%p_tend%sedflox          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflsi         , hamocc_state%p_tend%sedflsi          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflfe         , hamocc_state%p_tend%sedflfe          , cells)
    CALL add_fields(hamocc_state%p_acc%sedfln2         , hamocc_state%p_tend%sedfln2          , cells)
    CALL add_fields(hamocc_state%p_acc%sedflno3        , hamocc_state%p_tend%sedflno3         , cells)
    CALL add_fields(hamocc_state%p_acc%sedflic         , hamocc_state%p_tend%sedflic          , cells)
    CALL add_fields(hamocc_state%p_acc%sedrs           , hamocc_state%p_tend%sedrs            , cells,levels=ks)
    CALL add_fields(hamocc_state%p_acc%sedrn           , hamocc_state%p_tend%sedrn            , cells,levels=ks)
    CALL add_fields(hamocc_state%p_acc%sedro2          , hamocc_state%p_tend%sedro2           , cells,levels=ks)
    
    
  END SUBROUTINE update_hamocc_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE compute_mean_hamocc_statistics(p_acc,nsteps_since_last_output)
    TYPE(t_hamocc_acc), INTENT(inout) :: p_acc
    INTEGER,INTENT(in)                     :: nsteps_since_last_output
    
    
!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE
    p_acc%npp                         = p_acc%npp                        /REAL(nsteps_since_last_output,wp)
    p_acc%graz                        = p_acc%graz                       /REAL(nsteps_since_last_output,wp)
    p_acc%exudz                       = p_acc%exudz                      /REAL(nsteps_since_last_output,wp)
    p_acc%zoomor                      = p_acc%zoomor                     /REAL(nsteps_since_last_output,wp)
    p_acc%phymor                      = p_acc%phymor                     /REAL(nsteps_since_last_output,wp)
    p_acc%graton                      = p_acc%graton                     /REAL(nsteps_since_last_output,wp)
    p_acc%exud                        = p_acc%exud                       /REAL(nsteps_since_last_output,wp)
    p_acc%nfix                        = p_acc%nfix                       /REAL(nsteps_since_last_output,wp)
    p_acc%phoc                        = p_acc%phoc                       /REAL(nsteps_since_last_output,wp)
    p_acc%cyloss                      = p_acc%cyloss                     /REAL(nsteps_since_last_output,wp)
    p_acc%nfixd                       = p_acc%nfixd                      /REAL(nsteps_since_last_output,wp)
    p_acc%remina                      = p_acc%remina                     /REAL(nsteps_since_last_output,wp)
    p_acc%remins                      = p_acc%remins                     /REAL(nsteps_since_last_output,wp)
    p_acc%reminn                      = p_acc%reminn                     /REAL(nsteps_since_last_output,wp)
    p_acc%cflux                       = p_acc%cflux                      /REAL(nsteps_since_last_output,wp)
    p_acc%oflux                       = p_acc%oflux                      /REAL(nsteps_since_last_output,wp)
    p_acc%dmsflux                     = p_acc%dmsflux                    /REAL(nsteps_since_last_output,wp)
    p_acc%n2oflux                     = p_acc%n2oflux                    /REAL(nsteps_since_last_output,wp)
    p_acc%prcaca                      = p_acc%prcaca                     /REAL(nsteps_since_last_output,wp)
    p_acc%prorca                      = p_acc%prorca                     /REAL(nsteps_since_last_output,wp)
    p_acc%silpro                      = p_acc%silpro                     /REAL(nsteps_since_last_output,wp)
    p_acc%produs                      = p_acc%produs                     /REAL(nsteps_since_last_output,wp)
    p_acc%coex90                      = p_acc%coex90                     /REAL(nsteps_since_last_output,wp)
    p_acc%calex90                     = p_acc%calex90                    /REAL(nsteps_since_last_output,wp)
    p_acc%opex90                      = p_acc%opex90                     /REAL(nsteps_since_last_output,wp)
    p_acc%coex1000                    = p_acc%coex1000                   /REAL(nsteps_since_last_output,wp)
    p_acc%calex1000                   = p_acc%calex1000                  /REAL(nsteps_since_last_output,wp)
    p_acc%opex1000                    = p_acc%opex1000                   /REAL(nsteps_since_last_output,wp)
    p_acc%coex2000                    = p_acc%coex2000                   /REAL(nsteps_since_last_output,wp)
    p_acc%calex2000                   = p_acc%calex2000                  /REAL(nsteps_since_last_output,wp)
    p_acc%opex2000                    = p_acc%opex2000                   /REAL(nsteps_since_last_output,wp)
    p_acc%delsil                      = p_acc%delsil                     /REAL(nsteps_since_last_output,wp)
    p_acc%orginp                      = p_acc%orginp                     /1._wp
    p_acc%silinp                      = p_acc%silinp                     /1._wp
    p_acc%calinp                      = p_acc%calinp                     /1._wp
    p_acc%sedflic                     = p_acc%sedflic                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflal                     = p_acc%sedflal                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflph                     = p_acc%sedflph                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflox                     = p_acc%sedflox                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflsi                     = p_acc%sedflsi                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflfe                     = p_acc%sedflfe                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedfln2                     = p_acc%sedfln2                    /REAL(nsteps_since_last_output,wp)
    p_acc%sedflno3                    = p_acc%sedflno3                   /REAL(nsteps_since_last_output,wp)
    p_acc%sedrs                       = p_acc%sedrs                      /REAL(nsteps_since_last_output,wp)
    p_acc%sedrn                       = p_acc%sedrn                      /REAL(nsteps_since_last_output,wp)
    p_acc%sedro2                      = p_acc%sedro2                     /REAL(nsteps_since_last_output,wp)
    p_acc%dmsprod                     = p_acc%dmsprod                    /REAL(nsteps_since_last_output,wp)
    p_acc%dmsbac                      = p_acc%dmsbac                     /REAL(nsteps_since_last_output,wp)
    p_acc%dmsuv                       = p_acc%dmsuv                      /REAL(nsteps_since_last_output,wp)
    p_acc%euexp                       = p_acc%euexp                      /REAL(nsteps_since_last_output,wp)
    p_acc%flim                        = p_acc%flim                      /REAL(nsteps_since_last_output,wp)
    p_acc%nlim                        = p_acc%nlim                      /REAL(nsteps_since_last_output,wp)
    p_acc%plim                        = p_acc%plim                      /REAL(nsteps_since_last_output,wp)
    p_acc%aou                         = p_acc%aou                       /REAL(nsteps_since_last_output,wp)
    p_acc%cTlim                       = p_acc%cTlim                     /REAL(nsteps_since_last_output,wp)
    p_acc%cLlim                       = p_acc%cLlim                     /REAL(nsteps_since_last_output,wp)
    p_acc%cPlim                       = p_acc%cPlim                     /REAL(nsteps_since_last_output,wp)
    p_acc%cFlim                       = p_acc%cFlim                     /REAL(nsteps_since_last_output,wp)

!ICON_OMP_END_WORKSHARE
!ICON_OMP_END_PARALLEL
    
      
  END SUBROUTINE compute_mean_hamocc_statistics
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE reset_hamocc_statistics(p_acc,nsteps_since_last_output)
    TYPE(t_hamocc_acc), INTENT(inout)  :: p_acc
    INTEGER,OPTIONAL,        INTENT(inout)  :: nsteps_since_last_output
    
    IF (PRESENT(nsteps_since_last_output)) nsteps_since_last_output        = 0
    
    p_acc%npp                         = 0._wp
    p_acc%graz                        = 0._wp
    p_acc%graton                      = 0._wp
    p_acc%exud                        = 0._wp
    p_acc%exudz                       = 0._wp
    p_acc%zoomor                      = 0._wp
    p_acc%phymor                      = 0._wp
    p_acc%nfix                        = 0._wp
    p_acc%phoc                        = 0._wp
    p_acc%cyloss                      = 0._wp
    p_acc%nfixd                       = 0._wp
    p_acc%remina                      = 0._wp
    p_acc%remins                      = 0._wp
    p_acc%reminn                      = 0._wp
    p_acc%bacfra                      = 0._wp
    p_acc%cflux                       = 0._wp
    p_acc%oflux                       = 0._wp
    p_acc%dmsflux                     = 0._wp
    p_acc%nflux                       = 0._wp
    p_acc%n2oflux                     = 0._wp
    p_acc%prcaca                      = 0._wp
    p_acc%produs                      = 0._wp
    p_acc%prorca                      = 0._wp
    p_acc%silpro                      = 0._wp
    p_acc%coex90                      = 0._wp
    p_acc%calex90                     = 0._wp
    p_acc%opex90                      = 0._wp
    p_acc%coex1000                    = 0._wp
    p_acc%calex1000                   = 0._wp
    p_acc%opex1000                    = 0._wp
    p_acc%coex2000                    = 0._wp
    p_acc%calex2000                   = 0._wp
    p_acc%opex2000                    = 0._wp
    p_acc%delsil                      = 0._wp
    p_acc%delcar                      = 0._wp
    p_acc%sedflic                     = 0._wp
    p_acc%sedflal                     = 0._wp
    p_acc%sedflph                     = 0._wp
    p_acc%sedflox                     = 0._wp
    p_acc%sedflsi                     = 0._wp
    p_acc%sedflfe                     = 0._wp
    p_acc%sedfln2                     = 0._wp
    p_acc%sedflno3                    = 0._wp
    p_acc%sedrn                       = 0._wp
    p_acc%sedrs                       = 0._wp
    p_acc%sedro2                      = 0._wp
!    p_acc%orginp                      = 0._wp
!    p_acc%silinp                      = 0._wp
!    p_acc%calinp                      = 0._wp
    p_acc%dmsbac                      = 0._wp
    p_acc%dmsuv                       = 0._wp
    p_acc%dmsprod                      = 0._wp
    p_acc%euexp                      = 0._wp
    p_acc%flim                       = 0._wp
    p_acc%plim                       = 0._wp
    p_acc%nlim                       = 0._wp
    p_acc%aou                        = 0._wp
    p_acc%cTlim                      = 0._wp
    p_acc%cLlim                      = 0._wp
    p_acc%cPlim                      = 0._wp
    p_acc%cFlim                      = 0._wp
    

  END SUBROUTINE reset_hamocc_statistics
  !---------------------------------------------------------------------
  
  
END MODULE mo_hamocc_statistics
