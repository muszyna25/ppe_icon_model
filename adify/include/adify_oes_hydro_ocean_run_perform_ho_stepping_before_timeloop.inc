! Adjont related memory management {{{
#ifdef __CATCH_NANS__
    CALL catch_nans()
#endif /* __CATCH_NANS__ */


#ifdef __COMPAD_T1S__
    CALL DCO_T1S_SET( ocean_state(jg)%p_prog(nold(1))%h(:,:)         , 1.D0, 1 )
    CALL DCO_T1S_SET( ocean_state(jg)%p_prog(nold(1))%vn(:,:,:)      , 1.D0, 1 )
    CALL DCO_T1S_SET( ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1), 1.D0, 1 )
    CALL DCO_T1S_SET( ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2), 1.D0, 1 )
#endif /* __COMPAD_T1S__ */
#ifdef __COMPAD_A1S__
    ! allocate adjoint inputs
    szAB_H = SHAPE(ocean_state(jg)%p_prog(nnew(1))%h(:,:))
    ALLOCATE( ADJINP_H( szAB_H(1),szAB_H(2) ) )
    ALLOCATE( ADJINP_H0( szAB_H(1),szAB_H(2) ) )
    szAB_VN = SHAPE(ocean_state(jg)%p_prog(nnew(1))%vn(:,:,:))
    ALLOCATE( ADJINP_VN( szAB_VN(1),szAB_VN(2),szAB_VN(3) ) )
    ALLOCATE( ADJINP_VN0( szAB_VN(1),szAB_VN(2),szAB_VN(3) ) )
    szAB_TR = SHAPE(ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,1))
    ALLOCATE( ADJINP_TR( szAB_TR(1),szAB_TR(2),szAB_TR(3), 2 ) )
    ALLOCATE( ADJINP_TR0( szAB_TR(1),szAB_TR(2),szAB_TR(3), 2 ) )
    szAB_GNM1 = SHAPE(ocean_state(jg)%p_aux%g_n(:,:,:))
    ALLOCATE( ADJINP_GNM1( szAB_GNM1(1),szAB_GNM1(2),szAB_GNM1(3) ) )
    szAB_DW = SHAPE(ocean_state(jg)%p_diag%w(:,:,:))
    ALLOCATE( ADJINP_DW( szAB_DW(1),szAB_DW(2),szAB_DW(3) ))

#  ifdef __COMPAD_A1S_NOCP__
    ! assign values to adjoint inputs
    ADJINP_H(:,:)        = ocean_state(jg)%p_prog(nold(1))%h(:,:)
    ADJINP_VN(:,:,:)     = ocean_state(jg)%p_prog(nold(1))%vn(:,:,:)
    ADJINP_TR(:,:,:,:)   = ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1:2)
    ADJINP_H0(:,:)       = ocean_state(jg)%p_prog(nnew(1))%h(:,:)
    ADJINP_VN0(:,:,:)    = ocean_state(jg)%p_prog(nnew(1))%vn(:,:,:)
    ADJINP_TR0(:,:,:,:)  = ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,1:2)
    ADJINP_GNM1(:,:,:)   = ocean_state(jg)%p_aux%g_nm1(:,:,:)
    ADJINP_DW(:,:,:)     = ocean_state(jg)%p_diag%w(:,:,:)

    ! CREATE tape
    CALL DCO_A1S_TAPE_CREATE
    ! SWITCH ON tape
    CALL DCO_A1S_TAPE_SWITCH_TO_ACTIVE
    !    WRITE(0,*) " TAPE STAT ::", DCO_A1S_TAPE_GET_TAPE_MEMORY_SIZE()
    ! register inputs
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_h(:,:)          )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_vn(:,:,:)       )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_tr(:,:,:,:)     )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_h0(:,:)         )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_vn0(:,:,:)      )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_tr0(:,:,:,:)    )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_GNM1(:,:,:)     )
    CALL DCO_A1S_TAPE_REGISTER_VARIABLE( ADJINP_DW(:,:,:)       )
    !    WRITE(0,*) " TAPE STAT ::", DCO_A1S_TAPE_GET_TAPE_MEMORY_SIZE()
    ! assign inputs to variables
    ocean_state(jg)%p_prog(nold(1))%h(:,:)          = ADJINP_H(:,:)
    ocean_state(jg)%p_prog(nold(1))%vn(:,:,:)       = ADJINP_VN(:,:,:)
    ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,1) = ADJINP_TR(:,:,:,1)
    ocean_state(jg)%p_prog(nold(1))%tracer(:,:,:,2) = ADJINP_TR(:,:,:,2)
    ocean_state(jg)%p_prog(nnew(1))%h(:,:)          = ADJINP_H0(:,:)
    ocean_state(jg)%p_prog(nnew(1))%vn(:,:,:)       = ADJINP_VN0(:,:,:)
    ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,1) = ADJINP_TR0(:,:,:,1)
    ocean_state(jg)%p_prog(nnew(1))%tracer(:,:,:,2) = ADJINP_TR0(:,:,:,2)
    ocean_state(jg)%p_aux%g_nm1(:,:,:)              = ADJINP_GNM1(:,:,:)
    ocean_state(jg)%p_diag%w(:,:,:)                 = ADJINP_DW(:,:,:)

#  endif /* __COMPAD_A1S_NOCP__ */
#endif /* __COMPAD_A1S__ */
! }}}
