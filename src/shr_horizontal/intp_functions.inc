
!-------------------------------------------------------------------------   
!>
!! Gaussian kernel for RBF interpolation.
!! 
!! @f$\phi(r)=e^{-r^2}@f$
!! 
!! @par Revision History
!! Developed by L.Bonaventura  (2004).
!! 
FUNCTION gaussi (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = -1._wp * p_rbf_val * p_rbf_val
p_rbf_val = EXP(p_rbf_val)

END FUNCTION gaussi



!-------------------------------------------------------------------------
!
!  
!>
!! Multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{\frac{1}{2}}@f$
!! 
!! @par Revision History
!! Developed  by L.Bonaventura  (2004).
!! 
FUNCTION multiq (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = SQRT(1._wp + p_rbf_val)

END FUNCTION multiq

!-------------------------------------------------------------------------
!
!  
!>
!! Inverse multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{-\frac{1}{2}}@f$
!! 
!! @par Revision History
!! Developed by L.Bonaventura  (2004).
!! 
FUNCTION inv_multiq (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = SQRT(1._wp + p_rbf_val)
p_rbf_val = 1._wp / p_rbf_val

END FUNCTION inv_multiq

!-------------------------------------------------------------------------
!
!  
!>
!! Weighting function for IDW interpolation.
!! 
!! @f$\phi(r)=\frac{1}{(r)^p_{exp}}@f$
!! 
!! @par Revision History
!! Developed by G. Zaengl  (June 2008).
!! 
FUNCTION invdwgt (p_x, p_exp)  RESULT (wgtfac)
REAL(wp) , INTENT(in) :: p_x            ! radial distance
REAL(wp) , INTENT(in) :: p_exp          ! exponent

REAL(wp)              :: wgtfac       ! weighting factor

!-----------------------------------------------------------------------  

  wgtfac = 1._wp/max(1.e-10_wp , p_x**p_exp)

END FUNCTION invdwgt


!-------------------------------------------------------------------------  
!
!  
!>
!! Modified inverse multiquadric kernel for RBF interpolation.
!! 
!! @f$\phi(r)=(1+r^2)^{-1}@f$
!! 
!! @par Revision History
!! Implemented by G. Zaengl  (2008-09-10).
!! 
FUNCTION inv_multiq2 (p_x, p_scale)  RESULT (p_rbf_val)
REAL(wp) , INTENT(in) :: p_x             ! radial distance
REAL(wp) , INTENT(in) :: p_scale         ! scale parameter

REAL(wp)              :: p_rbf_val       ! RBF value

!-----------------------------------------------------------------------  

p_rbf_val = p_x / p_scale
p_rbf_val = p_rbf_val * p_rbf_val
p_rbf_val = 1._wp / (1._wp + p_rbf_val)

END FUNCTION inv_multiq2

