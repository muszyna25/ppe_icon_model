!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!option! -Nv -NO
#ifdef VPP
!OCL SCALAR
#endif
#ifdef __xlC__
@PROCESS NOOPTIMIZE
#endif
SUBROUTINE srtm_kgb26

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF

!     ------------------------------------------------------------------

USE mo_kind, ONLY : wp

USE mo_yoesrta26, ONLY : sfluxref, rayl

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
sfluxref = (/ &
 !  &     129.462_wp, 15*_ZERO_ /)
 & 29.0079_wp,  28.4088_wp,     20.3099_wp,  13.0283_wp &
 & ,  11.8619_wp,  9.95840_wp,     6.68696_wp,  5.38987_wp &
 & ,  3.49829_wp, 0.407693_wp,    0.299027_wp, 0.236827_wp &
 & , 0.188502_wp, 0.163489_wp, 4.64335e-02_wp, 2.72662e-03_wp /)

!     Rayleigh extinction coefficient at all v
rayl = (/ &
 & 1.21263e-06_wp,1.43428e-06_wp,1.67677e-06_wp,1.93255e-06_wp &
 & , 2.19177e-06_wp,2.44195e-06_wp,2.66926e-06_wp,2.85990e-06_wp &
 & , 3.00380e-06_wp,3.06996e-06_wp,3.08184e-06_wp,3.09172e-06_wp &
 & , 3.09938e-06_wp,3.10456e-06_wp,3.10727e-06_wp,3.10818e-06_wp /)

!     ------------------------------------------------------------------
END SUBROUTINE srtm_kgb26
