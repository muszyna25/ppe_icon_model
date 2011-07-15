!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_diffusion_config

  USE mo_kind,           ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish, print_value
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom
  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, apzero

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: configure_diffusion, t_diffusion_config, diffusion_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for diffusion
  !--------------------------------------------------------------------------
  TYPE t_diffusion_config

  INTEGER :: hdiff_order  ! order of horizontal diffusion
                          ! 2: 2nd order linear diffusion on all vertical levels 
                          ! 3: Smagorinsky diffusion for hexagonal model
                          ! 4: 4th order linear diffusion on all vertical levels 
                          ! 5: Smagorinsky diffusion for triangular model
                          ! 24 or 42: 2nd order linear diffusion for upper levels,
                          !           4th order for lower levels

  REAL(wp) :: k2_pres_max  ! (relevant only when hdiff_order = 24 or 42)
                           ! pressure (in Pa) specified by the user
                           ! to determine the lowest vertical level 
                           ! to which 2nd order linear diffusion is applied.
                           ! For the levels with pressure > k2_pres_max, 
                           ! 4th order linear diffusion is applied. 

  INTEGER  :: k2_klev_max  ! (relevant only when hdiff_order = 24 or 42)
                           ! vertical level index specified by the user
                           ! to determine the lowest vertical level 
                           ! to which 2nd order linear diffusion is applied.
                           ! For the levels with k > k2_klev_max, 
                           ! 4th order linear diffusion is applied. 

  REAL(wp) ::           &
    & hdiff_efdt_ratio, &! ratio of e-folding time to (2*)time step
    & hdiff_min_efdt_ratio, &! minimum value of hdiff_efdt_ratio (for upper sponge layer)
    & hdiff_tv_ratio,   &! the ratio of diffusion coefficient: temp:mom
    & hdiff_smag_fac,   &! scaling factor for Smagorinsky diffusion
    & hdiff_multfac      ! multiplication factor of normalized diffusion coefficient
                         ! for nested domains

  REAL(wp) :: &
    & k6, k4, k2       ! numerical diffusion coefficients
                       ! Values for these parameters are not directly
                       ! specified by the user, but derived from the ratio 
                       ! between the e-folding time and the model time step
                       ! (hdiff_efdt_ratio above), and the horizontal 
                       ! resolution of the model

  INTEGER k2s, k2e, k4s, k4e  ! indices defining to which vertical levels
                              ! 2nd and 4th linear diffusion are applied.
                              ! The values are not specified by the user via namelist,
                              ! but determined from k2_klev_max, k2_pres_max
                              ! and the configuration of the vertical coordinate

  LOGICAL ::          &
    & lhdiff_temp,    &! if .TRUE., apply horizontal diffusion to temp.
    & lhdiff_vn        ! if .TRUE., apply horizontal diffusion to momentum.


  END TYPE t_diffusion_config
  !>
  !!
  TYPE(t_diffusion_config) :: diffusion_config(max_dom)

CONTAINS

  SUBROUTINE configure_diffusion(i_ndom,parent_id,nlev)

    INTEGER, INTENT(IN) :: parent_id(max_dom-1) !< list of parent ID's
    INTEGER, INTENT(IN) :: i_ndom !< dimension for time level variables
    INTEGER, INTENT(IN) :: nlev

    INTEGER  :: istat, jg, ist, jk, funit
    REAL(wp) :: zpres(nlev+1)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'setup_diffusion_config:'


    DO jg= 1,i_ndom
      !-----------------------------------------------------------
      ! If using hybrid linear diffusion, set the starting and 
      ! ending vertical level indices for each diffusion order. 
      !-----------------------------------------------------------
      IF (  diffusion_config(jg)%hdiff_order==24 .OR. &
        &   diffusion_config(jg)%hdiff_order==42) THEN                                 
        
        CALL message('','')
        CALL message('----- horizontal diffusion','')
        
        IF ( diffusion_config(jg)%k2_pres_max >0._wp) THEN  ! User has specified a pressure value

          CALL print_value('hdiff:  k2_pres_max (Pa) = ',  diffusion_config(jg)%k2_pres_max)

        ! Calculate the pressure values at layer interfaces
        ! assuming surface pressure is apzero.

        zpres(:) = vct_a(:) + vct_b(:)*apzero

        IF (  diffusion_config(jg)%k2_pres_max <= zpres(1)) THEN
        ! Model does not include mass of the whole atmosphere; User
        ! specified a pressure value located above the model top.
        ! 2nd order diffusion will not be applied. Only 4th order.
        
            diffusion_config(jg)%k2_klev_max = 0
          CALL print_value('hdiff: ptop (Pa)        = ',vct_a(1))
          CALL message('--- hdiff',' k2_pres_max <= ptop')

        ELSE IF (  diffusion_config(jg)%k2_pres_max >= zpres(nlev+1)) THEN
        ! User specified a very high pressure. 2nd order diffusion 
        ! will be applied to all vertical levels.

            diffusion_config(jg)%k2_klev_max = nlev
          CALL print_value('hdiff: pres_sfc (Pa)    = ',zpres(nlev+1))
          CALL message('--- hdiff',' k2_pres_max >= pres_sfc')

        ELSE ! Search for the layer in which k2_pres_max is located.

          DO jk = 1,nlev
            IF ((  diffusion_config(jg)%k2_pres_max > zpres(jk)).AND.&
              & (  diffusion_config(jg)%k2_pres_max <= zpres(jk+1))) THEN
                   diffusion_config(jg)%k2_klev_max = jk 
              EXIT
            END IF
          END DO
          CALL print_value('hdiff: half level pressure (-) = ',&
            &                                      zpres( diffusion_config(jg)%k2_klev_max))
          CALL print_value('hdiff: half level pressure (+) = ',&
            &                                     zpres(diffusion_config(jg)%k2_klev_max+1))

        END IF ! k2_pres_max
      END IF   ! k2_pres_max >0
      ! If the user didn't specifiy a pressure value, then use the 
      ! default or user-specified level index k2_klev_max.

      diffusion_config(jg)%k2s = 1
      diffusion_config(jg)%k2e = diffusion_config(jg)%k2_klev_max
      diffusion_config(jg)%k4s = diffusion_config(jg)%k2_klev_max +1
      diffusion_config(jg)%k4e = nlev

      ! Inform the user about the configuration

      CALL print_value('hdiff: k2s =',diffusion_config(jg)%k2s)
      CALL print_value('hdiff: k2e =',diffusion_config(jg)%k2e)
      CALL print_value('hdiff: k4s =',diffusion_config(jg)%k4s)
      CALL print_value('hdiff: k4e =',diffusion_config(jg)%k4e)
      IF (diffusion_config(jg)%k2e >= diffusion_config(jg)%k2s) THEN
        WRITE(message_text,'(2(a,i4.4))') '2nd order from level ',&
          &                               diffusion_config(jg)%k2s,' to ',diffusion_config(jg)%k2e
        CALL message('--- hdiff',TRIM(message_text))
      END IF
      IF (diffusion_config(jg)%k4e >= diffusion_config(jg)%k4s) THEN
        WRITE(message_text,'(2(a,i4.4))') '4nd order from level ',&
          &                               diffusion_config(jg)%k4s,' to ',diffusion_config(jg)%k4e
        CALL message('--- hdiff',TRIM(message_text))
      END IF
      CALL message('--- hdiff','------')
      CALL message('','')

    END IF !  hdiff_order==24.OR. hdiff_order==42

  ENDDO ! n_dom

   IF ( diffusion_config(jg)%hdiff_efdt_ratio <= 0._wp) THEN
      diffusion_config(:)%k2 = 0._wp
      diffusion_config(:)%k4 = 0._wp
      diffusion_config(:)%k6 = 0._wp
    ELSE
      diffusion_config(1)%k2= 1._wp/( diffusion_config(jg)%hdiff_efdt_ratio*8._wp)
      diffusion_config(1)%k4= 1._wp/( diffusion_config(jg)%hdiff_efdt_ratio*64._wp)
      diffusion_config(1)%k6= 1._wp/( diffusion_config(jg)%hdiff_efdt_ratio*512._wp)

      DO jg = 2, i_ndom
         diffusion_config(jg)%k2 = diffusion_config(parent_id(jg-1))%k2 &
           &                     * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k4 = diffusion_config(parent_id(jg-1))%k4 &
           &                      * diffusion_config(jg)%hdiff_multfac
         diffusion_config(jg)%k6 = diffusion_config(parent_id(jg-1))%k6 &
           &                     * diffusion_config(jg)%hdiff_multfac
      ENDDO
    ENDIF

END SUBROUTINE configure_diffusion


END MODULE mo_diffusion_config
