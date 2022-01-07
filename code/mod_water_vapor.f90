
!=======================================================================
!  This code was adapted from WaterVaporMod.F90                        !
!  in the link below.                                                  !
!  https://github.com/gbonan/CLM-ml_v0                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Bonan G. B., Williams M., Fisher R. A., Oleson K. W. (2014),        !
!  Modeling stomatal conductance in the earth system: linking leaf     !
!  water-use efficiency and water transport along the                  !
!  soil–plant–atmosphere continuum, Geoscientific Model Development,   !
!  7, 2019-2222.                                                       !
!                                                                      !
!=======================================================================

module mod_water_vapor

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate saturation vapor pressure and latent heat of vaporization
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: sat_vap     ! Saturation vapor pressure and derivative
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine sat_vap (tfrz, t, es, desdt)
    !
    ! !DESCRIPTION:
    ! Compute saturation vapor pressure and change in saturation vapor pressure
    ! with respect to temperature. Polynomial approximations are from:
    ! Flatau et al (1992) Polynomial fits to saturation vapor pressure.
    ! Journal of Applied Meteorology 31:1507-1513
    !
    ! !USES:
    !
    !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: tfrz     ! Freezing point of water (K)
    real(8), intent(in)  :: t        ! Temperature (K)
    real(8), intent(out) :: es       ! Vapor pressure (Pa)
    real(8), intent(out) :: desdt    ! d(es)/d(t) (Pa/K)
    !
    ! !LOCAL VARIABLES:
    real(8) :: tc                    ! Temperature (C)
    !---------------------------------------------------------------------

    ! For water vapor (temperature range is 0C to 100C)

    real(8), parameter :: a0 =  6.11213476d0
    real(8), parameter :: a1 =  0.444007856d0
    real(8), parameter :: a2 =  0.143064234e-01
    real(8), parameter :: a3 =  0.264461437e-03
    real(8), parameter :: a4 =  0.305903558e-05
    real(8), parameter :: a5 =  0.196237241e-07
    real(8), parameter :: a6 =  0.892344772e-10
    real(8), parameter :: a7 = -0.373208410e-12
    real(8), parameter :: a8 =  0.209339997e-15

    ! and for derivative

    real(8), parameter :: b0 =  0.444017302d0
    real(8), parameter :: b1 =  0.286064092e-01
    real(8), parameter :: b2 =  0.794683137e-03
    real(8), parameter :: b3 =  0.121211669e-04
    real(8), parameter :: b4 =  0.103354611e-06
    real(8), parameter :: b5 =  0.404125005e-09
    real(8), parameter :: b6 = -0.788037859e-12
    real(8), parameter :: b7 = -0.114596802e-13
    real(8), parameter :: b8 =  0.381294516e-16

    tc = t - tfrz

    es    = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 &
          + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))))
    desdt = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 &
          + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))))

    es    = es    * 100.0d0            ! Convert from mb to Pa
    desdt = desdt * 100.0d0            ! Convert from mb to Pa

  end subroutine sat_vap

end module mod_water_vapor
