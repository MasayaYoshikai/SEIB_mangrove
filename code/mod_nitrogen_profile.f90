
!=======================================================================
!  This code was adapted from CanopyNitrogenProfileMod.F90             !
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

module mod_nitrogen_profile

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy profile of nitrogen and photosynthetic capacity
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: nitrogen_profile
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine nitrogen_profile ( &
                                       ! *** Input ***
    t_growth                      , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home                        , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    vcmaxpft                      , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
    universal                     , &  ! Universal variables
                                       !
                                       ! *** Input/Output ***
    layer                           &  ! Layer variables
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: t_growth, t_home
    real(8), intent(in) :: vcmaxpft
    type(universal_type), intent(in) :: universal
    type(layer_type), intent(inout)  :: layer
    !
    ! !LOCAL VARIABLES:
    real(8) :: vcmax25top       ! Canopy top - Maximum carboxylation rate at 25C (umol/m2/s)
    real(8) :: jmax25top        ! Canopy top - Maximum electron transport rate at 25C (umol/m2/s)
    real(8) :: rd25top          ! Canopy top - Leaf respiration rate at 25C (umol CO2/m2/s)
    real(8) :: kn               ! Leaf nitrogen decay coefficient
    real(8) :: nscale           ! Nitrogen scaling coefficient
    integer :: i                ! Layer index
    !---------------------------------------------------------------------

    associate ( &
                                          ! *** Input ***
    n_layer    => universal%n_layer  , &  ! Number of layers (-)
    leaf_layer => layer%leaf_layer   , &  ! 1: Leaf, 0: No leaf
    sumpai     => layer%sumpai       , &  ! Cumulative leaf area (m2 leaf/m2 ground)
                                          !
                                          ! *** Output ***
    vcmax25    => layer%vcmax25      , &  ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    jmax25     => layer%jmax25       , &  ! Maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    rd25       => layer%rd25           &  ! Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
    )

    ! vcmax and other parameters at 25C and top of canopy
    ! Temperature acclimation: Kumarathunge et al. (2019)

    vcmax25top = vcmaxpft
!    jmax25top = 1.67d0 * vcmax25top
    jmax25top = 1.54d0 * vcmax25top        ! <= Aspinwall et al. (2020)
!    rd25top = 0.015d0 * vcmax25top
    rd25top = -0.047d0 * t_growth + 2.30d0 ! <= Aspinwall et al. (2020)

    ! Leaf nitrogen decay coefficient

    kn = exp(0.00963d0 * vcmax25top - 2.43d0)

    ! Layer values

    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer

          nscale = exp (-kn * sumpai(i))
          vcmax25(i) = vcmax25top * nscale
          jmax25(i) = jmax25top * nscale
          rd25(i) = rd25top * nscale

       else ! non-leaf layer

          vcmax25(i) = 0.0d0
          jmax25(i) = 0.0d0
          rd25(i) = 0.0d0

       end if
    end do

    end associate
  end subroutine nitrogen_profile

end module mod_nitrogen_profile
