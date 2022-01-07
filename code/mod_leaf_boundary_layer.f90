
!=======================================================================
!  This code was adapted from LeafBoundaryLayerMod.F90                 !
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

module mod_leaf_boundary_layer

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf boundary layer conductances
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: leaf_boundary_layer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine leaf_boundary_layer ( &
                                       ! *** Input ***
    dleaf                         , &  ! Leaf dimension (m)
    rwind_profile                 , &  ! Relative wind speed profile to the canopy top (-)
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: dleaf
    real(8), intent(in) :: rwind_profile(:)
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    integer :: i                  ! Layer index
    real(8) :: wind_canopy        ! Wind speed within canopy at layer i (m/s)
    real(8) :: fac                ! Correction factor for temperature and pressure
    real(8) :: visc               ! Kinematic viscosity (m2/s)
    real(8) :: Dh                 ! Molecular diffusivity, heat (m2/s)
    real(8) :: Dv                 ! Molecular diffusivity, H2O (m2/s)
    real(8) :: Dc                 ! Molecular diffusivity, CO2 (m2/s)
    real(8) :: Re                 ! Reynolds number (dimensionless)
    real(8) :: Pr                 ! Prandtl number (dimensionless)
    real(8) :: Scv                ! Schmidt number, H2O (dimensionless)
    real(8) :: Scc                ! Schmidt number, CO2 (dimensionless)
    real(8) :: b1                 ! Empirical correction factor for Nu
    real(8) :: Nu_lam             ! Forced convection - laminar: Nusselt number (dimensionless)
    real(8) :: Shv_lam            ! Forced convection - laminar: Sherwood number, H2O (dimensionless)
    real(8) :: Shc_lam            ! Forced convection - laminar: Sherwood number, CO2 (dimensionless)
    real(8) :: Nu_turb            ! Forced convection - turbulent: Nusselt number (dimensionless)
    real(8) :: Shv_turb           ! Forced convection - turbulent: Sherwood number, H2O (dimensionless)
    real(8) :: Shc_turb           ! Forced convection - turbulent: Sherwood number, CO2 (dimensionless)
    real(8) :: Nu_forced          ! Forced convection: Nusselt number (dimensionless)
    real(8) :: Shv_forced         ! Forced convection: Sherwood number, H2O (dimensionless)
    real(8) :: Shc_forced         ! Forced convection: Sherwood number, CO2 (dimensionless)
    real(8) :: Nu                 ! Nusselt number (dimensionless)
    real(8) :: Shv                ! Sherwood number, H2O (dimensionless)
    real(8) :: Shc                ! Sherwood number, CO2 (dimensionless)
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    tfrz       => universal%tfrz    , &  ! Freezing point of water (K)
    n_layer    => universal%n_layer , &  ! Number of leaf layers
    visc0      => universal%visc0   , &  ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
    Dh0        => universal%Dh0     , &  ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
    Dv0        => universal%Dv0     , &  ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
    Dc0        => universal%Dc0     , &  ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
    pair       => atmos%pair        , &  ! Atmospheric pressure (Pa)
    rhomol     => atmos%rhomol      , &  ! Molar density (mol/m3)
    wind       => atmos%wind        , &  ! Wind speed (m/s)
    tair       => atmos%tair        , &  ! Air temperature (K)
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
                                         !
                                         ! *** Output ***
    gbh        => flux%gbh          , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv        => flux%gbv          , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gbc        => flux%gbc            &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    )

    ! Empirical correction factor for Nu and Sh

    b1 = 1.5d0

    ! Computation for leaf layer

    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer

          ! Wind speed at the layer (m/s)

          wind_canopy = wind * rwind_profile(i)

          ! Adjust diffusivity for temperature and pressure

          fac = 101325.0d0 / pair * (tair / tfrz)**1.81d0

          visc = visc0 * fac                               ! Kinematic viscosity (m2/s)
          Dh = Dh0 * fac                                   ! Molecular diffusivity, heat (m2/s)
          Dv = Dv0 * fac                                   ! Molecular diffusivity, H2O (m2/s)
          Dc = Dc0 * fac                                   ! Molecular diffusivity, CO2 (m2/s)

          ! Dimensionless numbers

          Re = wind_canopy * dleaf / visc                  ! Reynolds number
          Pr = visc / Dh                                   ! Prandtl number
          Scv = visc / Dv                                  ! Schmidt number for H2O
          Scc = visc / Dc                                  ! Schmidt number for CO2

          ! Nusselt number (Nu) and Sherwood numbers (H2O: Shv, CO2: Shc)

             ! Forced convection - laminar flow

             Nu_lam  = b1 * 0.66d0 *  Pr**0.33d0 * Re**0.5d0     ! Nusselt number
             Shv_lam = b1 * 0.66d0 * Scv**0.33d0 * Re**0.5d0     ! Sherwood number, H2O
             Shc_lam = b1 * 0.66d0 * Scc**0.33d0 * Re**0.5d0     ! Sherwood number, CO2

             ! Forced convection - turbulent flow

             Nu_turb  = b1 * 0.036d0 *  Pr**0.33d0 * Re**0.8d0   ! Nusselt number
             Shv_turb = b1 * 0.036d0 * Scv**0.33d0 * Re**0.8d0   ! Sherwood number, H2O
             Shc_turb = b1 * 0.036d0 * Scc**0.33d0 * Re**0.8d0   ! Sherwood number, CO2

             ! Choose correct flow regime for forced convection

             Nu_forced = max(Nu_lam, Nu_turb)
             Shv_forced = max(Shv_lam, Shv_turb)
             Shc_forced = max(Shc_lam, Shc_turb)

             ! Effects of free convection was removed

             Nu = Nu_forced
             Shv = Shv_forced
             Shc = Shc_forced

             ! Boundary layer conductances (m/s)

             gbh(i) = Dh *  Nu / dleaf
             gbv(i) = Dv * Shv / dleaf
             gbc(i) = Dc * Shc / dleaf

             ! Convert conductance (m/s) to (mol/m2/s)

             gbh(i) = gbh(i) * rhomol
             gbv(i) = gbv(i) * rhomol
             gbc(i) = gbc(i) * rhomol

       else ! No leaf layer

          gbh(i) = 0.0d0
          gbv(i) = 0.0d0
          gbc(i) = 0.0d0

       end if

    end do

    end associate
  end subroutine leaf_boundary_layer

end module mod_leaf_boundary_layer
