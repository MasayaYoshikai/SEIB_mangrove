
!=======================================================================
!  This code was adapted from LeafTemperatureMod.F90                   !
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

module mod_leaf_temperature

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf temperature and energy fluxes
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: leaf_temperature
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine leaf_temperature ( &
                                       ! *** Input ***
    i                             , &  ! Layer index
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Leaf temperature and energy fluxes
    !
    ! !USES:
    use mod_param
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: i
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: lambda                               ! Latent heat of vaporization (J/mol)
    real(8) :: esat                                 ! Saturation vapor pressure (Pa)
    real(8) :: desat                                ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(8) :: term1, term2, term3, term4, term5    ! For computation of leaf temerature
    real(8) :: energy_check                         ! Energy balance check
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    tfrz       => universal%tfrz    , &  ! Freezing point of water (K)
    tair       => atmos%tair        , &  ! Air temperature (K)
    pair       => atmos%pair        , &  ! Air pressure (Pa)
    cp         => atmos%cp          , &  ! Specific heat of air at constant pressure (J/mol/K)
    eair       => atmos%eair        , &  ! Vapor pressure in air (Pa)
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    rn         => flux%rn           , &  ! Leaf net radiation profile (W/m2 leaf)
    gbh        => flux%gbh          , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv        => flux%gbv          , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gs         => flux%gs           , &  ! Leaf stomatal conductance (mol H2O/m2/s)
                                         !
                                         ! *** Output ***
    tleaf      => flux%tleaf        , &  ! Leaf temperature (K)
    shleaf     => flux%shleaf       , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf     => flux%lhleaf       , &  ! Leaf latent heat flux (W/m2 leaf)
    etflux     => flux%etflux         &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    )

    if (leaf_layer(i) == 1) then ! Leaf layer

       ! Saturation vapor pressure

       call sat_vap (tfrz, tair, esat, desat)

       ! Latent heat of vaporization (J/mol)

       lambda = 56780.3d0 - 42.84d0 * tair

       ! Linearized leaf temperature calculation that balances the energy budget

!       term1 = (rn(i) + 2*cp*gbh*tair) * (gbv + gs(i))
!       term2 = (lambda * (esat - desat*tair - eair) * gbv * gs(i)) / pair
!       term3 = 2*cp*gbh*(gbv + gs(i))
!       term4 = (lambda * desat *gbv * gs(i)) / pair
!       tleaf(i) = (term1 - term2) / (term3 + term4)
       term1 = rn(i) * pair * (gs(i) + gbv(i))
       term2 = 2.0d0 * cp * gbh(i) * tair * pair * (gs(i) + gbv(i))
       term3 = lambda * gs(i) * gbv(i) * (esat - desat * tair - eair)
       term4 = 2.0d0 * cp * gbh(i) * pair * (gs(i) + gbv(i))
       term5 = desat * lambda * gs(i) * gbv(i)
       tleaf(i) = (term1 + term2 - term3) / (term4 + term5)

       ! Leaf transpiration rate (mol H2O/m2 leaf/s)

       etflux(i) = (esat + desat*(tleaf(i) - tair) - eair) * gs(i) * gbv(i)    &
                   / (pair * (gs(i) + gbv(i)))

       ! Sensible heat flux (W/m2 leaf)

       shleaf(i) = 2.0d0 * cp * (tleaf(i) - tair) * gbh(i)

       ! Latent heat flux (W/m2 leaf)

       lhleaf(i) = etflux(i) * lambda

       ! Energy balance check

       energy_check = rn(i) - shleaf(i) - lhleaf(i)

       if (abs(energy_check) > 0.001d0) then
          write(*,*) 'Error: Energy balance error', energy_check
       end if

    else ! non-leaf layer

       tleaf(i) = 0.0d0
       etflux(i) = 0.0d0
       shleaf(i) = 0.0d0
       lhleaf(i) = 0.0d0

    end if

    end associate
  end subroutine leaf_temperature

end module mod_leaf_temperature
