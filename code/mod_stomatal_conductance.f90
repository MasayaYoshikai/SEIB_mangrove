
!=======================================================================
!  This code was adapted from StomataOptimizationMod.F90               !
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
!  Xu X., Medvigy D., Powers J. S., Becknell J. M., Guan K. (2016),    !
!  Diversity in plant hydraulic traits explains seasonal and           !
!  inter-annual variations of vegetation dynamics in seasonally dry    !
!  tropical forests, New Phytologist, 212, 80-95.                      !
!                                                                      !
!=======================================================================

module mod_stomatal_conductance

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: stomatal_conductance
  public  :: dark_condition
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: stomata_efficiency        ! Water-use efficiency for maximum gs
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine stomatal_conductance ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota0                         , &  ! Reference marginal water use efficiency in well-watered condition (umol CO2/mol H2O)
    leaf_b0                       , &  ! Sensitivity of marginal water use efficiency to leaf water potential (MPa-1)
    leaf_pdwp                     , &  ! Pre-dawn leaf water potential of the last day (MPa)
    leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Photosynthesis and stomatal conductance
    ! with optimization of water use efficiency
    !
    ! !USES:
    use mod_param
    use mod_math_tools
    use mod_leaf_temperature
    use mod_photosynthesis
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    real(8), intent(in) :: iota0, leaf_b0, leaf_pdwp, leaf_wp
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    integer :: i                         ! Layer index
    integer :: skip_flag                 ! Calculation skip flag (0 for calculation)
    real(8) :: iota                      ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    real(8) :: gs1, gs2                  ! Initial guess for gs (mol H2O/m2/s)
    real(8) :: check1, check2            ! Water-use efficiency for gs1 and gs2
    real(8), parameter :: tol = 0.004d0  ! gs is updated to accuracy tol (mol H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    n_layer    => universal%n_layer , &  ! Number of leaf layers
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    top_layer  => layer%top_layer   , &  ! Top layer of the canopy
    skip_layer => layer%skip_layer  , &  ! Layer interval to skip the flux calculation
                                         ! *** Output ***
    tleaf      => flux%tleaf        , &  ! Leaf temperature (K)
    shleaf     => flux%shleaf       , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf     => flux%lhleaf       , &  ! Leaf latent heat flux (W/m2 leaf)
    etflux     => flux%etflux       , &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    rd         => flux%rd           , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci         => flux%ci           , &  ! Leaf intercellular CO2 (umol/mol)
    hs         => flux%hs           , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd        => flux%vpd          , &  ! Leaf vapor pressure deficit at surface (Pa)
    ac         => flux%ac           , &  ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj         => flux%aj           , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ag         => flux%ag           , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an         => flux%an           , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs         => flux%cs           , &  ! Leaf surface CO2 (umol/mol)
    gs         => flux%gs             &  ! Leaf stomatal conductance (mol H2O/m2/s)
    )

    ! Stomatal water-use efficiency at the given pre-dawn
    ! leaf water potential (umol CO2/ mol H2O)

    iota = iota0 * exp(leaf_b0 * leaf_pdwp)

    ! Low and high initial estimates for gs (mol H2O/m2/s)

    gs1 = 0.002d0
    gs2 = 2.0d0

    ! Calculate gs for each canopy layer

    do i = n_layer, 1, -1
       if (leaf_layer(i) == 1) then ! Leaf layer

          skip_flag = (top_layer-i) - ((top_layer-i)/skip_layer)*skip_layer

          if (skip_flag == 0) then

             ! Check for minimum stomatal conductance linked to low light
             ! or drought stress based on the water-use efficiency and
             ! cavitation checks for gs1 and gs2

             check1 =stomata_efficiency (    &
                                                ! *** Input ***
             o2air                         , &  ! Atmospheric O2 (mmol/mol)
             co2air                        , &  ! Atmospheric CO2 (umol/mol)
             iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
             i                             , &  ! Layer index
             gs1                           , &  ! Value for gs to use in calculations
             universal                     , &  ! Universal variables
             atmos                         , &  ! Atmospheric variables
             layer                         , &  ! Layer variables
                                                !
                                                ! *** Input/Output ***
             flux                            &  ! Flux variables
             )

             check2 = stomata_efficiency (    &
                                                ! *** Input ***
             o2air                         , &  ! Atmospheric O2 (mmol/mol)
             co2air                        , &  ! Atmospheric CO2 (umol/mol)
             iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
             i                             , &  ! Layer index
             gs2                           , &  ! Value for gs to use in calculations
             universal                     , &  ! Universal variables
             atmos                         , &  ! Atmospheric variables
             layer                         , &  ! Layer variables
                                                !
                                                ! *** Input/Output ***
             flux                            &  ! Flux variables
             )

             if (check1 * check2 < 0.0d0) then

                ! Calculate gs using the function stomata_efficiency to
                ! iterate gs to an accuracy of tol (mol H2O/m2/s)

                gs(i) = zbrent ('Stomatal optimaization', stomata_efficiency,  &
                               o2air, co2air, iota, i, universal, atmos,       &
                               layer, flux, gs1, gs2, tol)

             else

                ! Low light. Set gs to minimum conductance

                gs(i) = 0.002d0

             end if

             ! Leaf fluxes with the specified gs

             call leaf_temperature (           &
                                                  ! *** Input ***
             i                               , &  ! Layer index
             universal                       , &  ! Universal variables
             atmos                           , &  ! Atmospheric variables
             layer                           , &  ! Layer variables
                                                  !
                                                  ! *** Input/Output ***
             flux                              &  ! Flux variables
             )

             ! Leaf photosynthesis

             call leaf_photosynthesis (        &
                                                  ! *** Input ***
             o2air                           , &  ! Atmospheric O2 (mmol/mol)
             co2air                          , &  ! Atmospheric CO2 (umol/mol)
             i                               , &  ! Layer index
             universal                       , &  ! Universal variables
             atmos                           , &  ! Atmospheric variables
             layer                           , &  ! Layer variables
                                                  !
                                                  ! *** Input/Output ***
             flux                              &  ! Flux variables
             )

          else ! skip calculation

             gs(i)     = gs(i+1)
             tleaf(i)  = tleaf(i+1)
             shleaf(i) = shleaf(i+1)
             lhleaf(i) = lhleaf(i+1)
             etflux(i) = etflux(i+1)
             rd(i)     = rd(i+1)
             ci(i)     = ci(i+1)
             hs(i)     = hs(i+1)
             vpd(i)    = vpd(i+1)
             ac(i)     = ac(i+1)
             aj(i)     = aj(i+1)
             ag(i)     = ag(i+1)
             an(i)     = an(i+1)
             cs(i)     = cs(i+1)

          end if

       else ! non-leaf layer

          gs(i) = 0.0d0
          tleaf(i)  = 0.0d0
          shleaf(i) = 0.0d0
          lhleaf(i) = 0.0d0
          etflux(i) = 0.0d0
          rd(i)     = 0.0d0
          ci(i)     = 0.0d0
          hs(i)     = 0.0d0
          vpd(i)    = 0.0d0
          ac(i)     = 0.0d0
          aj(i)     = 0.0d0
          ag(i)     = 0.0d0
          an(i)     = 0.0d0
          cs(i)     = 0.0d0

       end if

    end do

    end associate
  end subroutine stomatal_conductance

  !-----------------------------------------------------------------------
  function stomata_efficiency ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    i                             , &  ! Layer index
    gs_val                        , &  ! Value for gs to use in calculations
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    ) result(val)
    !
    ! !DESCRIPTION:
    ! Stomata water-use efficiency check to determine maximum gs.
    ! For the stomatal conductance gs_val, calculate photosynthesis
    ! for an increase in stomatal conductance equal to "delta".
    ! The returned value is positive if this increase produces a change in
    ! photosynthesis > iota*vpd*delta.
    ! The returned value is negative if the increase produces a change in
    ! photosynthesis < iota*vpd*delta.
    !
    ! !USES:
    use mod_param
    use mod_leaf_temperature
    use mod_photosynthesis
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    real(8), intent(in) :: iota
    integer, intent(in) :: i
    real(8), intent(in) :: gs_val
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: delta        ! Small difference for gs (mol H2O/m2/s)
    real(8) :: gs2          ! Lower value for gs (mol H2O/m2/s)
    real(8) :: an2          ! Leaf photosynthesis at gs2 (umol CO2/m2/s)
    real(8) :: gs1          ! Higher value for gs (mol H2O/m2/s)
    real(8) :: an1          ! Leaf photosynthesis at gs1 (umol CO2/m2/s)
    real(8) :: wue          ! Water-use efficiency check
    real(8) :: val          ! Returned value
    !---------------------------------------------------------------------

    associate ( &
    pair       => atmos%pair         , &  ! Air pressure (Pa)
    gs         => flux%gs            , &  ! Leaf stomatal conductance (mol H2O/m2/s)
    an         => flux%an            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    vpd        => flux%vpd             &  ! Leaf vapor pressure deficit at surface (Pa)
    )

    ! Specify "delta" as a small difference in gs (mol H2O/m2/s)

    delta = 0.001d0

    ! --- Flux calculation with lower gs (gs_val - delta)

    gs2 = gs_val - delta
    gs(i)  = gs2

    ! Leaf temperature and energy fluxes

    call leaf_temperature (           &
                                         ! *** Input ***
    i                               , &  ! Layer index
    universal                       , &  ! Universal variables
    atmos                           , &  ! Atmospheric variables
    layer                           , &  ! Layer variables
                                         !
                                         ! *** Input/Output ***
    flux                              &  ! Flux variables
    )

    ! Leaf photosynthesis

    call leaf_photosynthesis (        &
                                         ! *** Input ***
    o2air                           , &  ! Atmospheric O2 (mmol/mol)
    co2air                          , &  ! Atmospheric CO2 (umol/mol)
    i                               , &  ! Layer index
    universal                       , &  ! Universal variables
    atmos                           , &  ! Atmospheric variables
    layer                           , &  ! Layer variables
                                         !
                                         ! *** Input/Output ***
    flux                              &  ! Flux variables
    )

    ! Photosynthesis at lower gs (gs_val - delta)

    an2 = an(i)

    ! --- Flux calculation with higher gs (gs_val)

    gs1 = gs_val
    gs(i)  = gs1

    ! Leaf temperature and energy fluxes

    call leaf_temperature (           &
                                         ! *** Input ***
    i                               , &  ! Layer index
    universal                       , &  ! Universal variables
    atmos                           , &  ! Atmospheric variables
    layer                           , &  ! Layer variables
                                         !
                                         ! *** Input/Output ***
    flux                              &  ! Flux variables
    )

    ! Leaf photosynthesis

    call leaf_photosynthesis (        &
                                         ! *** Input ***
    o2air                           , &  ! Atmospheric O2 (mmol/mol)
    co2air                          , &  ! Atmospheric CO2 (umol/mol)
    i                               , &  ! Layer index
    universal                       , &  ! Universal variables
    atmos                           , &  ! Atmospheric variables
    layer                           , &  ! Layer variables
                                         !
                                         ! *** Input/Output ***
    flux                              &  ! Flux variables
    )

    ! Photosynthesis at higher gs (gs_val)

    an1 = an(i)

    ! Efficiency check: wue < 0 when d(An) / d(gs) < iota * vpd

    wue = (an1 - an2) - iota * delta * (vpd(i) / pair)

    ! Return value

    val = wue

    end associate
  end function stomata_efficiency

  !-----------------------------------------------------------------------
  subroutine dark_condition ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Leaf fluxes for dark condition. gs is regulated to its minimum value.
    !
    ! !USES:
    use mod_param
    use mod_leaf_temperature
    use mod_photosynthesis
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    integer :: i
    integer :: skip_flag
    !---------------------------------------------------------------------

    associate ( &
                                           ! *** Input ***
    n_layer      => universal%n_layer , &  ! Number of leaf layers
    leaf_layer   => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    top_layer    => layer%top_layer   , &  ! Top layer of the canopy
    skip_layer   => layer%skip_layer  , &  ! Layer interval to skip the flux calculation
                                           ! *** Output ***
    tleaf        => flux%tleaf        , &  ! Leaf temperature (K)
    shleaf       => flux%shleaf       , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf       => flux%lhleaf       , &  ! Leaf latent heat flux (W/m2 leaf)
    etflux       => flux%etflux       , &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    rd           => flux%rd           , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci           => flux%ci           , &  ! Leaf intercellular CO2 (umol/mol)
    hs           => flux%hs           , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd          => flux%vpd          , &  ! Leaf vapor pressure deficit at surface (Pa)
    ac           => flux%ac           , &  ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj           => flux%aj           , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ag           => flux%ag           , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an           => flux%an           , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs           => flux%cs           , &  ! Leaf surface CO2 (umol/mol)
    gs           => flux%gs           , &  ! Leaf stomatal conductance (mol H2O/m2/s)
    dark_gs      => flux%dark_gs      , &  ! Flux variables for dark condition
    dark_tleaf   => flux%dark_tleaf   , &  ! Flux variables for dark condition
    dark_shleaf  => flux%dark_shleaf  , &  ! Flux variables for dark condition
    dark_lhleaf  => flux%dark_lhleaf  , &  ! Flux variables for dark condition
    dark_etflux  => flux%dark_etflux  , &  ! Flux variables for dark condition
    dark_rd      => flux% dark_rd     , &  ! Flux variables for dark condition
    dark_ci      => flux%dark_ci      , &  ! Flux variables for dark condition
    dark_hs      => flux%dark_hs      , &  ! Flux variables for dark condition
    dark_vpd     => flux%dark_vpd     , &  ! Flux variables for dark condition
    dark_ac      => flux%dark_ac      , &  ! Flux variables for dark condition
    dark_aj      => flux%dark_aj      , &  ! Flux variables for dark condition
    dark_ag      => flux%dark_ag      , &  ! Flux variables for dark condition
    dark_an      => flux%dark_an      , &  ! Flux variables for dark condition
    dark_cs      => flux%dark_cs        &  ! Flux variables for dark condition
    )

    do i = n_layer, 1, -1
       if (leaf_layer(i) == 1) then ! Leaf layer
          skip_flag = (top_layer-i)-((top_layer-i)/skip_layer)*skip_layer
          if (skip_flag == 0) then
             ! Set gs to minimum value.
             gs(i) = 0.002d0
             ! Leaf temperature and energy fluxes
             call leaf_temperature (           &
                                                  ! *** Input ***
             i                               , &  ! Layer index
             universal                       , &  ! Universal variables
             atmos                           , &  ! Atmospheric variables
             layer                           , &  ! Layer variables
                                                  !
                                                  ! *** Input/Output ***
             flux                              &  ! Flux variables
             )
             ! Leaf photosynthesis
             call leaf_photosynthesis (        &
                                                  ! *** Input ***
             o2air                           , &  ! Atmospheric O2 (mmol/mol)
             co2air                          , &  ! Atmospheric CO2 (umol/mol)
             i                               , &  ! Layer index
             universal                       , &  ! Universal variables
             atmos                           , &  ! Atmospheric variables
             layer                           , &  ! Layer variables
                                                  !
                                                  ! *** Input/Output ***
             flux                              &  ! Flux variables
             )
          else
             gs(i)     = gs(i+1)
             tleaf(i)  = tleaf(i+1)
             shleaf(i) = shleaf(i+1)
             lhleaf(i) = lhleaf(i+1)
             etflux(i) = etflux(i+1)
             rd(i)     = rd(i+1)
             ci(i)     = ci(i+1)
             hs(i)     = hs(i+1)
             vpd(i)    = vpd(i+1)
             ac(i)     = ac(i+1)
             aj(i)     = aj(i+1)
             ag(i)     = ag(i+1)
             an(i)     = an(i+1)
             cs(i)     = cs(i+1)
          end if
       else ! non-leaf layer
         gs(i)     = 0.0d0
         tleaf(i)  = 0.0d0
         shleaf(i) = 0.0d0
         lhleaf(i) = 0.0d0
         etflux(i) = 0.0d0
         rd(i)     = 0.0d0
         ci(i)     = 0.0d0
         hs(i)     = 0.0d0
         vpd(i)    = 0.0d0
         ac(i)     = 0.0d0
         aj(i)     = 0.0d0
         ag(i)     = 0.0d0
         an(i)     = 0.0d0
         cs(i)     = 0.0d0
       end if
    end do

    ! Save flux variables for the dark condition.

    dark_gs     = gs
    dark_tleaf  = tleaf
    dark_shleaf = shleaf
    dark_lhleaf = lhleaf
    dark_etflux = etflux
    dark_rd     = rd
    dark_ci     = ci
    dark_hs     = hs
    dark_vpd    = vpd
    dark_ac     = ac
    dark_aj     = aj
    dark_ag     = ag
    dark_an     = an
    dark_cs     = cs

    end associate
  end subroutine dark_condition

end module mod_stomatal_conductance
