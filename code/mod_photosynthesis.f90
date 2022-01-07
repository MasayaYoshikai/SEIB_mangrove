
!=======================================================================
!  This code was adapted from LeafPhotosynthesisMod.F90                !
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

module mod_photosynthesis

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate leaf photosynthesis and stomatal conductance
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  !
  ! Leaf-level photosynthesis variables
  !
  ! *** Input parameters ***
  !
  private
  real(8) :: kc25        ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
  real(8) :: ko25        ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
  real(8) :: cp25        ! CO2 compensation point at 25C (umol/mol)

  real(8) :: vcmaxha     ! Activation energy for vcmax (J/mol)
  real(8) :: jmaxha      ! Activation energy for jmax (J/mol)
  real(8) :: rdha        ! Activation energy for rd (J/mol)
  real(8) :: kcha        ! Activation energy for kc (J/mol)
  real(8) :: koha        ! Activation energy for ko (J/mol)
  real(8) :: cpha        ! Activation energy for cp (J/mol)

  real(8) :: vcmaxhd     ! Deactivation energy for vcmax (J/mol)
  real(8) :: jmaxhd      ! Deactivation energy for jmax (J/mol)
  real(8) :: rdhd        ! Deactivation energy for rd (J/mol)

  real(8) :: vcmaxse     ! Entropy term for vcmax (J/mol/K)
  real(8) :: jmaxse      ! Entropy term for jmax (J/mol/K)
  real(8) :: rdse        ! Entropy term for rd (J/mol/K)

  real(8) :: vcmaxc      ! Scaling factor for high temperature inhibition (25 C = 1.0)
  real(8) :: jmaxc       ! Scaling factor for high temperature inhibition (25 C = 1.0)
  real(8) :: rdc         ! Scaling factor for high temperature inhibition (25 C = 1.0)

  real(8) :: phi_psii    ! Quantum yield of PS II
  real(8) :: theta_j     ! Empirical curvature parameter for electron transport rate
  !
  ! *** Calculated variables ***
  !
  real(8) :: vcmax       ! Maximum carboxylation rate (umol/m2/s)
  real(8) :: jmax        ! Maximum electron transport rate (umol/m2/s)
  real(8) :: je          ! Electron transport rate (umol/m2/s)
  real(8) :: kc          ! Michaelis-Menten constant for CO2 (umol/mol)
  real(8) :: ko          ! Michaelis-Menten constant for O2 (mmol/mol)
  real(8) :: cp          ! CO2 compensation point (umol/mol)
  real(8) :: ceair       ! Vapor pressure of air, constrained (Pa)
  real(8) :: esat        ! Saturation vapor pressure (Pa)
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: photosynthesis_param  ! Leaf-level parameters for photosynthesis model
  public :: leaf_photosynthesis   ! Leaf photosynthesis and stomatal conductance
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ci_func_gs         ! Function to evaluate Ci, with gs specified as input
  private :: ft                 ! Photosynthesis temperature response
  private :: fth                ! Photosynthesis temperature inhibition
  private :: fth25              ! Scaling factor for photosynthesis temperature inhibition
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine photosynthesis_param ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    t_growth                      , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home                        , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    universal                       &  ! Universal variables
    )
    !
    ! !DESCRIPTION:
    ! Leaf-level parameters for photosynthesis model
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air
    real(8), intent(in) :: t_growth, t_home
    type(universal_type), intent(in) :: universal
    !
    ! !LOCAL VARIABLES:
    real(8) :: sco           ! Relative specificity of rubisco
    !---------------------------------------------------------------------

    associate ( &
    rgas        => universal%rgas      , &  ! Universal gas constant (J/K/mol)
    tfrz        => universal%tfrz        &  ! Freezing point of water (K)
    )

    !---------------------------------------------------------------------
    ! kc, ko, cp at 25C: Bonan et al. (2014)
    ! Derive sco from cp with o2=0.209 mol/mol and re-calculate cp to allow
    ! variation in o2
    !---------------------------------------------------------------------

    kc25 = 404.9d0                 ! umol/mol
    ko25 = 278.4d0                 ! mmol/mol
    cp25 = 42.75d0                 ! umol/mol

    sco = 0.5d0 * 0.209d0 / (cp25 * 1.e-06) ! cp25 (umol/mol) -> (mol/mol)
    cp25 = 0.5d0 * o2air / sco * 1000.0d0   ! O2 is mmol/mol. Multiply by 1000 for umol/mol

    !---------------------------------------------------------------------
    ! Activation energy: Bonan et al. (2014)
    ! Temperature acclimation: Kumarathunge et al. (2019)
    !---------------------------------------------------------------------

    kcha    = 79430.0d0
    koha    = 36380.0d0
    cpha    = 37830.0d0
!    vcmaxha = 65330.0d0
!    jmaxha  = 43540.0d0
    vcmaxha = 108245.0d0     ! <= Aspinwall et al. (2020)
    jmaxha  = 73079.1d0      ! <= Aspinwall et al. (2020)
    rdha    = 46390.0d0

    !---------------------------------------------------------------------
    ! High temperature deactivation: Bonan et al. (2014)
    ! Temperature acclimation: Kumarathunge et al. (2019)
    ! The factor "c" scales the deactivation to a value of 1.0 at 25C
    !---------------------------------------------------------------------

!    vcmaxhd = 150000.0d0
!    jmaxhd  = 150000.0d0
    vcmaxhd = 200000.0d0   ! <= Kumarathunge et al. (2019), Aspinwall et al. (2020)
    jmaxhd  = 200000.0d0   ! <= Kumarathunge et al. (2019), Aspinwall et al. (2020)
    rdhd    = 150000.0d0

!    vcmaxse = 490.0d0
!    jmaxse  = 490.0d0
    vcmaxse = 655.7d0      ! <= Aspinwall et al. (2020)
    jmaxse  = 655.0d0      ! <= Aspinwall et al. (2020)
    rdse    = 490.0d0

    vcmaxc = fth25 (tfrz, rgas, vcmaxhd, vcmaxse)
    jmaxc  = fth25 (tfrz, rgas, jmaxhd, jmaxse)
    rdc    = fth25 (tfrz, rgas, rdhd, rdse)

    !---------------------------------------------------------------------
    ! Miscellaneous parameters
    !---------------------------------------------------------------------

    phi_psii = 0.70d0
    theta_j = 0.90d0

    end associate
  end subroutine photosynthesis_param

  !-----------------------------------------------------------------------
  subroutine leaf_photosynthesis ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
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
    ! Leaf photosynthesis and stomatal conductance
    !
    ! !USES:
    use mod_param
    use mod_math_tools
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    integer, intent(in) :: i
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: qabs                ! PAR utilized by PS II (umol photons/m2/s)
    real(8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(8) :: r1,r2               ! Roots of quadratic equation
    real(8) :: an_err              ! An for error check
    real(8) :: desat               ! Derivative of saturation vapor pressure (Pa/K)
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    rgas       => universal%rgas    , &  ! Universal gas constant (J/K/mol)
    tfrz       => universal%tfrz    , &  ! Freezing point of water (K)
    eair       => atmos%eair        , &  ! Vapor pressure in air (Pa)
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    vcmax25    => layer%vcmax25     , &  ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    jmax25     => layer%jmax25      , &  ! Maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    rd25       => layer%rd25        , &  ! Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
    apar       => flux%apar         , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    tleaf      => flux%tleaf        , &  ! Leaf temperature profile (K)
    gbv        => flux%gbv          , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gbc        => flux%gbc          , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    gs         => flux%gs           , &  ! Leaf stomatal conductance (mol H2O/m2/s)
                                         !
                                         ! *** Output ***
    rd         => flux%rd           , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci         => flux%ci           , &  ! Leaf intercellular CO2 (umol/mol)
    hs         => flux%hs           , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd        => flux%vpd          , &  ! Leaf vapor pressure deficit at surface (Pa)
                                         !
                                         ! *** Output from calls to CiFuncGs ***
    ac         => flux%ac           , &  ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj         => flux%aj           , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ag         => flux%ag           , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an         => flux%an           , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs         => flux%cs             &  ! Leaf surface CO2 (umol/mol)
    )

    if (leaf_layer(i) == 1) then ! Leaf layer

       !------------------------------------------------------------------
       ! Adjust photosynthetic parameters for temperature
       !------------------------------------------------------------------

       kc     = kc25          * ft(tfrz, rgas, tleaf(i), kcha)
       ko     = ko25          * ft(tfrz, rgas, tleaf(i), koha)
       cp     = cp25          * ft(tfrz, rgas, tleaf(i), cpha)
       vcmax  = vcmax25(i) * ft(tfrz, rgas, tleaf(i), vcmaxha) * fth(rgas, tleaf(i), vcmaxhd, vcmaxse, vcmaxc)
       jmax   = jmax25(i)  * ft(tfrz, rgas, tleaf(i), jmaxha)  * fth(rgas, tleaf(i), jmaxhd, jmaxse, jmaxc)
       rd(i)  = rd25(i)    * ft(tfrz, rgas, tleaf(i), rdha)    * fth(rgas, tleaf(i), rdhd, rdse, rdc)

       !------------------------------------------------------------------
       ! Electron transport rate
       !------------------------------------------------------------------

       qabs = 0.5d0 * phi_psii * apar(i)
       aquad = theta_j
       bquad = -(qabs + jmax)
       cquad = qabs * jmax
       call quadratic (aquad, bquad, cquad, r1, r2)
       je = min(r1,r2)

       !------------------------------------------------------------------
       ! Ci calculation
       !------------------------------------------------------------------

       ! Calculate photosynthesis for a specified stomatal conductance

       ci(i) = ci_func_gs (i, o2air, co2air, layer, flux)

       !------------------------------------------------------------------
       ! Make sure iterative solution is correct
       !------------------------------------------------------------------

       if (gs(i) < 0.0d0) then
          write(*,*) 'Error: Leaf photosynthesis: negative stomatal conductance'
       end if

       ! Compare with diffusion equation: An = (ca - ci) * gleaf

       an_err = (co2air - ci(i)) / (1.0d0 / gbc(i) + 1.6d0 / gs(i))
       if (an(i) > 0.0d0 .and. abs(an(i)-an_err) > 0.01d0) then
          write(*,*) 'Error: Leaf photosynthesis: failed diffusion error check'
       end if

       !------------------------------------------------------------------
       ! Saturation vapor pressure at leaf temperature
       !------------------------------------------------------------------

       call sat_vap (tfrz, tleaf(i), esat, desat)

       !------------------------------------------------------------------
       ! Relative humidity and vapor pressure at leaf surface
       !------------------------------------------------------------------

       hs(i) = (gbv(i)*eair + gs(i)*esat) / ((gbv(i)+gs(i))*esat)
       vpd(i) = max(esat - hs(i)*esat, 0.1d0)

       else ! non-leaf layer

       rd(i) = 0.0d0
       ci(i) = ci_func_gs (i, o2air, co2air, layer, flux)
       hs(i) = 0.0d0
       vpd(i) = 0.0d0

    end if

    end associate
  end subroutine leaf_photosynthesis

  !-----------------------------------------------------------------------
  function ci_func_gs (i, o2air, co2air, layer, flux) result(ci_val)
    !
    ! !DESCRIPTION:
    ! Calculate leaf photosynthesis for a specified stomatal conductance.
    ! Then calculate Ci from the diffusion equation.
    !
    ! !USES:
    use mod_param
    use mod_math_tools, only : quadratic
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: i        ! Layer index
    real(8), intent(in) :: o2air    ! Atmospheric O2 at reference height (mmol/mol)
    real(8), intent(in) :: co2air   ! Atmospheric CO2 profile (umol/mol)
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: gleaf               ! Leaf CO2 conductance (mol CO2/m2/s)
    real(8) :: a0,e0,d0            ! Terms for quadratic photosynthesis calculation
    real(8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(8) :: r1,r2               ! Roots of quadratic equation
    real(8) :: ci_val              ! Calculated value for Ci (umol/mol)
    !---------------------------------------------------------------------

    associate ( &
                                        ! *** Input ***
    leaf_layer => layer%leaf_layer , &  ! 1: Leaf, 0: No leaf
    gbc        => flux%gbc         , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    gs         => flux%gs          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rd         => flux%rd          , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
                                        !
                                        ! *** Output ***
    ac         => flux%ac          , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj         => flux%aj          , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ag         => flux%ag          , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an         => flux%an          , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs         => flux%cs            &  ! Leaf surface CO2 (umol/mol)
    )

    !---------------------------------------------------------------------
    ! Calculate leaf photosynthesis for a specified stomatal conductance.
    ! Then calculate Ci from the diffusion equation.
    !
    ! This routine uses a quadratic equation to solve for net photosynthesis (An).
    ! A general equation for C3 photosynthesis is:
    !
    !      a (Ci - Cp)
    ! An = ----------- - Rd
    !        e Ci + d
    !
    ! where:
    !
    ! An = Net leaf photosynthesis (umol CO2/m2/s)
    ! Rd = Leaf respiration (umol CO2/m2/s)
    ! Ci = Intercellular CO2 concentration (umol/mol)
    ! Cp = CO2 compensation point (umol/mol)
    !
    ! Rubisco-limited photosynthesis (Ac)
    ! a  = Vcmax
    ! e  = 1
    ! d  = Kc (1 + Oi/Ko)
    !
    ! RuBP regeneration-limited photosynthesis (Aj)
    ! a = J
    ! e = 4
    ! d = 8 Cp
    !
    ! where:
    !
    ! Vcmax = Maximum carboxylation rate (umol/m2/s)
    ! Kc    = Michaelis-Menten constant for CO2 (umol/mol)
    ! Ko    = Michaelis-Menten constant for O2 (mmol/mol)
    ! Oi    = Intercellular O2 concentration (mmol/mol)
    ! J     = Electron transport rate (umol/m2/s)
    !
    ! Ci is calculated from the diffusion equation:
    !
    !                   1.4   1.6
    ! An = (Ca - Ci) / (--- + ---)
    !                   gb    gs
    !
    !            1.4   1.6
    ! Ci = Ca - (--- + ---) An
    !            gb    gs
    !
    ! where:
    !
    ! Ca  = Atmospheric CO2 concentration (umol/mol)
    ! gb  = Leaf boundary layer conductance (mol H2O/m2/s)
    ! gs  = Leaf stomatal conductance (mol H2O/m2/s)
    ! 1.4 = Corrects gb for the diffusivity of CO2 compared with H2O
    ! 1.6 = Corrects gs for the diffusivity of CO2 compared with H2O
    !
    ! The resulting quadratic equation is: a An**2 + b An + c = 0
    ! Correct solution is the smaller of the two roots.
    !---------------------------------------------------------------------

    if (leaf_layer(i) > 0.0d0) then ! leaf layer

       ! Leaf conductance: gbc has units mol CO2/m2/s, gs has units mol H2O/m2/s,
       ! gleaf has units mol CO2/m2/s

       gleaf = 1.0d0 / (1.0d0/gbc(i) + 1.6d0/gs(i))

       !------------------------------------------------------------------
       ! Gross assimilation rates
       !------------------------------------------------------------------

       ! C3: Rubisco-limited photosynthesis

       a0 = vcmax
       e0 = 1.0d0
       d0 = kc * (1.0d0 + o2air / ko)

       aquad = e0 / gleaf
       bquad = -(e0*co2air + d0) - (a0 - e0*rd(i)) / gleaf
       cquad = a0 * (co2air - cp) - rd(i) * (e0*co2air + d0)

       call quadratic (aquad, bquad, cquad, r1, r2)
       ac(i) = min(r1,r2) + rd(i)

       ! C3: RuBP regeneration-limited photosynthesis

       a0 = je
       e0 = 4.0d0
       d0 = 8.0d0 * cp

       aquad = e0 / gleaf
       bquad = -(e0*co2air + d0) - (a0 - e0*rd(i)) / gleaf
       cquad = a0 * (co2air - cp) - rd(i) * (e0*co2air + d0)

       call quadratic (aquad, bquad, cquad, r1, r2)
       aj(i) = min(r1,r2) + rd(i)

       !------------------------------------------------------------------
       ! Net assimilation as the minimum rate
       !------------------------------------------------------------------

       ag(i) = min(ac(i),aj(i))

       an(i) = ag(i) - rd(i)

       !------------------------------------------------------------------
       ! Leaf surface CO2
       !------------------------------------------------------------------

       cs(i) = co2air - an(i) / gbc(i)

       !------------------------------------------------------------------
       ! Intercelluar CO2
       !------------------------------------------------------------------

       ci_val = co2air - an(i) / gleaf

    else ! non-leaf layer

       ac(i) = 0.0d0
       aj(i) = 0.0d0
       ag(i) = 0.0d0
       an(i) = 0.0d0
       cs(i) = 0.0d0
       ci_val = 0.0d0

    end if

    end associate
  end function ci_func_gs

  !-----------------------------------------------------------------------
  function ft (tfrz, rgas, tl, ha) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature response
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: tfrz     ! Freezing point of water (K)
    real(8), intent(in) :: rgas     ! Universal gas constant (J/K/mol)
    real(8), intent(in) :: tl       ! Leaf temperature (K)
    real(8), intent(in) :: ha       ! Activation energy (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = exp( ha / (rgas * (tfrz+25.0d0)) * (1.0d0 - (tfrz+25.0d0) / tl) )

  end function ft

  !-----------------------------------------------------------------------
  function fth (rgas, tl, hd, se, c) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature inhibition
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: rgas     ! Universal gas constant (J/K/mol)
    real(8), intent(in) :: tl       ! Leaf temperature (K)
    real(8), intent(in) :: hd       ! Deactivation energy (J/mol)
    real(8), intent(in) :: se       ! Entropy term (J/mol/K)
    real(8), intent(in) :: c        ! Scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = c / ( 1.0d0 + exp( (-hd + se*tl) / (rgas*tl) ) )

  end function fth

  !-----------------------------------------------------------------------
  function fth25 (tfrz, rgas, hd, se) result(ans)
    !
    ! !DESCRIPTION:
    ! Scaling factor for photosynthesis temperature inhibition
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: tfrz     ! Freezing point of water (K)
    real(8), intent(in) :: rgas     ! Universal gas constant (J/K/mol)
    real(8), intent(in) :: hd       ! Deactivation energy (J/mol)
    real(8), intent(in) :: se       ! Entropy term (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = 1.0d0 + exp( (-hd + se * (tfrz+25.0d0)) / (rgas * (tfrz+25.0d0)) )

  end function fth25

end module mod_photosynthesis
