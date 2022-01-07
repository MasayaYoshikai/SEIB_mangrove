module mod_param

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Stand-level model data structure
  !
  ! !USES:
  use data_structure, only : Max_hgt
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_allocate
  !
  ! !PUBLIC DATA TYPES:
  !-----------------------------------------------------------------------

  type, public :: universal_type

    ! Universal variables in the model system (constant)

       ! physical constants

       real(8) :: rgas                              ! Universal gas constant (J/K/mol)
       real(8) :: iv                                ! Van't Hoff coefficient for NaCl (-)
       real(8) :: denh2o                            ! Density of liquid water (kg/m3)
       real(8) :: mmh2o                             ! Molecular mass of water (kg/mol)
       real(8) :: grav                              ! Gravitational acceleration (m/s2)
       real(8) :: tfrz                              ! Freezing point of water (K)
       real(8) :: mmdry                             ! Molecular mass of dry air (kg/mol)
       real(8) :: cpd                               ! Specific heat of dry air at constant pressure (J/kg/K)
       real(8) :: cpw                               ! Specific heat of water vapor at constant pressure (J/kg/K)
       real(8) :: visc0                             ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
       real(8) :: Dh0                               ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
       real(8) :: Dv0                               ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
       real(8) :: Dc0                               ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

       ! Miscellaneous variables

       integer :: n_layer                           ! Number of layers (-) (Max_hgt in SEIB-DGVM)

  end type universal_type
  !-----------------------------------------------------------------------
  type, public :: atmos_type

    ! Atmospheric input variables

    real(8) :: tair                                 ! Air temperature (K)
    real(8) :: pair                                 ! Air pressure (Pa)
    real(8) :: wind                                 ! Wind speed (m/s)
    real(8) :: eair                                 ! Vapor pressure in air (Pa)
    real(8) :: rad_dir                              ! Direct radiation at canopy top (W/m2)
    real(8) :: rad_dif                              ! Diffused radiation at canopy top (W/m2)

    ! Atmospheric derived variables

    real(8) :: rhomol                               ! Molar density (mol/m3)
    real(8) :: cp                                   ! Specific heat of air at constant pressure (J/mol/K)

  end type atmos_type
  !-----------------------------------------------------------------------
  type, public :: layer_type

    ! Layer variables

    integer, dimension(Max_hgt) :: leaf_layer    ! 1: Leaf, 0: No leaf
    integer                     :: top_layer     ! Top layer of the canopy
    integer                     :: bottom_layer  ! Bottom layer of the canopy
    integer                     :: skip_layer    ! Layer interval to skip the flux calculation
    real(8), dimension(Max_hgt) :: sumpai        ! Cumulative leaf area (m2 leaf/m2 ground)
    real(8) :: dpai                              ! Layer leaf area index (m2 leaf/m2 ground)
    real(8), dimension(Max_hgt) :: vcmax25       ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    real(8), dimension(Max_hgt) :: jmax25        ! Maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    real(8), dimension(Max_hgt) :: rd25          ! Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)

  end type layer_type
  !-----------------------------------------------------------------------
  type, public :: flux_type

    ! Flux variables

       ! Plant hydraulics related variables

       real(8) :: r_soil                            ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
       real(8) :: r_root                            ! Root-stem resistance (MPa.s.tree/mmol H2O)
       real(8) :: dr_sap                            ! Stem resistance per unit path length (MPa.s.tree/mmol H2O/m)
       real(8) :: r_sap                             ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
       real(8) :: r_whole                           ! Whole plant resistance (MPa.s.tree/mmol H2O)
       real(8) :: sapflow                           ! Stand-level sapflow rate (mol H2O/tree/s)

    real(8), dimension(Max_hgt) :: rn               ! Leaf net radiation profile (W/m2 leaf)
    real(8), dimension(Max_hgt) :: apar             ! Leaf absorbed PAR Profile (umol photon/m2 leaf/s)
    real(8), dimension(Max_hgt) :: tleaf            ! Leaf temperature profile (K)
    real(8), dimension(Max_hgt) :: shleaf           ! Leaf sensible heat flux (W/m2 leaf)
    real(8), dimension(Max_hgt) :: lhleaf           ! Leaf latent heat flux (W/m2 leaf)
    real(8), dimension(Max_hgt) :: etflux           ! Leaf transpiration rate profile (mol H2O/m2 leaf/s)
    real(8), dimension(Max_hgt) :: gbh              ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    real(8), dimension(Max_hgt) :: gbv              ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    real(8), dimension(Max_hgt) :: gbc              ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: gs               ! Leaf stomatal conductance (mol H2O/m2/s)
    real(8), dimension(Max_hgt) :: rd               ! Leaf respiration rate (umol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: ci               ! Leaf intercellular CO2 (umol/mol)
    real(8), dimension(Max_hgt) :: hs               ! Leaf fractional humidity at leaf surface (dimensionless)
    real(8), dimension(Max_hgt) :: vpd              ! Leaf vapor pressure deficit at surface (Pa)
    real(8), dimension(Max_hgt) :: ac               ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: aj               ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: ag               ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: an               ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(8), dimension(Max_hgt) :: cs               ! Leaf surface CO2 (umol/mol)

    real(8), dimension(Max_hgt) :: dark_gs          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_tleaf       ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_shleaf      ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_lhleaf      ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_etflux      ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_rd          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_ci          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_hs          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_vpd         ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_ac          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_aj          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_ag          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_an          ! Flux variables for dark condition
    real(8), dimension(Max_hgt) :: dark_cs          ! Flux variables for dark condition

  end type flux_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_allocate (universal, atmos, layer, flux)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate model data structure
    !
    ! !ARGUMENTS:
    type(universal_type), intent(inout) :: universal
    type(atmos_type)    , intent(inout) :: atmos
    type(layer_type)    , intent(inout) :: layer
    type(flux_type)     , intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    !
    !---------------------------------------------------------------------

    ! Universal variables

       ! Miscellaneous

       universal%n_layer = Max_hgt

       ! Physical constants

       universal%rgas   = 8.31446d0            ! Universal gas constant (J/K/mol)
       universal%iv     = 2.0d0                ! Van't Hoff coefficient for NaCl (-)
       universal%denh2o = 1000.0d0             ! Density of liquid water (kg/m3)
       universal%mmh2o  = 18.02d0 / 1000.0d0   ! Molecular mass of water (kg/mol)is:starred
       universal%grav   = 9.80665d0            ! Gravitational acceleration (m/s2)
       universal%tfrz   = 273.15d0             ! Freezing point of water (K)
       universal%mmdry  = 28.97d0 / 1000.0d0   ! Molecular mass of dry air (kg/mol)
       universal%cpd    = 1005.0d0             ! Specific heat of dry air at constant pressure (J/kg/K)
       universal%cpw    = 1846.0d0             ! Specific heat of water vapor at constant pressure (J/kg/K)
       universal%visc0  = 13.3e-06             ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
       universal%Dh0    = 18.9e-06             ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
       universal%Dv0    = 21.8e-06             ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
       universal%Dc0    = 13.8e-06             ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

    ! Atmospheric input variables

    atmos%tair = 0.0d0
    atmos%pair = 0.0d0
    atmos%wind = 0.0d0
    atmos%eair = 0.0d0
    atmos%rad_dir = 0.0d0
    atmos%rad_dif = 0.0d0

    ! Atmpspheric derived variables

    atmos%rhomol = 0.0d0
    atmos%cp = 0.0d0

    ! Layer derived variables

    layer%leaf_layer(:) = 0
    layer%top_layer = 0
    layer%bottom_layer = 0
    layer%skip_layer = 1
    layer%sumpai(:) = 0.0d0
    layer%dpai = 0.0d0
    layer%vcmax25(:) = 0.0d0
    layer%jmax25(:) = 0.0d0
    layer%rd25(:) = 0.0d0

    ! Flux variables

    flux%r_soil = 0.0d0
    flux%r_root = 0.0d0
    flux%dr_sap = 0.0d0
    flux%r_sap = 0.0d0
    flux%r_whole = 0.0d0
    flux%sapflow = 0.0d0

    flux%rn(:) = 0.0d0
    flux%apar(:) = 0.0d0
    flux%tleaf(:) = 0.0d0
    flux%shleaf(:) = 0.0d0
    flux%lhleaf(:) = 0.0d0
    flux%etflux(:) = 0.0d0
    flux%gbh(:) = 0.0d0
    flux%gbv(:) = 0.0d0
    flux%gbc(:) = 0.0d0
    flux%gs(:) = 0.0d0
    flux%rd(:) = 0.0d0
    flux%ci(:) = 0.0d0
    flux%hs(:) = 0.0d0
    flux%vpd(:) = 0.0d0
    flux%ac(:) = 0.0d0
    flux%aj(:) = 0.0d0
    flux%ag(:) = 0.0d0
    flux%an(:) = 0.0d0
    flux%cs(:) = 0.0d0

    flux%dark_gs(:) = 0.0d0
    flux%dark_tleaf(:) = 0.0d0
    flux%dark_shleaf(:) = 0.0d0
    flux%dark_lhleaf(:) = 0.0d0
    flux%dark_etflux(:) = 0.0d0
    flux%dark_rd(:) = 0.0d0
    flux%dark_ci(:) = 0.0d0
    flux%dark_hs(:) = 0.0d0
    flux%dark_vpd(:) = 0.0d0
    flux%dark_ac(:) = 0.0d0
    flux%dark_aj(:) = 0.0d0
    flux%dark_ag(:) = 0.0d0
    flux%dark_an(:) = 0.0d0
    flux%dark_cs(:) = 0.0d0

  end subroutine init_allocate

end module mod_param
