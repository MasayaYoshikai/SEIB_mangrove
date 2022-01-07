module mod_spac_photosynthesis

  !-----------------------------------------------------------------------
  ! !Description:
  ! SPAC-Photosynthesis model
  !
  ! !Uses:
  use mod_param
  use data_structure,       only : Max_hgt, File_no, Fn_diurnal,               &
                                   growth_calc, display_flag
  use time_counter,         only : year
  use vegi_status_current1, only : monitor, bottom_layer_monitor,              &
                                   top_layer_monitor, leaf_wp_monitor,         &
                                   leaf_wpmin_monitor, gs_top_monitor,         &
                                   an_top_monitor, an_bot_monitor,             &
                                   et_top_monitor, et_bot_monitor,             &
                                   plant_gpp_hour, plant_sap_hour,             &
                                   plant_et_hour, leaf_wp_hour
  !
  ! !Public types:
  implicit none
  !
  ! !Public member functions:
  public :: spac_photosynthesis
  !
  ! !Private member functions:
  !
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine spac_photosynthesis ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    hksat                         , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
    moist                         , &  ! Relative soil moisture (-)
    bsw                           , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat                        , &  ! Soil matric potential at saturation (MPa)
    soil_t                        , &  ! Soil water temperature (K)
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    root_filter                   , &  ! Filtering rate of salt at root surface (-)
    root_resist                   , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    root_density                  , &  ! Specific root density (fine root) (g root/m3 root)
    root_radius                   , &  ! Fine root radius (m)
    root_depth                    , &  ! Rooting depth (m)
    fine_root_ratio               , &  ! Fraction of fine root in below-ground biomass (-)
    k_sap                         , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
    p50_sap                       , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
    a2_sap                        , &  ! Conductivity vulnerability curve coefficient (-)
    wood_rho                      , &  ! Wood density (g/cm3)
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_root                      , &  ! C/N ratio in mol in root
    t_growth                      , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home                        , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    minlp                         , &  ! Minimum leaf water potential (MPa)
    dleaf                         , &  ! Leaf dimension (m)
    leaf_cp                       , &  ! Leaf water capacitance (mmol H2O/m2 leaf/MPa)
    vcmaxpft                      , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
    iota0                         , &  ! Reference marginal water use efficiency in well-watered condition (umol CO2/mol H2O)
    leaf_b0                       , &  ! Sensitivity of marginal water use efficiency to leaf water potential (MPa-1)
    grow_resp                     , &  ! Growth respiration rate (g DW/g DW)
    main_resp_stem                , &  ! Maintenance stem respiration rate at 15 degree (day-1)
    main_resp_root                , &  ! Maintenance stem respiration rate at 15 degree (day-1)
    coarse_root_turn              , &  ! Coarse root turnover rate (day-1)
    root_turn                     , &  ! Root turnover rate (day-1)
    leaf_turn                     , &  ! Leaf turnover rate (day-1)
    leaf_resorp                   , &  ! Nutrient resorption from Leaf (fraction)
    pr_s_a                        , &  ! Slope for the scaling factor of prop root system
    pr_s_b                        , &  ! Intercept for the scaling factor of prop root system
    pr_h_a                        , &  ! Slope for the maximum root height of prop root system
    pr_h_b                        , &  ! Intercept for the maximum root height of prop root system
    pr_d                          , &  ! Mean prop root diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    ALM5                          , &  ! Sapwood diameter proportion (m sapwood/m dbh)
    day_of_year                   , &  ! Day of the year (1 - 365)
    no                            , &  ! Tree index
    p                             , &  ! Species index (1: Rh, 2: Br)
    rPAR_dir                      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
    rPAR_dif                      , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
    rNIR_dir                      , &  ! Profile of relative intensity of direct NIR within canopy compared to canopy top (fraction)
    rNIR_dif                      , &  ! Profile of relative intensity of diffused NIR within canopy compared to canopy top (fraction)
    irleaf                        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
    kbm_vis                       , &  ! Scattering adjusted light extinction coefficient for VIS (-)
    kbm_nir                       , &  ! Scattering adjusted light extinction coefficient for NIR (-)
    rwind_profile                 , &  ! Relative wind speed profile to the canopy top (-)
    tree_h                        , &  ! Tree height (m)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    croot_biomass                 , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    lai                           , &  ! Leaf area index (m2 leaf/m2 ground)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g/tree)
    aroot_biomass                 , &  ! Stand-level above-ground root biomass (g/tree)
    seib_height                   , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
    time_series_tair              , &  ! Time-series air temperature (K)
    time_series_pair              , &  ! Time-series air pressure (Pa)
    time_series_wind              , &  ! Time-series wind speed (m/s)
    time_series_eair              , &  ! Time-series vapor pressure in air (Pa)
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** In/Output ***
    leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
    leaf_mdwp                     , &  ! Midday leaf water potential of the day (MPa)
    leaf_pdwp                     , &  ! Pre-dawn leaf water potential of the last day (MPa)
                                       !
                                       ! *** Output ***
    n_leaf_layer                  , &  ! Number of leaf layers
    gs_reg_flag                   , &  ! Flag of stomata regulation during the day (0: no, 1: yes)
    sapflow_day                   , &  ! Daily sapflow rate at stand-level (m3 H2O/tree/day)
    sh_day                        , &  ! Daily sensible heat flux at stand-level (MJ/tree/day)
    lh_day                        , &  ! Daily latent heat flux at stand-level (MJ/tree/day)
    leaf_wp_day_min               , &  ! Daily minimum leaf water potential (MPa)
    an_mean_day_max               , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_top_day_max                , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
    gpp_day                       , &  ! Daily gross photosynthesis at stand-level (mol C/tree/day)
    c_uptake_day                  , &  ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)
    c_uptake_bottom_day           , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day           , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    r_whole_root_increment        , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
    r_whole_d_increment           , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
    r_root_out                    , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    r_sap_out                     , &  ! Daily maximum of whole-plant stem resistance (MPa.s.tree/mmol H2O)
    tleaf_out                       &  ! Leaf temperature profile (K)
    )
    !
    ! !Description:
    ! Compute transpiration and photosynthesis for a tree
    !
    ! !Uses:
    use mod_water_vapor
    use mod_math_tools
    use mod_soil_char
    use mod_plant_hydraulics
    use mod_nitrogen_profile
    use mod_photosynthesis
    use mod_leaf_boundary_layer
    use mod_leaf_temperature
    use mod_leaf_water_potential
    use mod_stomatal_conductance
    use mod_tree_allometry
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: o2air, co2air
    real(8), intent(in)    :: hksat, moist, bsw, psisat, soil_t
    real(8), intent(in)    :: sal, din
    real(8), intent(in)    :: root_filter(:), root_resist(:)
    real(8), intent(in)    :: root_density(:), root_radius(:), root_depth
    real(8), intent(in)    :: fine_root_ratio(:)
    real(8), intent(in)    :: k_sap(:), p50_sap(:), a2_sap(:)
    real(8), intent(in)    :: wood_rho(:), c_n_leaf(:), c_n_stem(:), c_n_root(:)
    real(8), intent(in)    :: t_growth(:), t_home(:)
    real(8), intent(in)    :: minlp(:), dleaf(:), leaf_cp(:), vcmaxpft(:)
    real(8), intent(in)    :: iota0(:), leaf_b0(:)
    real(8), intent(in)    :: grow_resp(:), main_resp_stem(:), main_resp_root(:)
    real(8), intent(in)    :: coarse_root_turn(:), root_turn(:), leaf_turn(:)
    real(8), intent(in)    :: leaf_resorp(:), pr_s_a(:), pr_s_b(:)
    real(8), intent(in)    :: pr_h_a(:), pr_h_b(:), pr_d(:)
    real(8), intent(in)    :: C_in_drymass, ALM5(:)
    integer, intent(in)    :: day_of_year, no, p
    real(8), intent(in)    :: rPAR_dir(:), rPAR_dif(:), rNIR_dir(:)
    real(8), intent(in)    :: rNIR_dif(:), irleaf(:)
    real(8), intent(in)    :: kbm_vis(:), kbm_nir(:), rwind_profile(:)
    real(8), intent(in)    :: tree_h, root_biomass, croot_biomass
    real(8), intent(in)    :: leaf_biomass, lai, crown_area
    real(8), intent(in)    :: dbh_heartwood, dbh_sapwood
    real(8), intent(in)    :: trunk_biomass, aroot_biomass
    integer, intent(in)    :: seib_height, seib_bole
    real(8), intent(in)    :: time_series_tair(:), time_series_pair(:)
    real(8), intent(in)    :: time_series_wind(:), time_series_eair(:)
    real(8), intent(in)    :: time_series_rad_dir(:), time_series_rad_dif(:)
    real(8), intent(inout) :: leaf_wp, leaf_mdwp, leaf_pdwp
    integer, intent(out)   :: n_leaf_layer, gs_reg_flag
    real(8), intent(out)   :: sapflow_day, sh_day, lh_day
    real(8), intent(out)   :: leaf_wp_day_min, an_mean_day_max, an_top_day_max
    real(8), intent(out)   :: gpp_day, c_uptake_day
    real(8), intent(out)   :: c_uptake_bottom_day, n_uptake_bottom_day
    real(8), intent(out)   :: r_whole_root_increment, r_whole_d_increment
    real(8), intent(out)   :: r_root_out, r_sap_out
    real(8), intent(out)   :: tleaf_out(:)
    !
    ! !Local variables:
    type(universal_type) :: universal
    type(atmos_type)     :: atmos
    type(layer_type)     :: layer
    type(flux_type)      :: flux
    integer :: i
    integer :: count
    integer :: hour
    integer :: canopy_index
    integer :: hour_of_year
    real(8) :: root_area
    real(8) :: tot_psi
    real(8) :: hk
    real(8) :: swsky_vis_dir
    real(8) :: swsky_vis_dif
    real(8) :: swsky_nir_dir
    real(8) :: swsky_nir_dif
    real(8) :: par_dir
    real(8) :: par_dif
    real(8) :: swleaf
    real(8) :: qair
    real(8) :: rhoair
    real(8) :: mmair
    real(8) :: leaf_wp_new
    real(8) :: sapflow_m3
    real(8) :: stand_et
    real(8) :: stand_sh
    real(8) :: stand_lh
    real(8) :: stand_gpp
    real(8) :: an_mean
    real(8) :: an_top
    real(8) :: c_uptake
    real(8) :: c_uptake_bottom
    real(8) :: et_bottom
    real(8) :: n_uptake_day
    real(8) :: above_root_ratio
    real(8) :: d_trunk_dummy, d_aroot_dummy
    real(8) :: below_resource_day_stem
    real(8) :: below_resource_day_root
    real(8) :: new_dbh_heartwood
    real(8) :: new_dbh_sapwood
    real(8) :: new_tree_h
    real(8) :: new_root_biomass
    real(8) :: term1, term2, term3
    real(8) :: term4, term5, term6, term7
    real(8) :: term8, term9, term10, term11
    real(8) :: term12,term13, term14, term15
    real(8) :: term16
    !---------------------------------------------------------------------

    ! Allocate model parameters

    call init_allocate ( &
                              ! *** In/Output ***
    universal            , &  ! Universal variables
    atmos                , &  ! Atmospheric variables
    layer                , &  ! Layer variables
    flux                   &  ! Flux variables
    )

    associate ( &
                                                     ! *** Input ***
    n_layer      => universal%n_layer           , &  ! Number of layers (-) (Max_hgt in SEIB-DGVM)
    rgas         => universal%rgas              , &  ! Universal gas constant (J/K/mol)
    denh2o       => universal%denh2o            , &  ! Water density (kg/m3)
    mmh2o        => universal%mmh2o             , &  ! Molecular mass of water (kg/mol)
    tfrz         => universal%tfrz              , &  ! Freezing point of water (K)
    mmdry        => universal%mmdry             , &  ! Molecular mass of dry air (kg/mol)
    cpd          => universal%cpd               , &  ! Specific heat of dry air at constant pressure (J/kg/K)
    cpw          => universal%cpw               , &  ! Specific heat of water vapor at constant pressure (J/kg/K)
                                                     !
                                                     ! *** Derived variables
    tair         => atmos%tair                  , &  ! Air temperature (K)
    pair         => atmos%pair                  , &  ! Air pressure (Pa)
    wind         => atmos%wind                  , &  ! Wind speed (m/s)
    eair         => atmos%eair                  , &  ! Vapor pressure in air (Pa)
    rhomol       => atmos%rhomol                , &  ! Molar density (mol/m3)
    cp           => atmos%cp                    , &  ! Specific heat of air at constant pressure (J/mol/K)
    rad_dir      => atmos%rad_dir               , &  ! Direct radiation at canopy top (W/m2)
    rad_dif      => atmos%rad_dif               , &  ! Diffused radiation at canopy top (W/m2)
    leaf_layer   => layer%leaf_layer            , &  ! 1: Leaf, 0: No leaf
    top_layer    => layer%top_layer             , &  ! Top layer of the canopy
    bottom_layer => layer%bottom_layer          , &  ! Bottom layer of the canopy
    skip_layer   => layer%skip_layer            , &  ! Layer interval to skip the flux calculation
    sumpai       => layer%sumpai                , &  ! Cumulative leaf area (m2 leaf/m2 ground)
    dpai         => layer%dpai                  , &  ! Layer leaf area index (m2 leaf/m2 ground)
    rn           => flux%rn                     , &  ! Leaf net radiation profile (W/m2 leaf)
    apar         => flux%apar                   , &  ! Leaf absorbed PAR Profile (umol photon/m2 leaf/s)
    etflux       => flux%etflux                 , &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    tleaf        => flux%tleaf                  , &  ! Leaf temperature (K)
    shleaf       => flux%shleaf                 , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf       => flux%lhleaf                 , &  ! Leaf latent heat flux (W/m2 leaf)
    rd           => flux%rd                     , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci           => flux%ci                     , &  ! Leaf intercellular CO2 (umol/mol)
    hs           => flux%hs                     , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd          => flux%vpd                    , &  ! Leaf vapor pressure deficit at surface (Pa)
    ac           => flux%ac                     , &  ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj           => flux%aj                     , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ag           => flux%ag                     , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an           => flux%an                     , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs           => flux%cs                     , &  ! Leaf surface CO2 (umol/mol)
    gs           => flux%gs                     , &  ! Leaf stomatal conductance (mol H2O/m2/s)
    r_soil       => flux%r_soil                 , &  ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
    r_root       => flux%r_root                 , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    r_sap        => flux%r_sap                  , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    r_whole      => flux%r_whole                , &  ! Whole plant resistance (MPa.s.tree/mmol H2O)
    sapflow      => flux%sapflow                , &  ! Stand-level sapflow rate (mol H2O/tree/s)
    dark_gs      => flux%dark_gs                , &  ! Flux variables for dark condition
    dark_tleaf   => flux%dark_tleaf             , &  ! Flux variables for dark condition
    dark_shleaf  => flux%dark_shleaf            , &  ! Flux variables for dark condition
    dark_lhleaf  => flux%dark_lhleaf            , &  ! Flux variables for dark condition
    dark_etflux  => flux%dark_etflux            , &  ! Flux variables for dark condition
    dark_rd      => flux% dark_rd               , &  ! Flux variables for dark condition
    dark_ci      => flux%dark_ci                , &  ! Flux variables for dark condition
    dark_hs      => flux%dark_hs                , &  ! Flux variables for dark condition
    dark_vpd     => flux%dark_vpd               , &  ! Flux variables for dark condition
    dark_ac      => flux%dark_ac                , &  ! Flux variables for dark condition
    dark_aj      => flux%dark_aj                , &  ! Flux variables for dark condition
    dark_ag      => flux%dark_ag                , &  ! Flux variables for dark condition
    dark_an      => flux%dark_an                , &  ! Flux variables for dark condition
    dark_cs      => flux%dark_cs                  &  ! Flux variables for dark condition
    )

    ! Layer interval to skip the flux calculation

    skip_layer = 1

    ! Stand-level root coverage area (m2 ground/tree)
    ! assuming the root area is same as the crown area

    root_area = crown_area

    ! Leaf layer index (1: Leaf, 0: No leaf)

    n_leaf_layer = 0
    do i = 1, n_layer
       if (i>=(seib_bole + 1) .and. i<=seib_height) then
          leaf_layer(i) = 1
          n_leaf_layer = n_leaf_layer+1
       else
          leaf_layer(i) = 0
       end if
    end do

    ! Identifying canopy top layer

    canopy_index = 0
    do i = n_layer, 1, -1
       if (leaf_layer(i)==1 .and. canopy_index==0) then
          top_layer = i
          canopy_index = 1
       end if
    end do

    ! Identifying canopy bottom layer

    canopy_index = 0
    do i = 1, n_layer
       if (leaf_layer(i)==1 .and. canopy_index==0) then
          bottom_layer = i
          canopy_index = 1
       end if
    end do

    ! Layer leaf area index (m2 leaf/m2 ground)

    dpai = lai/real(n_leaf_layer)

    ! Cumulative leaf area from canopy top (m2 leaf/m2 ground)

    count = 0
    do i = n_layer, 1, -1
       if (leaf_layer(i) == 1) then
          count = count+1
          if (count == 1) then
             sumpai(i) = dpai
          else
             sumpai(i) = sumpai(i+1)+dpai
          end if
       else
          sumpai(i) = 0.0d0
       end if
    end do

    ! Set soil hydraulic parameters

    call soil_char (       &
                              ! *** Input ***
    hksat                , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
    moist                , &  ! Relative soil moisture (-)
    bsw                  , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat               , &  ! Soil matric potential at saturation (MPa)
    soil_t               , &  ! Soil water temperature (K)
    sal                  , &  ! Pore-water salinity (mol/m3)
    root_filter(p)       , &  ! Filtering rate of salt at root surface (-)
    universal            , &  ! Universal variables
                              !
                              ! *** Output ***
    tot_psi              , &  ! Total soil water potential (MPa)
    hk                     &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    )

    ! Calculate soil-to-root and root-to-stem hydraulic resistance

    call soil_root_hydraulics ( &
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    root_biomass         , &  ! Stand-level fine root biomass (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )

    ! Canopy profile of nitrogen and photosynthetic capacity

    call nitrogen_profile (&
                              ! *** Input ***
    t_growth(p)          , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home(p)            , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    vcmaxpft(p)          , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    layer                  &  ! Layer variables
    )

    ! Set leaf-level parameters for photosynthesis model

    call photosynthesis_param (&
                              ! *** Input ***
    o2air                , &  ! Atmospheric O2 (mmol/mol)
    t_growth(p)          , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home(p)            , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    universal              &  ! Universal variables
    )

    gs_reg_flag = 0
    sapflow_day = 0.0d0
    sh_day = 0.0d0
    lh_day = 0.0d0
    leaf_wp_day_min = -0.01d0
    an_mean_day_max = 0.0d0
    an_top_day_max = 0.0d0
    gpp_day = 0.0d0
    c_uptake_day = 0.0d0
    c_uptake_bottom_day = 0.0d0
    n_uptake_bottom_day = 0.0d0

    !---------------------------------------------------------------------
    ! Hourly computation
    !---------------------------------------------------------------------

    do hour = 0, 23

       hour_of_year = (day_of_year-1)*24+hour+1
       tair = time_series_tair(hour_of_year)          ! Air temperature (K)
       pair = time_series_pair(hour_of_year)          ! Air pressure (Pa)
       wind = time_series_wind(hour_of_year)          ! Wind speed (m/s)
       eair = time_series_eair(hour_of_year)          ! Vapor pressure in air (Pa)
       rad_dir = time_series_rad_dir(hour_of_year)    ! Direct radiation at canopy top (W/m2)
       rad_dif = time_series_rad_dif(hour_of_year)    ! Diffused radiation at canopy top (W/m2)

       ! Atmospheric VIS and NIR solar radiation (W/m2)
       ! Assuming 4% of solar radiation is UV.

       swsky_vis_dir = rad_dir*0.43d0
       swsky_vis_dif = rad_dif*0.57d0
       swsky_nir_dir = rad_dir*(1.0d0-0.43d0-0.04d0)
       swsky_nir_dif = rad_dif*(1.0d0-0.57d0-0.04d0)

       ! Direct and diffused PAR at canopy top (umol photon/m2 ground/s)

       par_dir = 4.6d0*swsky_vis_dir
       par_dif = 4.2d0*swsky_vis_dif

       qair = mmh2o/mmdry*eair/(pair-(1.0d0-mmh2o/mmdry)*eair)       ! Specific humidity (kg/kg)
       rhomol = pair/(rgas*tair)                                     ! Molar density (mol/m3)
       rhoair = rhomol*mmdry*(1.0d0-(1.0d0-mmh2o/mmdry)*eair/pair)   ! Air density (kg/m3)
       mmair = rhoair/rhomol                                         ! Molecular mass of air (kg/mol)
       cp = cpd*(1.0d0+(cpw/cpd-1.0d0)*qair)*mmair;                  ! Specific heat of air at constant pressure (J/mol/K)

       ! Leaf net radiation (W/m2 leaf) and absorbed PAR (umol photon/m2 leaf/s)

       do i = n_layer, 1, -1
          if (leaf_layer(i) == 1) then
             swleaf = (swsky_vis_dir*rPAR_dir(i)+swsky_vis_dif*rPAR_dif(i))    &
                      *kbm_vis(p)+(swsky_nir_dir*rNIR_dir(i)+swsky_nir_dif     &
                      *rNIR_dif(i))*kbm_nir(p)
             rn(i) = swleaf+irleaf(i)
             apar(i) = (par_dir*rPAR_dir(i)+par_dif*rPAR_dif(i))*kbm_vis(p)
          else
             rn(i) = 0.0d0
             apar(i) = 0.0d0
          end if
       end do

       !---------------------------------------------------------------------
       ! Plant hydraulic resistance and plant water uptake rate
       !---------------------------------------------------------------------

       if (leaf_wp > -0.0001d0) then ! Initial value of leaf_wp after the establishment
          leaf_wp = tot_psi - 0.1d0
       end if
       if (leaf_pdwp > -0.0001d0) then ! Initial value of leaf_pdwp after the establishment
          leaf_pdwp = tot_psi - 0.1d0
       end if

       ! Daily update of pre-dawn and midday leaf water potential (MPa)

       if ((rad_dir+rad_dif)<5.0d0 .and. hour_of_year<((365-1)*24+23+1)) then
          if (time_series_rad_dir(hour_of_year+1)                              &
             +time_series_rad_dif(hour_of_year+1)>=5.0d0) then
             leaf_pdwp = leaf_wp
          end if
       end if
       if (hour == 12) then
          leaf_mdwp = leaf_wp
       end if

       ! Calculate plant hydraulic resistance and water uptake rate

       call plant_hydraulics ( &
                                 ! *** Input ***
       dbh_heartwood        , &  ! Heartwood diameter (m)
       dbh_sapwood          , &  ! Sapwood diameter (m)
       k_sap(p)             , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
       p50_sap(p)           , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
       a2_sap(p)            , &  ! Conductivity vulnerability curve coefficient (-)
       tree_h               , &  ! Tree height (m)
       leaf_wp              , &  ! Leaf water potential at current time step (MPa)
       tot_psi              , &  ! Total soil water potential (MPa)
       universal            , &  ! Universal variables
                                 !
                                 ! *** Input/Output ***
       flux                   &  ! Flux variables
       )
       r_root_out = r_soil+r_root
       if (hour == 0) then
          r_sap_out = r_sap
       else
          r_sap_out = max(r_sap_out, r_sap)
       end if

       !---------------------------------------------------------------------
       ! Leaf boundary layer conductances
       !---------------------------------------------------------------------

       call leaf_boundary_layer (      &
                                          ! *** Input ***
       dleaf(p)                      , &  ! Leaf dimension (m)
       rwind_profile                 , &  ! Relative wind speed profile to the canopy top (-)
       universal                     , &  ! Universal variables
       atmos                         , &  ! Atmospheric variables
       layer                         , &  ! Layer variables
                                          !
                                          ! *** Input/Output ***
       flux                            &  ! Flux variables
       )

       !---------------------------------------------------------------------
       ! Leaf fluxes for dark condition.
       ! gs is regulated to its minimum value.
       !---------------------------------------------------------------------

       if (hour == 0) then
          call dark_condition ( &
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
       end if

       !---------------------------------------------------------------------
       ! Leaf temerature, energy fluxes, photosynthesis, and
       ! stomatal conductance using water-use efficiency optimization
       ! for each canopy layer.
       !---------------------------------------------------------------------

       if ((rad_dir+rad_dif) < 5.0d0) then
          gs     = dark_gs
          tleaf  = dark_tleaf
          shleaf = dark_shleaf
          lhleaf = dark_lhleaf
          etflux = dark_etflux
          rd     = dark_rd
          ci     = dark_ci
          hs     = dark_hs
          vpd    = dark_vpd
          ac     = dark_ac
          aj     = dark_aj
          ag     = dark_ag
          an     = dark_an
          cs     = dark_cs
       else
          call stomatal_conductance (      &
                                             ! *** Input ***
          o2air                         , &  ! Atmospheric O2 (mmol/mol)
          co2air                        , &  ! Atmospheric CO2 (umol/mol)
          iota0(p)                      , &  ! Reference marginal water use efficiency in well-watered condition (umol CO2/mol H2O)
          leaf_b0(p)                    , &  ! Sensitivity of marginal water use efficiency to leaf water potential (MPa-1)
          leaf_pdwp                     , &  ! Pre-dawn leaf water potential of the last day (MPa)
          leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
          universal                     , &  ! Universal variables
          atmos                         , &  ! Atmospheric variables
          layer                         , &  ! Layer variables
                                             !
                                             ! *** Input/Output ***
          flux                            &  ! Flux variables
          )
       end if

       !---------------------------------------------------------------------
       ! Calculate new leaf water potential based on the water balance
       ! detemined from the optimized gs.
       ! It is assumed that the leaf water potential is uniform over the canopy.
       !---------------------------------------------------------------------

       call leaf_water_potential (      &
                                          ! *** Input ***
       leaf_cp(p)                    , &  ! Leaf water capacitance (mmol H2O/m2 leaf/MPa)
       tree_h                        , &  ! Tree height (m)
       crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
       leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
       tot_psi                       , &  ! Total soil water potential (MPa)
       universal                     , &  ! Universal variables
       layer                         , &  ! Layer variables
       flux                          , &  ! Flux variables
                                          !
                                          ! *** Output ***
       leaf_wp_new                     &  ! New leaf water potential (MPa)
       )

       !---------------------------------------------------------------------
       ! Cavitation check.
       ! When the leaf water potential could fall below the acceptable
       ! value (minimum leaf water potential), stomata will close to
       ! reduce the risk of embolism.
       !---------------------------------------------------------------------

       if (leaf_wp_new > minlp(p)) then
          leaf_wp = leaf_wp_new
       else
          gs_reg_flag = 1   ! Flag of stomata regulation (0: no, 1: yes)
          do i = 1, n_layer
             if (leaf_layer(i) == 1) then

                ! Set gs to minimum conductance

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
                gs(i) = 0.0d0
             end if
          end do

          ! Leaf water potential with the regulated gs.

          call leaf_water_potential (      &
                                             ! *** Input ***
          leaf_cp(p)                    , &  ! Leaf water capacitance (mmol H2O/m2 leaf/MPa)
          tree_h                        , &  ! Tree height (m)
          crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
          leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
          tot_psi                       , &  ! Total soil water potential (MPa)
          universal                     , &  ! Universal variables
          layer                         , &  ! Layer variables
          flux                          , &  ! Flux variables
                                             !
                                             ! *** Output ***
          leaf_wp_new                     &  ! New leaf water potential (MPa)
          )
          leaf_wp = leaf_wp_new
       end if

       leaf_wp_day_min = min(leaf_wp_day_min, leaf_wp)
       if (leaf_wp < minlp(p)) then
          write(*,*) 'Error: Leaf water potential below the minimum value'
       end if
       if (leaf_wp > 0.0d0) then
          write(*,*) 'Error: Positive leaf water potential value.'
       end if

       ! Stand-level transpiration rate (mol H2O/tree/s),
       ! sensible heat and latent heat fluxes (W/tree)

       stand_et = 0.0d0
       stand_sh = 0.0d0
       stand_lh = 0.0d0
       do i = 1, n_layer
          if (leaf_layer(i) == 1) then
             stand_et = stand_et+etflux(i)*dpai*crown_area
             stand_sh = stand_sh+shleaf(i)*dpai*crown_area
             stand_lh = stand_lh+lhleaf(i)*dpai*crown_area
          end if
       end do

       ! Stand-level sapflow rate

       sapflow_m3 = stand_et*mmh2o/denh2o  ! (m3 H2O/tree/s): This is based on transpiration.
       ! sapflow_m3 = sapflow*mmh2o/denh2o     ! (m3 H2O/tree/s): This is based on the actual sapflow.

       ! Bottom layer transpiration rate (m3 H2O/m2 leaf/s)

       et_bottom = etflux(bottom_layer)*mmh2o/denh2o

       ! Carbon uptake rate by photosynthesis

       stand_gpp = 0.0d0
       c_uptake = 0.0d0
       an_mean = 0.0d0
       do i = 1, n_layer
          if (leaf_layer(i) == 1) then
             stand_gpp = stand_gpp+ag(i)*dpai*crown_area*1.e-06
             c_uptake = c_uptake+an(i)*dpai*crown_area*1.e-06
             an_mean = an_mean+an(i)
          end if
       end do
       an_mean = an_mean/real(n_leaf_layer)
       an_top = an(bottom_layer+n_leaf_layer-1)

       ! Daily maximum canopy mean and canopy top
       ! leaf net photosynthesis (umol CO2/m2 leaf/s)

       an_mean_day_max = max(an_mean_day_max, an_mean)
       an_top_day_max = max(an_top_day_max, an_top)

       ! For bottom layer (mol CO2/m2 leaf/s)

       c_uptake_bottom = an(bottom_layer)*1.e-06

       !---------------------------------------------------------------------
       ! Daily sum of fluxes
       !---------------------------------------------------------------------

       ! Daily sapflow rate at stand-level (m3 H2O/tree/day)

       sapflow_day = sapflow_day+sapflow_m3*60.0d0*60.0d0

       ! Daily sensible and latent heat flux at stand-level (MJ/tree/day)

       sh_day = sh_day+stand_sh*60.0d0*60.0d0*1.e-06
       lh_day = lh_day+stand_lh*60.0d0*60.0d0*1.e-06

       ! Daily gross photosynthesis at stand-level (mol C/tree/day)

       gpp_day = gpp_day+stand_gpp*60.0d0*60.0d0

       ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)

       c_uptake_day = c_uptake_day+c_uptake*60.0d0*60.0d0

       ! Daily bottom layer carbon uptake rate
       ! by photosynthesis (mol CO2/m2 leaf/day)

       c_uptake_bottom_day = c_uptake_bottom_day+c_uptake_bottom*60.0d0*60.0d0

       ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)

       n_uptake_bottom_day = n_uptake_bottom_day+et_bottom*din*60.0d0*60.0d0

       ! Stand-level diurnal variables

       plant_gpp_hour(no, hour+1) = stand_gpp*60.0d0*60.0d0
       plant_sap_hour(no, hour+1) = sapflow*(mmh2o/denh2o)*60.0d0*60.0d0
       plant_et_hour(no, hour+1) = stand_et*(mmh2o/denh2o)*60.0d0*60.0d0
       leaf_wp_hour(no, hour+1) = leaf_wp

       if (no == monitor) then
          if (hour == 12) then
             bottom_layer_monitor = bottom_layer
             top_layer_monitor = bottom_layer+n_leaf_layer-1
             leaf_wp_monitor = leaf_wp
             leaf_wpmin_monitor = leaf_wp_day_min
             gs_top_monitor = gs(top_layer_monitor)
             an_top_monitor = an(top_layer_monitor)
             an_bot_monitor = an(bottom_layer_monitor)
             et_top_monitor = etflux(top_layer_monitor)
             et_bot_monitor = etflux(bottom_layer_monitor)
          end if
       end if

       if (growth_calc .eqv. .false.) then
          if (no==monitor .and. year==3) then
             write(File_no(24), '(1(i4,a), 12(f12.5,a))') &
             hour_of_year                            , ',', &
             par_dir+par_dif                         , ',', &
             tair-tfrz                               , ',', &
             eair                                    , ',', &
             leaf_pdwp                               , ',', &
             leaf_wp                                 , ',', &
             (1.0d0/r_whole)*mmh2o                   , ',', &  ! (g H2O/s/tree/MPa)
             sapflow*mmh2o*1000.0d0                  , ',', &  ! (g H2O/tree/s)
             etflux(top_layer_monitor)               , ',', &
             vpd(top_layer_monitor)                  , ',', &
             tleaf(top_layer_monitor)-tfrz           , ',', &
             gs(top_layer_monitor)                   , ',', &
             an(top_layer_monitor)
          end if
       end if

    end do

    !---------------------------------------------------------------------
    ! Optimization of plant hydraulics
    !---------------------------------------------------------------------

    ! Required fraction of above-ground root biomass in above-ground

    d_trunk_dummy = 0.0d0
    call proot_allometry ( &
                                         ! *** Input ***
    wood_rho                        , &  ! Wood density (g/cm3)
    pr_s_a                          , &  ! Slope for the scaling factor of prop root system
    pr_s_b                          , &  ! Intercept for the scaling factor of prop root system
    pr_h_a                          , &  ! Slope for the maximum root height of prop root system
    pr_h_b                          , &  ! Intercept for the maximum root height of prop root system
    pr_d                            , &  ! Mean prop root diameter (m)
    p                               , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                   , &  ! Heartwood diameter (m)
    dbh_sapwood                     , &  ! Sapwood diameter (m)
    trunk_biomass                   , &  ! Stand-level trunk biomass (g trunk/tree)
    aroot_biomass                   , &  ! Stand-level above-ground root biomass (g/tree)
                                         !
                                         ! *** In/Output ***
    d_trunk_dummy                   , &  ! Dummy variable, dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                         !
                                         ! *** Output ***
    d_aroot_dummy                   , &  ! Dummy variable, dAGR/dt: Daily above-ground root biomass change (g/tree/day)
    above_root_ratio                  &  ! Required fraction of above-ground root biomass in above-ground
    )

    ! Daily nitrogen uptake rate (mol N/tree/day)

    n_uptake_day = sapflow_day * din
    ! n_uptake_day = 0.10d0

    ! Net biomass increment of stem and fine root from the daily
    ! nitrogen or carbon uptake rate (g DW/tree/day)

    tair = time_series_tair((day_of_year-1)*24+12+1)-tfrz
    term1 = fine_root_ratio(p)*(1.0d0-coarse_root_turn(p))                       ! Term for fine root fraction in below-ground
    term2 = (1.0d0-root_turn(p))*(1.0d0-fine_root_ratio(p))                      ! Term for fine root fraction in below-ground
    term3 = (12.0d0/C_in_drymass)*(1.0d0-root_turn(p))                           ! Accounting for fine root turnover
    term4 = root_biomass*root_turn(p)+croot_biomass*coarse_root_turn(p)          ! Fine and coarse root biomass to be lost by turnover
    term5 = leaf_resorp(p)*leaf_turn(p)*leaf_biomass                           &
            *C_in_drymass/12.0d0/c_n_leaf(p)                                     ! Nutrient resorption from leaves to fall
    term6 = (root_biomass*root_turn(p)/c_n_root(p)                             &
            +croot_biomass*coarse_root_turn(p)/c_n_stem(p))*C_in_drymass/12.0d0  ! Nitrogen demand for fine and coarse root turnover
    term7 = leaf_biomass*leaf_turn(p)*C_in_drymass/12.0d0/c_n_leaf(p)            ! Nitrogen demand for leaf turnover
    term8 = 2.0d0*exp(-0.009d0*(tair-15.0d0))                                    ! Q10
    term9 = exp((tair-15.0d0)*log(term8)/10.0d0)                                 ! Temperature sensitivity factor for respiration
    term10 = (main_resp_stem(p)*(trunk_biomass+aroot_biomass+croot_biomass)    &
             +main_resp_root(p)*root_biomass)*term9*C_in_drymass/12.0d0          ! Carbon demand for maintenance respiration
    term11 = 1.0d0/(1.0d0-grow_resp(p))                                          ! Term for growth respiration
    term12 = (term4+leaf_biomass*leaf_turn(p))*term11*C_in_drymass/12.0d0        ! Carbon demand for turnover
    term13 = (n_uptake_day+term5-term6-term7)*(c_n_stem(p)*term1               &
             /(c_n_stem(p)*term1+c_n_root(p)*term2))*c_n_root(p)*term3           ! Net root biomass increment from nitrogen
    term14 = (c_uptake_day-term10-term12)*(1.0d0-grow_resp(p))                 &
             *(term1/(term1+term2))*term3                                        ! Net root biomass increment from carbon
    term15 = (n_uptake_day+term5-term6-term7)*(1.0d0-above_root_ratio)         &
             *c_n_stem(p)*12.0d0/C_in_drymass                                    ! Net stem biomass increment from nitrogen
    term16 = (c_uptake_day-term10-term12)*(1.0d0-grow_resp(p))                 &
             *(1.0d0-above_root_ratio)*12.0d0/C_in_drymass                       ! Net stem biomass increment from carbon
    below_resource_day_root = min(term13, term14)
    below_resource_day_stem = min(term15, term16)

    ! New DBH after stem biomass increment

    call tree_allometry (  &
                                 ! *** Input ***
    dbh_heartwood           , &  ! Heartwood diameter (m)
    dbh_sapwood             , &  ! Sapwood diameter (m)
    tree_h                  , &  ! Tree height (m)
    trunk_biomass           , &  ! Stand-level trunk biomass (g/tree)
    wood_rho(p)             , &  ! Wood density (g/cm3)
    ALM5(p)                 , &  ! Sapwood diameter proportion (m sapwood/m dbh)
    below_resource_day_stem , &  ! Increment of trunk biomass (g/tree)
                                 !
                                 ! *** Output ***
    new_dbh_heartwood       , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood         , &  ! New sapwood diameter after biomass increment (m)
    new_tree_h                &  ! New tree height after biomass increment (m)
    )

    ! Hydraulic resistance after fine root biomass increment

    new_root_biomass = root_biomass + below_resource_day_root
    call soil_root_hydraulics ( &
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    new_root_biomass     , &  ! New fine root biomass after biomass increment (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    call plant_hydraulics ( &
                              ! *** Input ***
    dbh_heartwood        , &  ! Heartwood diameter (m)
    dbh_sapwood          , &  ! Sapwood diameter (m)
    k_sap(p)             , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
    p50_sap(p)           , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
    a2_sap(p)            , &  ! Conductivity vulnerability curve coefficient (-)
    tree_h               , &  ! Tree height (m)
    leaf_wp_day_min      , &  ! Daily minimum leaf water potential (MPa)
    tot_psi              , &  ! Total soil water potential (MPa)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    ! r_whole_root_increment = r_soil+r_root+r_sap
    r_whole_root_increment = r_root+r_sap

    ! Hydraulic resistance after stem biomass increment
    
    new_root_biomass = root_biomass
    call soil_root_hydraulics ( &
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    new_root_biomass     , &  ! New fine root biomass after biomass increment (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    call plant_hydraulics ( &
                              ! *** Input ***
    new_dbh_heartwood    , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood      , &  ! New sapwood diameter after biomass increment (m)
    k_sap(p)             , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
    p50_sap(p)           , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
    a2_sap(p)            , &  ! Conductivity vulnerability curve coefficient (-)
    tree_h               , &  ! Tree height (m)
    leaf_wp_day_min      , &  ! Daily minimum leaf water potential (MPa)
    tot_psi              , &  ! Total soil water potential (MPa)
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    ! r_whole_d_increment = r_soil+r_root+r_sap
    r_whole_d_increment = r_root+r_sap

    ! Save leaf tempearture profile
    tleaf_out = tleaf

    end associate
  end subroutine spac_photosynthesis

end module mod_spac_photosynthesis
