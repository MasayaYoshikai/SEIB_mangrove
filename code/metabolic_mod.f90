  !-----------------------------------------------------------------------
  subroutine daily_production ( &
    )
    !
    ! !Description:
    ! Compute daily production for all trees
    !
    ! !Uses:
    use data_structure
    use time_counter
    use grid_status_current2
    use vegi_status_current1
    use vegi_status_current2
    use mod_spac_photosynthesis
    use mod_metabolism
    use mod_crown_morphology
    use mod_growth
    use mod_tree_allometry
    !
    ! !Argumetns:
    implicit none
    !
    ! !Local variables:
    real(8), parameter :: crown_expand_rate = 1.2d0
    integer :: no
    integer :: hour
    integer :: hour_of_year
    real(8) :: seib_lai
    real(8) :: dpai
    real(8) :: tree_h_limit
    real(8) :: available_c
    real(8) :: available_n
    real(8) :: remaining_c
    real(8) :: remaining_n
    real(8) :: resp_c_cost
    real(8) :: d_trunk
    real(8) :: d_leaf
    real(8) :: d_root
    real(8) :: d_coarse_root
    real(8) :: d_above_root
    real(8) :: d_stock_c
    real(8) :: d_stock_n
    real(8) :: above_root_ratio
    integer :: stress_flag
    integer :: purge_flag
    integer :: n_limit_flag
    integer :: expand_flag
    real(8) :: new_dbh_heartwood
    real(8) :: new_dbh_sapwood
    real(8) :: new_tree_h
    real(8) :: dbh_incr
    real(8) :: d_la
    real(8) :: d_area
    real(8) :: crown_d_new
    real(8) :: crown_d_max
    real(8) :: crown_expand_allom
    real(8) :: x, y
    real(8), dimension(PFT_no) :: turn_dummy
    integer :: n_leaf_layer
    integer :: gs_reg_flag
    real(8) :: sapflow_day
    real(8) :: sh_day
    real(8) :: lh_day
    real(8) :: leaf_wp_day_min
    real(8) :: an_mean_day_max
    real(8) :: an_top_day_max
    real(8) :: gpp_day
    real(8) :: c_uptake_day
    real(8) :: c_uptake_bottom_day
    real(8) :: n_uptake_bottom_day
    real(8) :: r_whole_root_increment
    real(8) :: r_whole_d_increment
    real(8) :: r_root_out
    real(8) :: r_sap_out
    real(8), dimension(Max_hgt) :: tleaf_out
    !---------------------------------------------------------------------

    ! Computation of all trees
    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       if ( mass_leaf(no) <= 0.0 ) cycle

       seib_lai = mass_leaf(no)*SLA(pft(no))/crown_area(no)
       tree_h_limit = real(height_limit(no))*STEP+1.3d0
       plant_wnpp_day(no) = 0.0d0
       plant_fnpp_day(no) = 0.0d0
       dbh_incr = 0.0d0

       !---------------------------------------------------------------------
       ! SPAC-Photosynthesis model
       !---------------------------------------------------------------------

       call spac_photosynthesis ( &
                                            ! *** Input ***
       o2air                           , &  ! Atmospheric O2 (mmol/mol)
       co2air                          , &  ! Atmospheric CO2 (umol/mol)
       hksat                           , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
       moist                           , &  ! Relative soil moisture (-)
       bsw                             , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
       psisat                          , &  ! Soil matric potential at saturation (MPa)
       soil_t                          , &  ! Soil water temperature (K)
       sal                             , &  ! Pore-water salinity (mol/m3)
       din                             , &  ! DIN concentration in pore-water (mol N/m3) 1.e-03 for (umol/L) -> (mol/m3)
       root_filter                     , &  ! Filtering rate of salt at root surface (-)
       root_resist                     , &  ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
       root_density                    , &  ! Specific root density (fine root) (g biomass/m3 root)
       root_radius                     , &  ! Fine root radius (m)
       root_depth                      , &  ! Rooting depth (m)
       fine_root_ratio                 , &  ! Fraction of fine root in below-ground biomass (-)
       k_sap                           , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
       p50_sap                         , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
       a2_sap                          , &  ! Conductivity vulnerability curve coefficient (-)
       wood_rho                        , &  ! Wood density (g/cm3)
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_root                        , &  ! C/N ratio in mol in root
       t_growth                        , &  ! Growth temperature (degree), mean temperature for 30 days before
       t_home                          , &  ! Home temperature (degree), mean maximum temperature of the warmest month
       minlp                           , &  ! Minimum leaf water potential (MPa)
       dleaf                           , &  ! Leaf dimension (m)
       leaf_cp                         , &  ! Leaf water capacitance (mmol H2O/m2 leaf/MPa)
       vcmaxpft                        , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
       iota0                           , &  ! Reference marginal water use efficiency in well-watered condition (umol CO2/mol H2O)
       leaf_b0                         , &  ! Sensitivity of marginal water use efficiency to leaf water potential (MPa-1)
       grow_resp                       , &  ! Growth respiration rate (g DW/g DW)
       main_resp_stem                  , &  ! Maintenance stem respiration rate at 15 degree (day-1)
       main_resp_root                  , &  ! Maintenance stem respiration rate at 15 degree (day-1)
       coarse_root_turn                , &  ! Coarse root turnover rate (day-1)
       root_turn                       , &  ! Root turnover rate (day-1)
       leaf_turn                       , &  ! Leaf turnover rate (day-1)
       leaf_resorp                     , &  ! Nutrient resorption from Leaf (fraction)
       pr_s_a                          , &  ! Slope for the scaling factor of prop root system
       pr_s_b                          , &  ! Intercept for the scaling factor of prop root system
       pr_h_a                          , &  ! Slope for the maximum root height of prop root system
       pr_h_b                          , &  ! Intercept for the maximum root height of prop root system
       pr_d                            , &  ! Mean prop root diameter (m)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
       dble(ALM5)                      , &  ! Sapwood diameter proportion (m sapwood/m dbh)
       doy                             , &  ! Day of the year (1 - 365)
       no                              , &  ! Tree index
       pft(no)                         , &  ! Species index
       dble(par_direct_rel(no,:))      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
       dble(par_diffuse_rel)           , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
       dble(nir_direct_rel(no,:))      , &  ! Profile of relative intensity of direct NIR within canopy compared to canopy top (fraction)
       dble(nir_diffuse_rel)           , &  ! Profile of relative intensity of diffused NIR within canopy compared to canopy top (fraction)
       irleaf                          , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
       dble(kbm_vis)                   , &  ! Scattering adjusted light extinction coefficient for VIS (-)
       dble(kbm_nir)                   , &  ! Scattering adjusted light extinction coefficient for NIR (-)
       rwind_profile                   , &  ! Relative wind speed profile to the canopy top (-)
       tree_h(no)                      , &  ! Tree height (m)
       dble(mass_root(no))             , &  ! Stand-level fine root biomass (g/tree)
       dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
       dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
       dble(seib_lai)                  , &  ! Leaf area index (m2 leaf/m2 ground)
       dble(crown_area(no))            , &  ! Stand-level leaf area (m2/tree)
       dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
       dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
       dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g/tree)
       dble(mass_above_root(no))       , &  ! Stand-level above-ground root biomass (g/tree)
       height(no)                      , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
       bole(no)                        , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
       time_series_tair                , &  ! Time-series air temperature (K)
       time_series_pair                , &  ! Time-series air pressure (Pa)
       time_series_wind                , &  ! Time-series wind speed (m/s)
       time_series_eair                , &  ! Time-series vapor pressure in air (Pa)
       time_series_rad_dir             , &  ! Time-series direct radiation at canopy top (W/m2)
       time_series_rad_dif             , &  ! Time-series diffused radiation at canopy top (W/m2)
                                            !
                                            ! *** In/Output ***
       leaf_wp(no)                     , &  ! Leaf water potential at current time step (0:00 AM outside spac_photosynthesis) (MPa)
       leaf_mdwp(no)                   , &  ! Midday leaf water potential of the day (MPa)
       leaf_pdwp(no)                   , &  ! Pre-dawn leaf water potential of the last day (MPa)
                                            !
                                            ! *** Output ***
       n_leaf_layer                    , &  ! Number of leaf layers
       gs_reg_flag                     , &  ! Flag of stomata regulation during the day (0: no, 1: yes)
       sapflow_day                     , &  ! Daily sapflow rate at stand-level (m3 H2O/tree/day)
       sh_day                          , &  ! Daily sensible heat flux at stand-level (MJ/tree/day)
       lh_day                          , &  ! Daily latent heat flux at stand-level (MJ/tree/day)
       leaf_wp_day_min                 , &  ! Daily minimum leaf water potential (MPa)
       an_mean_day_max                 , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
       an_top_day_max                  , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
       gpp_day                         , &  ! Daily gross photosynthesis at stand-level (mol C/tree/day)
       c_uptake_day                    , &  ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)
       c_uptake_bottom_day             , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
       n_uptake_bottom_day             , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
       r_whole_root_increment          , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
       r_whole_d_increment             , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
       r_root_out                      , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
       r_sap_out                       , &  ! Daily maximum of whole-plant stem resistance (MPa.s.tree/mmol H2O)
       tleaf_out                         &  ! Leaf temperature profile (K)
       )

!       tleaf_all(no,:) = tleaf_out
       plant_water_uptake(no) = sapflow_day
       plant_sh_day(no) = sh_day
       plant_lh_day(no) = lh_day
       plant_gpp_day(no) = gpp_day
       plant_npp_day(no) = c_uptake_day
       if (leaf_wp_day_min < klp(pft(no)) .or. gs_reg_flag == 1) then
          stress_flag = 1
       else
          stress_flag = 0
       end if
       dpai = seib_lai / real(n_leaf_layer)

       !---------------------------------------------------------------------
       ! Woody respiration
       !---------------------------------------------------------------------

       available_c = c_uptake_day
       available_n = sapflow_day*din
       available_c = available_c+residual_c(no)
       available_n = available_n+residual_n(no)
       residual_c(no) = 0.0d0
       residual_n(no) = 0.0d0

       call woody_respiration ( &
                                            ! *** Input ***
       t_acclim                        , &  ! Flag of temperature acclimation (1: Use, 0: No use)
       main_resp_stem                  , &  ! Maintenance stem respiration rate at 15 degree (day-1)
       main_resp_root                  , &  ! Maintenance stem respiration rate at 15 degree (day-1)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
       t_growth                        , &  ! Growth temperature (degree), mean temperature for 30 days before
       available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
       doy                             , &  ! Day of the year (1 - 365)
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
       dble(mass_above_root(no))       , &  ! Stand-level above-ground root biomass (g/tree)
       dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
       dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
       dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
       mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
       time_series_tair                , &  ! Time-series air temperature (K)
                                            !
                                            ! *** Output ***
       resp_c_cost                     , &  ! Woody respiration cost (mol C/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
       d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
       remaining_c                       &  ! Remaining carbon for tree growth (mol C/tree/day)
       )

       available_c = remaining_c
       plant_npp_day(no) = plant_npp_day(no)-resp_c_cost
       plant_fnpp_day(no) = plant_fnpp_day(no)-d_leaf
       if (growth_calc) then
          mass_leaf(no) = max((mass_leaf(no)+d_leaf), 1.0d0)
          mass_root(no) = max((mass_root(no)+d_root), 1.0d0)
          mass_stock_c(no) = max((mass_stock_c(no)+d_stock_c), 1.0d0)
          available_n = available_n+((d_leaf*C_in_drymass/12.0d0)              &
                        /c_n_leaf(pft(no)))*leaf_resorp(pft(no))  ! Resorbed nitrogen (mol N/tree)
       end if
       mort_gpp(no) = mort_gpp(no)+available_c*12.0d0/C_in_drymass
       mort_regu1(no) = mort_regu1(no)+d_leaf+d_root+d_stock_c*12.0d0          &
                        /C_in_drymass+(remaining_c-residual_c(no))*12.0d0      &
                        /C_in_drymass
       mort_regu1(no) = mort_regu1(no)-mass_leaf(no)*leaf_turn(pft(no))        &
                        -mass_root(no)*root_turn(pft(no))-mass_coarse_root(no) &
                        *coarse_root_turn(pft(no))
       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock_c = 0.0d0
       ! Subtract growth respiration for carbon.
       available_c = available_c*(1.0d0-grow_resp(pft(no)))
       if (no == monitor) then
          c_uptake_bottom_day_monitor = c_uptake_bottom_day
          n_uptake_bottom_day_monitor = n_uptake_bottom_day
          water_uptake_day_monitor = sapflow_day
          available_c_monitor = available_c
          available_n_monitor = available_n
          lai_monitor = seib_lai
       end if

       !---------------------------------------------------------------------
       ! Decide whether to purge the bottom layer or not
       !---------------------------------------------------------------------

       call crown_bottom_purge ( &
                                            ! *** Input ***
       leaf_turn                       , &  ! Leaf turnover rate (day-1)
       grow_resp                       , &  ! Growth respiration rate (g DW/g DW)
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       leaf_resorp                     , &  ! Nutrient resorption from Leaf (fraction)
       dble(LA_max)                    , &  ! Maximum leaf area index (m2 leaf/m2 ground)
       dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
       c_uptake_bottom_day             , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
       n_uptake_bottom_day             , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
       seib_lai                        , &  ! Leaf area index (m2 leaf/m2 ground)
       counter                         , &  ! Days since the simulation start (days)
       no                              , &  ! Tree index
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(STEP)                      , &  ! Canopy layer thickness (m)
       height(no)                      , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
       bole(no)                        , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
                                            !
                                            ! *** In/Output ***
       npp_bottom_c                    , &  ! NPP of C from 1 m2 of leaves at crown bottom layer (mol C/m2 leaf)
       npp_bottom_n                    , &  ! NPP of N from 1 m2 of leaves at crown bottom layer (mol N/m2 leaf)
                                            !
                                            ! *** Output ***
       purge_flag                        &  ! Flag of perge crown bottom layer (0: no, 1: yes)
       )

       if (purge_flag == 1) then
          if (growth_calc) then
             bole(no) = bole(no)+1
             mass_leaf(no) = mass_leaf(no)-mass_leaf(no)/real(n_leaf_layer)
!             mort_regu1(no) = mort_regu1(no)-mass_leaf(no)/n_leaf_layer
             plant_fnpp_day(no) = plant_fnpp_day(no)+mass_leaf(no)/n_leaf_layer
             available_n = available_n+((mass_leaf(no)*C_in_drymass/12.0d0     &
             /real(n_leaf_layer))/c_n_leaf(pft(no)))*leaf_resorp(pft(no))  ! Resorbed nitrogen (mol N/tree)
             n_leaf_layer = n_leaf_layer-1
          end if
          if (n_leaf_layer < 1) then
             write(*,*) 'ERROR: no leaf layer after crown purge'
          end if
       end if

       !---------------------------------------------------------------------
       ! Coarse root, leaf and fine root turnover
       !---------------------------------------------------------------------

!       ! --- Attempt 1
!
!       if (age(no) >= 2) then
!          d_leaf = -1.0d0*mass_leaf(no)*leaf_turn(pft(no))
!          d_coarse_root = -1.0d0*mass_coarse_root(no)*coarse_root_turn(pft(no))
!          d_root = -1.0d0*mass_root(no)*root_turn(pft(no))
!          if (growth_calc) then
!             mass_leaf(no) = mass_leaf(no)+d_leaf
!             mass_coarse_root(no) = mass_coarse_root(no)+d_coarse_root
!             mass_root(no) = mass_root(no)+d_root
!             available_n = available_n+(-1.0d0*d_leaf*C_in_drymass/12.0d0      &
!                           /c_n_leaf(pft(no)))*leaf_resorp(pft(no))
!          end if
!       end if
!
!       ! --- Attempt 2
!
!       turn_dummy(:) = 0.0d0
!       if (age(no) >= 2) then
!          available_n = available_n+((mass_leaf(no)*leaf_turn(pft(no))         &
!                        *C_in_drymass/12.0d0)/c_n_leaf(pft(no)))               &
!                        *leaf_resorp(pft(no))
!          if (stress_flag == 1) then
!             d_leaf = -1.0d0*mass_leaf(no)*leaf_turn(pft(no))
!             if (growth_calc) then
!                mass_leaf(no) = mass_leaf(no)+d_leaf
!             end if
!             d_leaf = 0.0d0
!             call turnover ( &
!                                                  ! *** Input ***
!             coarse_root_turn                , &  ! Coarse root turnover rate (day-1)
!             turn_dummy                      , &  ! Leaf turnover rate (day-1)
!             root_turn                       , &  ! Root turnover rate (day-1)
!             c_n_stem                        , &  ! C/N ratio in mol in stem
!             c_n_leaf                        , &  ! C/N ratio in mol in leaf
!             c_n_root                        , &  ! C/N ratio in mol in root
!             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
!             stress_flag                     , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
!             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
!             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
!             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
!             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
!             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
!             dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
!             mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
!             mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
!                                                  !
!                                                  ! *** Output ***
!             d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
!             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
!             d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
!             d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
!             d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
!             remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
!             remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
!             )
!          else
!             d_coarse_root = -1.0d0*mass_coarse_root(no)                       &
!                             *coarse_root_turn(pft(no))
!             d_root = -1.0d0*mass_root(no)*root_turn(pft(no))
!             if (growth_calc) then
!                mass_coarse_root(no) = mass_coarse_root(no)+d_coarse_root
!                mass_root(no) = mass_root(no)+d_root
!             end if
!             d_coarse_root = 0.0d0
!             d_root = 0.0d0
!             call turnover ( &
!                                                  ! *** Input ***
!             turn_dummy                      , &  ! Coarse root turnover rate (day-1)
!             leaf_turn                       , &  ! Leaf turnover rate (day-1)
!             turn_dummy                      , &  ! Root turnover rate (day-1)
!             c_n_stem                        , &  ! C/N ratio in mol in stem
!             c_n_leaf                        , &  ! C/N ratio in mol in leaf
!             c_n_root                        , &  ! C/N ratio in mol in root
!             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
!             stress_flag                     , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
!             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
!             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
!             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
!             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
!             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
!             dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
!             mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
!             mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
!                                                  !
!                                                  ! *** Output ***
!             d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
!             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
!             d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
!             d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
!             d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
!             remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
!             remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
!             )
!          end if
!
!          available_c = remaining_c
!          available_n = remaining_n
!          plant_fnpp_day(no) = plant_fnpp_day(no)+mass_leaf(no)                &
!                               *leaf_turn(pft(no))
!          if (growth_calc) then
!             mass_coarse_root(no) = max((mass_coarse_root(no)+d_coarse_root),  &
!                                        1.0d0)
!             mass_leaf(no) = max((mass_leaf(no)+d_leaf), 1.0d0)
!             mass_root(no) = max((mass_root(no)+d_root), 1.0d0)
!             mass_stock_c(no) = max((mass_stock_c(no)+d_stock_c), 1.0d0)
!             mass_stock_n(no) = max((mass_stock_n(no)+d_stock_n), 1.0d0)
!          end if
!!          mort_regu1(no) = mort_regu1(no)+d_coarse_root+d_leaf+d_root          &
!!                           +min(d_stock_c*12.0d0/C_in_drymass,                 &
!!                           d_stock_n*c_n_leaf(pft(no))*12.0d0/C_in_drymass)
!       end if

       ! --- Attempt 3

       if (age(no) >= 2) then
          available_n = available_n+((mass_leaf(no)*leaf_turn(pft(no))         &
                        *C_in_drymass/12.0d0)/c_n_leaf(pft(no)))               &
                        *leaf_resorp(pft(no))
          call turnover ( &
                                               ! *** Input ***
          coarse_root_turn                , &  ! Coarse root turnover rate (day-1)
          leaf_turn                       , &  ! Leaf turnover rate (day-1)
          root_turn                       , &  ! Root turnover rate (day-1)
          c_n_stem                        , &  ! C/N ratio in mol in stem
          c_n_leaf                        , &  ! C/N ratio in mol in leaf
          c_n_root                        , &  ! C/N ratio in mol in root
          dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
          stress_flag                     , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
          available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
          available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
          pft(no)                         , &  ! Species index (1: Rh, 2: Br)
          dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
          dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
          dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
          mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
          mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
                                               !
                                               ! *** Output ***
          d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
          d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
          d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
          d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
          d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
          remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
          remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
          )

          available_c = remaining_c
          available_n = remaining_n
          plant_fnpp_day(no) = plant_fnpp_day(no)+mass_leaf(no)                &
                               *leaf_turn(pft(no))
          if (growth_calc) then
             mass_coarse_root(no) = max((mass_coarse_root(no)+d_coarse_root),  &
                                        1.0d0)
             mass_leaf(no) = max((mass_leaf(no)+d_leaf), 1.0d0)
             mass_root(no) = max((mass_root(no)+d_root), 1.0d0)
             mass_stock_c(no) = max((mass_stock_c(no)+d_stock_c), 1.0d0)
             mass_stock_n(no) = max((mass_stock_n(no)+d_stock_n), 1.0d0)
          end if
!         mort_regu1(no) = mort_regu1(no)+d_coarse_root+d_leaf+d_root           &
!                          +min(d_stock_c*12.0d0/C_in_drymass,                  &
!                          d_stock_n*c_n_leaf(pft(no))*12.0d0/C_in_drymass)
       end if

!       ! --- Attempt 4
!
!       turn_dummy(:) = 0.0d0
!       if (age(no) >= 2) then
!          available_n = available_n+((mass_leaf(no)*leaf_turn(pft(no))         &
!                        *C_in_drymass/12.0d0)/c_n_leaf(pft(no)))               &
!                        *leaf_resorp(pft(no))
!          d_coarse_root = -1.0d0*mass_coarse_root(no)*coarse_root_turn(pft(no))
!          d_root = -1.0d0*mass_root(no)*root_turn(pft(no))
!          if (growth_calc) then
!             mass_coarse_root(no) = mass_coarse_root(no)+d_coarse_root
!             mass_root(no) = mass_root(no)+d_root
!          end if
!          d_coarse_root = 0.0d0
!          d_root = 0.0d0
!          call turnover ( &
!                                               ! *** Input ***
!          turn_dummy                      , &  ! Coarse root turnover rate (day-1)
!          leaf_turn                       , &  ! Leaf turnover rate (day-1)
!          turn_dummy                      , &  ! Root turnover rate (day-1)
!          c_n_stem                        , &  ! C/N ratio in mol in stem
!          c_n_leaf                        , &  ! C/N ratio in mol in leaf
!          c_n_root                        , &  ! C/N ratio in mol in root
!          dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
!          stress_flag                     , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
!          available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
!          available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
!          pft(no)                         , &  ! Species index (1: Rh, 2: Br)
!          dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
!          dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
!          dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
!          mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
!          mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
!                                               !
!                                               ! *** Output ***
!          d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
!          d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
!          d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
!          d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
!          d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
!          remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
!          remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
!          )
!
!          available_c = remaining_c
!          available_n = remaining_n
!          plant_fnpp_day(no) = plant_fnpp_day(no)+mass_leaf(no)                &
!                               *leaf_turn(pft(no))
!          if (growth_calc) then
!             mass_coarse_root(no) = max((mass_coarse_root(no)+d_coarse_root),  &
!                                        1.0d0)
!             mass_leaf(no) = max((mass_leaf(no)+d_leaf), 1.0d0)
!             mass_root(no) = max((mass_root(no)+d_root), 1.0d0)
!             mass_stock_c(no) = max((mass_stock_c(no)+d_stock_c), 1.0d0)
!             mass_stock_n(no) = max((mass_stock_n(no)+d_stock_n), 1.0d0)
!          end if
!!          mort_regu1(no) = mort_regu1(no)+d_coarse_root+d_leaf+d_root          &
!!                           +min(d_stock_c*12.0d0/C_in_drymass,                 &
!!                           d_stock_n*c_n_leaf(pft(no))*12.0d0/C_in_drymass)
!       end if

       d_coarse_root = 0.0d0
       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock_c = 0.0d0
       d_stock_n = 0.0d0
       d_trunk = 0.0d0
       d_above_root = 0.0d0

       if (min(available_c, available_n) > 0.0d0) then

          !---------------------------------------------------------------------
          ! Reload to stock C and/or N
          !---------------------------------------------------------------------

          call stock_reload ( &
                                               ! *** Input ***
          stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
          c_n_stem                        , &  ! C/N ratio in mol in stem
          c_n_leaf                        , &  ! C/N ratio in mol in leaf
          dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
          available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
          available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
          pft(no)                         , &  ! Species index (1: Rh, 2: Br)
          dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
          dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
          mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
          mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
                                               !
                                               ! *** Output ***
          d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
          d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
          remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
          remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
          )

          available_c = remaining_c
          available_n = remaining_n
          if (growth_calc) then
             mass_stock_c(no) = mass_stock_c(no)+d_stock_c
             mass_stock_n(no) = mass_stock_n(no)+d_stock_n
          end if
!          mort_regu1(no) = mort_regu1(no)+min(d_stock_c*12.0d0/C_in_drymass,   &
!                           d_stock_n*c_n_leaf(pft(no))*12.0d0/C_in_drymass)
          d_stock_c = 0.0d0
          d_stock_n = 0.0d0
          if (available_c/available_n > c_n_leaf(pft(no))) then
             n_limit_flag = 1
          else
             n_limit_flag = 0
          end if
          if (n_limit_flag == 1 .or. stress_flag == 1) then

             !---------------------------------------------------------------------
             ! Below-ground limitation (N-limited) or
             ! stressed condition due to low water potential
             !---------------------------------------------------------------------

             call belowground_limit ( &
                                                  ! *** Input ***
             optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
             c_n_stem                        , &  ! C/N ratio in mol in stem
             c_n_leaf                        , &  ! C/N ratio in mol in leaf
             c_n_root                        , &  ! C/N ratio in mol in root
             fine_root_ratio                 , &  ! Fraction of fine root in below-ground biomass (-)
             stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
             dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
             dble(DBH_limit)                 , &  ! Limitation of trunk diameter (m)
             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
             stress_flag                     , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
             n_leaf_layer                    , &  ! Number of leaf layers
             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
             r_whole_root_increment          , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
             r_whole_d_increment             , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
!             r_sap_out                       , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
!             r_root_out                      , &  ! Daily maximum of Whole-plant stem resistance (MPa.s.tree/mmol H2O)
             tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
             tree_h(no)                      , &  ! Tree height (m)
             dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
             dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
             mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
             mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
             dble(crown_area(no))            , &  ! Stand-level leaf crown area (m2 ground/tree)
             dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
             dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                                                  !
                                                  ! *** Output ***
             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
             d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
             d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
             d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
             d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
             d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
             remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
             remaining_n                     , &  ! Remaining nitrogen for tree growth (mol N/tree/day)
             expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
             )

          else

             !---------------------------------------------------------------------
             ! Above-ground limitation (C-limited) and
             ! current leaf water potential is not stressful.
             !---------------------------------------------------------------------

             call aboveground_limit ( &
                                                  ! *** Input ***
             optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
             c_n_stem                        , &  ! C/N ratio in mol in stem
             c_n_leaf                        , &  ! C/N ratio in mol in leaf
             stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
             dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
             dble(DBH_limit)                 , &  ! Limitation of trunk diameter (m)
             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
             n_leaf_layer                    , &  ! Number of leaf layers
             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
             tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
             dble(par_direct)                , &  ! Direct PAR on canopy top of midday (umol photon/m2/s)
             dble(par_diffuse)               , &  ! Diffused PAR on canopy top of midday (umol photon/m2/s)
             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
             dble(par_direct_rel(no,:))      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
             dble(par_diffuse_rel)           , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
             height(no)                      , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
             tree_h(no)                      , &  ! Tree height (m)
             dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
             mass_stock_c(no)                , &  ! Carbon in stock (mol C/tree)
             mass_stock_n(no)                , &  ! Nitrogen in stock (mol N/tree)
             dble(crown_area(no))            , &  ! Stand-level leaf crown area (m2 ground/tree)
             dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
             dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                                                  !
                                                  ! *** Output ***
             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
             d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
             d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
             d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
             remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
             remaining_n                     , &  ! Remaining nitrogen for tree growth (mol N/tree/day)
             expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
             )
          end if
          if (d_trunk > 0.0d0) then
             if (expand_flag == 1) then
                call proot_allometry ( &
                                                     ! *** Input ***
                wood_rho                        , &  ! Wood density (g/cm3)
                pr_s_a                          , &  ! Slope for the scaling factor of prop root system
                pr_s_b                          , &  ! Intercept for the scaling factor of prop root system
                pr_h_a                          , &  ! Slope for the maximum root height of prop root system
                pr_h_b                          , &  ! Intercept for the maximum root height of prop root system
                pr_d                            , &  ! Mean prop root diameter (m)
                pft(no)                         , &  ! Species index (1: Rh, 2: Br)
                dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
                dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
                dble(mass_above_root(no))       , &  ! Stand-level above-ground root biomass (g/tree)
                                                     !
                                                     ! *** In/Output ***
                d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                                     !
                                                     ! *** Output ***
                d_above_root                    , &  ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
                above_root_ratio                  &  ! Required fraction of above-ground root biomass in above-ground
                )
             end if
             call tree_allometry (  &
                                                 ! *** Input ***
             dble(dbh_heartwood(no))        , &  ! Heartwood diameter (m)
             dble(dbh_sapwood(no))          , &  ! Sapwood diameter (m)
             tree_h(no)                     , &  ! Tree height (m)
             dble(mass_trunk(no))           , &  ! Stand-level trunk biomass (g/tree)
             wood_rho(pft(no))              , &  ! Wood density (g/cm3)
             dble(ALM5(pft(no)))            , &  ! Sapwood diameter proportion (m sapwood/m dbh)
             d_trunk                        , &  ! Increment of trunk biomass (g/tree)
                                                 !
                                                 ! *** Output ***
             new_dbh_heartwood              , &  ! New heartwood diameter after biomass increment (m)
             new_dbh_sapwood                , &  ! New sapwood diameter after biomass increment (m)
             new_tree_h                       &  ! New tree height after biomass increment (m)
             )
             dbh_incr = max(new_dbh_heartwood+new_dbh_sapwood                  &
                        -dble(dbh_heartwood(no))-dble(dbh_sapwood(no)), 0.0d0)
             if (growth_calc) then
                if (expand_flag == 0) then ! Tree heght grows.
                   tree_h(no) = new_tree_h
                   height(no) = int((tree_h(no)-1.3)/STEP)
                   height(no) = max(2, height(no))
                else ! Stem diameter grows.
                   dbh_heartwood(no) = new_dbh_heartwood
                   dbh_sapwood(no) = new_dbh_sapwood
                end if
             end if
          end if
          if (growth_calc) then
             mass_leaf(no) = mass_leaf(no)+d_leaf
             mass_root(no) = mass_root(no)+d_root
             mass_trunk(no) = mass_trunk(no)+d_trunk
             mass_coarse_root(no) = mass_coarse_root(no)+d_coarse_root
             mass_above_root(no) = mass_above_root(no)+d_above_root
             mass_stock_c(no) = mass_stock_c(no)+d_stock_c
             mass_stock_n(no) = mass_stock_n(no)+d_stock_n
          end if
          plant_wnpp_day(no) = plant_wnpp_day(no)+d_trunk+d_above_root
       end if
       residual_c(no) = remaining_c
       residual_n(no) = remaining_n

       !---------------------------------------------------------------------
       ! Tree crown expansion
       !---------------------------------------------------------------------

       ! Option 1-1: new crown diameter and area based on increase in leaf area
       crown_d_max = crown_allometry ( &
                                            ! *** Input ***
       crown_a                         , &  ! Tree crown allometric parameter
       crown_b                         , &  ! Tree crown allometric parameter
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
       dble(dbh_sapwood(no))             &  ! Sapwood diameter (m)
       )
       d_la = d_leaf*SLA(pft(no))
       d_area = d_la/seib_lai
       crown_d_new = sqrt(4*(crown_area(no)+d_area)/PI)
       crown_d_max = min(crown_d_max, radius_limit(no)*2.0d0)
       if (growth_calc) then
          crown_diameter(no) = min(crown_d_new, crown_d_max)
          crown_area(no) = (crown_diameter(no)*0.5d0)                          &
                           *(crown_diameter(no)*0.5d0)*PI
       end if

!       ! Option 1-2: new crown diameter and area based on increase in leaf area
!       !             <= expansion only when dpai satidfies optimum_dpai
!       if (dpai >= optimum_dpai(pft(no))*0.95d0) then
!          d_la = d_leaf*SLA(pft(no))
!          d_area = d_la/seib_lai
!          crown_d_new = sqrt(4*(crown_area(no)+d_area)/PI)
!       else
!          crown_d_new = crown_diameter(no)
!       end if
!       crown_d_max = min(crown_d_max, radius_limit(no)*2.0d0)
!       if (growth_calc) then
!          crown_diameter(no) = min(crown_d_new, crown_d_max)
!          crown_area(no) = (crown_diameter(no)*0.5d0)                          &
!                           *(crown_diameter(no)*0.5d0)*PI
!       end if
!
!       ! Option 2: new crown diameter and area besed on allometric curve
!       crown_expand_allom = dbh_incr*crown_a(pft(no))*crown_b(pft(no))         &
!                            *(dble(dbh_heartwood(no))+dble(dbh_sapwood(no)))   &
!                            **(crown_a(pft(no))-1.0d0)
!!       if (growth_calc .and. stress_flag==0) then
!       if (growth_calc) then
!!          x = min(crown_d_max, crown_diameter(no)*(1.0d0+(crown_expand_rate    &
!!              -1.0d0)/Day_in_Year))
!          x = min(crown_d_max, crown_diameter(no)+crown_expand_allom)
!          y = min(x, radius_limit(no)*2.0d0)
!          crown_diameter(no) = max(crown_diameter(no), y)
!          crown_area(no) = (crown_diameter(no)*0.5d0)                          &
!                           *(crown_diameter(no)*0.5d0)*PI
!       end if

       la(no) = mass_leaf(no)*SLA(pft(no))
       mort_regu2(no) = mort_regu2(no)+la(no)/Day_in_Year

    end do

  end subroutine daily_production
