module mod_monitoring

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: monitor_index
  public  :: monitor_tree
  public  :: monitor_plot
  public  :: monitor_diurnal
  private :: ind_pft_sample
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine monitor_index ( &
                                       ! *** Input ***
    flag_tree_presence            , &  ! Flag of tree presence
    tree_exist                    , &  ! Flag of tree presence
                                       !
                                       ! *** In/Output ***
    monitor                         &  ! Tree index for monitoring
    )
    !
    ! !DESCRIPTION:
    ! Decide a tree to monitor.
    !
    ! !USES:
    use data_structure, only : Max_no, randf
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: flag_tree_presence
    logical, intent(in)    :: tree_exist(:)
    integer, intent(inout) :: monitor
    !
    ! !LOCAL VARIABLES:
    integer, dimension(Max_no) :: tree_matrix
    integer :: tree_n
    integer :: no
    real :: rand_n
    integer :: loop
    !---------------------------------------------------------------------

    if (flag_tree_presence == 0) then
       return
    end if
    if (monitor > 0) then
       if (.not. tree_exist(monitor)) then
          monitor = 0
       end if
    end if
    if (monitor > 0) then
       return
    end if

    tree_matrix(:) = 0
    tree_n = 0
    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       tree_n = tree_n+1
       tree_matrix(tree_n) = no
    end do

    do loop = 1, 100
       do no = 1, tree_n
          rand_n = randf()
          if (rand_n < (1/real(tree_n))) exit
       end do
       if (rand_n < (1/real(tree_n))) exit
    end do
    monitor = tree_matrix(no)

  end subroutine monitor_index

  !-----------------------------------------------------------------------
  subroutine monitor_tree ( &
                                       ! *** Input ***
                                       !
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** From SEIB-DGVM ***
    Fn                            , &  ! File number for monitoring tree
    STEP                          , &  ! Canopy layer thickness (m)
    year                          , &  ! Year
    day_of_year                   , &  ! Day of the year (1 - 365)
    pft                           , &  ! Species index (1: Rh, 2: Br)
    bottom_layer_monitor          , &  ! Bottom layer of canopy
    top_layer_monitor             , &  ! Top layer of canopy
    par_direct_rel                , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
    par_diffuse_rel               , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
    leaf_wp_monitor               , &  ! Midday leaf water potential (MPa)
    leaf_wpmin_monitor            , &  ! Daily minimum leaf water potential (MPa)
    gs_top_monitor                , &  ! Midday top layer stomatal conductance (mol H2O/m2 leaf/s)
    an_top_monitor                , &  ! Midday top layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_bot_monitor                , &  ! Midday bottom layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    et_top_monitor                , &  ! Midday top layer leaf transpiration rate (mol H2O/m2 leaf/s)
    et_bot_monitor                , &  ! Midday bottom layer leaf transpiration rate (mol H2O/m2 leaf/s)
    c_uptake_bottom_day_monitor   , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day_monitor   , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    water_uptake_day_monitor      , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
    available_c_monitor           , &  ! Remaining carbon for tree growth after respiration and turnover (mol C/tree/day)
    available_n_monitor           , &  ! Remaining nitrogen for tree growth after turnover (mol N/tree/day)
    height_max_monitor            , &  ! Potential tree height based on allometric relation (m)
    height_limit                  , &  ! Tree height limitation based on proximate trees (unit is STEP!!)
    crown_d_max_monitor           , &  ! Potential crown diameter based on allometric relation (m)
    radius_limit                  , &  ! Crown radius limitation based on proximate trees (m)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    tree_h                        , &  ! Tree height (m)
    crown_diameter                , &  ! Crown diameter (m)
    seib_height                   , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
    lai_monitor                   , &  ! Leaf area index (m2 leaf/m2 ground)
    mass_trunk                    , &  ! Stand-level trunk biomass (g trunk/tree)
    mass_coarse_root              , &  ! Stand-level coarse root biomass (g/tree)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    mass_root                     , &  ! Stand-level fine root biomass (g root/tree)
    mass_stock                      &  ! Stand-level stock biomass (g stock/tree)
    )
    !
    ! !DESCRIPTION:
    ! Daily output for monitoring tree
    !
    ! !USES:
    use data_structure, only : display_flag
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: time_series_rad_dir(:), time_series_rad_dif(:)
    integer, intent(in) :: Fn
    real, intent(in)    :: STEP
    integer, intent(in) :: year, day_of_year, pft, bottom_layer_monitor
    integer, intent(in) :: top_layer_monitor
    real, intent(in)    :: par_direct_rel(:), par_diffuse_rel(:)
    real(8), intent(in) :: leaf_wp_monitor, leaf_wpmin_monitor
    real(8), intent(in) :: gs_top_monitor, an_top_monitor
    real(8), intent(in) :: an_bot_monitor, et_top_monitor, et_bot_monitor
    real(8), intent(in) :: c_uptake_bottom_day_monitor
    real(8), intent(in) :: n_uptake_bottom_day_monitor
    real(8), intent(in) :: water_uptake_day_monitor
    real(8), intent(in) :: available_c_monitor, available_n_monitor
    real(8), intent(in) :: height_max_monitor
    integer, intent(in) :: height_limit
    real(8), intent(in) :: crown_d_max_monitor
    real, intent(in)    :: radius_limit, dbh_heartwood, dbh_sapwood
    real(8), intent(in) :: tree_h
    real, intent(in)    :: crown_diameter
    integer, intent(in) :: seib_height, seib_bole
    real(8), intent(in) :: lai_monitor
    real, intent(in)    :: mass_trunk, mass_coarse_root, mass_above_root
    real, intent(in)    :: mass_leaf, mass_root, mass_stock
    !
    ! !LOCAL VARIABLES:
    real(8), parameter :: mmh2o = 18.02d0/1000.0d0
    real(8), parameter :: denh2o = 1000.0d0
    integer :: hour_of_year
    real(8) :: par_midday
    real(8) :: par_top
    real(8) :: par_bot
    real(8) :: water_uptake_day_monitor_kg
    real(8) :: height_limit_m
    real(8) :: crown_dep
    real(8) :: dpai
    !---------------------------------------------------------------------

    ! PAR at canopy top at midday (umol photon/m2 ground/s)
    hour_of_year = (day_of_year-1)*24+12+1
    par_midday = 4.6d0*0.43d0*time_series_rad_dir(hour_of_year)                &
                 +4.2d0*0.57d0*time_series_rad_dif(hour_of_year)

    ! PAR at crown top and bottom layer at midday (umol photon/m2 ground/s)
    par_top = par_direct_rel(top_layer_monitor)*4.6d0*0.43d0                   &
              *time_series_rad_dir(hour_of_year)                               &
              +par_diffuse_rel(top_layer_monitor)*4.2d0*0.57d0                 &
              *time_series_rad_dif(hour_of_year)
    par_bot = par_direct_rel(bottom_layer_monitor)*4.6d0*0.43d0                &
              *time_series_rad_dir(hour_of_year)                               &
              +par_diffuse_rel(bottom_layer_monitor)*4.2d0*0.57d0              &
              *time_series_rad_dif(hour_of_year)

    ! m3 H2O/tree/day -> kg H2O/tree/day
    water_uptake_day_monitor_kg = water_uptake_day_monitor*denh2o

    ! STEP -> m
    height_limit_m = real(height_limit)*STEP+1.3d0

    ! Crown depth (m)
    crown_dep = real(seib_height-seib_bole)*STEP

    ! Layer leaf area index (m2 leaf/m2 ground)
    dpai = lai_monitor/real(seib_height-seib_bole)

    ! Write monitoring tree output

    write (Fn, '(3(i4,a), 30(f12.5,a))') &
    year                            , ',', &
    day_of_year                     , ',', &
    pft                             , ',', &
    par_midday                      , ',', &
    par_top                         , ',', &
    par_bot                         , ',', &
    leaf_wp_monitor                 , ',', &
    gs_top_monitor                  , ',', &
    an_top_monitor                  , ',', &
    an_bot_monitor                  , ',', &
    et_top_monitor                  , ',', &
    et_bot_monitor                  , ',', &
    c_uptake_bottom_day_monitor     , ',', &
    n_uptake_bottom_day_monitor     , ',', &
    water_uptake_day_monitor_kg     , ',', &
    available_c_monitor             , ',', &
    available_n_monitor             , ',', &
    height_max_monitor              , ',', &
    height_limit_m                  , ',', &
    crown_d_max_monitor             , ',', &
    radius_limit*2.0d0              , ',', &
    dbh_heartwood+dbh_sapwood       , ',', &
    tree_h                          , ',', &
    crown_diameter                  , ',', &
    crown_dep                       , ',', &
    lai_monitor                     , ',', &
    dpai                            , ',', &
    mass_trunk/1000.0d0             , ',', &
    mass_coarse_root/1000.0d0       , ',', &
    mass_above_root/1000.0d0        , ',', &
    mass_leaf/1000.0d0              , ',', &
    mass_root/1000.0d0              , ',', &
    mass_stock/1000.0d0

  end subroutine monitor_tree

  !-----------------------------------------------------------------------
  subroutine monitor_plot ( &
                                       ! *** Input ***
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
    Fn1                           , &  ! File number for tree biomass output
    Fn2                           , &  ! File number for plot output
    year                          , &  ! Year
    day_of_year                   , &  ! Day of the year (1 - 365)
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    Max_hgt                       , &  ! Number of canopy layer
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    age                           , &  ! Tree age (year: 1~)
    tree_h                        , &  ! Tree height (m)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    mass_trunk                    , &  ! Stand-level trunk biomass (g trunk/tree)
    mass_coarse_root              , &  ! Stand-level coarse root biomass (g/tree)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    mass_root                     , &  ! Stand-level fine root biomass (g root/tree)
    dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    leaf_mdwp                     , &  ! Midday leaf water potential of the day (MPa)
    leaf_pdwp                     , &  ! Pre-dawn leaf water potential of the last day (MPa)
    plant_water_uptake            , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
    plant_sh_day                  , &  ! Daily sensible heat flux at stand-level (MJ/tree/day)
    plant_lh_day                  , &  ! Daily latent heat flux at stand-level (MJ/tree/day)
    plant_gpp_day                 , &  ! Daily gross photosynthesis at stand-level (mol C/tree/day)
    plant_npp_day                 , &  ! Daily net photosynthesis at stand-level (mol C/tree/day)
    plant_wnpp_day                , &  ! Daily woody net primary production of each tree (g DW/tree/day)
    plant_fnpp_day                , &  ! Daily foliage net primary production of each tree (g DW/tree/day)
    irveg                         , &  ! Absorbed longwave radiation, vegetation (W/m2)
    irsoi                         , &  ! Absorbed longwave radiation, ground (W/m2)
    ir_tree                       , &  ! VIS-interruption-coefficient by tree corwn
    rnsoi_day                     , &  ! Daily mean net radiation, ground (W/m2)
    rnveg_day                     , &  ! Daily mean net radiation, vegetation (W/m2)
    etsoi_day                     , &  ! Daily soil evaporation rate (mm/day)
    etveg_day                     , &  ! Daily transpiration rate (mm/day)
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    dip                           , &  ! DIP concentration in pore-water (mol P/m3)
    tree_n_est                    , &  ! Total number of established trees of each PFT in the entire simulation (trees)
    tree_n_die                      &  ! Total number of died trees of each PFT in the entire simulation (trees)
    )
    !
    ! !DESCRIPTION:
    ! Daily output for plot-scale variables
    !
    ! !USES:
    use data_structure, only : n_spe, Day_in_Year, display_flag
    use time_counter,   only : doy
    use Statistics
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: SLA(:), time_series_rad_dir(:)
    real(8), intent(in) :: time_series_rad_dif(:)
    integer, intent(in) :: Fn1, Fn2, year, day_of_year, Max_loc, Max_hgt, Max_no
    logical, intent(in) :: tree_exist(:)
    integer, intent(in) :: pft(:), age(:)
    real(8), intent(in) :: tree_h(:)
    real, intent(in)    :: dbh_heartwood(:), dbh_sapwood(:), mass_trunk(:)
    real, intent(in)    :: mass_coarse_root(:), mass_above_root(:)
    real, intent(in)    :: mass_leaf(:), mass_root(:)
    real(8), intent(in) :: dpai_layer_sum(:)
    real(8), intent(in) :: leaf_mdwp(:), leaf_pdwp(:)
    real(8), intent(in) :: plant_water_uptake(:), plant_sh_day(:)
    real(8), intent(in) :: plant_lh_day(:), plant_gpp_day(:)
    real(8), intent(in) :: plant_npp_day(:), plant_wnpp_day(:)
    real(8), intent(in) :: plant_fnpp_day(:), irveg, irsoi
    real, intent(in)    :: ir_tree
    real(8), intent(in) :: rnsoi_day, rnveg_day, etsoi_day, etveg_day
    real(8), intent(in) :: sal, din, dip
    integer, intent(in) :: tree_n_est(:), tree_n_die(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p
    integer :: no
    integer :: i
    real(8), dimension(3)      :: qua
    real(8), dimension(Max_no) :: ind_var_pft
    real,    dimension(n_spe)  :: tree_density_plot
    real(8), dimension(n_spe)  :: mean_tree_h
    real,    dimension(n_spe)  :: mean_dbh
    real,    dimension(n_spe)  :: above_biomass_plot
    real,    dimension(n_spe)  :: leaf_biomass_plot
    real,    dimension(n_spe)  :: root_biomass_plot
    real,    dimension(n_spe)  :: coarse_root_biomass_plot
    real,    dimension(n_spe)  :: above_root_biomass_plot
    real(8), dimension(n_spe)  :: water_uptake_plot
    real(8), dimension(n_spe)  :: plant_sh_plot
    real(8), dimension(n_spe)  :: plant_lh_plot
    real(8), dimension(n_spe)  :: gpp_plot
    real(8), dimension(n_spe)  :: npp_plot
    real(8), dimension(n_spe)  :: wnpp_plot
    real(8), dimension(n_spe)  :: fnpp_plot
    real(8), dimension(n_spe)  :: lai_species_plot
    real(8)                    :: lai_plot
    real(8)                    :: leaf_mdwp_median
    real(8), dimension(n_spe)  :: leaf_mdwp_pft_median
    real(8)                    :: leaf_mdwp_quart_low
    real(8), dimension(n_spe)  :: leaf_mdwp_pft_quart_low
    real(8)                    :: leaf_mdwp_quart_high
    real(8), dimension(n_spe)  :: leaf_mdwp_pft_quart_high
    real(8)                    :: leaf_pdwp_median
    real(8), dimension(n_spe)  :: leaf_pdwp_pft_median
    real(8), parameter         :: leaf_wp_blank = 0.0d0
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Tree biomass output
    !---------------------------------------------------------------------

    tree_density_plot(:) = 0.0d0
    mean_tree_h(:) = 0.0d0
    mean_dbh(:) = 0.0d0
    above_biomass_plot(:) = 0.0d0
    leaf_biomass_plot(:) = 0.0d0
    root_biomass_plot(:) = 0.0d0
    coarse_root_biomass_plot(:) = 0.0d0
    above_root_biomass_plot(:) = 0.0d0
    water_uptake_plot(:) = 0.0d0
    plant_sh_plot(:) = 0.0d0
    plant_lh_plot(:) = 0.0d0
    gpp_plot(:) = 0.0d0
    npp_plot(:) = 0.0d0
    wnpp_plot(:) = 0.0d0
    fnpp_plot(:) = 0.0d0
    lai_species_plot(:) = 0.0d0

    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       p = pft(no)
       tree_density_plot(p) = tree_density_plot(p)+1.0d0                          ! tree number
       mean_tree_h(p) = mean_tree_h(p)+tree_h(no)                                 ! m
       mean_dbh(p) = mean_dbh(p)+dbh_heartwood(no)+dbh_sapwood(no)                ! m
       above_biomass_plot(p) = above_biomass_plot(p)+mass_trunk(no)*1.e-03        ! kg
       leaf_biomass_plot(p) = leaf_biomass_plot(p)+mass_leaf(no)*1.e-03           ! kg
       root_biomass_plot(p) = root_biomass_plot(p)+mass_root(no)*1.e-03           ! kg
       coarse_root_biomass_plot(p) = coarse_root_biomass_plot(p)               &
                                     +mass_coarse_root(no)*1.e-03                 ! kg
       above_root_biomass_plot(p) = above_root_biomass_plot(p)                 &
                                    +mass_above_root(no)*1.e-03                   ! kg
       water_uptake_plot(p) = water_uptake_plot(p)+plant_water_uptake(no)         ! m3 H2O/day in the computational domain
       plant_sh_plot(p) = plant_sh_plot(p)+plant_sh_day(no)                       ! MJ/day in the computational domain
       plant_lh_plot(p) = plant_lh_plot(p)+plant_lh_day(no)                       ! MJ/day in the computational domain
       gpp_plot(p) = gpp_plot(p)+plant_gpp_day(no)*12.0d0                         ! g C/day in the computational domain
       npp_plot(p) = npp_plot(p)+plant_npp_day(no)*12.0d0                         ! g C/day in the computational domain
       wnpp_plot(p) = wnpp_plot(p)+plant_wnpp_day(no)*1.e-03                      ! kg/day in the computational domain
       fnpp_plot(p) = fnpp_plot(p)+plant_fnpp_day(no)*1.e-03                      ! kg/day in the computational domain
       lai_species_plot(p) = lai_species_plot(p)+mass_leaf(no)*SLA(p)             ! m2 leaf in the computational domain
    end do
    do p = 1, n_spe
       mean_tree_h(p) = mean_tree_h(p)/max(1.0d0, tree_density_plot(p))           ! m
       mean_dbh(p) = mean_dbh(p)/max(1.0d0, tree_density_plot(p))                 ! m
       above_biomass_plot(p) = (above_biomass_plot(p)/(real(Max_loc)**2.0d0))  &
                               *10.0d0                                            ! kg/m2 -> Mg/ha
       leaf_biomass_plot(p) = (leaf_biomass_plot(p)/(real(Max_loc)**2.0d0))    &
                               *10.0d0                                            ! kg/m2 -> Mg/ha
       root_biomass_plot(p) = (root_biomass_plot(p)/(real(Max_loc)**2.0d0))    &
                               *10.0d0                                            ! kg/m2 -> Mg/ha
       coarse_root_biomass_plot(p) = (coarse_root_biomass_plot(p)              &
                                     /(real(Max_loc)**2.0d0))*10.0d0              ! kg/m2 -> Mg/ha
       above_root_biomass_plot(p) = (above_root_biomass_plot(p)                &
                                    /(real(Max_loc)**2.0d0))*10.0d0               ! kg/m2 -> Mg/ha
       water_uptake_plot(p) = (water_uptake_plot(p)/(real(Max_loc)**2.0d0))    &
                              *1000.0d0                                           ! m/day -> mm/day
       plant_sh_plot(p) = plant_sh_plot(p)/(real(Max_loc)**2.0d0)                 ! MJ/m2 ground/day
       plant_lh_plot(p) = plant_lh_plot(p)/(real(Max_loc)**2.0d0)                 ! MJ/m2 ground/day
       gpp_plot(p) = gpp_plot(p)/(real(Max_loc)**2.0d0)                           ! g C/m2 ground/day
       npp_plot(p) = npp_plot(p)/(real(Max_loc)**2.0d0)                           ! g C/m2 ground/day
       wnpp_plot(p) = (wnpp_plot(p)/(real(Max_loc)**2.0d0))*10.0d0                ! kg/m2/day -> Mg/ha/day
       fnpp_plot(p) = (fnpp_plot(p)/(real(Max_loc)**2.0d0))*10.0d0                ! kg/m2/day -> Mg/ha/day
       lai_species_plot(p) = lai_species_plot(p)/(real(Max_loc)**2.0d0)           ! m2 leaf/m2 ground
    end do
    do p = 1, n_spe
       tree_density_plot(p) = tree_density_plot(p)/(real(Max_loc)**2.0d0)         ! tree/m2
    end do

    ! Quartile of leaf water potential for all species
    call Quartile_1d(leaf_mdwp, qua, leaf_wp_blank)
    leaf_mdwp_median = qua(2)
    leaf_mdwp_quart_low = qua(1)
    leaf_mdwp_quart_high = qua(3)
    call Quartile_1d(leaf_pdwp, qua, leaf_wp_blank)
    leaf_pdwp_median = qua(2)

    ! Quartile of leaf water potential per species
    do p = 1, n_spe
       ind_var_pft = ind_pft_sample ( &
                                            ! *** Input ***
       leaf_mdwp                       , &  ! Individual variable
       p                               , &  ! Specified pft to extract samples
       Max_no                          , &  ! Maximum number of individual stands
       tree_exist                      , &  ! Flag of tree presence
       pft                             , &  ! Species index (1: Rh, 2: Br)
       age                             , &  ! Tree age (year: 1~)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                       &  ! Sapwood diameter (m)
       )
       call Quartile_1d(ind_var_pft, qua, leaf_wp_blank)
       leaf_mdwp_pft_median(p)     = qua(2)
       leaf_mdwp_pft_quart_low(p)  = qua(1)
       leaf_mdwp_pft_quart_high(p) = qua(3)
       ind_var_pft = ind_pft_sample ( &
                                            ! *** Input ***
       leaf_pdwp                       , &  ! Individual variable
       p                               , &  ! Specified pft to extract samples
       Max_no                          , &  ! Maximum number of individual stands
       tree_exist                      , &  ! Flag of tree presence
       pft                             , &  ! Species index (1: Rh, 2: Br)
       age                             , &  ! Tree age (year: 1~)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                       &  ! Sapwood diameter (m)
       )
       call Quartile_1d(ind_var_pft, qua, leaf_wp_blank)
       leaf_pdwp_pft_median(p) = qua(2)
    end do

    ! Write output
    write (Fn1, '(2(i4,a))',advance='no') year,',', day_of_year, ','
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') tree_density_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') mean_tree_h(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') mean_dbh(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') above_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') coarse_root_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') above_root_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') root_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') water_uptake_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') plant_sh_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') plant_lh_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') gpp_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') npp_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') wnpp_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') fnpp_plot(p),','
    end do
    write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_median,','
    write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_quart_low,','
    write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_quart_high,','
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_pft_median(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_pft_quart_low(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_mdwp_pft_quart_high(p),','
    end do
    write (Fn1, '(f12.5,a)',advance='no') leaf_pdwp_median,','
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_pdwp_pft_median(p),','
    end do
    do p = 1, n_spe
       if (p < n_spe) then
       write (Fn1, '(f12.5,a)',advance='no') lai_species_plot(p),','
       else
       write (Fn1, '(f12.5,a)') lai_species_plot(p)
       end if
    end do

    if (display_flag > 0) then
       if (day_of_year <= 31) write(*,*) 'Jan'
       if (31 < day_of_year .and. day_of_year <= 59) write(*,*) 'Feb'
       if (59 < day_of_year .and. day_of_year <= 90) write(*,*) 'Mar'
       if (90 < day_of_year .and. day_of_year <= 120) write(*,*) 'Apr'
       if (120 < day_of_year .and. day_of_year <= 151) write(*,*) 'May'
       if (151 < day_of_year .and. day_of_year <= 181) write(*,*) 'Jun'
       if (181 < day_of_year .and. day_of_year <= 212) write(*,*) 'Jul'
       if (212 < day_of_year .and. day_of_year <= 243) write(*,*) 'Aug'
       if (243 < day_of_year .and. day_of_year <= 273) write(*,*) 'Sep'
       if (273 < day_of_year .and. day_of_year <= 304) write(*,*) 'Oct'
       if (304 < day_of_year .and. day_of_year <= 334) write(*,*) 'Nov'
       if (334 < day_of_year .and. day_of_year <= 365) write(*,*) 'Dec'
       write(*,*) 'Rh      year', '      day', '     tree number', '     height_mean', '     dbh_mean'
       write(*,*) year, day_of_year, tree_density_plot(1) * 7.0d0 * 7.0d0 * 3.1415d0, mean_tree_h(1), mean_dbh(1)
       write(*,*) 'Br      year', '      day', '     tree number', '     height_mean', '     dbh_mean'
       write(*,*) year, day_of_year, tree_density_plot(2) * 7.0d0 * 7.0d0 * 3.1415d0, mean_tree_h(2), mean_dbh(2)
       write(*,*) 'Rh      year', '      day', '     LAI', '                       AGB'
       write(*,*) year, day_of_year, lai_species_plot(1), above_biomass_plot(1)
       write(*,*) 'Br      year', '      day', '     LAI', '                       AGB'
       write(*,*) year, day_of_year, lai_species_plot(2), above_biomass_plot(2)
    end if

    !---------------------------------------------------------------------
    ! Plot output
    !---------------------------------------------------------------------

    lai_plot = 0.0d0
    do i = 1, Max_hgt
       lai_plot = lai_plot+dpai_layer_sum(i)
    end do
    ! Write plot output
    if (doy == Day_in_Year) then
    write (Fn2, '(2(i4,a))',advance='no') year,',', day_of_year, ','
    write (Fn2, '(5(f12.5,a))',advance='no') irveg,',', irsoi, ',', lai_plot,  &
          ',', ir_tree, ',', rnsoi_day, ','
    write (Fn2, '(3(f12.5,a))',advance='no') rnveg_day,',', etsoi_day,         &
          ',', etveg_day, ','
    write (Fn2, '(2(f12.5,a))',advance='no') sal * 58.44d0 / 1000.0d0,',',     &
          din, ','
    write (Fn2, '(f12.5,a)') dip
    end if

  end subroutine monitor_plot

  !-----------------------------------------------------------------------
  subroutine monitor_diurnal ( &
                                       ! *** Input ***
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    Fn                            , &  ! File number for diurnal output
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    age                           , &  ! Tree age (year: 1~)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    plant_gpp_hour                , &  ! Hourly gross photosynthesis at stand-level (mol C/tree/hour)
    plant_sap_hour                , &  ! Hourly sapflow at stand-level (m3 H2O/tree/day)
    plant_et_hour                 , &  ! Hourly transpiration at stand-level (m3 H2O/tree/day)
    leaf_wp_hour                    &  ! Hourly leaf water potential (MPa)
    )
    !
    ! !DESCRIPTION:
    ! Daily output for diurnal dynamics
    !
    ! !USES:
    use data_structure, only : n_spe
    use Statistics
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: SLA(:)
    integer, intent(in) :: Fn, Max_loc, Max_no
    logical, intent(in) :: tree_exist(:)
    integer, intent(in) :: pft(:), age(:)
    real, intent(in)    :: dbh_heartwood(:), dbh_sapwood(:)
    real, intent(in)    :: mass_leaf(:)
    real(8), intent(in) :: plant_gpp_hour(:,:), plant_sap_hour(:,:)
    real(8), intent(in) :: plant_et_hour(:,:), leaf_wp_hour(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: p
    integer :: no
    integer :: hour
    real(8) :: stat1
    real(8) :: stat2
    real(8), dimension(n_spe, 24) :: gpp_hour_plot
    real(8), dimension(n_spe, 24) :: sap_hour_plot
    real(8), dimension(n_spe, 24) :: et_hour_plot
    real(8), dimension(n_spe)     :: lai_species_plot
    real(8), dimension(24)        :: leaf_wp_hour_plot
    real(8), dimension(n_spe, 24) :: leaf_wp_hour_pft_plot
    real(8), dimension(Max_no)    :: ind_var_pft
    real(8), parameter            :: leaf_wp_blank = 0.0d0
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Diurnal dynamics output
    !---------------------------------------------------------------------

    gpp_hour_plot(:,:) = 0.0d0
    sap_hour_plot(:,:) = 0.0d0
    et_hour_plot(:,:) = 0.0d0
    lai_species_plot(:) = 0.0d0
    leaf_wp_hour_plot(:) = 0.0d0
    leaf_wp_hour_pft_plot(:,:) = 0.0d0

    do hour = 1, 24
       do no = 1, Max_no
          if ( .not. tree_exist(no) ) cycle
          p = pft(no)
          gpp_hour_plot(p, hour) = gpp_hour_plot(p, hour)                      &
                                   +plant_gpp_hour(no, hour)*12.0d0               ! g C/hour in the computational domain
          sap_hour_plot(p, hour) = sap_hour_plot(p, hour)                      &
                                   +plant_sap_hour(no, hour)                      ! m3 H2O/hour in the computational domain
          et_hour_plot (p, hour) = et_hour_plot (p, hour)                      &
                                   +plant_et_hour (no, hour)                      ! m3 H2O/hour in the computational domain
       end do
    end do

    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       p = pft(no)
       lai_species_plot(p) = lai_species_plot(p)+mass_leaf(no)*SLA(p)  ! m2 leaf in the computational domain
    end do

    ! Averaging for each species
    do p = 1, n_spe
       gpp_hour_plot(p, :) = gpp_hour_plot(p, :)/(real(Max_loc)**2.0d0)           ! g C/m2 ground/hour
       sap_hour_plot(p, :) = (sap_hour_plot(p, :)/(real(Max_loc)**2.0d0))      &
                             *1000.0d0                                            ! m/hour -> mm/hour
       et_hour_plot (p, :) = (et_hour_plot (p, :)/(real(Max_loc)**2.0d0))      &
                             *1000.0d0                                            ! m/hour -> mm/hour
       lai_species_plot(p) = lai_species_plot(p)/(real(Max_loc)**2.0d0)           ! m2 leaf/m2 ground
    end do

    ! Statistics of diurnal leaf water potential variation for all species (MPa)
    do hour = 1, 24
       ! Mean value
       call Mean_1d(leaf_wp_hour(:, hour), stat1, leaf_wp_blank)
       ! Median value
       call Median_1d(leaf_wp_hour(:, hour), stat2, leaf_wp_blank)
       ! leaf_wp_hour_plot(hour) = stat1
       leaf_wp_hour_plot(hour) = stat2
    end do

    ! Statistics of diurnal leaf water potential variation per all species (MPa)
    do hour = 1, 24
       do p = 1, n_spe
          ind_var_pft = ind_pft_sample ( &
                                               ! *** Input ***
          leaf_wp_hour(:,hour)            , &  ! Individual variable
          p                               , &  ! Specified pft to extract samples
          Max_no                          , &  ! Maximum number of individual stands
          tree_exist                      , &  ! Flag of tree presence
          pft                             , &  ! Species index (1: Rh, 2: Br)
          age                             , &  ! Tree age (year: 1~)
          dbh_heartwood                   , &  ! Heartwood diameter (m)
          dbh_sapwood                       &  ! Sapwood diameter (m)
          )
          ! Mean value
          call Mean_1d(ind_var_pft, stat1, leaf_wp_blank)
          ! Median value
          call Median_1d(ind_var_pft, stat2, leaf_wp_blank)
          ! leaf_wp_hour_pft_plot(p, hour) = stat1
          leaf_wp_hour_pft_plot(p, hour) = stat2
       end do
    end do

    ! Write diurnal output
    do p = 1, n_spe
       write (Fn, '(f12.5,a)',advance='no') lai_species_plot(p),','
    end do
    do p = 1, n_spe
       do hour = 1, 24
          write (Fn, '(f12.5,a)',advance='no') gpp_hour_plot(p, hour),','
       end do
    end do
    do p = 1, n_spe
       do hour = 1, 24
          write (Fn, '(f12.5,a)',advance='no') sap_hour_plot(p, hour),','
       end do
    end do
    do p = 1, n_spe
       do hour = 1, 24
          write (Fn, '(f12.5,a)',advance='no') et_hour_plot(p, hour),','
       end do
    end do
    do hour = 1, 24
       write (Fn, '(f12.5,a)',advance='no') leaf_wp_hour_plot(hour),','
    end do
    do p = 1, n_spe
       do hour = 1, 24
          if (p < n_spe .or. hour < 24) then
             write (Fn, '(f12.5,a)',advance='no') leaf_wp_hour_pft_plot(p, hour),','
          else
             write (Fn, '(f12.5,a)') leaf_wp_hour_pft_plot(p, hour)
          end if
       end do
    end do

  end subroutine monitor_diurnal

  !-----------------------------------------------------------------------
  function ind_pft_sample ( &
                                       ! *** Input ***
    ind_var                       , &  ! Individual variable
    pft_spe                       , &  ! Specified pft to extract samples
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    age                           , &  ! Tree age (year: 1~)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ind_var_pft)
    !
    ! !DESCRIPTION:
    ! Extract samples for a specified pft.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: ind_var(:)
    integer, intent(in) :: pft_spe
    integer, intent(in) :: Max_no
    logical, intent(in) :: tree_exist(:)
    integer, intent(in) :: pft(:), age(:)
    real, intent(in)    :: dbh_heartwood(:), dbh_sapwood(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p                                  ! Species index (1: Rh, 2: Br)
    integer :: no                                 ! Tree index
    integer :: count                              ! For counting
    real(8), dimension(Max_no) :: ind_var_pft     ! Extracted samples for the specified pft
    !---------------------------------------------------------------------

    ind_var_pft(:) = 0.0d0
    count = 0
    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       if (age(no) <= 2) cycle                                ! Ignore very small trees whose survival is not affected by tree density and production rate.
       if (dbh_heartwood(no) + dbh_sapwood(no) < 0.02) cycle  ! Ignore very small trees
       p = pft(no)
       if (p == pft_spe) then
          count = count + 1
          ind_var_pft(count) = ind_var(no)
       end if
    end do

  end function ind_pft_sample

end module mod_monitoring
