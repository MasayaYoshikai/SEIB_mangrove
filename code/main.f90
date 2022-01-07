!**************************************************************************
!* Spatially Explicit Individual-Based Dynamic Global Vegetation Model    *
!* SEIB-DGVM ver. 2.81                                                    *
!*                                                                        *
!*   All rights are reserved by Dr. Hisashi SATO (JAMSTEC)                *
!*                                                                        *
!*   Latest information and code using policy can be obtained at          *
!*   http://seib-dgvm.com/                                                *
!**************************************************************************

!*************************************************************************************************
! MAIN SIMULATION LOOP
!*************************************************************************************************
SUBROUTINE main_loop ( LAT, LON, GlobalZone, YearMaxClimate, YearMaxCO2, &
           tmp_air, prec, cloud, wind, humid, tmp_air_range, tmp_soil, &
           aco2_annual, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat )

!_____________ Global Variables
   USE data_structure
   USE time_counter
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current1
   USE grid_status_current2
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   USE mod_meteorology
   USE mod_radiation
   USE mod_soil_water_flux
   USE mod_monitoring
   USE Statistics
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
   implicit none

!_____________ Set Local Parameters
   !Interval for intensive computation of radiation distribution (day)
   integer,parameter::Days_LightComp_Int = 14

   !Default CO2 concentration in the atmosphere @ ppm
   real,parameter::CO2atm_default = 368.0

!_____________ Set Augments
!Coordination
   real   ,intent(IN):: LAT       !latitude  (degree)
   real   ,intent(IN):: LON       !longitude (degree)

!ID number of global zone
   integer,intent(IN):: GlobalZone

!Length of Inputted Climate Data (yr)
   integer,intent(IN):: YearMaxClimate

!Length of Inputted CO2 Data (yr)
   integer,intent(IN):: YearMaxCO2

!Climatic data
   real,dimension(Day_in_Year, YearMaxClimate),intent(IN) :: &
    tmp_air       , & !Surface air temperature (Celcius)
    prec          , & !Precipitation (mm day-1)
    cloud         , & !Cloudness (fraction)
    wind          , & !Wind velocity (m s-1)
    humid         , & !Specific humidity (kg kg-1)
    tmp_air_range     !Daily range of tmp_air (Celcius)

   real,dimension(Day_in_Year, YearMaxClimate, NumSoil),intent(IN):: &
    tmp_soil      !Soil temperature for each layers (Celcius)

!Atomospheric CO2 time-series @ ppm
   real,dimension(YearMaxCO2),intent(IN)::aco2_annual !Atomospheric co2 concentration (ppm)

!Location data
   real,intent(IN):: &
    ALT          ,& !Altitude  (m above MSL)
    Albedo_soil0 ,& !Soil albedo
    W_fi         ,& !Filed capacity   (m3/m3, 0.0 -> 1.0)
    W_wilt       ,& !Wilting point    (m3/m3, 0.0 -> 1.0)
    W_sat        ,& !Saturate point   (m3/m3, 0.0 -> 1.0)
    W_mat           !Matrix potential
!_____________ Set Local variables

!Daily mean meteological variables
   real :: &
    tmp_air_Today       , & !Surface air temperature (Celcius)
    prec_Today          , & !Precipitation           (mm day-1)
    cloud_Today         , & !Cloudness (fraction)
    wind_Today          , & !Wind velocity     (m s-1)
    humid_Today         , & !Wind velocity     (m s-1)
    tmp_air_range_Today     !Daily range of tmp_air  (Celcius)

   real,dimension(NumSoil):: &
    tmp_soil_Today=0.0   !Soil temperature for each layers (Celcius)

!Other variables
   integer year_climate !Current climate year
   integer i, p, no     !Loop counters
   real    x, y         !For General Usage

!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
!For standard output
   integer t1
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   integer flag_tree_presence

   ! Initial of flag of tree presence
   flag_tree_presence = 0
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!_____________ Intialize variables
!Initialize variables
   Call init_value (W_fi, tmp_air(:,1), tmp_soil(:,1,:), prec(:,1))
   Call radiation_seasonal_change (LAT)
!Read spinup files
   Spinup_year  = 0
   if (Flag_spinup_read) then
      open (File_no(1), file=Fn_spnin)
      Call spinup_in (File_no(1))
      close (File_no(1))
   endif

!_____________ Open output files
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   open ( File_no( 1), file = './output/log.txt'          )
!   open ( File_no( 2), file = './output/output.txt'       )
!   open ( File_no( 3), file = './output/climate.txt'      )
!   open ( File_no( 4), file = './output/air.txt'          )
!   open ( File_no( 5), file = './output/radiation.txt'    )
!   open ( File_no( 6), file = './output/water.txt'        )
!   open ( File_no( 7), file = './output/grass.txt'        )
!   open ( File_no( 8), file = './output/lai.txt'          )
!   open ( File_no( 9), file = './output/cflux.txt'        )
!   open ( File_no(10), file = './output/netradiation.txt' )
!   open ( File_no(11), file = './output/wflux.txt'        )
!   open ( File_no(12), file = './output/ld_vertical.txt'  )
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
   open ( File_no(13), file = './output/annual.txt'       )
   open ( File_no(14), file = Fn_forest_out )
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   open ( File_no(15), file = './output/biomass.txt'      )
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   open ( File_no(16), file = Fn_forest2_out )
!   open ( File_no(19), file = './output/forest_daily.csv' )
   open ( File_no(20), file = Fn_monitor_tree_out )
   open ( File_no(21), file = Fn_monitor_biomass_out )
   open ( File_no(23), file = Fn_monitor_plot_out )
   if (diurnal_out_flag == 1) then
      open ( File_no(18), file = Fn_monitor_diurnal_out )
   end if
   if (growth_calc .eqv. .false.) then  ! Only for no-growth calculation
      open ( File_no(24), file = Fn_diurnal )
   end if
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!_____________ Initialize various counters for the daily simulation loop
   doy          = 0
   year         = 1
   year_climate = 1

   if (Flag_spinup_read) then
      counter_begin = 1 + Spinup_year*Day_in_Year
      counter_end   = (Simulation_year+Spinup_year) * Day_in_Year
   else
      counter_begin = 1
      counter_end   = Simulation_year * Day_in_Year
   endif

!*****************************************************************************
!This loop corresponds to a simulation day, which calls subroutines sequentially.
!For readability of the code, I tried to avoid to call another subroutine
!from the subroutine that was called from this loop.
DO counter = counter_begin, counter_end

!_____________ Daily update of field statuses
   if (Logging==1) write (File_no(1),*) 'Daily update of field statuses'

!Time counters1
   doy = doy + 1
   if (doy == Day_in_Year + 1) then
      doy   = 1
      year  = year + 1

      year_climate = year_climate + 1
      if (year_climate==YearMaxClimate+1) year_climate=1
   endif

!Time counters2 (wild fire and phenology related)
   dfl_fire = dfl_fire + 1 !Day from the last fire

   Do p=1, PFT_no
      dfl_leaf_onset(p) = dfl_leaf_onset(p) + 1 !Day from the last leaf onset
      dfl_leaf_shed (p) = dfl_leaf_shed (p) + 1 !Day from the last leaf shedding
   End do

!_____________ Prepare Climatic data for this cycle
!Daily mean meteological properties
   tmp_air_Today       = tmp_air       (doy, year_climate)
   prec_Today          = prec          (doy, year_climate)
   cloud_Today         = cloud         (doy, year_climate)
   wind_Today          = wind          (doy, year_climate)
   humid_Today         = humid         (doy, year_climate)
   tmp_air_range_Today = tmp_air_range (doy, year_climate)

   do i=1, NumSoil
    tmp_soil_Today(i)  = tmp_soil (doy, year_climate, i)
   end do

!_____________ Daily update of metabolic status
!stat_water(1:PFT_no), a vegetation growth limitter due to soil water shortage
   Do p = 1, PFT_no
      x  = 0.0
      y  = 0.0
      no = 0

      do i=1, max(1,RootDepth(p))
      if (tmp_soil_Today(i)<=0.0) cycle
         no = no + 1
         x  =     x + (pool_w(i)/Depth - W_wilt) / max(W_fi-W_wilt, 0.001)    !Method1
         y  = max(y,  (pool_w(i)/Depth - W_wilt) / max(W_fi-W_wilt, 0.001) )  !Method2
      enddo
      x = x / max(1,no)

      x = max(min(x, 1.0), 0.0)
      y = max(min(y, 1.0), 0.0)

      if (Life_type(p)==2) then
         stat_water(p) = y !for larch
      else
         stat_water(p) = x !for other PFTs
      endif

   End do

!Growth suppresss regulator for each tree
!(Suppressed trees do not conduct tree growth)
   if (doy==Day_in_Year) then
      flag_suppress(:) = 0
      do no=1, Max_no
         if ( tree_exist(no) ) then
         if ( mort_regu1(no)<10.0 .and. age(no)>3 ) then
            flag_suppress(no)=1
         endif
         endif
      enddo
   end if

!_____________ Daily update of field records (running recorders, annual means, etc..)
   if (Logging==1) write (File_no(1),*) 'Daily update of field records'

!Daily update of Climate statuses
!Calculate climatic statistics, their running record, and their running mean
   Call stat_climate (prec_Today, tmp_air_Today, tmp_soil_Today)

!Daily update of Carbon pools and Carbon Fluxes statuses
!Reset variables, update their record in array variables, and calculate their running mean
   Call stat_carbon ()

!Daily update of vegetation statuses
!Reset variables, update their record in array variables, and calculate their running mean
   Call stat_vegetation ()

   !mortality related
   if (doy==1) then
      mort_regu1(:) = 0.0 !NPP annual (g / individual)
      mort_regu2(:) = 0.0 !average leaf area of last year (m2/day) {update on: growth_wood}
      mort_regu4(:) = 0.0 !stem diameter increament in last year (m year-1)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
      mort_gpp(:) = 0.0d0
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
   endif

!_____________ DAILY PSYSICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'DAILY PSYSIOLOGICAL PROCESSES'

   !Specify CO2 concentration in the atmosphere
    co2atm = CO2atm_default                                 !Fixed @ Default value
   !co2atm = aco2_annual(max(year_climate+51, YearMaxCO2))  !Same as year_climate

   !Calculate variables in relation to atmospheric physics (e.g. Air pressure, Vapor density)
   Call air (tmp_air_Today, humid_Today, ALT)

   Call radiation (LAT, cloud_Today)   !!!<<<<<<<<<<<<<<<<<<TN changed: cloud -> cloud_Today
   par_RunningRecord(2:Day_in_Year) = par_RunningRecord(1:Day_in_Year-1)
   par_RunningRecord(1)             = par

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   ! Decide a tree to monitor
   if (flag_tree_presence == 1 .and. year >= 2) then
      call monitor_index ( &
                                           ! *** Input ***
      flag_tree_presence              , &  ! Flag of tree presence
      tree_exist                      , &  ! Flag of tree presence
                                           !
                                           ! *** In/Output ***
      monitor                           &  ! Tree index for monitoring
      )
   end if
   ! Scattering adjusted light extinction coefficient for VIS/NIR
   call extinction_coefficient ( &
                                        ! *** Input ***
   n_spe                           , &  ! Number of mangrove species
   rhol_vis                        , &  ! Leaf relflectance to VIS (-)
   rhol_nir                        , &  ! Leaf relflectance to NIR (-)
   taul_vis                        , &  ! Leaf transmittance to VIS (-)
   taul_nir                        , &  ! Leaf transmittance to NIR (-)
   xl                              , &  ! Departure of leaf angle from spherical orientation (-)
   Max_loc                         , &  ! Dimension of the virtual forest (m)
   sl_hgt(doy)                     , &  ! Solar elevation angle at midday (degree)
   tree_exist                      , &  ! Flag of tree presence
   dble(la)                        , &  ! Leaf area per tree (m2 leaf/tree)
   height                          , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
   bole                            , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
   pft                             , &  ! Species index (1: Rh, 2: Br)
                                        !
                                        ! *** Output ***
   dpai_layer                      , &  ! Layer leaf area index for each PFT for each forest layer (m2 leaf/me ground)
   dpai_layer_sum                  , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
   kbm_vis                         , &  ! Scattering adjusted light extinction coefficient for VIS (-)
   kbm_nir                         , &  ! Scattering adjusted light extinction coefficient for NIR (-)
   albvegb_vis                     , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
   albvegb_nir                     , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
   nbot                            , &  ! Index for bottom leaf layer
   ntop                              &  ! Index for top leaf layer
   )
   ! Read growth temperature when temperature acclimation is active.
   call t_growth_read ( &
                                        ! *** Input ***
   n_spe                           , &  ! Number of mangrove species
   t_acclim                        , &  ! Flag of temperature acclimation
   time_series_tair                , &  ! Time-series air temperature (K)
   doy                             , &  ! Day of the year (1 - 365)
                                        !
                                        ! *** In/Output
   t_growth                          &  ! Growth temperature (degree)
   )
   ! Calculate wind profile within canopy.
   call wind_profile ( &
                                        ! *** Input ***
   n_spe                           , &  ! Number of mangrove species
   dleaf                           , &  ! Leaf dimension (m)
   dble(STEP)                      , &  ! Canopy layer thickness (m)
   dpai_layer_sum                  , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
   nbot                            , &  ! Index for bottom leaf layer
   ntop                            , &  ! Index for top leaf layer
                                        !
                                        ! *** In/Output
   rwind_profile                     &  ! Relative wind speed profile to the canopy top (-)
   )
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

   Call diffused_radiation ()

   if (mod(doy, Days_LightComp_Int)==1) then
      !Compute relative intensity of direct radiation for each crown disk of each tree.
!      write(*,*) 'Calc. direct radiation' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
      Call direct_radiation ()

      !Compute relative intensity of radiation for each cell of tree establishment.
!      write(*,*) 'Calc. floor radiation' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
!      Call floor_radiation ()!!!>>>>>>>>>>>>TN:rm
      Call floor_radiation2 () !!!<<<<<<<<<<<<TN:add
!      Call crown_coverage ()!!!>>>>>>>>>>>>TN:rm
!      write(*,*) 'fin' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   endif

   !Calculate coefficients ir_tree and ir_grass, which describe fractions of incoming radiation
   !absorbed by crown layer and grass layer, respectively.
   !This subroutine also calculates a variable par_grass,
   !intensity of photosynthesis active radiation on the top of grass layer.
   Call ir_index ()

!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!   !Wild fire subroutines
!   if     (GlobalZone==1) then
!     !African continent
!     Call fire_regime2 (wind_Today)   !!!<<<<<<<<<<<<<<<<<TN changed: wind -> wind_Today
!   else
!     !Default
!     Call fire_regime (W_fi)
!   endif
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm

!*****************************************************************************

!_____________ DAILY BIOLOGICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'DAILY BIOLOGICAL PROCESSES'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!write(*,*) 'photosynthesis_condition' !!!<<<<<<<<<<<<TN:add
!   !Calculate photosynthesis rate controlling variables
!   Call photosynthesis_condition (tmp_air_Today)
!
!write(*,*) 'photosynthesis' !!!<<<<<<<<<<<<TN:add
!   !Photosynthesis process
!   Call photosynthesis ()
!
!write(*,*) 'lai_optimum' !!!<<<<<<<<<<<<TN:add
!   !Calculate optimal leaf area index for grass layer
!   Call lai_optimum (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
!
!write(*,*) 'maintenance_resp' !!!<<<<<<<<<<<<TN:add
!   !Maintenance respiration
!   Call maintenance_resp (tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )
!
!write(*,*) 'turnover' !!!<<<<<<<<<<<<TN:add
!   !Turnover of plant organs
!   Call turnover ()
!
!write(*,*) 'leaf_season' !!!<<<<<<<<<<<<TN:add
!   !Penology controller, and biological processes before and after penology change
!   !(such as gradual release of stock biomass after onset of foliage phase).
!   Call leaf_season (LAT, tmp_soil_Today )
!
!write(*,*) 'growth_wood' !!!<<<<<<<<<<<<TN:add
!   !Daily growth procedure for woody PFTs (foliation, fine root growth, and reproduction)
!   Call growth_wood ()
!
!write(*,*) 'growth_grass' !!!<<<<<<<<<<<<TN:add
!   !Growth and reproduction procedure for grass PFTs.
!   Call growth_grass ()
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

!write(*,*) 'spatial_limitation' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   !Calculate free space around each individual tree.
   !This information will be employed by following subroutine growth_trunc,
   !where stem and crown expansion occurs.
   Call spatial_limitation ()

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:moved from below
   !Compute Albedo
   Call albedo_calc (Albedo_soil0)

   !Compute net radiation
   Call net_radiation (tmp_air_Today, cloud_Today)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:moved from below

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   ! Initialization of stand-level diurnal variables
   plant_gpp_hour(:,:) = 0.0d0
   plant_sap_hour(:,:) = 0.0d0
   plant_et_hour (:,:) = 0.0d0
   leaf_wp_hour  (:,:) = 0.0d0

   !---------------------------------------------------------------------
   ! Daily production
   !---------------------------------------------------------------------

   call daily_production ( &
   )

   !---------------------------------------------------------------------
   ! Daily soil evaporation rate
   !---------------------------------------------------------------------

   call soil_evap ( &
                                        ! *** Input ***
   ground_elev                     , &  ! Ground elevation (cm)
   soilresis                       , &  ! Soil evaporative resistance (s/m)
   albsoi_vis                      , &  ! Direct beam and diffuse albedo of ground for VIS (soil)
   albsoi_nir                      , &  ! Direct beam and diffuse albedo of ground for NIR (soil)
   moist                           , &  ! Relative soil moisture (-)
   bsw                             , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
   psisat                          , &  ! Soil matric potential at saturation (MPa)
   time_series_tide                , &  ! Dimension of the virtual forest (m)
   time_series_tair                , &  ! Time-series air temperature (K)
   time_series_pair                , &  ! Time-series air pressure (Pa)
   time_series_wind                , &  ! Time-series wind speed (m/s)
   time_series_eair                , &  ! Time-series vapor pressure in air (Pa)
   time_series_rad_dir             , &  ! Time-series direct radiation at canopy top (W/m2)
   time_series_rad_dif             , &  ! Time-series diffused radiation at canopy top (W/m2)
   dble(STEP)                      , &  ! Canopy layer thickness (m)
   doy                             , &  ! Day of the year (1 - 365)
   ntop                            , &  ! Index for top leaf layer
   dble(albvegb_vis)               , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
   dble(albvegb_nir)               , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
   dble(ir_tree)                   , &  ! VIS-interruption-coefficient by tree corwn
   dble(ir_tree_nir)               , &  ! NIR-interruption-coefficient by tree corwn
   irsoi                           , &  ! Absorbed longwave radiation, ground (W/m2)
   irveg                           , &  ! Absorbed longwave radiation, vegetation (W/m2)
                                        !
                                        ! *** Output ***
   rnsoi_day                       , &  ! Daily mean net radiation, ground (W/m2)
   rnveg_day                       , &  ! Daily mean net radiation, vegetation (W/m2)
   etsoi_day                         &  ! Daily soil evaporation rate (mm/day)
   )

   !---------------------------------------------------------------------
   ! Box model for soil water flux
   !---------------------------------------------------------------------

   ! Switch of soil water flux calculation (1: ON, 0: OFF)
   ! Only after tree establishment
   flux_calc_switch = 0
   if (flux_calculation == 1 .and. flag_tree_presence == 1) then
      flux_calc_switch = 1
   end if

   call soil_water_flux ( &
                                        ! *** Input ***
   sal_ini                         , &  ! Initial value of pore-water salinity (psu -> mol/m3)
   din_ini                         , &  ! Initial value of DIN concentration in pore-water (mol N/m3)
   dip_ini                         , &  ! Initial value of DIP concentration in pore-water (mol P/m3)
   flux_fw                         , &  ! Freshwater inputs to soil (mm/hr)
   din_fw                          , &  ! DIN concentration in freshwater (mol N/m3)
   dip_fw                          , &  ! DIP concentration in freshwater (mol P/m3)
   sal_sw                          , &  ! Salinity in seawater (mol/m3)
   din_sw                          , &  ! DIN concentration in seawater (mol N/m3)
   dip_sw                          , &  ! DIP concentration in seawater (mol P/m3)
   root_filter                     , &  ! Filtering rate of salt at root surface (-)
   root_depth                      , &  ! Rooting depth (m)
   flux_calc_switch                , &  ! Switch of soil water flux calculation (1: ON, 0: OFF)
   Max_loc                         , &  ! Dimension of the virtual forest (m)
   Max_no                          , &  ! Maximum number of individual stands
   tree_exist                      , &  ! Flag of tree presence
   pft                             , &  ! Species index (1: Rh, 2: Br)
   dble(mass_leaf)                 , &  ! Stand-level leaf biomass (g leaf/tree)
   plant_water_uptake              , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
   etsoi_day                       , &  ! Daily soil evaporation rate (mm/day)
                                        !
                                        ! *** In/Output ***
   etveg_day                       , &  ! Daily transpiration rate (mm/day)
   sal                             , &  ! Pore-water salinity (mol/m3)
   din                             , &  ! DIN concentration in pore-water (mol N/m3)
   dip                               &  ! DIP concentration in pore-water (mol P/m3)
   )

   !---------------------------------------------------------------------
   ! Daily output for forest structure
   !---------------------------------------------------------------------

!   call output_forest_daily (File_no(19))

   !---------------------------------------------------------------------
   ! Daily output for plot-scale variables
   !---------------------------------------------------------------------

   call monitor_plot ( &
                                        ! *** Input ***
   dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
   time_series_rad_dir             , &  ! Time-series direct radiation at canopy top (W/m2)
   time_series_rad_dif             , &  ! Time-series diffused radiation at canopy top (W/m2)
   File_no(21)                     , &  ! File number for tree biomass output
   File_no(23)                     , &  ! File number for plot output
   year                            , &  ! Year
   doy                             , &  ! Day of the year (1 - 365)
   Max_loc                         , &  ! Dimension of the virtual forest (m)
   Max_hgt                         , &  ! Number of canopy layer
   Max_no                          , &  ! Maximum number of individual stands
   tree_exist                      , &  ! Flag of tree presence
   pft                             , &  ! Species index (1: Rh, 2: Br)
   age                             , &  ! Tree age (year: 1~)
   tree_h                          , &  ! Tree height (m)
   dbh_heartwood                   , &  ! Heartwood diameter (m)
   dbh_sapwood                     , &  ! Sapwood diameter (m)
   mass_trunk                      , &  ! Stand-level trunk biomass (g trunk/tree)
   mass_coarse_root                , &  ! Stand-level coarse root biomass (g/tree)
   mass_above_root                 , &  ! Stand-level above-ground root biomass (g/tree)
   mass_leaf                       , &  ! Stand-level leaf biomass (g leaf/tree)
   mass_root                       , &  ! Stand-level fine root biomass (g root/tree)
   dpai_layer_sum                  , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
   leaf_mdwp                       , &  ! Midday leaf water potential of the day (MPa)
   leaf_pdwp                       , &  ! Pre-dawn leaf water potential of the last day (MPa)
   plant_water_uptake              , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
   plant_sh_day                    , &  ! Daily sensible heat flux at stand-level (MJ/tree/day)
   plant_lh_day                    , &  ! Daily latent heat flux at stand-level (MJ/tree/day)
   plant_gpp_day                   , &  ! Daily gross photosynthesis at stand-level (mol C/tree/day)
   plant_npp_day                   , &  ! Daily net photosynthesis at stand-level (mol C/tree/day)
   plant_wnpp_day                  , &  ! Daily woody net primary production of each tree (g DW/tree/day)
   plant_fnpp_day                  , &  ! Daily foliage net primary production of each tree (g DW/tree/day)
   irveg                           , &  ! Absorbed longwave radiation, vegetation (W/m2)
   irsoi                           , &  ! Absorbed longwave radiation, ground (W/m2)
   ir_tree                         , &  ! VIS-interruption-coefficient by tree corwn
   rnsoi_day                       , &  ! Daily mean net radiation, ground (W/m2)
   rnveg_day                       , &  ! Daily mean net radiation, vegetation (W/m2)
   etsoi_day                       , &  ! Daily soil evaporation rate (mm/day)
   etveg_day                       , &  ! Daily transpiration rate (mm/day)
   sal                             , &  ! Pore-water salinity (mol/m3)
   din                             , &  ! DIN concentration in pore-water (mol N/m3)
   dip                             , &  ! DIP concentration in pore-water (mol P/m3)
   tree_n_est                      , &  ! Total number of established trees of each PFT in the entire simulation (trees)
   tree_n_die                        &  ! Total number of died trees of each PFT in the entire simulation (trees)
   )

   !---------------------------------------------------------------------
   ! Daily output for diurnal dynamics
   !---------------------------------------------------------------------

   if (diurnal_out_flag == 1) then
      call monitor_diurnal ( &
                                           ! *** Input ***
      dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
      File_no(18)                     , &  ! File number for diurnal output
      Max_loc                         , &  ! Dimension of the virtual forest (m)
      Max_no                          , &  ! Maximum number of individual stands
      tree_exist                      , &  ! Flag of tree presence
      pft                             , &  ! Species index (1: Rh, 2: Br)
      age                             , &  ! Tree age (year: 1~)
      dbh_heartwood                   , &  ! Heartwood diameter (m)
      dbh_sapwood                     , &  ! Sapwood diameter (m)
      mass_leaf                       , &  ! Stand-level leaf biomass (g leaf/tree)
      plant_gpp_hour                  , &  ! Hourly gross photosynthesis at stand-level (mol C/tree/hour)
      plant_sap_hour                  , &  ! Hourly sapflow at stand-level (m3 H2O/tree/day)
      plant_et_hour                   , &  ! Hourly transpiration at stand-level (m3 H2O/tree/day)
      leaf_wp_hour                      &  ! Hourly leaf water potential (MPa)
      )
   end if

   !---------------------------------------------------------------------
   ! Daily output for monitoring tree
   !---------------------------------------------------------------------

   if (doy == Day_in_Year) then ! Output day filter
   if (monitor > 0) then
      call monitor_tree ( &
                                           ! *** Input ***
      time_series_rad_dir             , &  ! Time-series direct radiation at canopy top (W/m2)
      time_series_rad_dif             , &  ! Time-series diffused radiation at canopy top (W/m2)
      File_no(20)                     , &  ! File number for monitoring tree
      STEP                            , &  ! Canopy layer thickness (m)
      year                            , &  ! Year
      doy                             , &  ! Day of the year (1 - 365)
      pft (monitor)                   , &  ! Species index (1: Rh, 2: Br)
      bottom_layer_monitor            , &  ! Bottom layer of canopy
      top_layer_monitor               , &  ! Top layer of canopy
      par_direct_rel (monitor,:)      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
      par_diffuse_rel                 , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
      leaf_wp_monitor                 , &  ! Midday leaf water potential (MPa)
      leaf_wpmin_monitor              , &  ! Daily minimum leaf water potential (MPa)
      gs_top_monitor                  , &  ! Midday top layer stomatal conductance (mol H2O/m2 leaf/s)
      an_top_monitor                  , &  ! Midday top layer leaf net photosynthesis (umol CO2/m2 leaf/s)
      an_bot_monitor                  , &  ! Midday bottom layer leaf net photosynthesis (umol CO2/m2 leaf/s)
      et_top_monitor                  , &  ! Midday top layer leaf transpiration rate (mol H2O/m2 leaf/s)
      et_bot_monitor                  , &  ! Midday bottom layer leaf transpiration rate (mol H2O/m2 leaf/s)
      c_uptake_bottom_day_monitor     , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
      n_uptake_bottom_day_monitor     , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
      water_uptake_day_monitor        , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
      available_c_monitor             , &  ! Remaining carbon for tree growth after respiration and turnover (mol C/tree/day)
      available_n_monitor             , &  ! Remaining nitrogen for tree growth after turnover (mol N/tree/day)
      height_max_monitor              , &  ! Potential tree height based on allometric relation (m)
      height_limit (monitor)          , &  ! Tree height limitation based on proximate trees (unit is STEP!!)
      crown_d_max_monitor             , &  ! Potential crown diameter based on allometric relation (m)
      radius_limit (monitor)          , &  ! Crown radius limitation based on proximate trees (m)
      dbh_heartwood (monitor)         , &  ! Heartwood diameter (m)
      dbh_sapwood (monitor)           , &  ! Sapwood diameter (m)
      tree_h (monitor)                , &  ! Tree height (m)
      crown_diameter (monitor)        , &  ! Crown diameter (m)
      height (monitor)                , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
      bole (monitor)                  , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
      lai_monitor                     , &  ! Leaf area index (m2 leaf/m2 ground)
      mass_trunk (monitor)            , &  ! Stand-level trunk biomass (g trunk/tree)
      mass_coarse_root (monitor)      , &  ! Stand-level coarse root biomass (g/tree)
      mass_above_root (monitor)       , &  ! Stand-level above-ground root biomass (g/tree)
      mass_leaf (monitor)             , &  ! Stand-level leaf biomass (g leaf/tree)
      mass_root (monitor)             , &  ! Stand-level fine root biomass (g root/tree)
      mass_stock (monitor)              &  ! Stand-level stock biomass (g stock/tree)
      )
   end if
   end if ! Output day filter
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!write(*,*) 'growth_trunk' !!!<<<<<<<<<<<<TN:add
!   !Expands stem diameter, which cause stem height and crown diameter increases
!   Call growth_trunk ()
!
!write(*,*) 'decomposition' !!!<<<<<<<<<<<<TN:add
!   !Heterotrophic respiration (Litter and soil organic matter decomposition)
!   Call decomposition (W_fi)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

!_____________ ANNUALY BIOLOGICAL PROCESSES
   if (Logging==1) write (File_no(1),*) 'ANNUALY BIOLOGICAL PROCESSES'

IF ( doy==1 ) then
   !Listing all PFT that exist in this virtual forest
   Call pft_present ()
ENDIF

IF ( doy==Day_in_Year ) then

   !Determine which tree die
!   write(*,*) 'Calc. mortality' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   Call mortality ()

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   !Adjust crown depth by purging crown disks from the bottom of the crown layer
!   write(*,*) 'Calc. crown_adjust' !!!<<<<<<<<<<<<TN:add
!   Call crown_adjust ()
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

   !Cause horizontal movement of crown of each tree to the open direction
!   write(*,*) 'Calc. crown_shake' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   Call crown_shake ()

   !Establishment of woody PFTs
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   x = sum(prec_RunningRecord(1:Day_in_Year))                  !annual precipitation (mm/year)
!   y = sum(tmp_air_RunningRecord(1:Day_in_Year)) / Day_in_Year !annual mean temperature (Celcius)
!   if ( x > max(100.0, 20.0*y) ) then                          !from Koppen's criteria (1936)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
      !Determine ground mesh point, where newly saplings can establish.
      !This subroutine should be called before subroutine 'establish'.
!      write(*,*) 'Calc. ground_vacant' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
      Call ground_vacant ()

      !Recruitment of trees
!      write(*,*) 'Calc. establish' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
      Call establish (GlobalZone)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
      ! After first establishment
      flag_tree_presence = 1
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   end if
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

   !Compute fraction of tree crown coverage.
!   write(*,*) 'Calc. crown_coverage' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
!   Call crown_coverage ()

   !Determine dominant grass PFT in the next year
!   write(*,*) 'Calc. grass_priority' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   Call grass_priority ()

   !Determine current biome type of this grid cell
!   write(*,*) 'Calc. biome_determine' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
   Call biome_determine ()

END IF

!write(*,*) 'PASS:05' !!!<<<<<<<<<<<<TN:add <*********** rm
!_____________ DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS
   if (Logging==1) write (File_no(1),*) 'DAILY UPDATE NET-RADIATION & SOIL-WATER STATUS'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:moved above
!   !Compute Albedo
!   Call albedo_calc (Albedo_soil0)
!
!   !Compute net radiation
!   Call net_radiation (tmp_air_Today, cloud_Today)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:moved above

   !Compute water cycle
   Call waterbudget &
   (W_fi, W_wilt, prec_Today, wind_Today, tmp_air_Today, tmp_soil_Today(:) )

!_____________ DAILY INCREMENT of LITTER AND FUEL POOLS
   ! Increment of litters pools and fuel pools are basically conducted here
   ! by considering litter fluxes, which are computed in descendent subroutines.
   ! This way of convergent-computation is for handling
   ! a variety of litter pools and fuel pools in the model.

   pool_litter_leaf  = max(0.0, pool_litter_leaf  + flux_litter_leaf )
   pool_litter_trunk = max(0.0, pool_litter_trunk + flux_litter_trunk)
   pool_litter_root  = max(0.0, pool_litter_root  + flux_litter_root )
   pool_litter_ag    = max(0.0, pool_litter_ag    + flux_litter_ag   )
   pool_litter_bg    = max(0.0, pool_litter_bg    + flux_litter_bg   )

   pool_fuel_standT  = max(0.0, pool_fuel_standT  + flux_litter_leaf )
   pool_fuel_standG  = max(0.0, pool_fuel_standG  + flux_litter_ag   )

!write(*,*) 'PASS:06' !!!<<<<<<<<<<<<TN:add >>>>>>>>> MY:rm
!_____________ MAKE OUTPUT FILES (Write simulation results in output files)
IF (Flag_output_write) then

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   !Daily output (everyday)
!   Call output_for_viewer &
!   (File_no(2), LAT, LON, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat, &
!   cloud_Today, prec_Today, humid_Today, wind_Today, tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0 )   !!!<<<<<<<<<<<<<<<<<<TN changed: cloud -> cloud_Today
!
!   !Daily output (on each specifyed year)
!   Call output_climate &
!   (File_no(3), &
!   cloud_Today, prec_Today, humid_Today, wind_Today, tmp_air_Today, sum(tmp_soil_Today(1:5))/5.0)   !!!<<<<<<<<<<<<<<<<<<TN changed: cloud -> cloud_Today
!   Call output_air          (File_no(4))
!   Call output_radiation    (File_no(5))
!   Call output_water        (File_no(6), W_fi, W_wilt, W_sat)
!   Call output_grass        (File_no(7))
!   Call output_lai          (File_no(8))
!   Call output_cflux        (File_no(9))
!   Call output_netradiation (File_no(10))
!   Call output_wflux        (File_no(11), prec_Today)
!
!   !Monthly output
!   if ( Day_of_Month(doy) == Day_in_month(Month(doy)) ) then
!      Call output_ld_vertical (File_no(12))
!   endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

   !Annual output (At the end of the year)
   if (doy == Day_in_Year) then
!  if (doy == 213) then !@ 1 Aug
      Call output_annual (File_no(13))
      Call output_forest (File_no(14))
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
      Call output_forest2 (File_no(16))
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
     Call system_clock(t1)
!     write(*,*) t1 >>>>>>>>> MY:rm
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   !Annual output (At the middle of gwoeing season in Northern Hemisphere)
!   if (doy == 195) then
!      Call output_biomass(File_no(15))
!   endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm

ENDIF
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add >>>>>>>>> MY:rm
! write(*,*) year, doy
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add >>>>>>>>> MY:rm

END DO !The end of daily loop
!*****************************************************************************

!_____________ After Main-simulation-loop procedures
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!   close (File_no( 1))   !output_log.txt
!   close (File_no( 2))   !output.txt
!   close (File_no( 3))   !output_climate.txt
!   close (File_no( 4))   !output_air.txt
!   close (File_no( 5))   !output_radiation.txt
!   close (File_no( 6))   !output_water.txt
!   close (File_no( 7))   !output_grass.txt
!   close (File_no( 8))   !output_lai.txt
!   close (File_no( 9))   !output_cflux.txt
!   close (File_no(10))   !output_netradiation.txt
!   close (File_no(11))   !output_wflux.txt
!   close (File_no(12))   !output_ld_vertical.txt
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
   close (File_no(13))   !output_annual.txt
   close (File_no(14))   !output_forest.txt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!   close (File_no(15))   !output_biomass.txt
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   close (File_no(16))
!   close (File_no(19))
   close (File_no(20))
   close (File_no(21))
   close (File_no(23))
   if (diurnal_out_flag == 1) then
   close (File_no(18))
   end if
   if (growth_calc .eqv. .false.) then
   close (File_no(24))
   end if
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

IF ( Flag_spinup_write ) then
   !Write spinup data for restarting simualtion
   open(File_no(1), file = Fn_spnout)
   Call spinup_out (File_no(1))
   close (File_no(1))
END IF

END SUBROUTINE main_loop
