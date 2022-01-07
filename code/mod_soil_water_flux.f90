module mod_soil_water_flux

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
  public :: soil_evap
  public :: soil_water_flux
  !-----------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine soil_evap ( &
                                       ! *** Input ***
    ground_elev                   , &  ! Ground elevation (cm)
    soilresis                     , &  ! Soil evaporative resistance (s/m)
    albsoi_vis                    , &  ! Direct beam and diffuse albedo of ground for VIS (soil)
    albsoi_nir                    , &  ! Direct beam and diffuse albedo of ground for NIR (soil)
    moist                         , &  ! Relative soil moisture (-)
    bsw                           , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat                        , &  ! Soil matric potential at saturation (MPa)
    time_series_tide              , &  ! Dimension of the virtual forest (m)
    time_series_tair              , &  ! Time-series air temperature (K)
    time_series_pair              , &  ! Time-series air pressure (Pa)
    time_series_wind              , &  ! Time-series wind speed (m/s)
    time_series_eair              , &  ! Time-series vapor pressure in air (Pa)
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
    STEP                          , &  ! Canopy layer thickness (m)
    day_of_year                   , &  ! Day of the year (1 - 365)
    ntop                          , &  ! Index for top leaf layer
    albvegb_vis                   , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
    albvegb_nir                   , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
    ir_tree                       , &  ! VIS-interruption-coefficient by tree corwn
    ir_tree_nir                   , &  ! NIR-interruption-coefficient by tree corwn
    irsoi                         , &  ! Absorbed longwave radiation, ground (W/m2)
    irveg                         , &  ! Absorbed longwave radiation, vegetation (W/m2)
                                       !
                                       ! *** Output ***
    rnsoi_day                     , &  ! Daily mean net radiation, ground (W/m2)
    rnveg_day                     , &  ! Daily mean net radiation, vegetation (W/m2)
    etsoi_day                       &  ! Daily soil evaporation rate (mm/day)
    )
    !
    ! !DESCRIPTION:
    ! Soil evaporation rate.
    !
    ! !USES:
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: ground_elev, soilresis, albsoi_vis, albsoi_nir
    real(8), intent(in)  :: moist, bsw, psisat
    real(8), intent(in)  :: time_series_tide(:), time_series_tair(:)
    real(8), intent(in)  :: time_series_pair(:)
    real(8), intent(in)  :: time_series_wind(:), time_series_eair(:)
    real(8), intent(in)  :: time_series_rad_dir(:), time_series_rad_dif(:)
    real(8), intent(in)  :: STEP
    integer, intent(in)  :: day_of_year, ntop
    real(8), intent(in)  :: albvegb_vis, albvegb_nir, ir_tree, ir_tree_nir
    real(8), intent(in)  :: irsoi, irveg
    real(8), intent(out) :: rnsoi_day, rnveg_day, etsoi_day
    !
    ! !LOCAL VARIABLES:
    real(8), parameter :: denh2o = 1000.0d0            ! Density of liquid water (kg/m3)
    real(8), parameter :: grav = 9.80665d0             ! Gravitational acceleration (m/s2)
    real(8), parameter :: mmh2o = 18.02d0 / 1000.0d0   ! Molecular mass of water (kg/mol)is:starred
    real(8), parameter :: mmdry = 28.97d0 / 1000.0d0   ! Molecular mass of dry air (kg/mol)
    real(8), parameter :: rgas = 8.31446d0             ! Universal gas constant (J/K/mol)
    real(8), parameter :: cpd = 1005.0d0               ! Specific heat of dry air at constant pressure (J/kg/K)
    real(8), parameter :: cpw = 1846.0d0               ! Specific heat of water vapor at constant pressure (J/kg/K)
    real(8), parameter :: von_karman = 0.4d0           ! Von Karman constant
    real(8), parameter :: tfrz = 273.15d0              ! Freezing point of water (K)
    real(8) :: albcanb_vis
    real(8) :: albcanb_nir
    real(8) :: canopy_h
    real(8) :: dis_h
    real(8) :: mrough_h
    real(8) :: wrough_h
    real(8) :: head
    real(8) :: mat_psi
    integer :: hour
    integer :: hour_of_year
    real(8) :: tair
    real(8) :: pair
    real(8) :: wind
    real(8) :: eair
    real(8) :: rad_dir
    real(8) :: rad_dif
    real(8) :: swsky_vis
    real(8) :: swsky_nir
    real(8) :: tide
    real(8) :: rhomol
    real(8) :: cp
    real(8) :: qair
    real(8) :: rhoair
    real(8) :: mmair
    real(8) :: tg
    real(8) :: tg_old
    real(8) :: t_soisno
    real(8) :: ga
    real(8) :: gws
    real(8) :: gw
    real(8) :: swsoi
    real(8) :: rnsoi
    real(8) :: swveg
    real(8) :: rnveg
    real(8) :: lambda
    real(8) :: rhg
    real(8) :: gamma
    real(8) :: esat
    real(8) :: desat
    real(8) :: qsat
    real(8) :: dqsat
    real(8) :: num1, num2, num3, num4, den
    real(8) :: shsoi
    real(8) :: eg
    real(8) :: lhsoi
    real(8) :: gsoi
    real(8) :: err
    real(8) :: etsoi
    real(8), parameter :: thk = 1.58d0  ! Soil layer thermal conductivity (W/m/K): Clay soil (100% water content)
    real(8), parameter :: dz = 0.05d0   ! Thickness of soil layer (m)
    !---------------------------------------------------------------------

    ! Current ground temperature
    tg = 298.15d0
    tg_old = tg

    ! Effective canopy albedo including soil
    albcanb_vis = albvegb_vis+(albsoi_vis-albvegb_vis)*ir_tree**2
    albcanb_nir = albvegb_nir+(albsoi_nir-albvegb_nir)*ir_tree_nir**2

    ! Parameters for aerodynamic conductance (Perri et al. 2017)
    canopy_h = real(ntop)*STEP+1.3d0
    dis_h = 0.75d0*canopy_h
    mrough_h = 0.10d0*dis_h
    wrough_h = 0.20d0 * mrough_h

    ! Head of pressure (MPa/m)
    head = denh2o*grav*1.e-06

    ! Soil matrix potential (mm)
    ! Note: It is assumed that soil moisture is always saturated due to tidal water intrusion.
    mat_psi = psisat*moist**(-1.0d0*bsw)   ! Mpa
    mat_psi = mat_psi*head*1000.d0         ! MPa -> m -> mm

    ! Initialization
    rnsoi_day = 0.0d0
    rnveg_day = 0.0d0
    etsoi_day = 0.0d0

    !---------------------------------------------------------------------
    ! Hourly computation
    !---------------------------------------------------------------------

    do hour = 0, 23

       ! Soil temperature (set to the same as soil surface temperature.)
       t_soisno = tg

       ! Set atmospheric forcing parameters
       hour_of_year = (day_of_year-1)*24+hour+1
       tair = time_series_tair(hour_of_year)
       pair = time_series_pair(hour_of_year)
       wind = time_series_wind(hour_of_year)
       eair = time_series_eair(hour_of_year)
       rad_dir = time_series_rad_dir(hour_of_year)
       rad_dif = time_series_rad_dif(hour_of_year)

       ! Atmospheric VIS and NIR solar radiation (W/m2)
       ! Assuming 4% of solar radiation is UV.
       swsky_vis = rad_dir*0.43d0+rad_dif*0.57d0
       swsky_nir = rad_dir*(1.0d0-0.43d0-0.04d0)+rad_dif*(1.0d0-0.57d0-0.04d0)

       !---------------------------------------------------------------------
       ! Set atmospheric related parameters
       !---------------------------------------------------------------------

       qair = mmh2o/mmdry*eair/(pair-(1.0d0-mmh2o/mmdry)*eair)      ! Specific humidity (kg/kg)
       rhomol = pair/(rgas*tair)                                    ! Molar density (mol/m3)
       rhoair = rhomol*mmdry*(1.0d0-(1.0d0-mmh2o/mmdry)*eair/pair)  ! Air density (kg/m3)
       mmair = rhoair/rhomol                                        ! Molecular mass of air (kg/mol)
       cp = cpd*(1.0d0+(cpw/cpd-1.0d0)*qair)*mmair;                 ! Specific heat of air at constant pressure (J/mol/K)
       lambda = 56780.3d0-42.84d0*tair                              ! Latent heat of vaporization
       gamma = cp*pair/lambda                                       ! Psychrometric constant (Pa/K)

       ! Tide level (cm)
       tide = time_series_tide(hour_of_year)

       ! When the ground is exposed
       if (ground_elev >= tide) then

          !---------------------------------------------------------------------
          ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
          !---------------------------------------------------------------------

          ! Sato et al. (2007) forest biome
          ga = (wind*von_karman**2.0d0)/(log(17.4d0))**2.0d0  ! m/s
          ga = ga*rhomol                                      ! m/s -> mol H2O/m2/s

          ! Perri et al. (2017)
          ga = (wind*von_karman**2.0d0)/(log((canopy_h-dis_h)/mrough_h)        &
               *log((canopy_h-dis_h)/wrough_h))               ! m/s
          ga = ga*rhomol                                      ! m/s -> mol H2O/m2/s

          ! Soil conductance to water vapour diffusion
          gws = 1.0d0/soilresis    ! s/m -> m/s
          gws = gws*rhomol         ! m/s -> mol H2O/m2/s
          gw = ga*gws/(ga+gws)     ! Total conductance

          !---------------------------------------------------------------------
          ! Solar radiation absorbed by ground (W/m2)
          !---------------------------------------------------------------------

          ! Goudriaan radiative transfer model
          swsoi = swsky_vis*(1.0d0-albcanb_vis)*ir_tree                        &
                  +swsky_nir*(1.0d0-albcanb_nir)*ir_tree_nir

          ! Net radiation
          rnsoi = swsoi + irsoi

          !---------------------------------------------------------------------
          ! Soil evaporation rate based on energy balance
          ! <= Currently, the energy is not balanced.
          !---------------------------------------------------------------------

          ! Relative humidity in soil airspace
          rhg = exp(grav*mmh2o*mat_psi*1.e-03/(rgas*t_soisno))

          ! Saturation vapor pressure at ground temperature (Pa -> mol/mol)
          call sat_vap (tfrz, tg, esat, desat)
          qsat = esat/pair
          dqsat = desat/pair

          ! Calculate soil surface temperature
          ! Based on Rn = Hg + lambda*E
          num1 = (cp*pair/lambda)*rnsoi
          num2 = cp*rhg*(esat-eair)*gw
          num3 = (cp*pair/lambda)*ga*cp
          tg = (num1-num2)/num3+tair

          ! Calculate soil surface temperature
          ! Based on Rn = Hg + lambda*E + G
          num1 = cp*ga
          num2 = lambda*gw
          num3 = thk/dz
          num4 = rnsoi-num2*rhg*(qsat-dqsat*tg)+num3*t_soisno
          den = num1+num2*dqsat*rhg+num3
          tg = (num1*tair+num2*eair/pair+num4)/den

          ! Sensible heat flux
          shsoi = cp*(tg-tair)*ga

          ! Latent heat flux
          eg = rhg*(esat+desat*(tg-tg_old))
          lhsoi = lambda/pair*(eg-eair)*gw

          ! Soil heat flux
          gsoi = thk*(tg-t_soisno)/dz

          ! Error check

          err = rnsoi-shsoi-lhsoi-gsoi
          if (abs(err) > 0.001d0) then
!             write(*,*) 'ERROR: Soil water flux energy balance error', err
          end if

          ! Water vapor flux: W/m2 -> mol H2O/m2/s
          etsoi = lhsoi/lambda               ! (mol H2O/m2/s)
          etsoi = etsoi*mmh2o                ! (kg H2O/m2/s)
          etsoi = etsoi/denh2o               ! (m3 H2O/m2/s)
          etsoi = max(etsoi, 0.0d0)          ! Prevent negative value. (needed?)

          ! Save current ground temperature
          tg_old = tg

          !---------------------------------------------------------------------
          ! Soil evaporation rate based on Penman-Monteith approach
          !---------------------------------------------------------------------

          ! ga, gws. mol H2O/m2/s -> m/s
          ga = ga/rhomol
          gws = gws/rhomol

          ! cp. J/mol/K -> J/kg/K
          cp = cp/mmair

          ! Saturation vapor pressure at air temperature (Pa)
          call sat_vap (tfrz, tg, esat, desat)

          ! Calculate soil evaporation rate
          num1 = desat*rnsoi
          num2 = cp*rhoair*(esat-eair)*ga
          num3 = lambda*(desat+gamma*(1.0d0+ga/gws))
          etsoi = (num1+num2)/num3
          etsoi = etsoi*mmh2o                ! (kg H2O/m2/s)
          etsoi = etsoi/denh2o               ! (m3 H2O/m2/s)
          etsoi = max(etsoi, 0.0d0)          ! Prevent negative value. (needed?)

       !---------------------------------------------------------------------
       ! When the ground is submerged.
       ! Assumed Rn_soil = 0
       !---------------------------------------------------------------------
       else
          rnsoi = 0.0d0
          tg = 298.15d0
          tg_old = tg
          etsoi = 0.0d0
       end if

       !---------------------------------------------------------------------
       ! Solar radiation absorbed by vegetation (W/m2)
       ! Not affected by tide, but the changes in soil albedo by tidal
       ! water inundation is neglected.
       !---------------------------------------------------------------------

       swveg = swsky_vis*(1.0d0-albcanb_vis)*(1.0d0-ir_tree)                   &
               +swsky_nir*(1.0d0-albcanb_nir)*(1.0d0-ir_tree_nir)

       ! Net radiation
       rnveg = swveg+irveg

       ! Calculating daily mean net radiation, ground and vegetation (W/m2)
       rnsoi_day = rnsoi_day+rnsoi/24.0d0
       rnveg_day = rnveg_day+rnveg/24.0d0

       ! Calculating daily soil evaporation rate (mm/day)
       ! Multiplying 1000 is for (m3 H2O/m2/s) -> (mm H2O/s)
       ! Multiplying 60*60 is for (mm H2O/s) -> (mm H2O/hr)
       ! Summing every hour for (mm H2O/hr) -> (mm H2O/day)
       etsoi_day = etsoi_day+etsoi*1000.0d0*60.0d0*60.0d0

    end do

  end subroutine soil_evap

  !-----------------------------------------------------------------------
  subroutine soil_water_flux ( &
                                       ! *** Input ***
    sal_ini                       , &  ! Initial value of pore-water salinity (psu -> mol/m3)
    din_ini                       , &  ! Initial value of DIN concentration in pore-water (mol N/m3)
    dip_ini                       , &  ! Initial value of DIP concentration in pore-water (mol P/m3)
    flux_fw                       , &  ! Freshwater inputs to soil (mm/hr)
    din_fw                        , &  ! DIN concentration in freshwater (mol N/m3)
    dip_fw                        , &  ! DIP concentration in freshwater (mol P/m3)
    sal_sw                        , &  ! Salinity in seawater (mol/m3)
    din_sw                        , &  ! DIN concentration in seawater (mol N/m3)
    dip_sw                        , &  ! DIP concentration in seawater (mol P/m3)
    root_filter                   , &  ! Filtering rate of salt at root surface (-)
    root_depth                    , &  ! Rooting depth (m)
    flux_calc_switch              , &  ! Switch of soil water flux calculation (1: ON, 0: OFF)
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    plant_water_uptake            , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
    etsoi_day                     , &  ! Daily soil evaporation rate (mm/day)
                                       !
                                       ! *** In/Output ***
    etveg_day                     , &  ! Daily transpiration rate (mm/day) <= output only
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    dip                             &  ! DIP concentration in pore-water (mol P/m3)
    )
    !
    ! !DESCRIPTION:
    ! Box model for soil water flux.
    ! Evaporation from soil is not considered yet. Needs to be added in the future.
    !
    ! !USES:
    use data_structure, only : display_flag
    use time_counter,   only : counter
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: sal_ini, din_ini, dip_ini
    real(8), intent(in)  :: flux_fw, din_fw, dip_fw, sal_sw, din_sw, dip_sw
    real(8), intent(in)  :: root_filter(:), root_depth
    integer, intent(in)  :: flux_calc_switch, Max_loc, Max_no
    logical, intent(in)  :: tree_exist(:)
    integer, intent(in)  :: pft(:)
    real(8), intent(in)  :: mass_leaf(:), plant_water_uptake(:), etsoi_day
    real(8), intent(out) :: etveg_day
    real(8), intent(inout) :: sal, din, dip
    !
    ! !LOCAL VARIABLES:
    real(8) :: box_vol
    real(8) :: sal_box
    real(8) :: din_box
    real(8) :: dip_box
    real(8) :: wat_remove
    real(8) :: sal_remove
    real(8) :: din_remove
    real(8) :: dip_remove
    integer :: no
    real(8) :: wat_input
    real(8) :: din_input
    real(8) :: dip_input
    real(8) :: wat_vol
    real(8) :: wat_recharge
    real(8) :: sal_recharge
    real(8) :: din_recharge
    real(8) :: dip_recharge
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Amounts of salt/N/P in the box (mol)
    !---------------------------------------------------------------------

    ! Volume of the box (m3)
    box_vol=real(Max_loc)*real(Max_loc)*root_depth

    ! salt/N/P in the box (mol)
    sal_box = sal*box_vol
    din_box = din*box_vol
    dip_box = dip*box_vol

    !---------------------------------------------------------------------
    ! Daily water/salt/N/P removal by plants
    !---------------------------------------------------------------------

    ! Computation of all trees
    wat_remove = 0.0d0
    sal_remove = 0.0d0
    din_remove = 0.0d0
    dip_remove = 0.0d0
    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       if ( mass_leaf(no) <= 0.0d0 ) cycle
       ! Daily water uptake by plants (m3 H2O/day)
       wat_remove = wat_remove+plant_water_uptake(no)
       ! Daily salt/N/P uptake by plants (mol/day)
       sal_remove = sal_remove+plant_water_uptake(no)*sal                      &
                    *(1.0d0-root_filter(pft(no)))
       din_remove = din_remove+plant_water_uptake(no)*din
       dip_remove = dip_remove+plant_water_uptake(no)*dip
    end do

    ! Daily transpiration rate (mm/day)
    etveg_day = wat_remove*1000.0d0/(real(Max_loc)**2.0d0)

    !---------------------------------------------------------------------
    ! Add daily soil evaporation rate (m3 H2O/day)
    !---------------------------------------------------------------------

    wat_remove = wat_remove+etsoi_day*1.e-03*real(Max_loc)*real(Max_loc)

    !---------------------------------------------------------------------
    ! Inputs of freshwater/N/P to the box
    !---------------------------------------------------------------------

    ! Freshwater inputs (m3 H2O/day)
    wat_input = flux_fw*real(Max_loc)*real(Max_loc)*24.0d0*1.e-03

    ! N and P inputs (mol/day)
    din_input = din_fw*wat_input
    dip_input = dip_fw*wat_input

    !---------------------------------------------------------------------
    ! Case the water overflows by excess freshwater inputs
    ! 淡水が土の中の塩分・栄養塩濃度を薄めて、そのあと余剰分がボックスから出て行く。
    !---------------------------------------------------------------------

    if (wat_input>wat_remove .and. flux_calc_switch==1) then

       ! No sea water recharge
       wat_recharge = 0.0d0

       ! Water volume just before the overflow (m3 H2O)
       wat_vol = box_vol-wat_remove

       ! Amounts of salt/N/P in the box just before the overflow (mol)
       sal_box = sal_box-sal_remove
       din_box = din_box-din_remove
       dip_box = dip_box-dip_remove

       ! Amounts of salt/N/P in the box after the overflow (mol)
       sal_box = sal_box-sal_box*(wat_input-wat_remove)/box_vol
       din_box = din_box-din_box*(wat_input-wat_remove)/box_vol+din_input
       dip_box = dip_box-dip_box*(wat_input-wat_remove)/box_vol+dip_input

       ! Recalculation of pore-water salinity and nutrient concentrations (mol/m3)
       sal = sal_box/box_vol
       din = din_box/box_vol
       dip = dip_box/box_vol

    !---------------------------------------------------------------------
    ! Case the water is recharged by sea water intrusion
    !---------------------------------------------------------------------
    else if (flux_calc_switch == 1) then

       ! Recharge of water by sea water intrusion (m3 H2O/day)
       wat_recharge = wat_remove-wat_input

       ! Recharge of salt/N/P by sea water intrusion (mol/day)
       sal_recharge = sal_sw*wat_recharge
       din_recharge = din_sw*wat_recharge
       dip_recharge = dip_sw*wat_recharge

       ! Water volume after water recharge
       wat_vol = box_vol-wat_remove+wat_input+wat_recharge

       ! Amounts of salt/N/P in the box after water recharge (mol)
       sal_box = sal_box-sal_remove+sal_recharge
       din_box = din_box-din_remove+din_input+din_recharge
       dip_box = dip_box-dip_remove+dip_input+dip_recharge

       ! Recalculation of pore-water salinity and nutrient concentrations (mol/m3)
       sal = sal_box/wat_vol
       din = din_box/wat_vol
       dip = dip_box/wat_vol

    !---------------------------------------------------------------------
    ! When the flux calculation is off
    !---------------------------------------------------------------------
    else
       sal = sal_ini
       din = din_ini
       dip = dip_ini
    end if

  end subroutine soil_water_flux

end module mod_soil_water_flux
