module mod_metabolism

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
  ! !PUBLIC MEMBER FUNCTIONS
  public  :: woody_respiration
  public  :: turnover
  private :: turnover_x
  !
  ! !PRIVATE MEMBER FUNCTIONS
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine woody_respiration ( &
                                       ! *** Input ***
    t_acclim                      , &  ! Flag of temperature acclimation (1: Use, 0: No use)
    main_resp_stem                , &  ! Maintenance stem respiration rate at 15 degree (day-1)
    main_resp_root                , &  ! Maintenance stem respiration rate at 15 degree (day-1)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    t_growth                      , &  ! Growth temperature (degree), mean temperature for 30 days before
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    day_of_year                   , &  ! Day of the year (1 - 365)
    p                             , &  ! Species index (1: Rh, 2: Br)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    above_root_biomass            , &  ! Stand-level above-ground root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    time_series_tair              , &  ! Time-series air temperature (K)
                                       !
                                       ! *** Output ***
    resp_c_cost                   , &  ! Woody respiration cost (mol C/tree/day)
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    remaining_c                     &  ! Remaining carbon for tree growth (mol C/tree/day)
    )
    !
    ! !Description:
    ! Calculate carbon consumption by stem and root respiration
    ! Step 1: available_c is used
    ! Step 2: stock_biomass is used
    ! Step 3: A portion of leaf and root is dropped (as in SEIB-DGVM)
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    integer, intent(in)  :: t_acclim(:)
    real(8), intent(in)  :: main_resp_stem(:), main_resp_root(:)
    real(8), intent(in)  :: C_in_drymass, t_growth(:), available_c
    integer, intent(in)  :: day_of_year, p
    real(8), intent(in)  :: trunk_biomass, above_root_biomass, leaf_biomass
    real(8), intent(in)  :: root_biomass, coarse_root_biomass, stock_c
    real(8), intent(in)  :: time_series_tair(:)
    real(8), intent(out) :: resp_c_cost, d_leaf, d_root, d_stock_c, remaining_c
    !
    ! !Local variables:
    integer :: hour_of_year
    real(8) :: tmp_air_ref
    real(8) :: tmp_sensitiv_air
    real(8) :: qt
    real(8) :: required_c
    real(8) :: required_c_remaining
    real(8), parameter :: tfrz = 273.15d0
    real(8), parameter :: frac_organ_remove = 0.01d0
    !---------------------------------------------------------------------

    if (t_acclim(p) == 1) then
       tmp_air_ref = t_growth(p)
    else
       hour_of_year = (day_of_year-1)*24+12+1
       tmp_air_ref = time_series_tair(hour_of_year)-tfrz
    end if
    qt = 2.0d0*exp(-0.009d0*(tmp_air_ref-15.0d0))
    tmp_sensitiv_air = exp((tmp_air_ref-15.0d0)*log(qt)/10.0d0)
    remaining_c = available_c*12.0d0
    required_c = main_resp_stem(p)*(trunk_biomass+above_root_biomass           &
                 +coarse_root_biomass)*tmp_sensitiv_air                        &
                 +main_resp_root(p)*root_biomass*tmp_sensitiv_air
    required_c = required_c*C_in_drymass
    resp_c_cost = required_c/12.0d0
    if (remaining_c < 0.0d0) then
       required_c = required_c+abs(remaining_c)
       remaining_c = 0.0d0
    end if

    if (remaining_c > required_c) then
       remaining_c = remaining_c-required_c
       required_c_remaining = 0.0d0
       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock_c = 0.0d0
    else
       required_c_remaining = required_c-remaining_c
       remaining_c = 0.0d0
       if (stock_c*12.0d0 > required_c_remaining) then
          d_stock_c = -1.0d0*required_c_remaining/12.0d0
          required_c_remaining = 0.0d0
          d_leaf = 0.0d0
          d_root = 0.0d0
       else
          d_stock_c = -1.0d0*stock_c
          required_c_remaining = required_c_remaining-stock_c*12.0d0
          d_leaf = -1.0d0*frac_organ_remove*leaf_biomass
          d_root = -1.0d0*frac_organ_remove*root_biomass
          required_c_remaining = 0.0d0
       end if
    end if

    ! Unit conversion (g C/tree/day) -> (mol C/tree/day)
    remaining_c = remaining_c / 12.0d0

    if (required_c_remaining > 0.0d0) then
       write(*,*) 'ERROR: required carbon for respiration is still remaining'
    end if

  end subroutine woody_respiration

  !-----------------------------------------------------------------------
  subroutine turnover ( &
                                       ! *** Input ***
    coarse_root_turn              , &  ! Coarse root turnover rate (day-1)
    leaf_turn                     , &  ! Leaf turnover rate (day-1)
    root_turn                     , &  ! Root turnover rate (day-1)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    c_n_root                      , &  ! C/N ratio in mol in root
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    stress_flag                   , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    p                             , &  ! Species index (1: Rh, 2: Br)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    stock_n                       , &  ! Nitrogen in stock (mol N/tree)
                                       !
                                       ! *** Output ***
    d_coarse_root                 , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    d_stock_n                     , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
    remaining_c                   , &  ! Remaining carbon for tree growth (mol C/tree/day)
    remaining_n                     &  ! Remaining nitrogen for tree growth (mol N/tree/day)
    )
    !
    ! !Description:
    ! Coarse root, leaf, and root turnover.
    ! Firstly, turnover of coarse root is compensated.
    ! Then, leaf or root is compensated depending on nutrient limitation
    ! and water uptake ability.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: coarse_root_turn(:), leaf_turn(:), root_turn(:)
    real(8), intent(in)  :: c_n_stem(:), c_n_leaf(:), c_n_root(:), C_in_drymass
    integer, intent(in)  :: stress_flag
    real(8), intent(in)  :: available_c, available_n
    integer, intent(in)  :: p
    real(8), intent(in)  :: coarse_root_biomass, leaf_biomass, root_biomass
    real(8), intent(in)  :: stock_c, stock_n
    real(8), intent(out) :: d_coarse_root, d_leaf, d_root, d_stock_c, d_stock_n
    real(8), intent(out) :: remaining_c, remaining_n
    !
    ! !Local variables:
    integer :: n_limit_flag
    integer :: root_flag
    real(8) :: coarse_root_c_required
    real(8) :: root_c_required
    real(8) :: leaf_c_required
    real(8) :: stock_c_remaining
    real(8) :: stock_n_remaining
    !---------------------------------------------------------------------

    remaining_c = available_c
    remaining_n = available_n
    coarse_root_c_required = coarse_root_turn(p)*coarse_root_biomass           &
                             *C_in_drymass/12.0d0
    root_c_required = root_turn(p)*root_biomass*C_in_drymass/12.0d0
    leaf_c_required = leaf_turn(p)*leaf_biomass*C_in_drymass/12.0d0
    coarse_root_c_required = max(coarse_root_c_required, 0.0d0)
    root_c_required = max(root_c_required, 0.0d0)
    leaf_c_required = max(leaf_c_required, 0.0d0)
    d_coarse_root = 0.0d0
    d_root = 0.0d0
    d_leaf = 0.0d0
    stock_c_remaining = stock_c
    stock_n_remaining = stock_n

    !---------------------------------------------------------------------
    ! Coarse root turnover
    !---------------------------------------------------------------------

    call turnover_x ( &
                                         ! *** Input ***
    c_n_stem(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
    C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
    coarse_root_c_required          , &  ! Required carbon for turnover (mol C/tree/day)
                                         !
                                         ! *** In/Output ***
    remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_coarse_root                   , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
    stock_c_remaining               , &  ! Remaining carbon in stock (mol C stock/tree)
    stock_n_remaining                 &  ! Remaining nitrogen in stock (mol N stock/tree)
    )

    ! Flag of nitrogen limitation (1: N-limited, 0: C-limited)

    if (remaining_c > remaining_n*c_n_leaf(p)) then
       n_limit_flag = 1
    else
       n_limit_flag = 0
    end if

    ! Flag of priority of turnover compensation (1: root, 0: leaf)
    
    if (n_limit_flag==1 .or. stress_flag==1) then
       root_flag = 1
    else
       root_flag = 0
    end if

    if (root_flag == 1) then

       !---------------------------------------------------------------------
       ! Fine root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
       c_n_root(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       root_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_root                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_c_remaining               , &  ! Remaining carbon in stock (mol C stock/tree)
       stock_n_remaining                 &  ! Remaining nitrogen in stock (mol N stock/tree)
       )

       !---------------------------------------------------------------------
       ! Leaf root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
       c_n_leaf(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       leaf_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_c_remaining               , &  ! Remaining carbon in stock (mol C stock/tree)
       stock_n_remaining                 &  ! Remaining nitrogen in stock (mol N stock/tree)
       )

    else

       !---------------------------------------------------------------------
       ! Leaf root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
       c_n_leaf(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       leaf_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_c_remaining               , &  ! Remaining carbon in stock (mol C stock/tree)
       stock_n_remaining                 &  ! Remaining nitrogen in stock (mol N stock/tree)
       )

       !---------------------------------------------------------------------
       ! Fine root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
       c_n_root(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       root_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_root                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_c_remaining               , &  ! Remaining carbon in stock (mol C stock/tree)
       stock_n_remaining                 &  ! Remaining nitrogen in stock (mol N stock/tree)
       )

    end if

    d_stock_c = stock_c_remaining-stock_c
    d_stock_n = stock_n_remaining-stock_n

  end subroutine turnover

  !-----------------------------------------------------------------------
  subroutine turnover_x ( &
                                       ! *** Input ***
    c_n_ratio                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    c_required                    , &  ! Required carbon for turnover (mol C/tree/day)
                                       !
                                       ! *** In/Output ***
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_x                           , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
    stock_c_remaining             , &  ! Remaining carbon in stock (mol C stock/tree)
    stock_n_remaining               &  ! Remaining nitrogen in stock (mol N stock/tree)
    )
    !
    ! !Description:
    ! Turnover of component x (x = coarse root, leaf, root)
    ! Step 1: Gained resource is used.
    ! Step 2: stock_biomass is used.
    ! Step 3: Portion of remaining required resource is dropped.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: c_n_ratio, C_in_drymass
    real(8), intent(in)    :: c_required
    real(8), intent(inout) :: remaining_c, remaining_n
    real(8), intent(inout) :: d_x, stock_c_remaining, stock_n_remaining
    !
    ! !Local variables:
    real(8) :: c_required_remaining
    !---------------------------------------------------------------------

    if (min(remaining_c, remaining_n*c_n_ratio)>c_required) then
       c_required_remaining = 0.0d0
       remaining_c = remaining_c-c_required
       remaining_n = remaining_n-c_required/c_n_ratio
    else

       ! Step 1: Gained resource is used for turnover.

       c_required_remaining = c_required-min(remaining_c, remaining_n*c_n_ratio)
       remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_ratio)
       remaining_n = remaining_n-min(remaining_c/c_n_ratio, remaining_n)

       ! Load remaining C and N to stock

       stock_c_remaining = stock_c_remaining+max(remaining_c, 0.0d0)
       stock_n_remaining = stock_n_remaining+max(remaining_n, 0.0d0)

       remaining_c = 0.0d0
       remaining_n = 0.0d0

       if (min(stock_c_remaining, stock_n_remaining*c_n_ratio)                 &
           > c_required_remaining) then

          ! Step 2: remaining stock C and N is used for turnover.

          c_required_remaining = 0.0d0
          stock_c_remaining = stock_c_remaining-c_required_remaining
          stock_n_remaining = stock_n_remaining-c_required_remaining/c_n_ratio

       else

          ! Step 2: remaining stock C and N is used for turnover.

          c_required_remaining = c_required_remaining-min(stock_c_remaining,   &
                                 stock_n_remaining*c_n_ratio)
          stock_c_remaining = stock_c_remaining-min(stock_c_remaining,         &
                              stock_n_remaining*c_n_ratio)
          stock_n_remaining = stock_n_remaining-min(stock_c_remaining          &
                              /c_n_ratio, stock_n_remaining)

          ! Step 3: Portion of remaining required resource is dropped.

          c_required_remaining = 0.0d0
          d_x = d_x-c_required_remaining*12.0d0/C_in_drymass

       end if
    end if

  end subroutine turnover_x

end module mod_metabolism
