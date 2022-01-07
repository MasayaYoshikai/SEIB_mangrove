module mod_growth

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Biomass investment based on tree growth optimization
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS
  public  :: stock_reload
  public  :: belowground_limit
  public  :: aboveground_limit
  private :: invest_pattern1
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !LOCAL VARIABLES
  private
  real(8), parameter :: invest_ratio = 0.3d0
  real(8), parameter :: dpai_max = 1.0d0
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine stock_reload ( &
                                       ! *** Input ***
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    p                             , &  ! Species index (1: Rh, 2: Br)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    stock_n                       , &  ! Nitrogen in stock (mol N/tree)
                                       !
                                       ! *** Output ***
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    d_stock_n                     , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
    remaining_c                   , &  ! Remaining carbon for tree growth (mol C/tree/day)
    remaining_n                     &  ! Remaining nitrogen for tree growth (mol N/tree/day)
    )
    !
    ! !Description:
    ! Reload available resource to stock C and N to realize stock_trunk_ratio
    ! If stock C and/or N is lacking for stock_trunk_ratio, a portion of available
    ! resource (= invest_ratio) is invested to stock_c and/or n.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: stock_trunk_ratio(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)  :: C_in_drymass, available_c, available_n
    integer, intent(in)  :: p
    real(8), intent(in)  :: trunk_biomass, coarse_root_biomass, stock_c, stock_n
    real(8), intent(out) :: d_stock_c, d_stock_n, remaining_c, remaining_n
    !
    ! !Local variables:
    real(8) :: optimum_stock_c
    !---------------------------------------------------------------------

    optimum_stock_c = (trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p) &
                      *C_in_drymass/12.0d0
    if (min(stock_c, stock_n*c_n_stem(p))>=optimum_stock_c) then
       d_stock_c = 0.0d0
       d_stock_n = 0.0d0
       remaining_c = available_c
       remaining_n = available_n
    else
       if (stock_c > stock_n*c_n_stem(p)) then
          d_stock_c = 0.0d0
          d_stock_n = available_n*invest_ratio
          remaining_c = available_c
          remaining_n =available_n*(1.0d0-invest_ratio)
       else
          d_stock_c = available_c*invest_ratio
          d_stock_n = 0.0d0
          remaining_c = available_c*(1.0d0-invest_ratio)
          remaining_n = available_n
       end if
    end if

  end subroutine stock_reload

  !-----------------------------------------------------------------------
  subroutine belowground_limit ( &
                                       ! *** Input ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    c_n_root                      , &  ! C/N ratio in mol in root
    fine_root_ratio               , &  ! Fraction of fine root in below-ground biomass (-)
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    stress_flag                   , &  ! Flag of plant stress due to low leaf water potential (0: not stressed, 1: stressed)
    n_leaf_layer                  , &  ! Number of leaf layers
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    r_whole_root_increment        , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
    r_whole_d_increment           , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
    p                             , &  ! Species index (1: Rh, 2: Br)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    stock_n                       , &  ! Nitrogen in stock (mol N/tree)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
                                       !
                                       ! *** Output ***
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_coarse_root                 , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    d_stock_n                     , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment in below-ground limitation case
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)  :: c_n_root(:), fine_root_ratio(:)
    real(8), intent(in)  :: stock_trunk_ratio(:)
    real(8), intent(in)  :: SLA(:), DBH_limit(:), C_in_drymass
    integer, intent(in)  :: stress_flag, n_leaf_layer
    real(8), intent(in)  :: available_c, available_n
    real(8), intent(in)  :: r_whole_root_increment, r_whole_d_increment
    real(8), intent(in)  :: tree_h_limit
    integer, intent(in)  :: p
    real(8), intent(in)  :: tree_h, trunk_biomass, coarse_root_biomass
    real(8), intent(in)  :: leaf_biomass, root_biomass, stock_c
    real(8), intent(in)  :: stock_n, crown_area, dbh_heartwood, dbh_sapwood
    real(8), intent(out) :: d_leaf, d_root, d_trunk, d_coarse_root, d_stock_c
    real(8), intent(out) :: d_stock_n, remaining_c, remaining_n
    integer, intent(out) :: expand_flag
    !
    ! !Local variables:
    integer :: trunk_invest
    real(8) :: dpai
    real(8) :: leaf_required
    real(8) :: optimum_coarse_root
    real(8) :: coarse_root_required
    integer :: case_check
    !---------------------------------------------------------------------

    trunk_invest = 1
    if (min(stock_c, stock_n*c_n_stem(p))*12.0d0/C_in_drymass                  &
       < (trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p)*0.5d0) then
       trunk_invest = 0
    end if
    dpai = (leaf_biomass*SLA(p)/crown_area)/real(n_leaf_layer)
    if (dpai < optimum_dpai(p)) then
       leaf_required = ((optimum_dpai(p)*crown_area*real(n_leaf_layer)         &
                       /SLA(p))-leaf_biomass)*C_in_drymass/12.0d0
    else
       leaf_required = 0.0d0
    end if
    optimum_coarse_root = root_biomass*(1.0d0-fine_root_ratio(p))              &
                          /fine_root_ratio(p)
    if (coarse_root_biomass < optimum_coarse_root) then
       coarse_root_required = (optimum_coarse_root-coarse_root_biomass)        &
                              *C_in_drymass/12.0d0
    else
       coarse_root_required = 0.0d0
    end if
    remaining_c = available_c
    remaining_n = available_n
    d_trunk = 0.0d0
    d_coarse_root = 0.0d0
    d_leaf = 0.0d0
    d_root = 0.0d0
    d_stock_c = 0.0d0
    d_stock_n = 0.0d0
    case_check = 0
    expand_flag = 1

    !---------------------------------------------------------------------
    ! When current transpiration rate is not very stressful for the plant,
    ! plant tries to increase transpiration rate by investing
    ! biomass to leaf until it reaches the optimum_dpai.
    ! After that, if resource is still remaining, plant invests
    ! residual resoueces to stem to be higher.
    !---------------------------------------------------------------------

    if (stress_flag == 0) then
       call invest_pattern1 ( &
                                            ! *** Input ***
       optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
       DBH_limit                       , &  ! Limitation of trunk diameter (m)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       leaf_required                   , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
       dpai                            , &  ! Layer leaf area index (m2 leaf/m2 ground)
       trunk_invest                    , &  ! Flag of trunk investment (1: okay, 2: not okay)
       tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
       p                               , &  ! Species index (1: Rh, 2: Br)
       tree_h                          , &  ! Tree height (m)
       trunk_biomass                   , &  ! Stand-level trunk biomass (g trunk/tree)
       coarse_root_biomass             , &  ! Stand-level coarse root biomass (g/tree)
       stock_c                         , &  ! Carbon in stock (mol C/tree)
       stock_n                         , &  ! Nitrogen in stock (mol N/tree)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                     , &  ! Sapwood diameter (m)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
       d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
       d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
       case_check                      , &  ! Error check
       expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
       )

    !---------------------------------------------------------------------
    ! When current transpiration rate is stressful for the plant,
    ! plant tries to increase root water uptake rate by
    ! investing biomass based on plant hydraulics optimization.
    !---------------------------------------------------------------------

    else

       ! Case 1: Investing to root is most effective.

       if (r_whole_root_increment < r_whole_d_increment) then
          if (coarse_root_required>0.0d0 .and. trunk_invest==1) then
             d_coarse_root = d_coarse_root+min(remaining_c,                    &
                             remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
             remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
             remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
             case_check = case_check+1
          else
             d_root = d_root+min(remaining_c,                                  &
                      remaining_n*c_n_root(p))*12.0d0/C_in_drymass
             remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_root(p))
             remaining_n = remaining_n-min(remaining_c/c_n_root(p), remaining_n)
             case_check = case_check+1
          end if

       ! Case 2: Investing to trunk is most effective.

       else
          if (trunk_invest==1 .and. (dbh_heartwood+dbh_sapwood)                &
             <=DBH_limit(p)) then
             d_trunk = d_trunk+min(remaining_c,                                &
                       remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
             remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
             remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
             case_check = case_check+1
          elseif (coarse_root_required > 0.0d0) then
             d_coarse_root = d_coarse_root+min(remaining_c,                    &
                             remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
             remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
             remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
             case_check = case_check+1
          elseif (min(stock_c, stock_n*c_n_stem(p))*12.0d0/C_in_drymass        &
                 <(trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p)) then
             d_stock_c = d_stock_c+remaining_c
             d_stock_n = d_stock_n+remaining_n
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check+1
          else
             d_root = d_root+min(remaining_c,                                  &
                      remaining_n*c_n_root(p))*12.0d0/C_in_drymass
             remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_root(p))
             remaining_n = remaining_n-min(remaining_c/c_n_root(p), remaining_n)
             case_check = case_check+1
          end if
       end if
    end if

    remaining_c = max(remaining_c, 0.0d0)
    remaining_n = max(remaining_n, 0.0d0)

    if (case_check > 1) then
       write(*,*) 'ERROR: case_check is more than one in belowground_limit.'
    elseif (case_check < 1) then
       write(*,*) 'ERROR: case_check is still zero in belowground_limit'
    end if

  end subroutine belowground_limit

  !-----------------------------------------------------------------------
  subroutine aboveground_limit ( &
                                       ! *** Input ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    n_leaf_layer                  , &  ! Number of leaf layers
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
    par_direct                    , &  ! Direct PAR on canopy top of midday (umol photon/m2/s)
    par_diffuse                   , &  ! Diffused PAR on canopy top of midday (umol photon/m2/s)
    p                             , &  ! Species index (1: Rh, 2: Br)
    rPAR_dir                      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
    rPAR_dif                      , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
    seib_height                   , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    stock_n                       , &  ! Nitrogen in stock (mol N/tree)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
                                       !
                                       ! *** Output ***
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    d_stock_n                     , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment in above-ground limitation case
    !
    ! !Uses:
    use data_structure, only : Max_hgt
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)  :: stock_trunk_ratio(:), SLA(:), DBH_limit(:)
    real(8), intent(in)  :: C_in_drymass
    integer, intent(in)  :: n_leaf_layer
    real(8), intent(in)  :: available_c, available_n, tree_h_limit
    real(8), intent(in)  :: par_direct, par_diffuse
    integer, intent(in)  :: p
    real(8), intent(in)  :: rPAR_dir(:), rPAR_dif(:)
    integer, intent(in)  :: seib_height
    real(8), intent(in)  :: tree_h, trunk_biomass, coarse_root_biomass
    real(8), intent(in)  :: leaf_biomass, stock_c
    real(8), intent(in)  :: stock_n, crown_area, dbh_heartwood, dbh_sapwood
    real(8), intent(out) :: d_leaf, d_trunk, d_stock_c, d_stock_n
    real(8), intent(out) :: remaining_c, remaining_n
    integer, intent(out) :: expand_flag
    !
    ! !Local variables:
    real(8), parameter :: critical_ratio = 0.15d0  ! Critical ratio for determining low photosynthesis efficiency (fraction)
    real(8) :: rPAR
    integer :: efficiency_flag
    integer :: trunk_invest
    real(8) :: dpai
    real(8) :: leaf_required
    integer :: case_check
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Relative PAR intensity at canopy top
    !---------------------------------------------------------------------

    rPAR = (rPAR_dir(seib_height)*par_direct+rPAR_dif(seib_height)*par_diffuse)&
           /(par_direct+par_diffuse)

    !---------------------------------------------------------------------
    ! Flag of photosynthesis efficiency.
    ! 1: Still efficient => Low photosynthesis is due to small amount of leaves
    ! 0: Inefficient     => Low photosynthesis is due to low light
    !---------------------------------------------------------------------

    if (rPAR > critical_ratio) then
       efficiency_flag = 1
    else
       efficiency_flag = 0
    end if

    trunk_invest = 1
    if (min(stock_c, stock_n*c_n_stem(p))*12.0d0/C_in_drymass                  &
       <(trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p)*0.5d0) then
       trunk_invest = 0
    end if
    dpai = (leaf_biomass*SLA(p)/crown_area)/real(n_leaf_layer)
    if (dpai < optimum_dpai(p)) then
       leaf_required = ((optimum_dpai(p)*crown_area*real(n_leaf_layer)/SLA(p)) &
                       -leaf_biomass)*C_in_drymass/12.0d0
    else
       leaf_required = 0.0d0
    end if
    remaining_c = available_c
    remaining_n = available_n
    d_trunk = 0.0d0
    d_leaf = 0.0d0
    d_stock_c = 0.0d0
    d_stock_n = 0.0d0
    case_check = 0
    expand_flag = 0

    !---------------------------------------------------------------------
    ! When low productivity is due to small amount of leaves,
    ! plant invests biomass to leaf to realize optimum_dpai.
    ! Then if resource is still remaining, plant invests residual resource
    ! to stem to be higher.
    !---------------------------------------------------------------------

    if (efficiency_flag == 1) then
       call invest_pattern1 ( &
                                            ! *** Input ***
       optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
       DBH_limit                       , &  ! Limitation of trunk diameter (m)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
       leaf_required                   , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
       dpai                            , &  ! Layer leaf area index (m2 leaf/m2 ground)
       trunk_invest                    , &  ! Flag of trunk investment (1: okay, 2: not okay)
       tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
       p                               , &  ! Species index (1: Rh, 2: Br)
       tree_h                          , &  ! Tree height (m)
       trunk_biomass                   , &  ! Stand-level trunk biomass (g trunk/tree)
       coarse_root_biomass             , &  ! Stand-level coarse root biomass (g/tree)
       stock_c                         , &  ! Carbon in stock (mol C/tree)
       stock_n                         , &  ! Nitrogen in stock (mol N/tree)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                     , &  ! Sapwood diameter (m)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
       d_stock_c                       , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
       d_stock_n                       , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
       case_check                      , &  ! Error check
       expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
       )

    !---------------------------------------------------------------------
    ! When low productivity is due to lack of light,
    ! plant invests biomass to stem to be higher.
    !---------------------------------------------------------------------

    else
       if (trunk_invest==1 .and. tree_h<tree_h_limit) then
          d_trunk = d_trunk+min(remaining_c,                                   &
                    remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
          remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
          case_check = case_check+1
          expand_flag = 0
       elseif (trunk_invest==1 .and. (dbh_heartwood+dbh_sapwood)               &
               <=DBH_limit(p)) then
          d_trunk = d_trunk+min(remaining_c,                                   &
                    remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
          remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
          case_check = case_check+1
          expand_flag = 1
       elseif (min(stock_c, stock_n*c_n_stem(p))*12.0d0/C_in_drymass           &
              < (trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p)) then
          d_stock_c = d_stock_c+remaining_c
          d_stock_n = d_stock_n+remaining_n
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check+1
       elseif (dpai > optimum_dpai(p)*dpai_max) then
          d_stock_c = d_stock_c+remaining_c
          d_stock_n = d_stock_n+remaining_n
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check+1
       else
          d_leaf = d_leaf+min(remaining_c,                                     &
                   remaining_n*c_n_leaf(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_leaf(p))
          remaining_n = remaining_n-min(remaining_c/c_n_leaf(p), remaining_n)
          case_check = case_check+1
       end if
    end if

    remaining_c = max(remaining_c, 0.0d0)
    remaining_n = max(remaining_n, 0.0d0)

    if (case_check > 1) then
       write(*,*) 'ERROR: case_check is more than one in aboveground_limit.'
    elseif (case_check < 1) then
       write(*,*) 'ERROR: case_check is still zero in aboveground_limit'
    end if

  end subroutine aboveground_limit

  !-----------------------------------------------------------------------
  subroutine invest_pattern1 ( &
                                       ! *** Input ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    leaf_required                 , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
    dpai                          , &  ! Layer leaf area index (m2 leaf/m2 ground)
    trunk_invest                  , &  ! Flag of trunk investment (1: okay, 2: not okay)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
    p                             , &  ! Species index (1: Rh, 2: Br)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_c                       , &  ! Carbon in stock (mol C/tree)
    stock_n                       , &  ! Nitrogen in stock (mol N/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
                                       !
                                       ! *** In/Output ***
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_stock_c                     , &  ! dSTOCK_C/dt: Daily stock carbon change (mol C/tree/day)
    d_stock_n                     , &  ! dSTOCK_N/dt: Daily stock nitrogen change (mol N/tree/day)
    case_check                    , &  ! Error check
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment pattern 1
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)    :: stock_trunk_ratio(:), DBH_limit(:), C_in_drymass
    real(8), intent(in)    :: leaf_required, dpai
    integer, intent(in)    :: trunk_invest
    real(8), intent(in)    :: tree_h_limit
    integer, intent(in)    :: p
    real(8), intent(in)    :: tree_h, trunk_biomass, coarse_root_biomass, stock_c, stock_n
    real(8), intent(in)    :: dbh_heartwood, dbh_sapwood
    real(8), intent(inout) :: remaining_c, remaining_n
    real(8), intent(inout) :: d_leaf, d_trunk, d_stock_c, d_stock_n
    integer, intent(inout) :: case_check
    integer, intent(inout) :: expand_flag
    !
    ! !Local variables:
    !
    !
    !---------------------------------------------------------------------

    if (min(remaining_c, remaining_n*c_n_leaf(p)) > leaf_required) then
       d_leaf = d_leaf+leaf_required*12.0d0/C_in_drymass
       remaining_c = remaining_c-leaf_required
       remaining_n = remaining_n-leaf_required/c_n_leaf(p)
       if (trunk_invest==1 .and. tree_h<tree_h_limit) then
          d_trunk = d_trunk+min(remaining_c,                                   &
                    remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
          remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
          case_check = case_check+1
          expand_flag = 0
       elseif (trunk_invest==1 .and. (dbh_heartwood+dbh_sapwood)               &
               <=DBH_limit(p)) then
          d_trunk = d_trunk+min(remaining_c,                                   &
                    remaining_n*c_n_stem(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_stem(p))
          remaining_n = remaining_n-min(remaining_c/c_n_stem(p), remaining_n)
          case_check = case_check+1
          expand_flag = 1
       elseif (min(stock_c, stock_n*c_n_stem(p))*12.0d0/C_in_drymass           &
              < (trunk_biomass+coarse_root_biomass)*stock_trunk_ratio(p)) then
          d_stock_c = d_stock_c+remaining_c
          d_stock_n = d_stock_n+remaining_n
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check+1
       elseif (dpai > optimum_dpai(p)*dpai_max) then
          d_stock_c = d_stock_c+remaining_c
          d_stock_n = d_stock_n+remaining_n
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check+1
       else
          d_leaf = d_leaf+min(remaining_c,                                     &
                   remaining_n*c_n_leaf(p))*12.0d0/C_in_drymass
          remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_leaf(p))
          remaining_n = remaining_n-min(remaining_c/c_n_leaf(p), remaining_n)
          case_check = case_check+1
       end if
    else
       d_leaf = d_leaf+min(remaining_c,                                        &
                remaining_n*c_n_leaf(p))*12.0d0/C_in_drymass
       remaining_c = remaining_c-min(remaining_c, remaining_n*c_n_leaf(p))
       remaining_n = remaining_n-min(remaining_c/c_n_leaf(p), remaining_n)
       case_check = case_check+1
    end if

  end subroutine invest_pattern1

end module mod_growth
