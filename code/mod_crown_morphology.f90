module mod_crown_morphology

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
  public :: crown_bottom_purge
  public :: crown_allometry
  !
  ! !PRIVATE MEMBER FUNCTIONS
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine crown_bottom_purge ( &
                                       ! *** Input ***
    leaf_turn                     , &  ! Leaf turnover rate (day-1)
    grow_resp                     , &  ! Growth respiration rate (g DW/g DW)
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    leaf_resorp                   , &  ! Nutrient resorption from Leaf (fraction)
    LA_max                        , &  ! Maximum leaf area index (m2 leaf/m2 ground)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    c_uptake_bottom_day           , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day           , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    seib_lai                      , &  ! Leaf area index (m2 leaf/m2 ground)
    day_count                     , &  ! Days since the simulation start (days)
    no                            , &  ! Tree index
    p                             , &  ! Species index (1: Rh, 2: Br)
    STEP                          , &  ! Canopy layer thickness (m)
    seib_height                   , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
                                       !
                                       ! *** In/Output ***
    npp_bottom_c                  , &  ! NPP of C from 1 m2 of leaves at crown bottom layer (mol C/m2 leaf)
    npp_bottom_n                  , &  ! NPP of N from 1 m2 of leaves at crown bottom layer (mol N/m2 leaf)
                                       !
                                       ! *** Output ***
    purge_flag                      &  ! Flag of perge crown bottom layer (0: no, 1: yes)
    )
    !
    ! !Description:
    ! Calculate net production of 1 m2 of leaves at crown bottom layer.
    ! Then decide whether purge the bottom layer or not.
    ! Also, if LAI is already reached at the maximum,
    ! the bottom layer will be purged.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: leaf_turn(:), grow_resp(:), c_n_leaf(:)
    real(8), intent(in)    :: leaf_resorp(:), LA_max(:), SLA(:), C_in_drymass
    real(8), intent(in)    :: c_uptake_bottom_day, n_uptake_bottom_day, seib_lai
    integer, intent(in)    :: day_count, no, p
    real(8), intent(in)    :: STEP
    integer, intent(in)    :: seib_height, seib_bole
    real(8), intent(inout) :: npp_bottom_c(:), npp_bottom_n(:)
    integer, intent(out)   :: purge_flag
    !
    ! !Local variables:
    integer, parameter :: crown_day = 183
    real(8) :: turn_cost_c
    real(8) :: turn_cost_n
    real(8) :: leaf_production_c
    real(8) :: leaf_production_n
    real(8) :: crown_depth
    real(8) :: crown_depth_min = 0.5d0
    real(8), parameter :: crit_ratio = 1.05d0
    !---------------------------------------------------------------------

    ! Turnover cost of 1 m2 of leaves (mol C or N/m2 leaf)

    turn_cost_c = (1.0d0/SLA(p))*leaf_turn(p)*real(crown_day)                  &
                  *C_in_drymass/12.0d0
    turn_cost_n = (1.0d0/SLA(p))*leaf_turn(p)*real(crown_day)                  &
                  *C_in_drymass/12.0d0/c_n_leaf(p)
    turn_cost_n = turn_cost_n*(1.0d0-leaf_resorp(p))

    ! Leaf biomass production from C and N uptake rates
    ! per unit leaf area (mol C or N/m2 leaf)

    npp_bottom_c(no) = npp_bottom_c(no)+(1.0d0-grow_resp(p))                   &
                       *c_uptake_bottom_day
    npp_bottom_n(no) = npp_bottom_n(no)+n_uptake_bottom_day

    crown_depth = (real(seib_height)-real(seib_bole))*STEP
    purge_flag = 0

    ! Day to purge crown bottom layer

    if (mod(day_count-2, crown_day)==(crown_day-1)) then

       ! Case bottom layer production is not effective

       if (min(npp_bottom_c(no)/turn_cost_c, npp_bottom_n(no)/turn_cost_n)     &
           <crit_ratio) then
          purge_flag = 1
          if (crown_depth <= crown_depth_min) then
             purge_flag = 0
          end if

       ! Case bottom layer production is effective

       else
          purge_flag = 0
       end if

       ! When LAI is already reached at the maximum
       
       if (seib_lai>LA_max(p)) then
          purge_flag = 1
          if (crown_depth<=crown_depth_min) then
             purge_flag = 0
          end if
       end if
       npp_bottom_c(no) = 0.0d0
       npp_bottom_n(no) = 0.0d0
    end if

  end subroutine crown_bottom_purge

  !-----------------------------------------------------------------------
  function crown_allometry ( &
                                       ! *** Input ***
    crown_a                       , &  ! Tree crown allometric parameter
    crown_b                       , &  ! Tree crown allometric parameter
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ans)
    !
    ! !DESCRIPTION:
    ! Calculate potential crown diameter based on allometric relation.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: crown_a(:), crown_b(:)
    integer, intent(in) :: p
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: ans
    !---------------------------------------------------------------------

    dbh = dbh_heartwood+dbh_sapwood
!     ans = crown_a(p)*log(dbh)+crown_b(p)  ! Old equation
    ans = crown_b(p)*dbh**crown_a(p)
    ans = max(ans, 1.0d0)

  end function crown_allometry

end module mod_crown_morphology
