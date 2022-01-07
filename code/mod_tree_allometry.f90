module mod_tree_allometry

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate increment of D and H after stem biomass increment
  ! based on allometric equation.
  !
  ! Reference:
  ! Chave et al. (2005) Oecologia 145: 87-99
  !
  ! The allometric model is expresses as:
  ! M = 0.051 * rho * D^2 * H
  ! where
  ! M: above-ground biomass (kg)
  ! rho: wood density (g/cm3)
  ! D: stem diameter (cm)
  ! H: tree height (m)
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: tree_allometry
  public :: height_allometry
  public :: min_height_allometry
  public :: proot_allometry
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine tree_allometry ( &
                                       ! *** Input ***
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g/tree)
    wood_rho                      , &  ! Wood density (g/cm3)
    ALM5                          , &  ! Sapwood diameter proportion (m sapwood/m dbh)
    biomass_increment             , &  ! Increment of trunk biomass (g/tree)
                                       !
                                       ! *** Output ***
    new_dbh_heartwood             , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood               , &  ! New sapwood diameter after biomass increment (m)
    new_tree_h                      &  ! New tree height after biomass increment (m)
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: dbh_heartwood, dbh_sapwood, tree_h
    real(8), intent(in)  :: trunk_biomass, wood_rho, ALM5, biomass_increment
    real(8), intent(out) :: new_dbh_heartwood, new_dbh_sapwood, new_tree_h
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: new_dbh
    real(8) :: new_trunk_biomass
    !
    ! -------
    ! Use Komiyama
    real(8) :: term1
    !---------------------------------------------------------------------

    new_trunk_biomass = (trunk_biomass+biomass_increment)/1000.0d0
    dbh = (dbh_heartwood+dbh_sapwood)*100.0d0

    ! New DBH after biomass increment (m)
    new_dbh = (sqrt(new_trunk_biomass/(0.051d0*wood_rho*tree_h)))*0.01d0
    new_dbh = max(new_dbh, dbh/100.0d0)
    new_dbh_heartwood = (1.0d0-ALM5)*new_dbh
    new_dbh_sapwood = ALM5*new_dbh

    ! New tree height after biomass increment (m)
    new_tree_h = new_trunk_biomass/(0.051d0*wood_rho*dbh**2.0d0)
    new_tree_h = max(new_tree_h, tree_h)

    ! ------------------
    ! Use Komiyama

    term1 = (new_trunk_biomass/(0.0696d0*wood_rho))**(1.0d0/0.931d0)

    ! New DBH after biomass increment (m)
    new_dbh = sqrt(term1/tree_h)*0.01d0
    new_dbh = max(new_dbh, dbh/100.0d0)
    new_dbh_heartwood = (1.0d0-ALM5)*new_dbh
    new_dbh_sapwood = ALM5*new_dbh

    ! New tree height after biomass increment (m)
    new_tree_h = term1/(dbh**2.0d0)
    new_tree_h = max(new_tree_h, tree_h)
    ! ------------------

  end subroutine tree_allometry

  !-----------------------------------------------------------------------
  function height_allometry ( &
                                       ! *** Input ***
    tree_h_a                      , &  ! Tree height allometric parameter
    tree_h_b                      , &  ! Tree height allometric parameter
    HGT_max                       , &  ! Maximum tree height (m)
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ans)
    !
    ! !DESCRIPTION:
    ! Calculate potential tree height based on allometric relation.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: tree_h_a(:), tree_h_b(:), HGT_max(:)
    integer, intent(in) :: p
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: ans
    !---------------------------------------------------------------------

    dbh = dbh_heartwood+dbh_sapwood
    ans = tree_h_a(p)*(dbh**tree_h_b(p))
    ans = max(ans, 2.0d0)
    ans = min(ans, HGT_max(p))

  end function height_allometry

  !-----------------------------------------------------------------------
  function min_height_allometry ( &
                                       ! *** Input ***
    tree_h_min_a                  , &  ! Minimum tree height allometric parameter
    tree_h_min_b                  , &  ! Minimum tree height allometric parameter
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ans)
    !
    ! !DESCRIPTION:
    ! Calculate potential minimum tree height given DBH,
    ! which will be related to motrality rate.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: tree_h_min_a(:), tree_h_min_b(:)
    integer, intent(in) :: p
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: ans
    !---------------------------------------------------------------------

    dbh = dbh_heartwood+dbh_sapwood
    ans = tree_h_min_a(p)*(dbh**tree_h_min_b(p))
!    ans = max(ans, 3.0d0)

  end function min_height_allometry

  !-----------------------------------------------------------------------
  subroutine proot_allometry ( &
                                       ! *** Input ***
    wood_rho                      , &  ! Wood density (g/cm3)
    pr_s_a                        , &  ! Slope for the scaling factor of prop root system
    pr_s_b                        , &  ! Intercept for the scaling factor of prop root system
    pr_h_a                        , &  ! Slope for the maximum root height of prop root system
    pr_h_b                        , &  ! Intercept for the maximum root height of prop root system
    pr_d                          , &  ! Mean prop root diameter (m)
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    mass_trunk                    , &  ! Stand-level trunk biomass (g trunk/tree)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
                                       !
                                       ! *** In/Output ***
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                       !
                                       ! *** Output ***
    d_above_root                  , &  ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
    above_root_ratio                &  ! Required fraction of above-ground root biomass in above-ground
    )
    !
    ! !DESCRIPTION:
    ! Calculate needed prop root biomass corresponding to the DBH
    ! based on allometric relationship.
    !
    ! !USES:
    use data_structure, only : PI
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)    :: wood_rho(:), pr_s_a(:), pr_s_b(:)
    real(8), intent(in)    :: pr_h_a(:), pr_h_b(:), pr_d(:)
    integer, intent(in)    :: p
    real(8), intent(in)    :: dbh_heartwood, dbh_sapwood
    real(8), intent(in)    :: mass_trunk, mass_above_root
    real(8), intent(inout) :: d_trunk
    real(8), intent(out)   :: d_above_root, above_root_ratio
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: s_factor
    real(8) :: hr_max
    real(8) :: hr
    integer :: z
    integer :: r
    real(8), parameter :: hr_min   = 5.0d0  ! Minimum prop root height (cm)
    integer, parameter :: n_layer  = 20     ! Number of vertical layers
    real(8), parameter :: dz       = 10.0d0 ! Layer thickness (cm)
    real(8), parameter :: l_factor = 14.5d0 ! Factor for converting prop root number to total length in cm
    real(8), dimension(n_layer) :: h_layer
    real(8), dimension(n_layer) :: n_root
    real(8), dimension(n_layer) :: v_root
    real(8) :: total_v_root
    real(8) :: proot_biomass
    !---------------------------------------------------------------------

    if (pr_s_a(p) < 0.0d0) then ! For Rhizophora genus

       dbh = dbh_heartwood+dbh_sapwood
       dbh = max(dbh, 0.03d0)

       ! Scaling factor for prop root system
       s_factor = 1.0d0-10.0d0**(pr_s_a(p)*log10(dbh)+pr_s_b(p))
       s_factor = max(s_factor, 0.3d0)  ! Prevent very low values

       ! Maximum prop root height (cm)
       hr_max = (pr_h_a(p)*dbh+pr_h_b(p))*100.0d0

       ! Height of upper boundary of each vertical layer (cm)
       h_layer(:) = 0.0d0
       do z = 1, n_layer
          h_layer(z) = real(n_layer-z+1)*dz
       end do

       ! Computation of number of prop roots in each vertical layer
       hr = hr_max
       n_root(:) = 0.0d0
       v_root(:) = 0.0d0
       do r = 1, 1000
          do z = 1, n_layer
             if (hr > h_layer(z)) then
                n_root(z) = n_root(z)+1.0d0
             else
                if (z == n_layer) then
                   n_root(z) = n_root(z)+hr/dz
                elseif (hr > h_layer(z+1)) then
                   n_root(z) = n_root(z)+(hr-h_layer(z+1))/dz
                end if
             end if
          end do

          ! Next prop root height
          hr = hr*s_factor

          ! Exit the loop when the prop root height is lower than minimum value.
          if (hr < hr_min) exit

       end do

       ! Total volume of prop roots in each vertical layer (cm3)
       do z = 1, n_layer
          v_root(z) = l_factor*n_root(z)*PI*(pr_d(p)*100.0d0*0.5d0)**2.0d0
       end do

       ! Integrate for the vertical layers (cm3)
       total_v_root = sum(v_root)

       ! Prop root biomass needed based on allometric relationship (g/tree)
       proot_biomass = total_v_root*wood_rho(p)

       ! Required fraction of above-ground root biomass in above-ground
       above_root_ratio = proot_biomass/(proot_biomass+mass_trunk)

       ! Regulate very large fraction of above-ground roots
       above_root_ratio = min(above_root_ratio, 0.50d0)

       ! Recalculation of proot_biomass from above_root_ratio
       proot_biomass = mass_trunk*above_root_ratio/(1.0d0-above_root_ratio)

!       ! Biomass allocation to stem and above-ground root
!
!       d_above_root = d_trunk*above_root_ratio
!       d_trunk = d_trunk*(1.0d0-above_root_ratio)
!
       if (mass_above_root < proot_biomass) then
          if (d_trunk >= (proot_biomass-mass_above_root)) then
             d_above_root = proot_biomass-mass_above_root
             d_trunk = d_trunk-d_above_root
          else
             d_above_root = 0.50d0*d_trunk
             d_trunk = d_trunk-d_above_root
          end if
       else
          d_above_root = 0.0d0
       end if

    else ! For other genus

       d_above_root = 0.0d0
       above_root_ratio = 0.0d0

    end if

  end subroutine proot_allometry

end module mod_tree_allometry
