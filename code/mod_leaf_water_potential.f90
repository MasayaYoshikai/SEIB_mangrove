
!=======================================================================
!  This code was adapted from StomataOptimizationMod.F90               !
!  in the link below.                                                  !
!  https://github.com/gbonan/CLM-ml_v0                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Bonan G. B., Williams M., Fisher R. A., Oleson K. W. (2014),        !
!  Modeling stomatal conductance in the earth system: linking leaf     !
!  water-use efficiency and water transport along the                  !
!  soil–plant–atmosphere continuum, Geoscientific Model Development,   !
!  7, 2019-2222.                                                       !
!                                                                      !
!=======================================================================

module mod_leaf_water_potential

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf water potential based on plant water flux
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: leaf_water_potential
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine leaf_water_potential ( &
                                       ! *** Input ***
    leaf_cp                       , &  ! Leaf water capacitance (mmol H2O/m2 leaf/MPa)
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
    !
    ! !DESCRIPTION:
    ! New leaf water potential based on the stand-level water balance.
    ! It is assumed that the leaf water potential is uniform over the canopy.
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: leaf_cp, tree_h, crown_area
    real(8), intent(in) :: leaf_wp, tot_psi
    type(universal_type), intent(in)    :: universal
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(in)    :: flux
    real(8), intent(out) :: leaf_wp_new
    !
    ! !LOCAL VARIABLES:
    integer :: i             ! Layer index
    real(8) :: stand_et      ! Total transpiration rate at stand-level (mol H2O/tree/s)
    real(8) :: total_la      ! Total leaf area (m2 leaf/tree)
    real(8) :: total_cp      ! Total Leaf water capacitance (mmol H2O/tree/MPa)
    real(8) :: height_psi    ! Height potential (MPa)
    real(8) :: dtime         ! Model time step (s)
    real(8) :: y0            ! Leaf water potential at beginning of timestep (MPa)
    real(8) :: dy            ! Change in leaf water potential (MPa)
    real(8) :: a, b          ! Intermediate calculation
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    denh2o     => universal%denh2o  , &  ! Water density (kg/m3)
    grav       => universal%grav    , &  ! Gravitational acceleration (m/s2)
    n_layer    => universal%n_layer , &  ! Number of leaf layers
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    dpai       => layer%dpai        , &  ! Layer leaf area index (m2 leaf/m2 ground)
    etflux     => flux%etflux       , &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    r_whole    => flux%r_whole        &  ! Whole plant resistance (MPa.s.tree/mmol H2O)
    )

    ! Height potential (MPa)

    height_psi = denh2o * grav * tree_h * 1.e-06

    ! Stand-level transpiration rate (mol H2O/tree/s)

    stand_et = 0.0d0
    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer
          stand_et = stand_et + etflux(i) * dpai * crown_area
       end if
    end do

    ! Total leaf area (m2 leaf/tree)

    total_la =0.0d0
    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer
          total_la = total_la + dpai * crown_area
       end if
    end do

    ! Total Leaf water capacitance (mmol H2O/tree/MPa)

    total_cp = leaf_cp * total_la

    ! Time step (s): 1 hour

    dtime = 3600.0d0

    ! Change in leaf water potential is: dy / dt = (a -y) / b. The integrated change
    ! over a full model time step is: dy = (a - y0) * (1 - exp(-dt/b))

    y0 = leaf_wp
    a = tot_psi - height_psi - 1000.0d0 * stand_et * r_whole
    b = r_whole * total_cp
    dy = (a - y0) * (1.0d0 - exp(-dtime/b))
    leaf_wp_new = y0 + dy

    end associate
  end subroutine leaf_water_potential

end module mod_leaf_water_potential
