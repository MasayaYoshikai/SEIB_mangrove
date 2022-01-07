
!=======================================================================
!  This code was adapted from PlantHydraulicsMod.F90                   !
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
!  Xu X., Medvigy D., Powers J. S., Becknell J. M., Guan K. (2016),    !
!  Diversity in plant hydraulic traits explains seasonal and           !
!  inter-annual variations of vegetation dynamics in seasonally dry    !
!  tropical forests, New Phytologist, 212, 80-95.                      !
!                                                                      !
!=======================================================================

module mod_plant_hydraulics

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate plant hydraulic resistance
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: soil_root_hydraulics
  public :: plant_hydraulics
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !LOCAL VARIABLES
  private
  real(8), parameter :: pi = 3.14159265d0
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine soil_root_hydraulics ( &
                                       ! *** Input ***
    root_radius                   , &  ! Fine root radius (m)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    root_area                     , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth                    , &  ! Rooting depth (m)
    root_density                  , &  ! Specific root density (fine root) (g root/m3 root)
    hk                            , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist                   , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    universal                     , &  ! Universal variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Calculate soil-to-root and root-to-stem hydraulic resistance
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: root_radius, root_biomass, root_area, root_depth
    real(8), intent(in) :: root_density, hk, root_resist
    type(universal_type), intent(in) :: universal
    type(flux_type), intent(inout)   :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: root_cross_sec_area
    real(8) :: root_length_density
    real(8) :: root_dist
    !---------------------------------------------------------------------

    associate ( &
                                      ! *** Output ***
    r_soil  => flux%r_soil       , &  ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
    r_root  => flux%r_root         &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    )

    !---------------------------------------------------------------------
    ! Soil-root resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    ! Root cross-sectional area (m2 root)
    root_cross_sec_area = pi*root_radius**2.0d0

    ! Root length density per unit volume of soil (m root/m3 ground)
    root_length_density = root_biomass/(root_area*root_depth*root_density      &
                          *root_cross_sec_area)

    ! One-half distance between roots (m)
    root_dist = sqrt(1.0d0/(pi*root_length_density))

    ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
    r_soil = log(root_dist/root_radius)/(2.0d0*pi*root_length_density          &
             *root_depth*hk*root_area)

    !---------------------------------------------------------------------
    ! Root-stem resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    r_root = root_resist/root_biomass

    end associate
  end subroutine soil_root_hydraulics

  !-----------------------------------------------------------------------
  subroutine plant_hydraulics ( &
                                       ! *** Input ***
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    k_sap                         , &  ! Stem hydraulic conductivity at saturation (kg H2O.m/m2 sapwood/s/MPa)
    p50_sap                       , &  ! Stem water potential at which 50% of conductivity is lost (MPa)
    a2_sap                        , &  ! Conductivity vulnerability curve coefficient (-)
    tree_h                        , &  ! Tree height (m)
    leaf_wp                       , &  ! Leaf water potential at current time step (MPa)
    tot_psi                       , &  ! Total soil water potential (MPa)
    universal                     , &  ! Universal variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Calculate plant hydraulic resistance based on
    ! the current leaf water potential
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood, k_sap, p50_sap, a2_sap
    real(8), intent(in) :: tree_h, leaf_wp, tot_psi
    type(universal_type), intent(in) :: universal
    type(flux_type), intent(inout)   :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh
    real(8) :: sap_area
    real(8) :: k_sap_val
    real(8) :: k_sap_molar
    real(8) :: height_psi
    real(8), parameter :: coef_a1 = 1.20d0  ! Multiplier for converting tree height to water path length
    !---------------------------------------------------------------------

    associate ( &
                                      ! *** Input ***
    denh2o  => universal%denh2o  , &  ! Water density (kg/m3)
    grav    => universal%grav    , &  ! Gravitational acceleration (m/s2)
    mmh2o   => universal%mmh2o   , &  ! Molecular mass of water (kg/mol)
    r_soil  => flux%r_soil       , &  ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
    r_root  => flux%r_root       , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
                                      !
                                      ! *** Output ***
    dr_sap  => flux%dr_sap       , &  ! Stem resistance per unit path length (MPa.s.tree/mmol H2O/m)
    r_sap   => flux%r_sap        , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    r_whole => flux%r_whole      , &  ! Whole plant resistance (MPa.s.tree/mmol H2O)
    sapflow => flux%sapflow        &  ! Stand-level sapflow rate (mol H2O/tree/s)
    )

    !---------------------------------------------------------------------
    ! Stem resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    ! Stem diameter (= DBH, m)
    dbh = dbh_heartwood+dbh_sapwood

    ! Sapwood area (m2/tree)
    sap_area = pi*(dbh/2.0d0)**2.0d0-pi*(dbh_heartwood/2.0d0)**2.0d0

    ! Stem hydraulic conductivity at the given water potential
    ! (kg H2O.m/m2 sapwood/s/MPa)
    k_sap_val = k_sap*(1.0d0+(leaf_wp/p50_sap)**a2_sap)**(-1.0d0)

    ! Stem hydraulic conductivity unit conversion
    ! (kg H2O.m/m2 sapwood/s/MPa) -> (mmol H2O.m/m2 sapwood/s/MPa)
    k_sap_molar = (k_sap_val/mmh2o)*1000.0d0

    ! Stem resistance per unit path length (MPa.s.tree/mmol H2O/m)
    dr_sap = 1.0d0/(k_sap_molar*sap_area)

    ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    r_sap = dr_sap*tree_h*coef_a1

    !---------------------------------------------------------------------
    ! Whole plant resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    r_whole = r_soil+r_root+r_sap

    !---------------------------------------------------------------------
    ! Plant water uptake rate (mol H2O/tree/s)
    !---------------------------------------------------------------------

    ! Tree height potential (MPa)
    height_psi = denh2o*grav*tree_h*1.e-06

    ! Plant water uptake rate (mol H2O/tree/s)
    sapflow = ((tot_psi-height_psi-leaf_wp)/r_whole)*1.e-03

    end associate
  end subroutine plant_hydraulics

end module mod_plant_hydraulics
