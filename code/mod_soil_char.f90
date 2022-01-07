
!=======================================================================
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Perri S., Viola F., Noto L. V., Molini S. (2017), Salinity and      !
!  periodic inundation controls on the soil‐plant‐atmosphere continuum !
!  of gray mangroves, Hydrological Processes, 31, 1271–1282.           !
!                                                                      !
!=======================================================================

module mod_soil_char

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set soil hydraulic parameters
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: soil_char
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine soil_char ( &
                                       ! *** Input ***
    hksat                         , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
    moist                         , &  ! Relative soil moisture (-)
    bsw                           , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat                        , &  ! Soil matric potential at saturation (MPa)
    soil_t                        , &  ! Soil water temperature (K)
    sal                           , &  ! Pore-water salinity (mol/m3)
    root_filter                   , &  ! Filtering rate of salt at root surface (-)
    universal                     , &  ! Universal variables
                                       !
                                       ! *** Output ***
    tot_psi                       , &  ! Total soil water potential (MPa)
    hk                              &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: hksat, moist, bsw, psisat, soil_t, sal, root_filter
    type(universal_type), intent(in) :: universal
    real(8), intent(out) :: tot_psi, hk
    !
    ! !LOCAL VARIABLES:
    real(8) :: osmo_psi
    real(8) :: mat_psi
    real(8) :: head
    !---------------------------------------------------------------------

    associate ( &
                                      ! *** Input ***
    rgas    => universal%rgas    , &  ! Universal gas constant (J/K/mol)
    iv      => universal%iv      , &  ! Van't Hoff coefficient for NaCl
    denh2o  => universal%denh2o  , &  ! Water density (kg/m3)
    mmh2o   => universal%mmh2o   , &  ! Molecular mass of water (kg/mol)
    grav    => universal%grav      &  ! Gravitational acceleration (m/s2)
    )

    ! Head of pressure (MPa/m)

    head = denh2o * grav * 1.e-06

    ! Effective osmotic potential (MPa)

    osmo_psi = -1.0d0 * root_filter * sal * rgas * iv * soil_t * 1.e-06

    ! Soil matrix potential (MPa)

    mat_psi = psisat * moist ** (-1.0d0 * bsw)

    ! Total soil water potential (MPa)

    tot_psi = osmo_psi + mat_psi

    ! Soil hydraulic conductance (mm H2O/s2)

    hk = hksat * moist ** (2.0d0*bsw + 3.0d0)

    ! (mm H2O/s2) -> (m H2O/s) -> (m2 H2O/s/MPa)

    hk = hk * 1.e-03 / head

    ! (m2 H2O/s/MPa) -> (mmol H2O/m/s/MPa)

    hk = hk * denh2o / mmh2o * 1000.0d0

    end associate
  end subroutine soil_char

end module mod_soil_char
