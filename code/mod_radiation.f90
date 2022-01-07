module mod_radiation

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
  public :: extinction_coefficient
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine extinction_coefficient ( &
                                       ! *** Input ***
    n_spe                         , &  ! Number of mangrove species
    rhol_vis                      , &  ! Leaf relflectance to VIS (-)
    rhol_nir                      , &  ! Leaf relflectance to NIR (-)
    taul_vis                      , &  ! Leaf transmittance to VIS (-)
    taul_nir                      , &  ! Leaf transmittance to NIR (-)
    xl                            , &  ! Departure of leaf angle from spherical orientation (-)
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    sl_hgt                        , &  ! Solar elevation angle at midday (degree)
    tree_exist                    , &  ! Flag of tree presence
    la                            , &  ! Leaf area per tree (m2 leaf/tree)
    height                        , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
    bole                          , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
    pft                           , &  ! Species index (1: Rh, 2: Br)
                                       !
                                       ! *** Output ***
    dpai_layer                    , &  ! Layer leaf area index for each PFT for each forest layer (m2 leaf/me ground)
    dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    kbm_vis                       , &  ! Scattering adjusted light extinction coefficient for VIS (-)
    kbm_nir                       , &  ! Scattering adjusted light extinction coefficient for NIR (-)
    albvegb_vis                   , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
    albvegb_nir                   , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
    nbot                          , &  ! Index for bottom leaf layer
    ntop                            &  ! Index for top leaf layer
    )
    !
    ! !DESCRIPTION:
    ! Scattering adjusted light extinction coefficient for VIS/NIR
    !
    ! !USES:
    use data_structure, only : Max_no, Max_hgt, PFT_no
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: n_spe
    real, intent(in)     :: rhol_vis(:), rhol_nir(:), taul_vis(:)
    real, intent(in)     :: taul_nir(:), xl(:)
    integer, intent(in)  :: Max_loc
    real, intent(in)     :: sl_hgt
    logical, intent(in)  :: tree_exist(:)
    real(8), intent(in)  :: la(:)
    integer, intent(in)  :: height(:), bole(:), pft(:)
    real(8), intent(out) :: dpai_layer(:,:), dpai_layer_sum(:)
    real, intent(out)    :: kbm_vis(:), kbm_nir(:), albvegb_vis, albvegb_nir
    integer, intent(out) :: nbot, ntop
    !
    ! !LOCAL VARIABLES:
    real    :: solar_zen                     ! Solar zenith angle (radians)
    integer :: no                            ! Tree index
    real(8) :: dpai_tree                     ! Layer leaf area index for a tree (m2 leaf/m2 ground)
    integer :: p                             ! Species index
    integer :: i                             ! Layer index
    real, dimension(PFT_no) :: omega_vis     ! Leaf scattering coefficient for VIS
    real, dimension(PFT_no) :: omega_nir     ! Leaf scattering coefficient for NIR
    real, dimension(PFT_no) :: chil          ! Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
    real, dimension(PFT_no) :: phi1          ! Term in Ross-Goudriaan function for gdir
    real, dimension(PFT_no) :: phi2          ! Term in Ross-Goudriaan function for gdir
    real    :: gdir                          ! Relative projected area of leaf elements in the direction of solar beam
    real    :: kb                            ! Direct beam extinction coefficient
    real    :: albvegh_vis                   ! Vegetation albedo for VIS, horizontal leaves
    real    :: albvegh_nir                   ! Vegetation albedo for NIR, horizontal leaves
    integer :: canopy_index                  ! Canopy index
    integer :: j                             ! Sky angle index
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Solar zenith angle (radians)
    !---------------------------------------------------------------------

    ! Note: The effect of solar zenith angle on light attenuation was
    !       removed as in SEIB-DGVM.

!    solar_zen = (90.0d0 - sl_hgt) * (pi / 180.0d0)
    solar_zen = 0.0d0

    !---------------------------------------------------------------------
    ! Layer leaf area index for each PFT for each forest layer (m2 leaf/me ground)
    !---------------------------------------------------------------------

    dpai_layer(:,:) = 0.0d0
    do no = 1, Max_no
       if ( .not. tree_exist(no) ) cycle
       dpai_tree = la(no) / real(height(no) - bole(no))
       do i = bole(no)+1, height(no)
          dpai_layer(pft(no),i) = dpai_layer(pft(no),i) + dpai_tree
       end do
    end do
    dpai_layer(:,:) = dpai_layer(:,:) / real(Max_loc) / real(Max_loc)

    !------------------------------------------------------------------
    ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    !------------------------------------------------------------------

    do i = 1, Max_hgt
       dpai_layer_sum(i) = 0.0d0
       do p = 1, n_spe
          dpai_layer_sum(i) = dpai_layer_sum(i) + dpai_layer(p,i)
       end do
    end do

    !---------------------------------------------------------------------
    ! Leaf scattering coefficient (-)
    !---------------------------------------------------------------------

    do p = 1, n_spe
       omega_vis(p) = rhol_vis(p) + taul_vis(p)
       omega_nir(p) = rhol_nir(p) + taul_nir(p)
    end do

    !---------------------------------------------------------------------
    ! Light extinction coefficient
    !---------------------------------------------------------------------

    do p = 1, n_spe

       chil(p) = min(max(xl(p), -0.4d0), 0.6d0)
       if (abs(chil(p)) <= 0.01d0) chil(p) = 0.01d0

       phi1(p) = 0.5d0 - 0.633d0*chil(p) - 0.330d0*chil(p)*chil(p)
       phi2(p) = 0.877d0 * (1.0d0 - 2.0d0*phi1(p))

       gdir = phi1(p) + phi2(p) * cos(solar_zen)
       kb = gdir / cos(solar_zen)
       kb = min(kb, 40.0d0)

       ! Adjust for scattering

       kbm_vis(p) = kb * sqrt(1.0d0 - omega_vis(p))
       kbm_nir(p) = kb * sqrt(1.0d0 - omega_nir(p))

    end do

    !---------------------------------------------------------------------
    ! Vegetation albedo
    !---------------------------------------------------------------------

    ! Vegetation albedo, horizontal leaves
    ! Note: This equation is valid only for canopy with the kb (same leaf angle xl).

    albvegh_vis = (1.0d0 - sqrt(1.0d0 - omega_vis(1))) / (1.0d0 + sqrt(1.0d0   &
                  - omega_vis(1)))
    albvegh_nir = (1.0d0 - sqrt(1.0d0 - omega_nir(1))) / (1.0d0 + sqrt(1.0d0   &
                  - omega_nir(1)))

    ! (Direct beam) vegetation albedo, non-horizontal leaves
    ! Note: This equation is valid only for canopy with the kb (same leaf angle xl).

!   albvegb = 2.0d0 * kb / (kb + kd) * albvegh
    albvegb_vis = albvegh_vis
    albvegb_nir = albvegh_nir

    !---------------------------------------------------------------------
    ! Identifying top and bottom layer
    !---------------------------------------------------------------------

    ! Identifying canopy bottom layer

    canopy_index = 0
    do i = 1, Max_hgt
       if (dpai_layer_sum(i) > 0.0d0 .and. canopy_index == 0) then
          nbot = i
          canopy_index = 1
       end if
    end do

    ! Identifying canopy top layer

    canopy_index = 0
    do i = Max_hgt, 1, -1
       if (dpai_layer_sum(i) > 0.0d0 .and. canopy_index == 0) then
          ntop = i
          canopy_index = 1
       end if
    end do

  end subroutine extinction_coefficient

end module mod_radiation
