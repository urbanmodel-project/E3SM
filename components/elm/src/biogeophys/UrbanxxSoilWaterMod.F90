module UrbanxxSoilWaterMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set soil water boundary conditions and compute hydrology
  ! for urban areas in the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use ColumnDataType       , only : col_ws, col_wf
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  public :: urbanxx_soilWater

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilWater(num_urbanl, num_urbanc, filter_urbanc, soilhydrology_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Set soil water boundary conditions for urban areas
    !
    use WaterFluxType, only : waterflux_type
    use WaterStateType, only : waterstate_type
    use ColumnType, only : col_pp
    use column_varcon, only : icol_road_perv
    use elm_varpar, only : nlevgrnd
    use SoilHydrologyType          , only : soilhydrology_type
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)  ! urban column filter
    type(soilhydrology_type) , intent(in) :: soilhydrology_vars
    real(r8)      , intent(in) :: dtime                ! time step (s)
    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status
    integer                              :: fc, c, j, idx, idx_perv, nlevbed
    integer(c_int)                       :: totalSize
    integer(c_int), dimension(2)         :: size2D
    logical(c_bool)                      :: isLayoutLeft
    real(c_double), allocatable, target  :: qflxInfl(:)
    real(c_double), allocatable, target  :: zwt(:)
    real(c_double), allocatable, target  :: qflxTran(:)
    real(c_double), allocatable, target  :: h2oLiq(:)
    real(c_double), allocatable, target  :: h2oIce(:)
    real(c_double), allocatable, target  :: h2oVol(:)

    associate(                             &
         qflx_infl    => col_wf%qflx_infl    , & ! Input: [real(r8) (:)] infiltration (mm H2O /s)
         qflx_rootsoi => col_wf%qflx_rootsoi , & ! Input: [real(r8) (:,:)] vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         nlev2bed     => col_pp%nlevbed      , & ! Input: [integer (:)] number of layers to bedrock
         h2osoi_ice   => col_ws%h2osoi_ice   , & ! Input: [real(r8) (:,:)] ice lens (kg/m2)
         h2osoi_vol   => col_ws%h2osoi_vol   , & ! Input: [real(r8) (:,:)] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liq   => col_ws%h2osoi_liq   , & ! Input: [real(r8) (:,:)] liquid water (kg/m2)
         zwt_col      => soilhydrology_vars%zwt_col & ! Input: [real(r8) (:)] water table depth (m)
         )

      ! Set infiltration flux (1D: per landunit)
      allocate(qflxInfl(num_urbanl))
      allocate(zwt(num_urbanl))

      ! Loop through urban columns and extract infiltration flux for pervious road
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)

        if (col_pp%itype(c) == icol_road_perv) then
           idx_perv = idx_perv + 1
           qflxInfl(idx_perv) = qflx_infl(c)
           zwt(idx_perv) = zwt_col(c)
        end if
      end do

      call UrbanSetInfiltrationFlux(urbanxx, c_loc(qflxInfl), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetWaterTableDepth(urbanxx, c_loc(zwt), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      deallocate(qflxInfl)
      deallocate(zwt)

      ! Set soil water content and transpiration flux (2D: per landunit x nlevgrnd)
      totalSize = num_urbanl * nlevgrnd
      size2D(1) = num_urbanl
      size2D(2) = nlevgrnd

      allocate(h2oLiq(totalSize))
      allocate(h2oIce(totalSize))
      allocate(h2oVol(totalSize))
      allocate(qflxTran(totalSize))

      ! Check Kokkos memory layout
      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      if (isLayoutLeft) then
        ! LayoutLeft: First dimension (landunits) varies fastest
        ! Iterate: layer (outer), landunits (inner)
        idx = 0
        do j = 1, nlevgrnd
          idx_perv = 0
          do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            if (col_pp%itype(c) == icol_road_perv) then
              idx_perv = idx_perv + 1
              idx = idx + 1
              nlevbed = nlev2bed(c)
              if (j <= nlevbed) then
                h2oLiq(idx) = h2osoi_liq(c, j)
                h2oIce(idx) = h2osoi_ice(c, j)
                h2oVol(idx) = h2osoi_vol(c, j)
                qflxTran(idx) = qflx_rootsoi(c, j)
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
                h2oVol(idx) = 0.0_r8
                qflxTran(idx) = 0.0_r8
              end if
            end if
          end do
        end do
      else
        ! LayoutRight: Last dimension (layers) varies fastest
        ! Iterate: landunits (outer), layer (inner)
        idx = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            nlevbed = nlev2bed(c)
            do j = 1, nlevgrnd
              idx = idx + 1
              if (j <= nlevbed) then
                h2oLiq(idx) = h2osoi_liq(c, j)
                h2oIce(idx) = h2osoi_ice(c, j)
                h2oVol(idx) = h2osoi_vol(c, j)
                qflxTran(idx) = qflx_rootsoi(c, j)
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
                h2oVol(idx) = 0.0_r8
                qflxTran(idx) = 0.0_r8
              end if
            end do
          end if
        end do
      end if

      call UrbanSetSoilLiquidWater(urbanxx, c_loc(h2oLiq), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilIceContent(urbanxx, c_loc(h2oIce), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilVolumetricWater(urbanxx, c_loc(h2oVol), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      deallocate(h2oLiq)
      deallocate(h2oIce)
      deallocate(h2oVol)

      call UrbanSetTranspirationFlux(urbanxx, c_loc(qflxTran), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      deallocate(qflxTran)

      call UrbanComputeHydrology(urbanxx, dtime, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

     end associate

  end subroutine urbanxx_soilWater

end module UrbanxxSoilWaterMod
