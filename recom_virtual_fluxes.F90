!===============================================================================
! MODULE: recom_virtual_fluxes
! 30.09.2026 OG
!
! Purpose: virtual surface fluxes of REcoM tracers (DIN, DIC, Alk, DSi, DFe, O2)
!          caused by the surface freshwater flux, analogous to FESOM's
!          virtual_salt (src/ice_oce_coupling.F90, subroutine oce_fluxes).
!
! Why: with which_ALE='linfs' the volume of the surface layer is fixed, so
!      evaporation, precipitation, runoff and sea-ice growth/melt do not
!      concentrate/dilute tracers. For salinity FESOM adds
!          virtual_salt = S_ref * water_flux
!      as a surface boundary flux. BGC tracers need the same term, otherwise
!      e.g. DIC and Alk do not follow salinity (wrong salinity-normalised
!      carbonate chemistry, spurious pCO2 patterns near rivers/ice edge).
!
!      With zstar/zlevel the freshwater flux changes the layer thickness and the
!      dilution is done by the ALE volume change -> virtual fluxes must be ZERO
!      there (otherwise double counting), EXCEPT under ice-shelf cavities, where
!      FESOM keeps the surface fixed (linfs-like) also in zstar.
!
! Sign convention (same as FESOM):
!      water_flux > 0 : water leaves the ocean (evaporation, ice formation)
!      vflux      > 0 : tracer flux INTO the ocean   [tracer unit * m/s]
!      -> evaporation concentrates, precipitation/runoff/melt dilutes.
!
! Global balancing (virt_balance_mode):
!      0 : none (not conservative, only for testing)
!      1 : subtract the global net uniformly over the open ocean (exactly what
!          FESOM does for virtual_salt)
!      2 : subtract the global net weighted with the local surface
!          concentration (keeps the correction small where the tracer is
!          depleted, e.g. DIN/DFe in oligotrophic gyres, avoids pushing
!          concentrations negative)
!      Cavity nodes are never modified by the balancing; their net flux is
!      compensated over the open ocean (as for virtual_salt).
!===============================================================================
module recom_virtual_fluxes
    use recom_declarations, only: wp
    implicit none
    private

    public :: recom_virtual_flux_single
    public :: recom_virtual_fluxes_all

contains

    !---------------------------------------------------------------------------
    ! Virtual flux for ONE tracer
    !---------------------------------------------------------------------------
    subroutine recom_virtual_flux_single(tr, water_flux, ref_val, ref_local, &
            use_virt_salt, use_cavity, balance_mode, ulevels_nod2D, areasvol, &
            ocean_area, myDim_nod2D, eDim_nod2D, MPI_COMM_FESOM, vflux)

        use recom_extra, only: integrate_nod_2d_recom

        implicit none

        real(kind=wp), intent(in),    dimension(:, :) :: tr          ! [nl-1, nod2D] tracer
        real(kind=wp), intent(in),    dimension(:)    :: water_flux  ! [m/s] balanced FESOM water_flux
        real(kind=wp), intent(in)                     :: ref_val     ! constant reference concentration
        logical,       intent(in)                     :: ref_local   ! use local surface concentration
        logical,       intent(in)                     :: use_virt_salt ! .true. for linfs
        logical,       intent(in)                     :: use_cavity
        integer,       intent(in)                     :: balance_mode
        integer,       intent(in),    dimension(:)    :: ulevels_nod2D
        real(kind=wp), intent(in),    dimension(:, :) :: areasvol
        real(kind=wp), intent(in)                     :: ocean_area  ! open-ocean area (excl. cavities)
        integer,       intent(in)                     :: myDim_nod2D, eDim_nod2D, MPI_COMM_FESOM
        real(kind=wp), intent(inout), dimension(:)    :: vflux       ! [tracer unit * m/s]

        integer :: n, nn
        real(kind=wp) :: cref, net, wsum
        real(kind=wp), allocatable, dimension(:) :: wgt

        nn = myDim_nod2D + eDim_nod2D
        vflux(1:nn) = 0.0_wp

        ! zstar/zlevel without cavities: dilution is done by the volume change
        if ((.not. use_virt_salt) .and. (.not. use_cavity)) return

        !_______________________________________________________________________
        ! 1. local virtual flux  C_ref * water_flux
        !    linfs : every node
        !    zstar : only cavity nodes (ulevels_nod2D > 1), which are linfs-like
        do n = 1, nn
            if ((.not. use_virt_salt) .and. (ulevels_nod2D(n) == 1)) cycle
            cref = ref_val
            if (ref_local) cref = tr(ulevels_nod2D(n), n)
            vflux(n) = cref * water_flux(n)
        end do

        if (balance_mode == 0) return

        !_______________________________________________________________________
        ! 2. global net [tracer unit * m3/s]
        call integrate_nod_2d_recom(vflux, net, MPI_COMM_FESOM, myDim_nod2D, &
                ulevels_nod2D, areasvol)

        !_______________________________________________________________________
        ! 3. remove it over the open ocean (cavity nodes untouched)
        select case (balance_mode)

        case (2)
            allocate(wgt(size(vflux)))
            wgt = 0.0_wp
            do n = 1, nn
                if (ulevels_nod2D(n) > 1) cycle
                wgt(n) = max(tr(1, n), 0.0_wp)
            end do
            call integrate_nod_2d_recom(wgt, wsum, MPI_COMM_FESOM, myDim_nod2D, &
                    ulevels_nod2D, areasvol)
            if (wsum > tiny(1.0_wp)) then
                do n = 1, nn
                    if (ulevels_nod2D(n) > 1) cycle
                    vflux(n) = vflux(n) - net * wgt(n) / wsum
                end do
            else
                ! tracer (nearly) absent everywhere -> fall back to uniform
                net = net / ocean_area
                do n = 1, nn
                    if (ulevels_nod2D(n) > 1) cycle
                    vflux(n) = vflux(n) - net
                end do
            end if
            deallocate(wgt)

        case default ! 1: uniform, identical to FESOM virtual_salt
            net = net / ocean_area
            do n = 1, nn
                if (ulevels_nod2D(n) > 1) cycle
                vflux(n) = vflux(n) - net
            end do

        end select

    end subroutine recom_virtual_flux_single

    !---------------------------------------------------------------------------
    ! Driver: loop over the tracers and fill virtual_din ... virtual_oxy
    !---------------------------------------------------------------------------
    subroutine recom_virtual_fluxes_all(tracers_info, num_tracers, water_flux, &
            use_virt_salt, use_cavity, ulevels_nod2D, areasvol, ocean_area, &
            myDim_nod2D, eDim_nod2D, MPI_COMM_FESOM, mype)

        use recom_declarations, only: tracer_ids
        use recom_glovar, only: tracers_info_type, virtual_din, virtual_dic, &
                virtual_alk, virtual_dsi, virtual_dfe, virtual_oxy
        use recom_config, only: use_virt_bgc, virt_ref_local, virt_balance_mode, &
                virt_ref_din, virt_ref_dic, virt_ref_alk, virt_ref_dsi, &
                virt_ref_dfe, virt_ref_oxy, recom_debug

        implicit none

        type(tracers_info_type), intent(in)          :: tracers_info
        integer,       intent(in)                    :: num_tracers
        real(kind=wp), intent(in), dimension(:)      :: water_flux
        logical,       intent(in)                    :: use_virt_salt, use_cavity
        integer,       intent(in), dimension(:)      :: ulevels_nod2D
        real(kind=wp), intent(in), dimension(:, :)   :: areasvol
        real(kind=wp), intent(in)                    :: ocean_area
        integer,       intent(in)                    :: myDim_nod2D, eDim_nod2D
        integer,       intent(in)                    :: MPI_COMM_FESOM, mype

        integer :: i, id

        if (.not. use_virt_bgc) return   ! arrays stay zero (allocated with source=0)

        do i = 1, num_tracers
            id = tracers_info%ids(i)

            if (id == tracer_ids%dissolved_inorganic_nitrogen) then
                call one(virt_ref_din, virtual_din)
            else if (id == tracer_ids%dissolved_inorganic_carbon) then
                call one(virt_ref_dic, virtual_dic)
            else if (id == tracer_ids%alkalinity) then
                call one(virt_ref_alk, virtual_alk)
            else if (id == tracer_ids%silica) then
                call one(virt_ref_dsi, virtual_dsi)
            else if (id == tracer_ids%iron) then
                call one(virt_ref_dfe, virtual_dfe)
            else if (id == tracer_ids%oxygen) then
                call one(virt_ref_oxy, virtual_oxy)
            end if
        end do

        if (recom_debug .and. mype == 0) then
            print *, achar(27) // '[36m' // '     --> recom_virtual_fluxes' // achar(27) // '[0m'
        end if

    contains

        subroutine one(ref_val, vflux)
            real(kind=wp), intent(in)                 :: ref_val
            real(kind=wp), intent(inout), dimension(:) :: vflux
            call recom_virtual_flux_single(tracers_info%data_pointers(i)%tracer_data, &
                    water_flux, ref_val, virt_ref_local, use_virt_salt, use_cavity, &
                    virt_balance_mode, ulevels_nod2D, areasvol, ocean_area, &
                    myDim_nod2D, eDim_nod2D, MPI_COMM_FESOM, vflux)
        end subroutine one

    end subroutine recom_virtual_fluxes_all

end module recom_virtual_fluxes
