module recom_sinking
    !===============================================================================
    ! recom_sinking
    !   Vertical sinking / gravitational settling of REcoM particulate tracers
    !   (detritus, phytoplankton, diatoms, coccolithophores, phaeocystis, and
    !   their calcite/silica/chlorophyll analogues), including:
    !     - explicit vertical advection through the water column
    !       (ver_sinking_recom, 3rd-order DST flux-limited scheme)
    !     - deposition into the benthic layer (ver_sinking_recom_benthos)
    !     - remineralization flux back out of the benthos (diff_ver_recom_expl)
    !     - the Cram et al. (2018) ballasting parameterisation, which scales
    !       sinking speed by particle density and seawater viscosity
    !       (ballast, get_particle_density, get_seawater_viscosity)
    !===============================================================================
    implicit none
    private

    public :: diff_ver_recom_expl
    public :: ver_sinking_recom
    public :: ver_sinking_recom_benthos
    public :: ballast
    public :: get_particle_density
    public :: get_seawater_viscosity

contains

    !===============================================================================
    ! ver_sinking_recom_benthos
    !   Computes vertical sinking of particulate tracers into the benthos layer.
    !
    !   For each locally-owned open-ocean node the sinking velocity is selected
    !   by tracer ID, the downward flux through each layer interface is accumulated
    !   into str_bf, and the net material reaching the seafloor is added to the
    !   appropriate Benthos bucket (N, C, Si, Cal) and optionally to the MEDUSA
    !   sediment flux arrays.
    !
    !   Cavity nodes (ulevels_nod2D > 1) are skipped - there is no overlying
    !   water column above the ice-shelf cavity base that could deposit material
    !   onto a benthic boundary.
    !
    !   Note: second-detritus class (zoo2det) uses a separate sinking speed
    !   (VDet_zoo2) and is identified by hardcoded IDs (1025-1028) because no
    !   recom_det2_tracer_id array is defined yet.
    !   TODO: define recom_det2_tracer_id in recom_declarations for consistency
    !         with recom_det_tracer_id / recom_phy_tracer_id / recom_dia_tracer_id.
    !   TODO: replace all hardcoded tracer IDs with named constants from
    !         recom_declarations throughout this file.
    !===============================================================================
    subroutine ver_sinking_recom_benthos(tr_num, nl, ulevels_nod2D, nlevels_nod2D, zbar_3d_n, &
            nod_in_elem2D_num, nod_in_elem2D, nlevels, area, areasvol, tracer_id, tracer_data_values, &
            myDim_nod2d, str_bf, mype, MPI_COMM_FESOM, npes, sn, rn, s_mpitype_nod2D, &
            r_mpitype_nod2D, s_mpitype_nod3D, r_mpitype_nod3D, sPE, rPE, requests, nreq, dt)

        use recom_g_comm_auto, only: recom_exchange_nod

        use recom_declarations, only: wp, tracer_ids
        use recom_glovar, only: Benthos, Benthos_tr, SinkFlx, SinkFlx_tr

        use recom_config, only: allow_var_sinking, benthos_num, bottflx_num, ciso, Vdet, VPhy, &
            VDia, VDet_zoo2, enable_3zoo2det, use_MEDUSA, sedflx_num, recom_det_tracer_id, &
            recom_phy_tracer_id, recom_dia_tracer_id, SecondsPerDay, vdet_a

        use recom_ciso, only: ciso_organic_14

        implicit none

        integer, intent(in) :: tr_num, nl, tracer_id, myDim_nod2D
        integer, intent(in) :: mype, MPI_COMM_FESOM
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D
        integer, intent(in), dimension(:) :: nod_in_elem2D_num, nlevels
        integer, intent(in), dimension(:, :) :: nod_in_elem2D
        real(kind=WP), intent(in) :: dt
        real(kind=WP), intent(in), dimension(:, :) :: zbar_3d_n, area, areasvol, tracer_data_values
        real(kind=WP), intent(inout), dimension(:, :) :: str_bf

        ! These should all go into a dedicated REcoM type
        integer, intent(in) :: sn, rn, npes
        integer, intent(inout) :: nreq
        integer, intent(in), dimension(:) :: sPE, rPE
        integer, intent(inout), dimension(:) :: requests
        integer, intent(in), dimension(:), pointer :: s_mpitype_nod2D, r_mpitype_nod2D
        integer, intent(in), dimension(:, :, :), pointer :: s_mpitype_nod3D, r_mpitype_nod3D

        integer :: k
        integer :: nl1, ul1, nz, n
        real(kind=WP) :: Vben(nl), aux(nl - 1), add_benthos_2d(myDim_nod2D)
        integer :: nlevels_nod2D_minimum
        real(kind=WP) :: tv

        do n = 1, myDim_nod2D ! needs exchange_nod in the end

            if (ulevels_nod2D(n) > 1) cycle ! Cavity guard: no water column above cavity base

            nl1 = nlevels_nod2D(n) - 1   ! index of deepest layer centre
            ul1 = ulevels_nod2D(n)       ! index of shallowest (surface) layer centre

            aux = 0._WP
            Vben = 0._WP
            add_benthos_2d = 0._WP

            ! ----------------------------------------------------------------
            ! Select sinking velocity by tracer functional group
            ! ----------------------------------------------------------------
            if (any(recom_det_tracer_id == tracer_id)) Vben = Vdet
            if (any(recom_phy_tracer_id == tracer_id)) Vben = VPhy
            if (any(recom_dia_tracer_id == tracer_id)) Vben = VDia
            if (allow_var_sinking) then
                Vben = Vdet_a * abs(zbar_3d_n(:, n)) + Vben
            end if

            ! ----------------------------------------------------------------
            ! Second detritus class (macrozooplankton faecal pellets) uses a
            ! constant, separate sinking speed, overriding the above.
            ! ----------------------------------------------------------------
            if (enable_3zoo2det) then
                if (tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. & !idetz2n
                    tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. & !idetz2c
                    tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. & !idetz2si
                    tracer_id == tracer_ids%macrozooplankton_detrital_calcite) then !idetz2calc
                    Vben = VDet_zoo2
                end if
            end if

            ! conversion [m/d] --> [m/s] (vertical velocity, note that it is positive here)
            Vben = Vben / SecondsPerDay

            k = nod_in_elem2D_num(n)

            ! Screen minimum depth among neighbouring nodes around node n, so that
            ! the flux loop below does not run below the shallowest shared bottom.
            nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

            ! ----------------------------------------------------------------
            ! Downward mass flux through each layer, weighted by the change in
            ! cross-sectional area between adjacent layer interfaces.
            ! ----------------------------------------------------------------
            do nz = nlevels_nod2D_minimum, nl1
                tv = tracer_data_values(nz, n) * Vben(nz)
                aux(nz) = -tv * (area(nz, n) - area(nz + 1, n))
            end do

            do nz = ul1, nl1
                str_bf(nz, n) = str_bf(nz, n) &
                        + (aux(nz)) * dt / areasvol(nz, n) / (zbar_3d_n(nz, n) - zbar_3d_n(nz + 1, n))
                !!!!!!!!CHECK Maybe /area(nz,n) -> [mmol/m2]
                add_benthos_2d(n) = add_benthos_2d(n) - (aux(nz)) * dt
            end do

            ! ----------------------------------------------------------------
            ! Route accumulated deposition into the appropriate Benthos /
            ! SinkFlx bucket depending on tracer element (N, C, Si, Cal),
            ! plus isotope-resolved buckets (13C, 14C) when ciso is active.
            !
            ! Buckets: 1=N, 2=C, 3=Si, 4=Cal, 5=13C(DIC/POC), 6=13C(Cal),
            !          7=14C(DIC/POC), 8=14C(Cal)
            !
            ! __usetp branches accumulate per-tracer-index into *_tr buffers to
            ! keep global sums bit-reproducible under threaded tracer loops;
            ! these get summed into Benthos/SinkFlx across tr_num afterwards
            ! (in oce_ale_tracer.F90).
            ! ----------------------------------------------------------------

            !! * Particulate Organic Nitrogen *
            if (tracer_id == tracer_ids%phytoplankton_nitrogen .or. & !iphyn
                tracer_id == tracer_ids%detrital_nitrogen .or. & !idetn
                tracer_id == tracer_ids%diatom_nitrogen .or. & !idian
                tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen) then !idetz2n
#if defined(__usetp)
                Benthos_tr(n, 1, tr_num) = Benthos_tr(n, 1, tr_num) + add_benthos_2d(n) ![mmol]

                if (use_MEDUSA) then
                    SinkFlx_tr(n, 1, tr_num) = SinkFlx_tr(n, 1, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt ![mmol/m2]
                    ! SinkFlx has units mmol/timestep; MEDUSA needs mmol/m2/timestep,
                    ! hence the division by area here.
                end if
#else
                Benthos(n, 1) = Benthos(n, 1) + add_benthos_2d(n) ![mmol]

                if (use_MEDUSA) then
                    SinkFlx(n, 1) = SinkFlx(n, 1) + add_benthos_2d(n) / area(1, n) / dt ![mmol/m2]
                end if
#endif

            end if

            !! * Particulate Organic Carbon *
            if (tracer_id == tracer_ids%phytoplankton_carbon .or. & !iphyc
                tracer_id == tracer_ids%detrital_carbon .or. & !idetc
                tracer_id == tracer_ids%diatom_carbon .or. & !idiac
                tracer_id == tracer_ids%macrozooplankton_detrital_carbon) then !idetz2c
#if defined(__usetp)
                Benthos_tr(n, 2, tr_num) = Benthos_tr(n, 2, tr_num) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx_tr(n, 2, tr_num) = SinkFlx_tr(n, 2, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
#else
                Benthos(n, 2) = Benthos(n, 2) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx(n, 2) = SinkFlx(n, 2) + add_benthos_2d(n) / area(1, n) / dt ![mmol/m2]
                end if
#endif

            end if

            !! *Particulate Organic Silicon *
            if (tracer_id == tracer_ids%diatom_silica .or. & !idiasi
                    tracer_id == tracer_ids%detrital_silica .or. & !idetsi
                    tracer_id == tracer_ids%macrozooplankton_detrital_silica) then !idetz2si
#if defined(__usetp)
                Benthos_tr(n, 3, tr_num) = Benthos_tr(n, 3, tr_num) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx_tr(n, 3, tr_num) = SinkFlx_tr(n, 3, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
#else
                Benthos(n, 3) = Benthos(n, 3) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx(n, 3) = SinkFlx(n, 3) + add_benthos_2d(n) / area(1, n) / dt
                end if
#endif

            end if

            !! * Calcite *
            if (tracer_id == tracer_ids%phytoplankton_calcite .or. & !iphycal
                    tracer_id == tracer_ids%detrital_calcite .or. & !idetcal
                    tracer_id == tracer_ids%macrozooplankton_detrital_calcite) then !idetz2cal
#if defined(__usetp)
                Benthos_tr(n, 4, tr_num) = Benthos_tr(n, 4, tr_num) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx_tr(n, 4, tr_num) = SinkFlx_tr(n, 4, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
#else
                Benthos(n, 4) = Benthos(n, 4) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    SinkFlx(n, 4) = SinkFlx(n, 4) + add_benthos_2d(n) / area(1, n) / dt
                end if
#endif

            end if

            ! flux of 13C into the sediment
            if (ciso) then
                if (tracer_id == 1305 .or. & !iphyc_13
                        tracer_id == 1308 .or. & !idetc_13
                        tracer_id == 1314) then !idiac_13

#if defined(__usetp)
                    Benthos_tr(n, 5, tr_num) = Benthos_tr(n, 5, tr_num) + add_benthos_2d(n)

                    if (use_MEDUSA) then
                        SinkFlx_tr(n, 5, tr_num) = SinkFlx_tr(n, 5, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
#else
                    Benthos(n, 5) = Benthos(n, 5) + add_benthos_2d(n)
                    if (use_MEDUSA) then
                        SinkFlx(n, 5) = SinkFlx(n, 5) + add_benthos_2d(n) / area(1, n) / dt
                    end if
#endif

                end if

                if (tracer_id == 1320 .or. & !iphycal_13
                        tracer_id == 1321) then !idetcal_13

#if defined(__usetp)
                    Benthos_tr(n, 6, tr_num) = Benthos_tr(n, 6, tr_num) + add_benthos_2d(n)

                    if (use_MEDUSA) then
                        SinkFlx_tr(n, 6, tr_num) = SinkFlx_tr(n, 6, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
#else
                    Benthos(n, 6) = Benthos(n, 6) + add_benthos_2d(n)
                    if (use_MEDUSA) then
                        SinkFlx(n, 6) = SinkFlx(n, 6) + add_benthos_2d(n) / area(1, n) / dt
                    end if
#endif

                end if

            end if

            ! flux of 14C into the sediment
            if (ciso .and. ciso_organic_14) then
                if (tracer_id == 1405 .or. & !iphyc_14
                        tracer_id == 1408 .or. & !idetc_14
                        tracer_id == 1414) then !idiac_14

#if defined(__usetp)
                    Benthos_tr(n, 7, tr_num) = Benthos_tr(n, 7, tr_num) + add_benthos_2d(n)

                    if (use_MEDUSA) then
                        SinkFlx_tr(n, 7, tr_num) = SinkFlx_tr(n, 7, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
#else
                    Benthos(n, 7) = Benthos(n, 7) + add_benthos_2d(n)
                    if (use_MEDUSA) then
                        SinkFlx(n, 7) = SinkFlx(n, 7) + add_benthos_2d(n) / area(1, n) / dt
                    end if
#endif

                end if

                if (tracer_id == 1420 .or. & !iphycal_14
                        tracer_id == 1421) then !idetcal_14
#if defined(__usetp)
                    Benthos_tr(n, 8, tr_num) = Benthos_tr(n, 8, tr_num) + add_benthos_2d(n)

                    if (use_MEDUSA) then
                        SinkFlx_tr(n, 8, tr_num) = SinkFlx_tr(n, 8, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
#else
                    Benthos(n, 8) = Benthos(n, 8) + add_benthos_2d(n)
                    if (use_MEDUSA) then
                        SinkFlx(n, 8) = SinkFlx(n, 8) + add_benthos_2d(n) / area(1, n) / dt
                    end if
#endif
                end if

            end if

        end do

        ! ----------------------------------------------------------------
        ! Halo exchange of accumulated fluxes/benthos buckets across MPI
        ! partitions so that every rank has consistent boundary values.
        ! ----------------------------------------------------------------
        if (use_MEDUSA) then
            do n = 1, bottflx_num
#if defined(__usetp)
                call recom_exchange_nod(SinkFlx_tr(:, n, tr_num), npes, sn, rn, MPI_COMM_FESOM, &
                        mype, s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
#else
                call recom_exchange_nod(SinkFlx(:, n), npes, sn, rn, MPI_COMM_FESOM, mype, &
                        s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
#endif
            end do
        end if ! use_MEDUSA

        do n = 1, benthos_num
#if defined(__usetp)
            call recom_exchange_nod(Benthos_tr(:, n, tr_num), npes, sn, rn, MPI_COMM_FESOM, &
                    mype, s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
#else
            call recom_exchange_nod(Benthos(:, n), npes, sn, rn, MPI_COMM_FESOM, mype, &
                    s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
#endif
        end do

    end subroutine ver_sinking_recom_benthos


    !===============================================================================
    ! diff_ver_recom_expl
    !   Distributes benthic remineralization fluxes (from GloSed / GlodecayBenthos)
    !   back up into the bottom-most water-column cells as a source term
    !   (dtr_bf), for dissolved tracers (DIN, DIC, Alk, Si, DFe, O2, and their
    !   isotopologues).
    !
    !   Two flux sources depending on configuration:
    !     - use_MEDUSA & sedflx_num /= 0 : flux comes from GloSed (MEDUSA sediment
    !       model), converted from mol/time/area to mol/time by multiplying by
    !       area(1,:).
    !     - otherwise: flux comes from GlodecayBenthos (REcoM's own simple
    !       benthic remineralization), already combining terms as needed
    !       per tracer (e.g. DIC = organic decay + calcite dissolution;
    !       Alk = 2*calcite - Redfield N contribution; O2 = -redO2C*organic decay).
    !
    !   Area / cavity handling
    !   -----------------------
    !   bottom_flux(n) is built from GloSed(:,k) * area(1,n) (or, in the
    !   non-MEDUSA branch, is already a flux referenced to the surface),
    !   i.e. it is a *total* flux expressed relative to the surface
    !   cross-sectional area area(1,n) -- a quantity that is well-defined
    !   and cavity-independent (it is simply the node's surface footprint).
    !
    !   The redistribution across layers must therefore divide by that same
    !   reference area, area(1,n), to first recover a physically meaningful
    !   per-unit-area flux before spreading it across the topography-
    !   following interface-area differences (area(nz,n) - area(nz+1,n)).
    !   Only the final step -- converting the redistributed flux into a
    !   tracer tendency for a specific cell nz -- should divide by that
    !   cell's own volume-area areasvol(nz,n), since that is what turns a
    !   flux into a per-volume concentration tendency for cell nz.
    !===============================================================================
    subroutine diff_ver_recom_expl(nl, ulevels_nod2D, nlevels_nod2D, nod_in_elem2D_num, &
            nod_in_elem2D, nlevels, area, areasvol, hnode_new, tracer_id, myDim_nod2d, &
            eDim_nod2D, mype, MPI_COMM_FESOM, dtr_bf, dt)

        use recom_declarations, only: wp, tracer_ids
        use recom_glovar, only: GloSed, glodecayBenthos
        use recom_config, only: ciso, use_MEDUSA, sedflx_num, redO2C, Fe2N_benthos

        implicit none

        integer, intent(in) :: myDim_nod2d, eDim_nod2D, mype, MPI_COMM_FESOM
        integer, intent(in) :: nl, tracer_id
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D
        integer, intent(in), dimension(:) :: nod_in_elem2D_num, nlevels
        real(kind=WP), intent(in) :: dt
        integer, intent(in), dimension(:, :) :: nod_in_elem2D
        real(kind=WP), intent(in), dimension(:, :) :: area, areasvol, hnode_new
        real(kind=WP), intent(inout), dimension(:, :) :: dtr_bf

        integer :: k
        integer :: nl1, nz, n, ul1
        real(kind=WP) :: vd_flux(nl)
        integer :: nlevels_nod2D_minimum
        real(kind=WP) :: bottom_flux(myDim_nod2D + eDim_nod2D)

        bottom_flux = 0._WP

#if defined(__recom)
        ! ----------------------------------------------------------------
        ! Select the per-node bottom flux [mol/time] for this dissolved
        ! tracer, from either the MEDUSA sediment model or REcoM's simple
        ! benthic decay bucket. Both are total fluxes referenced to the
        ! node's surface area area(1,:) -- see module docstring.
        ! ----------------------------------------------------------------
        if (use_MEDUSA .and. (sedflx_num /= 0)) then
            ! Note (CV/OG): GloSed is a flux per area (mol/time/area); multiply
            ! by area(1,:) to match the mol/time convention used below.
            select case (tracer_id)
            case (1001)
                bottom_flux = GloSed(:, 1) * area(1, :) ! DIN     !CHECK: (areasvol(ul1,:) OG: 02.06.2026
            case (1002)
                bottom_flux = GloSed(:, 2) * area(1, :) ! DIC
            case (1003)
                bottom_flux = GloSed(:, 3) * area(1, :) ! Alk
            case (1018)
                bottom_flux = GloSed(:, 4) * area(1, :) ! Si
            case (1019)
                bottom_flux = GloSed(:, 1) * Fe2N_benthos * area(1, :) ! DFe, scaled from N flux
            case (1022)
                bottom_flux = GloSed(:, 5) * area(1, :) ! O2
            case (1302)
                if (ciso) then
                    bottom_flux = GloSed(:, 6) * area(1, :) ! DIC_13 (organic + calcite)
                end if
            case (1402)
                if (ciso) then
                    bottom_flux = GloSed(:, 7) * area(1, :) ! DIC_14 (organic + calcite)
                end if
            case default
                if (mype == 0) then
                    write(*, *) 'check specified in boundary conditions'
                    write(*, *) 'the model will stop!'
                end if
                call MPI_ABORT(MPI_COMM_FESOM, 1)
                stop
            end select
        else
            select case (tracer_id)
            case (1001)
                bottom_flux = GlodecayBenthos(:, 1) !*** DIN [mmolN/m^2/s] ***
            case (1002)
                !*** DIC = organic-matter decay + CaCO3 dissolution ***
                bottom_flux = GlodecayBenthos(:, 2) + GlodecayBenthos(:, 4)
            case (1003)
                !*** Alk: +2 per mol CaCO3 dissolved, -1.0625 per mol N remineralized
                !*** (Redfield-derived alkalinity correction for nitrate uptake) ***
                bottom_flux = GlodecayBenthos(:, 4) * 2.0_WP - 1.0625_WP * GlodecayBenthos(:, 1)
            case (1018)
                bottom_flux = GlodecayBenthos(:, 3) !*** Si ***
            case (1019)
                bottom_flux = GlodecayBenthos(:, 1) * Fe2N_benthos !*** DFe, scaled from N ***
            case (1022)
                bottom_flux = -GlodecayBenthos(:, 2) * redO2C !*** O2 consumed by C remin ***
            case (1302)
                if (ciso) then
                    !*** DIC_13: organic decay + calcite dissolution ***
                    bottom_flux = GlodecayBenthos(:, 5) + GlodecayBenthos(:, 6)
                end if
            case (1402)
                if (ciso) then
                    !*** DIC_14: organic decay + calcite dissolution ***
                    bottom_flux = GlodecayBenthos(:, 7) + GlodecayBenthos(:, 8)
                end if
            case default
                if (mype == 0) then
                    write(*, *) 'check specified in boundary conditions'
                    write(*, *) 'the model will stop!'
                end if
                call MPI_ABORT(MPI_COMM_FESOM, 1)
                stop
            end select
        end if ! (use_MEDUSA .and. (sedflux_num .gt. 0))
#endif

        do n = 1, myDim_nod2D

            nl1 = nlevels_nod2D(n) - 1
            ul1 = ulevels_nod2D(n)

            vd_flux = 0._WP

            k = nod_in_elem2D_num(n)
            ! Screen minimum depth among neighbouring nodes around node n
            nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

            !_______________________________________________________________________
            ! Bottom flux, distributed across layers proportional to the change
            ! in cross-sectional area between interfaces, then normalized by the
            ! surface-layer volume-area.
            !
            ! NOTE / NEEDS CONFIRMATION: every layer nz here divides by
            ! areasvol(ul1, n) -- the *surface* (ul1) layer's volume-area -- not
            ! areasvol(nz, n), the volume-area of the layer the flux is actually
            ! being assigned to. This was already flagged by OG
            ! (see "CHECK" comment above on the GloSed DIN line) as needing review.
            ! Left unchanged; verify against the intended discretization before
            ! trusting sub-surface benthic flux distribution.
            do nz = nlevels_nod2D_minimum, nl1
                vd_flux(nz) = (area(nz, n) - area(nz + 1, n)) * bottom_flux(n) / (area(1, n))
            end do
            nz = nl1
            vd_flux(nz + 1) = (area(nz + 1, n)) * bottom_flux(n) / (area(1, n))

            !_______________________________________________________________________
            ! Add bottom flux into the tracer tendency (rhs). Each cell nz only
            ! receives flux through its lower (bottom) face, vd_flux(nz+1).
            do nz = ul1, nl1
                dtr_bf(nz, n) = dtr_bf(nz, n) &
                        + vd_flux(nz + 1) * dt / areasvol(nz, n) / hnode_new(nz, n)
            end do
        end do
    end subroutine diff_ver_recom_expl


    !===============================================================================
    ! ver_sinking_recom
    !   Explicit vertical advection ("sinking") of a particulate tracer through
    !   the water column using a 3rd-order Direct Space-Time (DST3) scheme with
    !   flux limiting (a TVD-type limiter), so sinking fronts stay sharp without
    !   introducing spurious oscillations.
    !
    !   Sinking speed depends on tracer functional type (constant base speeds
    !   VDet/VPhy/VDia/VCocco/VPhaeo/VDet_zoo2), optionally combined with:
    !     - allow_var_sinking : depth-proportional acceleration (detritus only)
    !     - use_ballasting    : density/viscosity-scaled speed following
    !                           Cram et al. (2018), computed in `ballast`
    !
    !   Only tracers with Vsink > 0.1 m/day are advected; slower tracers are
    !   left untouched by this routine (their sinking contribution is ~0).
    !===============================================================================
    subroutine ver_sinking_recom(tr_num, nl, ulevels_nod2D, nlevels_nod2D, zbar_3d_n, z_3d_n, &
            nod_in_elem2D_num, nod_in_elem2D, nlevels, area, areasvol, hnode, hnode_new, &
            tracer_id, tracer_data_values, myDim_nod2d, vert_sink, dt)

        use REcoM_declarations, only: wp, tracer_ids

        use REcoM_GloVar, only: sinkvel1_tr, sinkvel2_tr, scaling_visc_3D, scaling_density1_3D, &
            scaling_density2_3D

        use recom_config, only: allow_var_sinking, depth_scaling1, depth_scaling2, &
            max_sinking_velocity, use_ballasting, enable_3zoo2det, VDet_zoo2, enable_coccos, &
            secondsPerDay, vcocco, vdet, vdet_a, vdia, vphaeo, vphy, w_ref1, w_ref2

        implicit none

        integer, intent(in) :: tr_num, myDim_nod2D
        integer, intent(in) :: nl, tracer_id
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D
        integer, intent(in), dimension(:) :: nod_in_elem2D_num, nlevels
        integer, intent(in), dimension(:, :) :: nod_in_elem2D
        real(kind=WP), intent(in) :: dt
        real(kind=WP), intent(in), dimension(:, :) :: zbar_3d_n, z_3d_n, area, areasvol, hnode, &
                hnode_new
        real(kind=WP), intent(in), dimension(:, :) :: tracer_data_values
        real(kind=WP), intent(inout), dimension(:, :) :: vert_sink

        integer :: nz, nzmin, nzmax, n, k, nlevels_nod2D_minimum
        real(kind=wp) :: wM, wPs
        real(kind=wp) :: Rjp, Rj, Rjm

        real(kind=wp) :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
        real(kind=wp) :: dt_sink
        real(kind=wp) :: Vsink, tv
        real(kind=wp), save :: onesixth = 1.d0 / 6.d0

        real(kind=wp), dimension(nl) :: Wvel_flux, vd_flux, dz_trr

        ! ----------------------------------------------------------------
        ! Select base sinking speed [m/day] by tracer functional group.
        ! Any tracer not matched below keeps Vsink = 0 and is skipped by
        ! the "Vsink > 0.1" guard further down.
        ! ----------------------------------------------------------------
        Vsink = 0.0_WP

        if (tracer_id == tracer_ids%detrital_nitrogen .or. & ! idetn
                tracer_id == tracer_ids%detrital_carbon .or. & ! idetc
                tracer_id == tracer_ids%detrital_silica .or. & ! idetsi
                tracer_id == tracer_ids%detrital_calcite) then ! idetcal
            Vsink = VDet

        elseif (tracer_id == tracer_ids%phytoplankton_nitrogen .or. & ! iphyn
                    tracer_id == tracer_ids%phytoplankton_carbon .or. & ! iphyc
                    tracer_id == tracer_ids%phytoplankton_chlorophyll) then ! ipchl
            Vsink = VPhy

        elseif (tracer_id == tracer_ids%diatom_nitrogen .or. & ! idian
                    tracer_id == tracer_ids%diatom_carbon .or. & ! idiac
                    tracer_id == tracer_ids%diatom_silica .or. & ! idiasi
                    tracer_id == tracer_ids%diatom_chlorophyll) then ! idchl
            Vsink = VDia

        elseif (enable_coccos .and. &
                    (tracer_id == tracer_ids%coccolithophore_nitrogen .or. & ! icocn
                    tracer_id == tracer_ids%coccolithophore_carbon .or. & ! icocc
                    tracer_id == tracer_ids%coccolithophore_chlorophyll)) then ! icchl
            Vsink = VCocco

        elseif (enable_coccos .and. &
                    (tracer_id == tracer_ids%phaeocystis_nitrogen .or. & ! iphan
                    tracer_id == tracer_ids%phaeocystis_carbon .or. & ! iphac
                    tracer_id == tracer_ids%phaeocystis_chlorophyll)) then ! iphachl
            Vsink = VPhaeo

        elseif (tracer_id == tracer_ids%phytoplankton_calcite) then ! iphycal
            ! Calcite is produced by either coccolithophores or (in simpler
            ! configs without coccos) generic phytoplankton, so it inherits
            ! whichever group's sinking speed applies.
            if (enable_coccos) then
                Vsink = VCocco
            else
                Vsink = VPhy
            end if

        elseif (enable_3zoo2det .and. &
                (tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. & !idetz2n
                 tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. & !idetz2c
                 tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. & !idetz2si
                 tracer_id == tracer_ids%macrozooplankton_detrital_calcite)) then !idetz2calc
            Vsink = VDet_zoo2

        end if

        !! ---- Skip advection entirely for negligibly slow-sinking tracers
        if (Vsink > 0.1) then

            do n = 1, myDim_nod2D
                if (ulevels_nod2D(n) > 1) cycle   ! skip cavity nodes
                nzmin = ulevels_nod2D(n)
                nzmax = nlevels_nod2D(n) - 1

                ! Distance between tracer points; surface and bottom cells use
                ! half the local layer thickness (ghost-cell-like treatment).
                dz_trr = 0.0d0
                dz_trr(nzmin + 1:nzmax) = abs(Z_3d_n(nzmin:nzmax - 1, n) &
                        - Z_3d_n(nzmin + 1:nzmax, n))
                dz_trr(nzmin) = hnode(nzmin, n) / 2.0d0
                dz_trr(nzmax + 1) = hnode(nzmax, n) / 2.0d0

                Wvel_flux(nzmin:nzmax + 1) = 0.d0 ! Vertical sinking velocity at interfaces

                ! ----------------------------------------------------------------
                ! Build per-interface sinking velocity Wvel_flux(nz) [m/s],
                ! negative = downward (surface -> bottom).
                ! ----------------------------------------------------------------
                do nz = nzmin, nzmax + 1

                    Wvel_flux(nz) = -Vsink / SecondsPerDay ! base case: constant speed

                    if (allow_var_sinking) then
                        ! Depth-proportional speed-up (applies the same
                        ! detritus-derived Vdet_a term to whichever tracer is
                        ! currently being processed -- consistent with the
                        ! main detritus branch below, since ballasting further
                        ! overwrites this for detritus specifically).
                        Wvel_flux(nz) = -((Vdet_a * abs(zbar_3d_n(nz, n)) / SecondsPerDay) &
                                + Vsink / SecondsPerDay)

                        if (use_ballasting) then
                            ! Ballasting (Cram et al. 2018) overrides the speed
                            ! for the *first* detritus class only, using
                            ! density/viscosity scaling factors precomputed by
                            ! `ballast`.
                            if (tracer_id == tracer_ids%detrital_nitrogen .or. & !idetn
                                    tracer_id == tracer_ids%detrital_carbon .or. & !idetc
                                    tracer_id == tracer_ids%detrital_silica .or. & !idetsi
                                    tracer_id == tracer_ids%detrital_calcite) then !idetcal
                                Wvel_flux(nz) = w_ref1 * scaling_density1_3D(nz, n) &
                                        * scaling_visc_3D(nz, n)

                                if (depth_scaling1 > 0.0) Wvel_flux(nz) = Wvel_flux(nz) &
                                        + (depth_scaling1 * abs(zbar_3d_n(nz, n)))

                                ! Cap at max_sinking_velocity for numerical stability
                                if (abs(Wvel_flux(nz)) > max_sinking_velocity) Wvel_flux(nz) = &
                                        max_sinking_velocity

                                ! [m/d] -> [m/s], sign flipped to point downward (negative)
                                Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz) / SecondsPerDay
                            end if
                        end if
                    end if

                    !! ---- Second detritus class always sinks at a constant speed
                    ! Parenthesized explicitly: .and. binds tighter than .or. in Fortran, so
                    ! without these parens this read as (enable_3zoo2det .and. ...nitrogen) .or.
                    ! ...carbon .or. ...silica .or. ...calcite instead of the intended
                    ! enable_3zoo2det .and. (...nitrogen .or. ...carbon .or. ...silica .or.
                    ! ...calcite). Harmless today only because disabled tracer_ids default to a
                    ! -1 sentinel that a real tracer_id never matches, but keep the parens.
                    if (enable_3zoo2det .and. &
                        (tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. & !idetz2n
                         tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. & !idetz2c
                         tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. & !idetz2si
                         tracer_id == tracer_ids%macrozooplankton_detrital_calcite)) then !idetz2calc
                        Wvel_flux(nz) = -VDet_zoo2 / SecondsPerDay

                        if (use_ballasting) then
                            Wvel_flux(nz) = w_ref2 * scaling_density2_3D(nz, n) &
                                    * scaling_visc_3D(nz, n)

                            if (depth_scaling2 > 0.0) Wvel_flux(nz) = Wvel_flux(nz) &
                                    + (depth_scaling2 * abs(zbar_3d_n(nz, n)))

                            if (abs(Wvel_flux(nz)) > max_sinking_velocity) Wvel_flux(nz) = &
                                    max_sinking_velocity

                            Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz) / SecondsPerDay ! [m/d] -> [m/s]
                        end if

                    end if

                    ! Diagnostic storage of final sinking velocity for the
                    ! calcite tracers of each detritus class (used elsewhere,
                    ! e.g. for output / other modules that need the actual
                    ! ballasted speed rather than recomputing it).
                    if (tracer_id == tracer_ids%detrital_calcite) then
                        Sinkvel1_tr(nz, n, tr_num) = Wvel_flux(nz)
                    end if

                    if (enable_3zoo2det .and. &
                        tracer_id == tracer_ids%macrozooplankton_detrital_calcite) then
                        Sinkvel2_tr(nz, n, tr_num) = Wvel_flux(nz) !idetz2calc
                    end if

                end do

                dt_sink = dt
                vd_flux = 0.0d0

                ! ----------------------------------------------------------------
                ! 3rd-order DST (Direct Space Time) advection scheme with a TVD
                ! flux limiter (ported from the original REcoM code). Computes
                ! the advective flux vd_flux(nz) through each interface nz.
                ! ----------------------------------------------------------------
                k = nod_in_elem2D_num(n)
                ! Screen minimum depth among neighbouring nodes around node n.
                ! NOTE: computed but not currently used to bound the loop below
                ! (loop runs nzmin+1..nzmax); kept for parity with other
                ! routines in this file and potential future use -- harmless
                ! as-is, but flagging as dead computation.
                nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

                vd_flux(nzmin:nzmax + 1) = 0.0_WP

                do nz = nzmax, nzmin + 1, -1

                    Rjp = tracer_data_values(nz, n) - tracer_data_values(min(nz + 1, nzmax), n)
                    Rj = tracer_data_values(max(nzmin, nz - 1), n) - tracer_data_values(nz, n)
                    Rjm = tracer_data_values(max(nzmin, nz - 2), n) &
                            - tracer_data_values(max(nzmin, nz - 1), n)

                    ! Courant number: [m/s] * [s] / [m]
                    cfl = abs(Wvel_flux(nz) * dt_sink / dz_trr(nz))

                    wPs = Wvel_flux(nz) + abs(Wvel_flux(nz)) ! positive part of velocity
                    wM = Wvel_flux(nz) - abs(Wvel_flux(nz))  ! negative part of velocity

                    d0 = (2.d0 - cfl) * (1.d0 - cfl) * onesixth
                    d1 = (1.d0 - cfl * cfl) * onesixth

                    ! Upwind-side limiter ratio/coefficient (psiP) and
                    ! downwind-side limiter ratio/coefficient (psiM). Both are
                    ! clamped into [0,1] and further bounded by a CFL-dependent
                    ! TVD envelope to guarantee monotonicity.
                    thetaP = Rjm / (1.d-20 + Rj)
                    psiP = d0 + d1 * thetaP
                    psiP = max(0.d0, min(min(1.d0, psiP), &
                            (1.d0 - cfl) / (1.d-20 + cfl) * thetaP))

                    thetaM = Rjp / (1.d-20 + Rj)
                    psiM = d0 + d1 * thetaM
                    ! BUG FIX: denominator was "1.d-20 - cfl" in the original code,
                    ! asymmetric with psiP's "1.d-20 + cfl" for no apparent physical
                    ! reason. Since cfl >= 0, "1.d-20 - cfl" shrinks (and can even
                    ! flip the sign of) the denominator instead of safely guarding
                    ! against division by zero, which would corrupt the limiter
                    ! envelope for psiM. Corrected to "+ cfl" to match psiP.
                    psiM = max(0.d0, min(min(1.d0, psiM), &
                            (1.d0 - cfl) / (1.d-20 + cfl) * thetaM))
!                            (1.d0 - cfl) / (1.d-20 - cfl) * thetaM))
                    tv = (0.5 * wPs * (tracer_data_values(nz, n) + psiM * Rj) + &
                            0.5 * wM * (tracer_data_values(max(nzmin, nz - 1), n) + psiP * Rj))
                    vd_flux(nz) = -tv * area(nz, n)
                end do

                ! ----------------------------------------------------------------
                ! Legacy simple first-order upwind scheme, kept for reference /
                ! debugging but permanently disabled (.false.). Switch this on
                ! manually (or wire up a proper namelist flag, see FIXME) to
                ! fall back to upwind advection instead of DST3.
                ! ----------------------------------------------------------------
                if (.false.) then ! simple upwind FIXME: use a flag here later

                    vd_flux(nzmin) = 0.0_WP      ! no flux through the free surface
                    vd_flux(nzmax + 1) = 0.0_WP  ! no flux through the seafloor here
                                                  ! (handled separately by the benthos routine)

                    k = nod_in_elem2D_num(n)
                    nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

                    do nz = nzmin + 1, nzmax !nlevels_nod2D_minimum-1
                        !         tv = trarr(nz,n)                                ! simple scheme
                        !    - test1
                        !         tv = 0.5_WP*(trarr(nz-1,n)+trarr(nz,n))         ! consider both
                        ! layers - test2
                        !         tv = tv*Wvel_flux(nz) ! Wvel_flux is negative
                        tv = -0.5 * & ! - test3
                                (tracer_data_values(nz - 1, n) &
                                * (Wvel_flux(nz) - abs(Wvel_flux(nz))) &
                                + &
                                tracer_data_values(nz, n) * (Wvel_flux(nz) + abs(Wvel_flux(nz))))
                        vd_flux(nz) = tv * area(nz, n)
                    end do
                end if ! simple upwind

                ! ----------------------------------------------------------------
                ! Accumulate the net flux divergence into the tracer tendency.
                ! ----------------------------------------------------------------
                do nz = nzmin, nzmax
                    vert_sink(nz, n) = vert_sink(nz, n) &
                            + (vd_flux(nz) - vd_flux(nz + 1)) * dt / areasvol(nz, n) &
                            / hnode_new(nz, n)
                end do
            end do
        end if ! Vsink > 0.1

    end subroutine ver_sinking_recom


    !===============================================================================
    ! ballast
    !   Computes per-node, per-level scaling factors for the ballasted sinking
    !   velocity parameterisation (Cram et al. 2018).
    !
    !   Scaling factors:
    !     scaling_density1/2_3D - ratio of (particle - seawater) excess density
    !                              to a reference excess density; set to 1.0 if
    !                              density scaling is disabled or particle is less
    !                              dense than seawater.
    !     scaling_visc_3D       - ratio of reference viscosity to local seawater
    !                              viscosity; set to 1.0 if viscosity scaling is
    !                              disabled or viscosity is zero.
    !
    !   The actual sinking velocities are assembled in ver_sinking_recom using
    !   w_ref1/2 * scaling_density * scaling_visc (+ optional depth term).
    !===============================================================================
    subroutine ballast(myDim_nod2D, ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, Z_3d_n, &
            tracer_data_values_1, tracer_data_values_2, rad)

        use recom_config, only: enable_3zoo2det, use_density_scaling, use_viscosity_scaling, &
            rho_ref_part, rho_ref_water, visc_ref_water, tiny

        use recom_glovar, only: scaling_density1_3D, scaling_density2_3D, rho_particle1, &
            rho_particle2, seawater_visc_3D, scaling_visc_3D

        use recom_declarations, only: wp
        use mdepth2press, only: depth2press
        use gsw_mod_toolbox, only: gsw_sa_from_sp, gsw_ct_from_pt, gsw_rho

        implicit none

        integer, intent(in) :: myDim_nod2D
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D
        real(kind=WP), intent(in) :: rad
        real(kind=WP), intent(in), dimension(:, :) :: geo_coord_nod2D, Z_3d_n
        real(kind=WP), intent(in), dimension(:, :) :: tracer_data_values_1, tracer_data_values_2

        integer       :: row, k, nzmin, nzmax
        real(kind=wp) :: depth_pos(1)    ! Layer mid-depth [m, positive]
        real(kind=wp) :: pres(1)         ! Pressure at layer mid-depth [dbar]
        real(kind=wp) :: sa(1)           ! Absolute salinity [g/kg]
        real(kind=wp) :: ct(1)           ! Conservative temperature [deg C]
        real(kind=wp) :: rho_seawater(1) ! In-situ seawater density [kg/m3]
        real(kind=wp) :: Lon_degree(1)
        real(kind=wp) :: Lat_degree(1)

        ! For ballasting, calculate scaling factors here and pass them to FESOM,
        ! where sinking velocities are computed from these factors plus a reference speed.
        !
        ! Sinking velocity depends on:
        !   a) particle composition (= density)
        !   b) sea-water viscosity
        !   c) depth (currently small detritus only)
        !   d) a constant reference sinking speed
        do row = 1, myDim_nod2D

            nzmin = ulevels_nod2D(row)

            ! BUG FIX: nzmax was set to nlevels_nod2D(row) (the number of layer
            ! *interfaces*), one too many compared to every other routine in this
            ! file, which consistently uses nlevels_nod2D(n) - 1 as the index of
            ! the deepest layer *centre* (see ver_sinking_recom, get_particle_density,
            ! diff_ver_recom_expl). With the old value, the "extend one level below
            ! bottom" step further down wrote to scaling_*_3D(nzmax+1, row), i.e.
            ! two levels past the last valid layer centre, and the main k-loop
            ! computed seawater density/pressure one level too deep as well.
            ! Corrected to be consistent with the rest of the module.
            nzmax = nlevels_nod2D(row) - 1

            Lon_degree(1) = geo_coord_nod2D(1, row) / rad
            Lat_degree(1) = geo_coord_nod2D(2, row) / rad

            do k = nzmin, nzmax
                ! --- in-situ seawater density at this node/level (TEOS-10 / GSW) ---
                depth_pos(1) = abs(Z_3d_n(k, row))
                ! Pass length-1 arrays (not scalar elements) to depth2press so the
                ! interface matches its expected array dummy arguments.
                call depth2press(depth_pos(1), Lat_degree(1), pres, 1)
                ! Pass length-1 array literals for scalar tracer values so all
                ! GSW arguments are consistently rank-1 arrays of length 1.
                sa = gsw_sa_from_sp(tracer_data_values_2(k, row), pres, Lon_degree(1), Lat_degree(1))
                ct = gsw_ct_from_pt(sa, tracer_data_values_1(k, row))
                rho_seawater = gsw_rho(sa, ct, pres)

                ! Default: no density effect (neutral scaling)
                scaling_density1_3D(k, row) = 1.0_WP
                scaling_density2_3D(k, row) = 1.0_WP

                ! MERGE-REVIEW: int_recom's ballast only recomputed scaling_density1_3D/
                ! scaling_density2_3D from their 1.0 default when the detrital-carbon tracer's
                ! own concentration exceeded a 0.001 guard (tracers%data(tr_num)%values(k,row)
                ! > 0.001), avoiding wild scaling ratios from near-zero detrital carbon -- see
                ! the commented-out remnants of that guard directly below. This subroutine no
                ! longer receives per-tracer identity/concentration (it's now called once per
                ! timestep with T/S, from oce_ale_tracer.F90, not per-tracer with a tr_num as
                ! in int_recom), so restoring the guard's intent would mean threading DetC/
                ! DetZ2C concentration into this subroutine's argument list -- left unresolved
                ! rather than guessing at that redesign.
                if (use_density_scaling) then
                    ! Ratio of the particle's excess density (relative to local
                    ! seawater) to the reference excess density used to define
                    ! w_ref1/2 -- scales the reference sinking speed up/down.
                    !if (tracers%data(tr_num)%ID ==1008)then !idetc
                    !if (tracers%data(tr_num)%values(k,row)>0.001) then ! only apply ballasting
                    scaling_density1_3D(k, row) = &
                        (rho_particle1(k, row) - rho_seawater(1)) / (rho_ref_part - rho_ref_water)

                    !endif
                    !endif
                    if (enable_3zoo2det) &
                    scaling_density2_3D(k, row) = &
                        (rho_particle2(k, row) - rho_seawater(1)) / (rho_ref_part - rho_ref_water)
                end if

                scaling_visc_3D(k, row) = 1.0
                if (use_viscosity_scaling .and. seawater_visc_3D(k, row) /= 0.0_WP) then
                    scaling_visc_3D(k, row) = visc_ref_water / seawater_visc_3D(k, row)
                end if

            end do

            ! Extend bottom values one level down for flux-point access in
            ! ver_sinking_recom (which indexes up to nzmax+1 at the seafloor
            ! interface).
            scaling_density1_3D(nzmax + 1, row) = scaling_density1_3D(nzmax, row)
            scaling_visc_3D(nzmax + 1, row) = scaling_visc_3D(nzmax, row)
            if (enable_3zoo2det) then
                scaling_density2_3D(nzmax + 1, row) = scaling_density2_3D(nzmax, row)
            end if
        end do

    end subroutine ballast


    !===============================================================================
    ! get_particle_density
    !   Computes the bulk density of each detritus class from its elemental
    !   composition, following Cram et al. (2018):
    !
    !     rho_particle = rho_CaCO3 * f_cal + rho_opal * f_si
    !                  + rho_POC   * f_c   + rho_PON  * f_n
    !
    !   where fractions f_x = component_x / sum(all components).
    !
    !   Computed separately for the primary detritus class (rho_particle1,
    !   from det{N,C,Si,Cal}) and, if enabled, the second/macrozooplankton
    !   detritus class (rho_particle2, from detz2{N,C,Si,Cal}).
    !===============================================================================
    subroutine get_particle_density(num_tracers, myDim_nod2d, eDim_nod2D, nl, ulevels_nod2D, &
            nlevels_nod2D, tracers_info)

        use recom_config, only: enable_3zoo2det, rho_CaCO3, rho_opal, rho_POC, rho_PON, tiny
        use recom_glovar, only: tracers_info_type, rho_particle1, rho_particle2
        use recom_declarations, only: wp, tracer_ids

        implicit none

        integer, intent(in) :: myDim_nod2d, eDim_nod2D, nl, num_tracers
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D
        type(tracers_info_type), intent(in) :: tracers_info

        integer :: row, nzmin, nzmax, tr_num

        real(kind=wp) :: a1(nl - 1, myDim_nod2D + eDim_nod2D) ! [n.d.] carbon fraction
        real(kind=wp) :: a2(nl - 1, myDim_nod2D + eDim_nod2D) ! [n.d.] nitrogen fraction
        real(kind=wp) :: a3(nl - 1, myDim_nod2D + eDim_nod2D) ! [n.d.] opal (Si) fraction
        real(kind=wp) :: a4(nl - 1, myDim_nod2D + eDim_nod2D) ! [n.d.] CaCO3 fraction

        real(kind=wp) :: b1(nl - 1, myDim_nod2D + eDim_nod2D) ! detritus carbon [mmol m-3]
        real(kind=wp) :: b2(nl - 1, myDim_nod2D + eDim_nod2D) ! detritus nitrogen [mmol m-3]
        real(kind=wp) :: b3(nl - 1, myDim_nod2D + eDim_nod2D) ! detritus Si [mmol m-3]
        real(kind=wp) :: b4(nl - 1, myDim_nod2D + eDim_nod2D) ! detritus CaCO3 [mmol m-3]
        real(kind=wp) :: aux(nl - 1, myDim_nod2D + eDim_nod2D) ! sum of components (normalizer)

        rho_particle1 = 0.0
        b1 = 0.0
        b2 = 0.0
        b3 = 0.0
        b4 = 0.0
        aux = 0.0

        ! ----------------------------------------------------------------
        ! Gather primary detritus class components. `max(tiny, ...)`
        ! guarantees strictly positive, non-zero tracer fields so the
        ! fraction division below never hits 0/0.
        ! ----------------------------------------------------------------
        do tr_num = 1, num_tracers
            if (tracers_info%ids(tr_num) == tracer_ids%detrital_carbon) then      !idetc
                b1 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_nitrogen) then !idetn
                b2 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_silica) then   !idetsi
                b3 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_calcite) then  !idetcal
                b4 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            end if
        end do

        do row = 1, myDim_nod2d
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row) - 1
            aux(nzmin:nzmax, row) = b1(nzmin:nzmax, row) + b2(nzmin:nzmax, row) &
                    + b3(nzmin:nzmax, row) + b4(nzmin:nzmax, row)
            a1(nzmin:nzmax, row) = b1(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
            a2(nzmin:nzmax, row) = b2(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
            a3(nzmin:nzmax, row) = b3(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
            a4(nzmin:nzmax, row) = b4(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
            rho_particle1(nzmin:nzmax, row) = rho_CaCO3 * a4(nzmin:nzmax, row) &
                    + rho_opal * a3(nzmin:nzmax, row) &
                    + rho_POC * a1(nzmin:nzmax, row) + rho_PON * a2(nzmin:nzmax, row)
        end do

        ! ----------------------------------------------------------------
        ! Repeat for second (macrozooplankton) detritus class, if enabled.
        ! ----------------------------------------------------------------
        if (enable_3zoo2det) then
            rho_particle2 = 0.0
            b1 = 0.0
            b2 = 0.0
            b3 = 0.0
            b4 = 0.0
            aux = 0.0
            do tr_num = 1, num_tracers
                if (tracers_info%ids(tr_num) == tracer_ids%macrozooplankton_detrital_carbon) then
                    !idetz2c
                    b1 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
                else if (tracers_info%ids(tr_num) == &
                         tracer_ids%macrozooplankton_detrital_nitrogen) then
                    !idetz2n
                    b2 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
                else if (tracers_info%ids(tr_num) == &
                         tracer_ids%macrozooplankton_detrital_silica) then
                    !idetz2si
                    b3 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
                else if (tracers_info%ids(tr_num) == &
                         tracer_ids%macrozooplankton_detrital_calcite) then
                    !idetz2calc
                    b4 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
                end if
            end do

            do row = 1, myDim_nod2d ! + eDim_nod2D ! myDim is sufficient
                nzmin = ulevels_nod2D(row)
                nzmax = nlevels_nod2D(row) - 1
                aux(nzmin:nzmax, row) = b1(nzmin:nzmax, row) + b2(nzmin:nzmax, row) &
                        + b3(nzmin:nzmax, row) + b4(nzmin:nzmax, row)
                a1(nzmin:nzmax, row) = b1(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
                a2(nzmin:nzmax, row) = b2(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
                a3(nzmin:nzmax, row) = b3(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
                a4(nzmin:nzmax, row) = b4(nzmin:nzmax, row) / aux(nzmin:nzmax, row)
                rho_particle2(nzmin:nzmax, row) = rho_CaCO3 * a4(nzmin:nzmax, row) &
                        + rho_opal * a3(nzmin:nzmax, row) &
                        + rho_POC * a1(nzmin:nzmax, row) + rho_PON * a2(nzmin:nzmax, row)
            end do
        end if

    end subroutine get_particle_density

    !===============================================================================
    ! get_seawater_viscosity
    !   Approximates dynamic seawater viscosity [kg/(m*s)] as a function of
    !   temperature and salinity, following Sharqawy et al. (2010):
    !   https://bitbucket.org/ohnoplus/ballasted-sinking/src/master/tools/waterviscosity.m
    !
    !   Validity: T in [0, 180] deg C, S in [0, 0.15] kg/kg.
    !   Salinity effects are secondary to temperature but included via the
    !   A/B correction terms; salinity is converted from practical units
    !   (assumed g/kg / psu) to kg/kg by multiplying by 0.001.
    !
    !   Cavity nodes are intentionally included here (unlike
    !   ver_sinking_recom_benthos) -- viscosity is still needed to scale
    !   sinking speeds within sub-shelf cavity waters, even though those
    !   nodes don't deposit into the benthos.
    !===============================================================================
    subroutine get_seawater_viscosity(tr_num, myDim_nod2d, ulevels_nod2D, nlevels_nod2D, &
            tracer_data_values_1, tracer_data_values_2)

        use recom_glovar, only: seawater_visc_3D
        use recom_declarations, only: wp

        implicit none

        !!  tracer_data_values_1 [degrees C] Ocean temperature
        !!  tracer_data_values_2 [g/kg or n.d.] Ocean salinity
        !!  seawater_visc_3D     [kg m-1 s-1] Ocean dynamic viscosity (output)

        integer, intent(in) :: myDim_nod2d
        integer, intent(in), target :: tr_num
        real(kind=wp), intent(in), dimension(:, :) :: tracer_data_values_1, tracer_data_values_2
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D

        integer :: row, k, nzmin, nzmax
        real(kind=wp), dimension(1) :: A, B, mu_w

        seawater_visc_3D(:, :) = 0.0
        do row = 1, myDim_nod2d
            ! Cavity nodes (ulevels_nod2D(row) > 1) are deliberately NOT skipped
            ! here -- see subroutine docstring above.
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row) - 1

            do k = nzmin, nzmax
                ! Sharqawy (2010) fit. Validity: 0 < T < 180 degC, 0 < S < 0.15 kg/kg.
                ! Salinity must be in kg/kg here, hence the 0.001 conversion factor
                ! applied to tracer_data_values_2 (assumed practical salinity / g kg-1).
                A(1) = 1.541 + 1.998 * 0.01 * tracer_data_values_1(k, row) &
                        - 9.52 * 1e-5 * tracer_data_values_1(k, row) * tracer_data_values_1(k, row)
                B(1) = 7.974 - 7.561 * 0.01 * tracer_data_values_1(k, row) &
                        + 4.724 * 1e-4 * tracer_data_values_1(k, row) * tracer_data_values_1(k, row)
                mu_w(1) = 4.2844 * 1.0e-5 &
                        + (1.0 &
                        / (0.157 * (tracer_data_values_1(k, row) + 64.993) &
                        * (tracer_data_values_1(k, row) + 64.993) &
                        - 91.296))
                seawater_visc_3D(k, row) = mu_w(1) &
                        * (1.0 + A(1) * tracer_data_values_2(k, row) * 0.001 &
                        + B(1) * tracer_data_values_2(k, row) * 0.001 &
                        * tracer_data_values_2(k, row) * 0.001)
            end do
        end do

    end subroutine get_seawater_viscosity

end module recom_sinking
