module recom_sinking
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
    ! YY: sinking of second detritus adapted from Ozgur's code
    ! but not using recom_det_tracer_id, since
    ! second detritus has a different sinking speed than the first
    ! define recom_det2_tracer_id to make it consistent???
    !===============================================================================
    subroutine ver_sinking_recom_benthos(tr_num, nl, ulevels_nod2D, nlevels_nod2D, zbar_3d_n, &
            nod_in_elem2D_num, nod_in_elem2D, nlevels, area, tracer_id, tracer_data_values, &
            myDim_nod2d, str_bf, mype, MPI_COMM_FESOM, npes, sn, rn, s_mpitype_nod2D, &
            r_mpitype_nod2D, s_mpitype_nod3D, r_mpitype_nod3D, sPE, rPE, requests, nreq, dt)

        use recom_g_comm_auto, only: recom_exchange_nod

        use recom_declarations, only: wp, tracer_ids
        use recom_glovar, only: Benthos, Benthos_tr, SinkFlx_tr

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
        real(kind=WP), intent(in), dimension(:, :) :: zbar_3d_n, area, tracer_data_values
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
            nl1 = nlevels_nod2D(n) - 1
            ul1 = ulevels_nod2D(n)

            aux = 0._WP
            Vben = 0._WP
            add_benthos_2d = 0._WP

            ! Calculate sinking velociy for vertical sinking case
            ! ******************************************************
            if (any(recom_det_tracer_id == tracer_id)) Vben = Vdet
            if (any(recom_phy_tracer_id == tracer_id)) Vben = VPhy
            if (any(recom_dia_tracer_id == tracer_id)) Vben = VDia
            if (allow_var_sinking) then
                Vben = Vdet_a * abs(zbar_3d_n(:, n)) + Vben
            end if

            ! Constant vertical sinking for the second detritus class
            ! *******************************************************

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

            !! * Screening minimum depth in neigbouring nodes around node n*
            nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

            do nz = nlevels_nod2D_minimum, nl1
                tv = tracer_data_values(nz, n) * Vben(nz)
                aux(nz) = -tv * (area(nz, n) - area(nz + 1, n))
            end do

            do nz = ul1, nl1
                str_bf(nz, n) = str_bf(nz, n) &
                        + (aux(nz)) * dt / area(nz, n) / (zbar_3d_n(nz, n) - zbar_3d_n(nz + 1, n))
                !!!!!!!!CHECK Maybe /area(nz,n) -> [mmol/m2]
                add_benthos_2d(n) = add_benthos_2d(n) - (aux(nz)) * dt
            end do

            !! * Particulate Organic Nitrogen *
            if (tracer_id == tracer_ids%phytoplankton_nitrogen .or. & !iphyn
                tracer_id == tracer_ids%detrital_nitrogen .or. & !idetn
                tracer_id == tracer_ids%diatom_nitrogen .or. & !idian
                tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen) then !idetz2n
                Benthos(n, 1) = Benthos(n, 1) + add_benthos_2d(n) ![mmol]

                if (use_MEDUSA) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    SinkFlx_tr(n, 1, tr_num) = SinkFlx_tr(n, 1, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt ![mmol/m2]
                    ! now SinkFlx hat the unit mmol/time step
                    ! but mmol/m2/time is needed for MEDUSA: thus /area
                end if
                if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    Benthos_tr(n, 1, tr_num) = Benthos_tr(n, 1, tr_num) + add_benthos_2d(n) ![mmol]
                end if

            end if

            !! * Particulate Organic Carbon *
            if (tracer_id == tracer_ids%phytoplankton_carbon .or. & !iphyc
                tracer_id == tracer_ids%detrital_carbon .or. & !idetc
                tracer_id == tracer_ids%diatom_carbon .or. & !idiac
                tracer_id == tracer_ids%macrozooplankton_detrital_carbon) then !idetz2c
                Benthos(n, 2) = Benthos(n, 2) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    SinkFlx_tr(n, 2, tr_num) = SinkFlx_tr(n, 2, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
                if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    Benthos_tr(n, 2, tr_num) = Benthos_tr(n, 2, tr_num) + add_benthos_2d(n)
                end if

            end if

            !! *Particulate Organic Silicon *
            if (tracer_id == tracer_ids%diatom_silica .or. & !idiasi
                    tracer_id == tracer_ids%detrital_silica .or. & !idetsi
                    tracer_id == tracer_ids%macrozooplankton_detrital_silica) then !idetz2si
                Benthos(n, 3) = Benthos(n, 3) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    SinkFlx_tr(n, 3, tr_num) = SinkFlx_tr(n, 3, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
                if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    Benthos_tr(n, 3, tr_num) = Benthos_tr(n, 3, tr_num) + add_benthos_2d(n)
                end if

            end if

            !! * Cal *
            if (tracer_id == tracer_ids%phytoplankton_calcite .or. & !iphycal
                    tracer_id == tracer_ids%detrital_calcite .or. & !idetcal
                    tracer_id == tracer_ids%macrozooplankton_detrital_calcite) then !idetz2cal
                Benthos(n, 4) = Benthos(n, 4) + add_benthos_2d(n)

                if (use_MEDUSA) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    SinkFlx_tr(n, 4, tr_num) = SinkFlx_tr(n, 4, tr_num) &
                            + add_benthos_2d(n) / area(1, n) / dt
                end if
                if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                    ! regarding global sums when running the tracer loop in parallel
                    Benthos_tr(n, 4, tr_num) = Benthos_tr(n, 4, tr_num) + add_benthos_2d(n)
                end if

            end if

            ! flux of 13C into the sediment
            if (ciso) then
                if (tracer_id == 1305 .or. & !iphyc_13
                        tracer_id == 1308 .or. & !idetc_13
                        tracer_id == 1314) then !idiac_14

                    if (use_MEDUSA) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        SinkFlx_tr(n, 5, tr_num) = SinkFlx_tr(n, 5, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
                    if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        Benthos_tr(n, 5, tr_num) = Benthos_tr(n, 5, tr_num) + add_benthos_2d(n)
                    end if

                end if

                if (tracer_id == 1320 .or. & !iphycal
                        tracer_id == 1321) then !idetcal

                    if (use_MEDUSA) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        SinkFlx_tr(n, 6, tr_num) = SinkFlx_tr(n, 6, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
                    if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        Benthos_tr(n, 6, tr_num) = Benthos_tr(n, 6, tr_num) + add_benthos_2d(n)
                    end if

                end if

            end if

            ! flux of 14C into the sediment
            if (ciso .and. ciso_organic_14) then
                if (tracer_id == 1405 .or. & !iphyc_13
                        tracer_id == 1408 .or. & !idetc_13
                        tracer_id == 1414) then !idiac_14

                    if (use_MEDUSA) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        SinkFlx_tr(n, 7, tr_num) = SinkFlx_tr(n, 7, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
                    if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        Benthos_tr(n, 7, tr_num) = Benthos_tr(n, 7, tr_num) + add_benthos_2d(n)
                    end if

                end if

                if (tracer_id == 1420 .or. & !iphycal
                        tracer_id == 1421) then !idetcal
                    if (use_MEDUSA) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        SinkFlx_tr(n, 8, tr_num) = SinkFlx_tr(n, 8, tr_num) &
                                + add_benthos_2d(n) / area(1, n) / dt
                    end if
                    if ((.not.use_MEDUSA) .or. (sedflx_num == 0)) then
                        ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical
                        ! results regarding global sums when running the tracer loop in parallel
                        Benthos_tr(n, 8, tr_num) = Benthos_tr(n, 8, tr_num) + add_benthos_2d(n)
                    end if
                end if

            end if

        end do

        if (use_MEDUSA) then
            do n = 1, bottflx_num
                !           SinkFlx(:,n) = Sinkflx(:,n)/dt
                ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results
                ! regarding global sums when running the tracer loop in parallel
                call recom_exchange_nod(SinkFlx_tr(:, n, tr_num), npes, sn, rn, MPI_COMM_FESOM, &
                        mype, s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
            end do
        end if ! use_MEDUSA

        do n = 1, benthos_num
            ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding
            ! global sums when running the tracer loop in parallel
            call recom_exchange_nod(Benthos_tr(:, n, tr_num), npes, sn, rn, MPI_COMM_FESOM, &
                    mype, s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)

            call recom_exchange_nod(Benthos(:, n), npes, sn, rn, MPI_COMM_FESOM, mype, &
                    s_mpitype_nod2D, r_mpitype_nod2D, sPE, rPE, requests, nreq)
        end do

    end subroutine ver_sinking_recom_benthos

    !
    !
    !===============================================================================
    subroutine diff_ver_recom_expl(nl, ulevels_nod2D, nlevels_nod2D, nod_in_elem2D_num, &
            nod_in_elem2D, nlevels, area, areasvol, hnode_new, tracer_id, myDim_nod2d, &
            eDim_nod2D, mype, MPI_COMM_FESOM, dtr_bf, dt)

        ! Remineralization from benthos
        ! bottom_flux

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
        if (use_MEDUSA .and. (sedflx_num /= 0)) then
            !CV update: the calculation later has been changed by Ozgur in such
            !a way  that now the  variable bottom_flux is in  (mol/time) units,
            !rather than  a flux in  (mol/time/area). I therefore  multiply the
            !Medusa fluxes by the area to get the same unit.

            select case (tracer_id)
            case (1001)
                bottom_flux = GloSed(:, 1) * area(1, :) ! DIN
            case (1002)
                bottom_flux = GloSed(:, 2) * area(1, :) ! DIC
            case (1003)
                bottom_flux = GloSed(:, 3) * area(1, :) ! Alk
            case (1018)
                bottom_flux = GloSed(:, 4) * area(1, :) ! Si
            case (1019)
                bottom_flux = GloSed(:, 1) * Fe2N_benthos * area(1, :)
            case (1022)
                bottom_flux = GloSed(:, 5) * area(1, :) ! Oxy
            case (1302)
                if (ciso) then
                    bottom_flux = GloSed(:, 6) * area(1, :) ! DIC_13 and Calc: DIC_13
                end if
            case (1402)
                if (ciso) then
                    bottom_flux = GloSed(:, 7) * area(1, :) ! DIC_14 and Calc: DIC_14
                end if
            case default
                if (mype == 0) then
                    write(*, *) 'check specified in boundary conditions'
                    write(*, *) 'the model will stop!'
                end if
                ! This can be improved later on
                call MPI_ABORT(MPI_COMM_FESOM, 1)
                stop
            end select
        else
            select case (tracer_id)
            case (1001)
                bottom_flux = GlodecayBenthos(:, 1) !*** DIN [mmolN/m^2/s] ***
            case (1002)
                !*** DIC + calcification ***
                bottom_flux = GlodecayBenthos(:, 2) + GlodecayBenthos(:, 4)
            case (1003)
                !*** Alk ***
                bottom_flux = GlodecayBenthos(:, 4) * 2.0_WP - 1.0625_WP * GlodecayBenthos(:, 1)
            case (1018)
                bottom_flux = GlodecayBenthos(:, 3) !*** Si ***
            case (1019)
                bottom_flux = GlodecayBenthos(:, 1) * Fe2N_benthos !*** DFe ***
            case (1022)
                bottom_flux = -GlodecayBenthos(:, 2) * redO2C !*** O2 ***
            case (1302)
                if (ciso) then
                    !*** DIC_13 and Calc: DIC_13 ***
                    bottom_flux = GlodecayBenthos(:, 5) + GlodecayBenthos(:, 6)
                end if
            case (1402)
                if (ciso) then
                    !*** DIC_14 and Calc: DIC_14 ***
                    bottom_flux = GlodecayBenthos(:, 7) + GlodecayBenthos(:, 8)
                end if
            case default
                if (mype == 0) then
                    write(*, *) 'check specified in boundary conditions'
                    write(*, *) 'the model will stop!'
                end if
                ! This can be improved later on
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
            ! Screening minimum depth in neigbouring nodes around node n
            nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

            !_______________________________________________________________________
            ! Bottom flux
            do nz = nlevels_nod2D_minimum, nl1
                vd_flux(nz) = (area(nz, n) - area(nz + 1, n)) * bottom_flux(n) / (area(1, n))
            end do
            nz = nl1
            vd_flux(nz + 1) = (area(nz + 1, n)) * bottom_flux(n) / (area(1, n))
            !_______________________________________________________________________
            ! writing flux into rhs
            do nz = ul1, nl1
                ! flux contribute only the cell through its bottom !!!
                !            dtr_bf(nz,n) = dtr_bf(nz,n) +
                ! vd_flux(nz+1)*dt/area(nz,n)/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
                dtr_bf(nz, n) = dtr_bf(nz, n) &
                        + vd_flux(nz + 1) * dt / areasvol(nz, n) / hnode_new(nz, n)
            end do
        end do
    end subroutine diff_ver_recom_expl

    subroutine ver_sinking_recom(tr_num, nl, ulevels_nod2D, nlevels_nod2D, zbar_3d_n, z_3d_n, &
            nod_in_elem2D_num, nod_in_elem2D, nlevels, area, areasvol, hnode, hnode_new, &
            tracer_id, tracer_data_values, myDim_nod2d, vert_sink, dt)
        ! Sinking in water column

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

        ! calculate scaling factors
        ! scaling_density1_3D, scaling_density2_3D
        ! scaling_visc_3D
        ! .OG. 04.11.2022

        ! Constant sinking velocities (we prescribe them under namelist recom)
        ! This hardcoded part is temporary
        ! .OG. 07.07.2021

        Vsink = 0.0_WP

        ! Assign sinking velocities based on tracer ID
        ! Groups tracers by functional type and assigns corresponding velocity

        ! Detritus tracers (nitrogen, carbon, silicate, calcite)
        if (tracer_id == tracer_ids%detrital_nitrogen .or. & ! idetn
                tracer_id == tracer_ids%detrital_carbon .or. & ! idetc
                tracer_id == tracer_ids%detrital_silica .or. & ! idetsi
                tracer_id == tracer_ids%detrital_calcite) then ! idetcal
            Vsink = VDet

            ! Phytoplankton tracers (nitrogen, carbon, chlorophyll)
        elseif (tracer_id == tracer_ids%phytoplankton_nitrogen .or. & ! iphyn
                    tracer_id == tracer_ids%phytoplankton_carbon .or. & ! iphyc
                    tracer_id == tracer_ids%phytoplankton_chlorophyll) then ! ipchl
            Vsink = VPhy

            ! Diatom tracers (nitrogen, carbon, silicate, chlorophyll)
        elseif (tracer_id == tracer_ids%diatom_nitrogen .or. & ! idian
                    tracer_id == tracer_ids%diatom_carbon .or. & ! idiac
                    tracer_id == tracer_ids%diatom_silica .or. & ! idiasi
                    tracer_id == tracer_ids%diatom_chlorophyll) then ! idchl
            Vsink = VDia

            ! Coccolithophore tracers (nitrogen, carbon, chlorophyll)
        elseif (enable_coccos .and. &
                    (tracer_id == tracer_ids%coccolithophore_nitrogen .or. & ! icocn
                    tracer_id == tracer_ids%coccolithophore_carbon .or. & ! icocc
                    tracer_id == tracer_ids%coccolithophore_chlorophyll)) then ! icchl
            Vsink = VCocco

            ! Phaeocystis tracers (nitrogen, carbon, chlorophyll)
        elseif (enable_coccos .and. &
                    (tracer_id == tracer_ids%phaeocystis_nitrogen .or. & ! iphan
                    tracer_id == tracer_ids%phaeocystis_carbon .or. & ! iphac
                    tracer_id == tracer_ids%phaeocystis_chlorophyll)) then ! iphachl
            Vsink = VPhaeo

            ! Phytoplankton calcite tracer (special case)
        elseif (tracer_id == tracer_ids%phytoplankton_calcite) then ! iphycal
            if (enable_coccos) then
                Vsink = VCocco
            else
                Vsink = VPhy
            end if

            ! Zooplankton-2 detritus tracers (nitrogen, carbon, silicate, calcite)
        elseif (enable_3zoo2det .and. &
                (tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. & !idetz2n
                 tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. & !idetz2c
                 tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. & !idetz2si
                 tracer_id == tracer_ids%macrozooplankton_detrital_calcite)) then !idetz2calc
            Vsink = VDet_zoo2

        end if

        !! ---- No sinking if Vsink < 0.1 m/day
        if (Vsink > 0.1) then

            do n = 1, myDim_nod2D
                if (ulevels_nod2D(n) > 1) cycle
                nzmin = ulevels_nod2D(n)
                nzmax = nlevels_nod2D(n) - 1

                !! distance between tracer points, surface and bottom dz_trr is half
                !! the layer thickness
                dz_trr = 0.0d0
                dz_trr(nzmin + 1:nzmax) = abs(Z_3d_n(nzmin:nzmax - 1, n) &
                        - Z_3d_n(nzmin + 1:nzmax, n))
                dz_trr(nzmin) = hnode(nzmin, n) / 2.0d0
                dz_trr(nzmax + 1) = hnode(nzmax, n) / 2.0d0

                Wvel_flux(nzmin:nzmax + 1) = 0.d0 ! Vertical velocity for BCG tracers

                do nz = nzmin, nzmax + 1

                    Wvel_flux(nz) = -Vsink / SecondsPerDay ! allow_var_sinking = .false.

                    if (allow_var_sinking) then
                        Wvel_flux(nz) = -((Vdet_a * abs(zbar_3d_n(nz, n)) / SecondsPerDay) &
                                + Vsink / SecondsPerDay)
                        if (use_ballasting) then
                            ! Apply ballasting on slow sinking detritus
                            !if (any(recom_sinking_tracer_id == tracer_id(tr_num))) then

                            if (tracer_id == tracer_ids%detrital_nitrogen .or. & !idetn
                                    tracer_id == tracer_ids%detrital_carbon .or. & !idetc
                                    tracer_id == tracer_ids%detrital_silica .or. & !idetsi
                                    tracer_id == tracer_ids%detrital_calcite) then !idetcal
                                Wvel_flux(nz) = w_ref1 * scaling_density1_3D(nz, n) &
                                        * scaling_visc_3D(nz, n)

                                if (depth_scaling1 > 0.0) Wvel_flux(nz) = Wvel_flux(nz) &
                                        + (depth_scaling1 * abs(zbar_3d_n(nz, n)))

                                if (abs(Wvel_flux(nz)) > max_sinking_velocity) Wvel_flux(nz) = &
                                        max_sinking_velocity

                                ! sinking velocity [m d-1] surface --> bottom (negative)
                                ! now in [m s-1]
                                Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz) / SecondsPerDay
                            end if
                        end if
                    end if

                    !! ---- We assume constant sinking for second detritus
                    ! Parenthesized explicitly: .and. binds tighter than .or. in Fortran, so
                    ! without these parens this read as (enable_3zoo2det .and. ...nitrogen) .or.
                    ! ...carbon .or. ...silica .or. ...calcite instead of the intended
                    ! enable_3zoo2det .and. (...nitrogen .or. ...carbon .or. ...silica .or.
                    ! ...calcite). Harmless today only because disabled tracer_ids default to a
                    ! -1 sentinel that a real tracer_id never matches.
                    if (enable_3zoo2det .and. &
                        (tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. & !idetz2n
                         tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. & !idetz2c
                         tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. & !idetz2si
                         tracer_id == tracer_ids%macrozooplankton_detrital_calcite)) &  !idetz2calc
                        then
                        Wvel_flux(nz) = -VDet_zoo2 / SecondsPerDay

                        if (use_ballasting) then

                            Wvel_flux(nz) = w_ref2 * scaling_density2_3D(nz, n) &
                                    * scaling_visc_3D(nz, n)

                            if (depth_scaling2 > 0.0) Wvel_flux(nz) = Wvel_flux(nz) &
                                    + (depth_scaling2 * abs(zbar_3d_n(nz, n)))

                            if (abs(Wvel_flux(nz)) > max_sinking_velocity) Wvel_flux(nz) = &
                                    max_sinking_velocity

                            ! sinking velocity [m d-1] surface --> bottom (negative)
                            Wvel_flux(nz) = -1.0d0 * Wvel_flux(nz) / SecondsPerDay ! now in [m s-1]
                        end if

                    end if

                    !-1.0d0/SecondsPerDay  !idetcal
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

                !FIXME: Having IF True and IF False is bad practice. Either throw away the old code,
                !or make a namelist switch...
                ! 3rd Order DST Sceheme with flux limiting. This code comes from old recom
                if (.true.) then

                    k = nod_in_elem2D_num(n)
                    ! Screening minimum depth in neigbouring nodes around node n
                    nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)

                    vd_flux(nzmin:nzmax + 1) = 0.0_WP

                    do nz = nzmax, nzmin + 1, -1

                        Rjp = tracer_data_values(nz, n) - tracer_data_values(min(nz + 1, nzmax), n)
                        Rj = tracer_data_values(max(nzmin, nz - 1), n) - tracer_data_values(nz, n)
                        Rjm = tracer_data_values(max(nzmin, nz - 2), n) &
                                - tracer_data_values(max(nzmin, nz - 1), n)

                        !(Z_n(nz-1)-Z_n(nz)))
                        ! [m/day] * [day] * [1/m]  ! NEW BALL changed dt to dt_sink
                        cfl = abs(Wvel_flux(nz) * dt_sink / dz_trr(nz))

                        wPs = Wvel_flux(nz) + abs(Wvel_flux(nz)) ! --> Positive vertical velocity
                        wM = Wvel_flux(nz) - abs(Wvel_flux(nz)) ! --> Negative vertical velocity

                        d0 = (2.d0 - cfl) * (1.d0 - cfl) * onesixth
                        d1 = (1.d0 - cfl * cfl) * onesixth

                        thetaP = Rjm / (1.d-20 + Rj)
                        psiP = d0 + d1 * thetaP
                        psiP = max(0.d0, min(min(1.d0, psiP), &
                                (1.d0 - cfl) / (1.d-20 + cfl) * thetaP))

                        thetaM = Rjp / (1.d-20 + Rj)
                        psiM = d0 + d1 * thetaM
                        psiM = max(0.d0, min(min(1.d0, psiM), &
                                (1.d0 - cfl) / (1.d-20 - cfl) * thetaM))

                        tv = (0.5 * wPs * (tracer_data_values(nz, n) + psiM * Rj) + &
                                0.5 * wM * (tracer_data_values(max(nzmin, nz - 1), n) + psiP * Rj))
                        vd_flux(nz) = -tv * area(nz, n)
                    end do
                end if ! 3rd Order DST Sceheme with flux limiting

                if (.false.) then ! simple upwind

                    ! Surface flux
                    vd_flux(nzmin) = 0.0_WP

                    ! Bottom flux
                    vd_flux(nzmax + 1) = 0.0_WP

                    k = nod_in_elem2D_num(n)
                    ! Screening minimum depth in neigbouring nodes around node n
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
                do nz = nzmin, nzmax
                    vert_sink(nz, n) = vert_sink(nz, n) &
                            + (vd_flux(nz) - vd_flux(nz + 1)) * dt / areasvol(nz, n) &
                            / hnode_new(nz, n) !/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
                end do
            end do
        end if ! Vsink .gt. 0.1

    end subroutine ver_sinking_recom

    !-------------------------------------------------------------------------------
    ! Subroutine calculate ballasting
    !-------------------------------------------------------------------------------
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

        integer :: row, k, nzmin, nzmax
        real(kind=wp) :: depth_pos(1)
        real(kind=wp) :: pres(1)
        real(kind=wp) :: sa(1)
        real(kind=wp) :: ct(1)
        real(kind=wp) :: rho_seawater(1)
        real(kind=wp) :: Lon_degree(1)
        real(kind=wp) :: Lat_degree(1)

        ! For ballasting, calculate scaling factors here and pass them to FESOM, where sinking
        ! velocities are calculated
        ! -----
        ! If ballasting is used, sinking velocities are a function of a) particle composition
        ! (=density),
        ! b) sea water viscosity, c) depth (currently for small detritus only), and d) a constant
        ! reference sinking speed
        ! -----

        ! check oce_ale_tracer.F90
        !     call get_seawater_viscosity(mesh) ! seawater_visc_3D
        !     call get_particle_density(mesh) ! rho_particle = density of particle class 1 and 2
        !___________________________________________________________________________
        ! loop over local nodes
        do row = 1, myDim_nod2D
            ! max. number of levels at node n
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row)
            !! lon
            Lon_degree(1) = geo_coord_nod2D(1, row) / rad !! convert from rad to degree
            !! lat
            Lat_degree(1) = geo_coord_nod2D(2, row) / rad !! convert from rad to degree

            ! get scaling vectors -> these need to be passed to FESOM to get sinking velocities
            ! get local seawater density
            do k = nzmin, nzmax

                !! level depth
                ! take depth of tracers instead of levels abs(zbar_3d_n(k,row))
                depth_pos(1) = abs(Z_3d_n(k, row))

                ! pres is output of function,1=number of records
                call depth2press(depth_pos(1), Lat_degree(1), pres, 1)
                sa = gsw_sa_from_sp(tracer_data_values_2(k, row), pres, Lon_degree(1), Lat_degree(1&
                        ))
                ct = gsw_ct_from_pt(sa, tracer_data_values_1(k, row))
                rho_seawater = gsw_rho(sa, ct, pres)

                ! (i.e. no density scaling)
                scaling_density1_3D(k, row) = 1.0
                scaling_density2_3D(k, row) = 1.0

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
                    !if (tracers%data(tr_num)%ID ==1008)then !idetc
                    !if (tracers%data(tr_num)%values(k,row)>0.001) then ! only apply ballasting
                    !above a certain biomass (OG Todo: remove)
                    scaling_density1_3D(k, row) = (rho_particle1(k, row) - rho_seawater(1)) &
                            / (rho_ref_part - rho_ref_water)
                    !endif
                    !endif
                    if (enable_3zoo2det) then
                    !if (enable_3zoo2det .and. &
                    !tracers%data(tr_num)%ID ==1026)then ! idetz2c
                    !if (tracers%data(tr_num)%values(k,row)>0.001) then ! only apply ballasting
                    !above a certain biomass (OG Todo: remove)
                    scaling_density2_3D(k, row) = (rho_particle2(k, row) - rho_seawater(1)) &
                            / (rho_ref_part - rho_ref_water)
                    !endif
                    endif
                end if

                scaling_visc_3D(k, row) = 1.0

                if (use_viscosity_scaling) then
                    if (seawater_visc_3D(k, row) < tiny) then
                        scaling_visc_3D(k, row) = 1.0
                    else
                        scaling_visc_3D(k, row) = visc_ref_water / seawater_visc_3D(k, row)
                    end if
                end if

            end do
            if (use_density_scaling) then
                rho_particle1(nzmax + 1, row) = rho_particle1(nzmax, row)
                if (enable_3zoo2det) then
                    rho_particle2(nzmax + 1, row) = rho_particle2(nzmax, row)
                end if
            end if
            scaling_visc_3D(nzmax + 1, row) = scaling_visc_3D(nzmax, row)
        end do
        ! in the unlikely (if possible at all...) case that rho_particle(k)-rho_seawater(1)<0,
        ! prevent the scaling factor from being negative

        ! tiny = 2.23D-16
        if (use_density_scaling) then
            if (any(scaling_density1_3D(:, :) <= tiny)) scaling_density1_3D(:, :) = 1.0_WP
            if (enable_3zoo2det) then
                ! tiny = 2.23D-16
                if (any(scaling_density2_3D(:, :) <= tiny)) scaling_density2_3D(:, :) = 1.0_WP
            end if
        end if

    end subroutine ballast

    !-------------------------------------------------------------------------------
    ! Subroutine calculate density of particle
    ! depending on composition (detC, detOpal, detCaCO3) based on Cram et al. (2018)
    !-------------------------------------------------------------------------------
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

        ! [n.d.] fraction of carbon in detritus class
        real(kind=wp) :: a1(nl - 1, myDim_nod2D + eDim_nod2D)

        ! [n.d.] fraction of nitrogen in detritus class
        real(kind=wp) :: a2(nl - 1, myDim_nod2D + eDim_nod2D)

        ! [n.d.] fraction of Opal in detritus class
        real(kind=wp) :: a3(nl - 1, myDim_nod2D + eDim_nod2D)

        ! [n.d.] fraction of CaCO3 in detritus class
        real(kind=wp) :: a4(nl - 1, myDim_nod2D + eDim_nod2D)

        real(kind=wp) :: b1(nl - 1, myDim_nod2D + eDim_nod2D)
        real(kind=wp) :: b2(nl - 1, myDim_nod2D + eDim_nod2D)
        real(kind=wp) :: b3(nl - 1, myDim_nod2D + eDim_nod2D)
        real(kind=wp) :: b4(nl - 1, myDim_nod2D + eDim_nod2D)
        real(kind=wp) :: aux(nl - 1, myDim_nod2D + eDim_nod2D)

        rho_particle1 = 0.0
        b1 = 0.0
        b2 = 0.0
        b3 = 0.0
        b4 = 0.0
        aux = 0.0

        ! Below guarantees non-negative tracer field
        do tr_num = 1, num_tracers
            if (tracers_info%ids(tr_num) == tracer_ids%detrital_carbon) then
                !idetc      ! [mmol m-3] detritus carbon
                b1 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_nitrogen) then
                !idetn      ! [mmol m-3] detritus nitrogen
                b2 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_silica) then
                !idetsi     ! [mmol m-3] detritus Si
                b3 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            else if (tracers_info%ids(tr_num) == tracer_ids%detrital_calcite) then
                !idetcal    ! [mmol m-3] detritus CaCO3
                b4 = max(tiny, tracers_info%data_pointers(tr_num)%tracer_data(:, :))
            end if
        end do

        do row = 1, myDim_nod2d
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row)
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
                !if (ulevels_nod2D(row)>1) cycle
                nzmin = ulevels_nod2D(row)
                nzmax = nlevels_nod2D(row)
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

    !-------------------------------------------------------------------------------
    ! Subroutine to approximate seawater viscosity with current temperature
    ! based on Cram et al. (2018)
    !-------------------------------------------------------------------------------

    ! neglecting salinity effects, which are much smaller than those of temperature
    ! https://bitbucket.org/ohnoplus/ballasted-sinking/src/master/tools/waterviscosity.m

    subroutine get_seawater_viscosity(tr_num, myDim_nod2d, ulevels_nod2D, nlevels_nod2D, &
            tracer_data_values_1, tracer_data_values_2)

        use recom_glovar, only: seawater_visc_3D
        use recom_declarations, only: wp

        implicit none

        !!  temp [degrees C] Ocean temperature
        !!  salt [g/kg or n.d.] Ocean salinity
        !!  seawater_visc_3D [kg m-1 s-1] Ocean viscosity

        integer, intent(in) :: myDim_nod2d
        integer, intent(in), target :: tr_num
        real(kind=wp), intent(in), dimension(:, :) :: tracer_data_values_1, tracer_data_values_2
        integer, intent(in), dimension(:) :: ulevels_nod2D, nlevels_nod2D

        integer :: row, k, nzmin, nzmax
        real(kind=wp), dimension(1) :: A, B, mu_w

        seawater_visc_3D(:, :) = 0.0
        do row = 1, myDim_nod2d
            !if (ulevels_nod2D(row)>1) cycle
            ! OG Do we need any limitation here?
            ! i.e., if (seawater_visc_3D(row)<=0.0_WP) cycle
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row)

            do k = nzmin, nzmax
                ! Eq from Sharaway 2010
                ! validity:
                !  0<temp<180 degC
                !  0<salt<0.15 kg/kg
                ! Note: because salinity is expected to be in kg/kg, use conversion factor 0.001
                ! below!
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
