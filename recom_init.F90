! ------------
! 23.03.2023
! OG
!===============================================================================
! allocate & initialise arrays for REcoM
module recom_init_interface
    implicit none
    private

    public :: recom_init
    public :: initialize_tracer_ids
    public :: get_tracer_init_value

contains
    !
    !===============================================================================
    ! Model Configuration Summary
    !===============================================================================
    ! Configuration 1: Base model (enable_3zoo2det=F, enable_coccos=F)
    !   - 2 Phytoplankton: General Phy, Diatoms
    !   - 1 Zooplankton: Heterotrophs
    !   - 1 Detritus pool
    !
    ! Configuration 2: 3Zoo2Det (enable_3zoo2det=T, enable_coccos=F)
    !   - 2 Phytoplankton: General Phy, Diatoms
    !   - 3 Zooplankton: Het, Zoo2, Zoo3
    !   - 2 Detritus pools
    !
    ! Configuration 3: Coccos (enable_3zoo2det=F, enable_coccos=T)
    !   - 4 Phytoplankton: General Phy, Diatoms, Coccos, Phaeo
    !   - 1 Zooplankton: Heterotrophs
    !   - 1 Detritus pool
    !
    ! Configuration 4: Full model (enable_3zoo2det=T, enable_coccos=T)
    !   - 4 Phytoplankton: General Phy, Diatoms, Coccos, Phaeo
    !   - 3 Zooplankton: Het, Zoo2, Zoo3
    !   - 2 Detritus pools
    !===============================================================================
    !
    ! Tracer layout (fixed order):
    !   T, S  |  BGC (bgc_num tracers)  |  age (optional, ID=100)  |  transit (optional)
    !
    ! BGC tracers always start immediately at slot 3, regardless of whether
    ! transit tracers are active. Transit tracers (SF6, CFC-11, CFC-12, R14C,
    ! R39Ar) are appended at the very end, after the optional age tracer.
    !===============================================================================
    subroutine recom_init(nl, ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, Z_3d_n, myDim_nod2d, &
            eDim_nod2D, mype, MPI_COMM_FESOM, myDim_elem2D, eDim_elem2D, tracers_info, &
            num_tracers, rad, use_age_tracer, use_transit, l_sf6, l_f11, l_f12, l_r14c, l_r39ar)

        use REcoM_declarations, only: wp, tracer_ids
        use REcoM_GloVar, only: tracers_info_type
        use recom_config, only: validate_recom_tracers, initialize_tracer_indices, &
                validate_tracer_id_sequence, bgc_num

        implicit none

        integer, intent(in) :: nl, mydim_nod2d, edim_nod2d, mype, num_tracers
        integer, intent(in) :: mpi_comm_fesom, mydim_elem2d, edim_elem2d
        real(kind=WP), intent(in) :: rad
        logical, intent(in) :: use_age_tracer, use_transit, l_sf6, l_f11, l_f12, l_r14c, l_r39ar
        integer, intent(in), dimension(:) :: ulevels_nod2d, nlevels_nod2d
        real(kind=wp), intent(in), dimension(:, :) :: geo_coord_nod2d, z_3d_n
        type(tracers_info_type), intent(in) :: tracers_info

        integer :: i, tracer_id
        integer :: actual_bgc_num
        integer :: num_physical_tracers
        integer :: n_transit_tracers   ! number of active transit tracers
        integer :: bgc_start, bgc_end  ! first/last slot of the BGC-only block

        call initialize_memory(myDim_nod2D + eDim_nod2D, nl, num_tracers)

        call initialize_tracer_ids

        ! After reading parecomsetup namelist
        call initialize_tracer_indices

        ! Validation check here
        call validate_recom_tracers(num_tracers, use_age_tracer, use_transit, l_sf6, l_f11, l_f12, l_r14c, l_r39ar, mype)

        ! After reading tracer namelist - validate actual IDs
        call validate_tracer_id_sequence(tracers_info%ids(1:num_tracers), num_tracers, use_age_tracer, use_transit, &
                                            l_sf6, l_f11, l_f12, l_r14c, l_r39ar, mype)

        ! T,S | BGC | [age] | [transit] num_physical_tracers
        ! is always just T,S (=2): BGC tracers start right after them, regardless
        ! of whether transit tracers are active. Transit tracers are appended at
        ! the tail instead and are not initialized here.
        num_physical_tracers = 2

        n_transit_tracers = 0
        if (use_transit) then
            if (l_sf6)   n_transit_tracers = n_transit_tracers + 1
            if (l_f11)   n_transit_tracers = n_transit_tracers + 1
            if (l_f12)   n_transit_tracers = n_transit_tracers + 1
            if (l_r14c)  n_transit_tracers = n_transit_tracers + 1
            if (l_r39ar) n_transit_tracers = n_transit_tracers + 1
        end if

        actual_bgc_num = num_tracers - num_physical_tracers
        if (use_age_tracer) actual_bgc_num = actual_bgc_num - 1
        if (use_transit) actual_bgc_num = actual_bgc_num - n_transit_tracers

        bgc_start = num_physical_tracers + 1
        bgc_end   = num_physical_tracers + actual_bgc_num

        do i = bgc_start, bgc_end
            tracer_id = tracers_info%ids(i)

            !Iron (unit conversion: mol/L => umol/m3)
            if (tracer_id == tracer_ids%iron) then

                tracers_info%data_pointers(i)%tracer_data(:, :) = &
                        tracers_info%data_pointers(i)%tracer_data(:, :) * 1.e9

                ! Avoids tracers 1001, 1002, 1003, 1018, 1022 and DICremin (read from gen_ic3d)
            else if (tracer_id > 1003 .and. tracer_id /= 1018 .and. tracer_id /= 1022 .and. &
                    tracer_id /= tracer_ids%dic_remineralization) then  ! allowed since call initialize_tracer_ids happens earlier
                tracers_info%data_pointers(i)%tracer_data(:, :) = get_tracer_init_value(tracer_id)
            end if
        end do

        call mask_hydrothermal_vents(tracers_info, myDim_nod2D, eDim_nod2D, ulevels_nod2D, &
                nlevels_nod2D, geo_coord_nod2D, Z_3d_n, rad)
        call initialization_diagnostics(tracers_info, myDim_nod2D, ulevels_nod2D, nlevels_nod2D, &
                MPI_COMM_FESOM, mype)
    end subroutine recom_init

    subroutine initialize_memory(node_size, nl, num_tracers)
        use recom_declarations, only: decayBenthos

        use recom_glovar, only: GloFeDust, AtmFeInput, GloNDust, AtmNInput, RiverDIN2D, &
                RiverDON2D, RiverDOC2D, RiverDSi2D, RiverDIC2D, RiverAlk2D, RiverFe, &
                ErosionTON2D, ErosionTOC2D, ErosionTSi2D, relax_alk, virtual_alk, cosAI, &
                GloPCO2surf, GloCO2flux, GloO2flux, GloCO2flux_seaicemask, GloO2flux_seaicemask, &
                GlodPCO2surf, DenitBen, PistonVelocity, alphaCO2, GlodecayBenthos, Benthos, &
                Benthos_tr, GloHplus, PAR3D, NPPn, NPPd, NPPc, NPPp, GPPn, GPPd, GPPc, GPPp, &
                NNAn, NNAd, NNAc, NNAp, Chldegn, Chldegd, Chldegc, Chldegp, grazmeso_tot, &
                grazmeso_n, grazmeso_d, grazmeso_c, grazmeso_p, grazmeso_det, grazmeso_mic, &
                grazmeso_det2, grazmacro_tot, grazmacro_n, grazmacro_d, grazmacro_c, &
                grazmacro_p, grazmacro_mes, grazmacro_det, grazmacro_mic, grazmacro_det2, &
                grazmicro_tot, grazmicro_n, grazmicro_d, grazmicro_c, grazmicro_p, respmeso, &
                respmacro, respmicro, calcdiss, calcif, aggn, aggd, aggc, aggp, docexn, &
                docexd, docexc, docexp, respn, respd, respc, respp, NPPn3D, NPPd3D, NPPc3D, &
                NPPp3D, TTemp_diatoms, TTemp_phyto, TTemp_cocco, TTemp_phaeo, TPhyCO2, &
                TDiaCO2, TCoccoCO2, TPhaeoCO2, TqlimitFac_phyto, TqlimitFac_diatoms, &
                TqlimitFac_cocco, TqlimitFac_phaeo, TCphotLigLim_diatoms, TCphotLigLim_phyto, &
                TCphotLigLim_cocco, TCphotLigLim_phaeo, TCphot_diatoms, TCphot_phyto, &
                TCphot_cocco, TCphot_phaeo, TSi_assimDia, CO23D, pH3D, pCO23D, HCO33D, &
                CO33D, OmegaC3D, kspc3D, rhoSW3D, rho_particle1, rho_particle2, &
                scaling_density1_3D, scaling_density2_3D, scaling_visc_3D, seawater_visc_3D, &
                Sinkingvel1, Sinkingvel2, Sinkvel1_tr, Sinkvel2_tr, GloSed, SinkFlx, &
                SinkFlx_tr, lb_flux

        use recom_locvar, only: LocBenthos
        use recom_config, only: Diags, benthos_num, use_MEDUSA, bottflx_num, sedflx_num

        implicit none

        integer, intent(in) :: node_size, nl, num_tracers

        !! *** Allocate and initialize ***

        !! * Fe and N deposition as surface boundary condition *
        allocate(GloFeDust(node_size), source=0.d0)
        allocate(AtmFeInput(node_size), source=0.d0)
        allocate(GloNDust(node_size), source=0.d0)
        allocate(AtmNInput(node_size), source=0.d0)

        !! * River nutrients as surface boundary condition *
        allocate(RiverDIN2D(node_size), source=0.d0)
        allocate(RiverDON2D(node_size), source=0.d0)
        allocate(RiverDOC2D(node_size), source=0.d0)
        allocate(RiverDSi2D(node_size), source=0.d0)
        allocate(RiverDIC2D(node_size), source=0.d0)
        allocate(RiverAlk2D(node_size), source=0.d0)
        allocate(RiverFe(node_size), source=0.d0)

        !! * Erosion nutrients as surface boundary condition *
        allocate(ErosionTON2D(node_size), source=0.d0)
        allocate(ErosionTOC2D(node_size), source=0.d0)
        allocate(ErosionTSi2D(node_size), source=0.d0)

        !! * Alkalinity restoring to climatology *
        allocate(relax_alk(node_size), source=0.d0)
        allocate(virtual_alk(node_size), source=0.d0)

        allocate(cosAI(node_size), source=0.d0)
        allocate(GloPCO2surf(node_size), source=0.d0)
        allocate(GloCO2flux(node_size), source=0.d0)
        allocate(GloO2flux(node_size), source=0.d0)
        allocate(GloCO2flux_seaicemask(node_size), source=0.d0)
        allocate(GloO2flux_seaicemask(node_size), source=0.d0)
        allocate(GlodPCO2surf(node_size), source=0.d0)
        allocate(DenitBen(node_size), source=0.d0)
        allocate(PistonVelocity(node_size), source=0.d0)
        allocate(alphaCO2(node_size), source=0.d0)
        allocate(GlodecayBenthos(node_size, benthos_num), source=0.d0)
        allocate(Benthos(node_size, benthos_num), source=0.d0)

        ! kh 25.03.22 buffer per tracer index
        allocate(Benthos_tr(node_size, benthos_num, num_tracers), source=0.d0)
        allocate(GloHplus(node_size), source=exp(-8.d0 * log(10.d0)))

        allocate(LocBenthos(benthos_num), source=0.d0)

        ! [1/day] Decay rate of detritus in the benthic layer
        allocate(decayBenthos(benthos_num), source=0.d0)
        allocate(PAR3D(nl - 1, node_size), source=0.d0)

        if (Diags) then
            !! *** Allocate 2D diagnostics ***
            allocate(NPPn(node_size), source=0.d0)
            allocate(NPPd(node_size), source=0.d0)
            allocate(NPPc(node_size), source=0.d0)
            allocate(NPPp(node_size), source=0.d0)
            allocate(GPPn(node_size), source=0.d0)
            allocate(GPPd(node_size), source=0.d0)
            allocate(GPPc(node_size), source=0.d0)
            allocate(GPPp(node_size), source=0.d0)
            allocate(NNAn(node_size), source=0.d0)
            allocate(NNAd(node_size), source=0.d0)
            allocate(NNAc(node_size), source=0.d0)
            allocate(NNAp(node_size), source=0.d0)
            allocate(Chldegn(node_size), source=0.d0)
            allocate(Chldegd(node_size), source=0.d0)
            allocate(Chldegc(node_size), source=0.d0)
            allocate(Chldegp(node_size), source=0.d0)

            allocate(grazmeso_tot(node_size), source=0.d0)
            allocate(grazmeso_n(node_size), source=0.d0)
            allocate(grazmeso_d(node_size), source=0.d0)
            allocate(grazmeso_c(node_size), source=0.d0)
            allocate(grazmeso_p(node_size), source=0.d0)
            allocate(grazmeso_det(node_size), source=0.d0)
            allocate(grazmeso_mic(node_size), source=0.d0)
            allocate(grazmeso_det2(node_size), source=0.d0)

            allocate(grazmacro_tot(node_size), source=0.d0)
            allocate(grazmacro_n(node_size), source=0.d0)
            allocate(grazmacro_d(node_size), source=0.d0)
            allocate(grazmacro_c(node_size), source=0.d0)
            allocate(grazmacro_p(node_size), source=0.d0)
            allocate(grazmacro_mes(node_size), source=0.d0)
            allocate(grazmacro_det(node_size), source=0.d0)
            allocate(grazmacro_mic(node_size), source=0.d0)
            allocate(grazmacro_det2(node_size), source=0.d0)

            allocate(grazmicro_tot(node_size), source=0.d0)
            allocate(grazmicro_n(node_size), source=0.d0)
            allocate(grazmicro_d(node_size), source=0.d0)
            allocate(grazmicro_c(node_size), source=0.d0)
            allocate(grazmicro_p(node_size), source=0.d0)

            !! *** Allocate 3D diagnostics ***
            allocate(respmeso(nl - 1, node_size), source=0.d0)
            allocate(respmacro(nl - 1, node_size), source=0.d0)
            allocate(respmicro(nl - 1, node_size), source=0.d0)
            allocate(calcdiss(nl - 1, node_size), source=0.d0)
            allocate(calcif(nl - 1, node_size), source=0.d0)
            allocate(aggn(nl - 1, node_size), source=0.d0)
            allocate(aggd(nl - 1, node_size), source=0.d0)
            allocate(aggc(nl - 1, node_size), source=0.d0)
            allocate(aggp(nl - 1, node_size), source=0.d0)
            allocate(docexn(nl - 1, node_size), source=0.d0)
            allocate(docexd(nl - 1, node_size), source=0.d0)
            allocate(docexc(nl - 1, node_size), source=0.d0)
            allocate(docexp(nl - 1, node_size), source=0.d0)
            allocate(respn(nl - 1, node_size), source=0.d0)
            allocate(respd(nl - 1, node_size), source=0.d0)
            allocate(respc(nl - 1, node_size), source=0.d0)
            allocate(respp(nl - 1, node_size), source=0.d0)
            allocate(NPPn3D(nl - 1, node_size), source=0.d0)
            allocate(NPPd3D(nl - 1, node_size), source=0.d0)
            allocate(NPPc3D(nl - 1, node_size), source=0.d0)
            allocate(NPPp3D(nl - 1, node_size), source=0.d0)

            !! From Hannahs new temperature function (not sure if needed as diagnostic):
            allocate(TTemp_diatoms(nl - 1, node_size), source=0.d0)
            allocate(TTemp_phyto(nl - 1, node_size), source=0.d0)
            allocate(TTemp_cocco(nl - 1, node_size), source=0.d0)
            allocate(TTemp_phaeo(nl - 1, node_size), source=0.d0)

            allocate(TPhyCO2(nl - 1, node_size), source=0.d0)
            allocate(TDiaCO2(nl - 1, node_size), source=0.d0)
            allocate(TCoccoCO2(nl - 1, node_size), source=0.d0)
            allocate(TPhaeoCO2(nl - 1, node_size), source=0.d0)

            allocate(TqlimitFac_phyto(nl - 1, node_size), source=0.d0)
            allocate(TqlimitFac_diatoms(nl - 1, node_size), source=0.d0)
            allocate(TqlimitFac_cocco(nl - 1, node_size), source=0.d0)
            allocate(TqlimitFac_phaeo(nl - 1, node_size), source=0.d0)

            allocate(TCphotLigLim_diatoms(nl - 1, node_size), source=0.d0)
            allocate(TCphotLigLim_phyto(nl - 1, node_size), source=0.d0)
            allocate(TCphotLigLim_cocco(nl - 1, node_size), source=0.d0)
            allocate(TCphotLigLim_phaeo(nl - 1, node_size), source=0.d0)

            allocate(TCphot_diatoms(nl - 1, node_size), source=0.d0)
            allocate(TCphot_phyto(nl - 1, node_size), source=0.d0)
            allocate(TCphot_cocco(nl - 1, node_size), source=0.d0)
            allocate(TCphot_phaeo(nl - 1, node_size), source=0.d0)

            allocate(TSi_assimDia(nl - 1, node_size), source=0.d0)
        end if

        !! *** Allocate 3D mocsy ***
        allocate(CO23D(nl - 1, node_size), source=0.d0)
        allocate(pH3D(nl - 1, node_size), source=0.d0)
        allocate(pCO23D(nl - 1, node_size), source=0.d0)
        allocate(HCO33D(nl - 1, node_size), source=0.d0)
        allocate(CO33D(nl - 1, node_size), source=0.d0)
        allocate(OmegaC3D(nl - 1, node_size), source=0.d0)
        allocate(kspc3D(nl - 1, node_size), source=0.d0)
        allocate(rhoSW3D(nl - 1, node_size), source=0.d0)

        !! *** Allocate ballasting ***
        allocate(rho_particle1(nl - 1, node_size), source=0.d0)
        allocate(rho_particle2(nl - 1, node_size), source=0.d0)
        allocate(scaling_density1_3D(nl, node_size), source=0.d0)
        allocate(scaling_density2_3D(nl, node_size), source=0.d0)
        allocate(scaling_visc_3D(nl, node_size), source=0.d0)
        allocate(seawater_visc_3D(nl - 1, node_size), source=0.d0)

        allocate(Sinkingvel1(nl, node_size), source=0.d0)
        allocate(Sinkingvel2(nl, node_size), source=0.d0)

        allocate(Sinkvel1_tr(nl, node_size, num_tracers), source=0.d0) ! OG 16.03.23
        allocate(Sinkvel2_tr(nl, node_size, num_tracers), source=0.d0) ! OG 16.03.23

        if (use_MEDUSA) then
            allocate(GloSed(node_size, sedflx_num), source=0.d0)
            allocate(SinkFlx(node_size, bottflx_num), source=0.d0)
            ! kh 25.03.22 buffer sums per tracer index
            allocate(SinkFlx_tr(node_size, bottflx_num, num_tracers), source=0.d0)
            allocate(lb_flux(node_size, 9), source=0.d0)
        end if
    end subroutine initialize_memory

    subroutine initialize_tracer_ids
        use recom_declarations, only: tracer_ids
        use REcoM_config, only: enable_3zoo2det, enable_coccos

        integer :: current_tracer_id, total_tracers

        tracer_ids%dissolved_inorganic_nitrogen = 1001
        tracer_ids%dissolved_inorganic_carbon = 1002
        tracer_ids%alkalinity = 1003
        tracer_ids%phytoplankton_nitrogen = 1004
        tracer_ids%phytoplankton_carbon = 1005
        tracer_ids%phytoplankton_chlorophyll = 1006
        tracer_ids%detrital_nitrogen = 1007
        tracer_ids%detrital_carbon = 1008
        tracer_ids%heterotroph_nitrogen = 1009
        tracer_ids%heterotroph_carbon = 1010
        tracer_ids%dissolved_organic_nitrogen = 1011
        tracer_ids%dissolved_organic_carbon = 1012
        tracer_ids%diatom_nitrogen = 1013
        tracer_ids%diatom_carbon = 1014
        tracer_ids%diatom_chlorophyll = 1015
        tracer_ids%diatom_silica = 1016
        tracer_ids%detrital_silica = 1017
        tracer_ids%silica = 1018
        tracer_ids%iron = 1019
        tracer_ids%phytoplankton_calcite = 1020
        tracer_ids%detrital_calcite = 1021
        tracer_ids%oxygen = 1022

        ! DIC remineralization diagnostic tracer (added by Sina) -- always present, always last
        tracer_ids%dic_remineralization = 1037

        current_tracer_id = 1023

        if (enable_3zoo2det) then
            tracer_ids%macrozooplankton_nitrogen = current_tracer_id
            tracer_ids%macrozooplankton_carbon = current_tracer_id + 1
            tracer_ids%macrozooplankton_detrital_nitrogen = current_tracer_id + 2
            tracer_ids%macrozooplankton_detrital_carbon = current_tracer_id + 3
            tracer_ids%macrozooplankton_detrital_silica = current_tracer_id + 4
            tracer_ids%macrozooplankton_detrital_calcite = current_tracer_id + 5

            current_tracer_id = current_tracer_id + 6
        end if

        if (enable_coccos) then
            tracer_ids%coccolithophore_nitrogen = current_tracer_id
            tracer_ids%coccolithophore_carbon = current_tracer_id + 1
            tracer_ids%coccolithophore_chlorophyll = current_tracer_id + 2
            tracer_ids%phaeocystis_nitrogen = current_tracer_id + 3
            tracer_ids%phaeocystis_carbon = current_tracer_id + 4
            tracer_ids%phaeocystis_chlorophyll = current_tracer_id + 5

            current_tracer_id = current_tracer_id + 6
        end if

        if (enable_3zoo2det) then
            tracer_ids%microzooplankton_nitrogen = current_tracer_id
            tracer_ids%microzooplankton_carbon = current_tracer_id + 1
            current_tracer_id = current_tracer_id + 2   ! << required so DICremin doesn't overwrite microzoo_carbon's slot
        end if

    end subroutine initialize_tracer_ids

    function get_tracer_init_value(tracer_id) result(init_value)
        use recom_declarations, only: tracer_ids, wp
        use REcoM_config, only: tiny, tiny_chl, chl2N_max, NCmax, chl2N_max_d, NCmax_d, SiCmax, &
                Redfield

        integer, intent(in) :: tracer_id
        real(kind=wp) :: init_value

        if (tracer_id == tracer_ids%phytoplankton_nitrogen .or. &
                tracer_id == tracer_ids%diatom_nitrogen .or. &
                tracer_id == tracer_ids%coccolithophore_nitrogen .or. &
                tracer_id == tracer_ids%phaeocystis_nitrogen) then

            init_value = tiny_chl / chl2N_max

        else if (tracer_id == tracer_ids%phytoplankton_carbon .or. &
                    tracer_id == tracer_ids%diatom_carbon .or. &
                    tracer_id == tracer_ids%coccolithophore_carbon .or. &
                    tracer_id == tracer_ids%phaeocystis_carbon) then

            init_value = tiny_chl / chl2N_max / NCmax

        else if (tracer_id == tracer_ids%phytoplankton_chlorophyll .or. &
                    tracer_id == tracer_ids%diatom_chlorophyll .or. &
                    tracer_id == tracer_ids%coccolithophore_chlorophyll .or. &
                    tracer_id == tracer_ids%phaeocystis_chlorophyll) then

            init_value = tiny_chl

        else if (tracer_id == tracer_ids%detrital_nitrogen .or. &
                    tracer_id == tracer_ids%detrital_carbon .or. &
                    tracer_id == tracer_ids%heterotroph_nitrogen .or. &
                    tracer_id == tracer_ids%dissolved_organic_nitrogen .or. &
                    tracer_id == tracer_ids%dissolved_organic_carbon .or. &
                    tracer_id == tracer_ids%detrital_silica .or. &
                    tracer_id == tracer_ids%detrital_calcite .or. &
                    tracer_id == tracer_ids%macrozooplankton_nitrogen .or. &
                    tracer_id == tracer_ids%macrozooplankton_detrital_nitrogen .or. &
                    tracer_id == tracer_ids%macrozooplankton_detrital_carbon .or. &
                    tracer_id == tracer_ids%macrozooplankton_detrital_silica .or. &
                    tracer_id == tracer_ids%macrozooplankton_detrital_calcite .or. &
                    tracer_id == tracer_ids%microzooplankton_nitrogen) then

            init_value = tiny

        else if (tracer_id == tracer_ids%heterotroph_carbon .or. &
                    tracer_id == tracer_ids%phytoplankton_calcite .or. &
                    tracer_id == tracer_ids%macrozooplankton_carbon .or. &
                    tracer_id == tracer_ids%microzooplankton_carbon) then

            init_value = tiny * Redfield

        else if (tracer_id == tracer_ids%diatom_silica) then

            init_value = tiny_chl / chl2N_max_d / NCmax_d / SiCmax

        else
            init_value = 0.0_wp
            write(*, *) 'Warning: No initial value defined for tracer ID ', tracer_id, '.' // &
                    ' Setting to 0'
        end if
    end function get_tracer_init_value

    subroutine mask_hydrothermal_vents(tracers_info, myDim_nod2D, eDim_nod2D, ulevels_nod2D, &
            nlevels_nod2D, geo_coord_nod2D, Z_3d_n, rad)

        use REcoM_glovar, only: tracers_info_type
        use REcoM_config, only: tiny
        use recom_declarations, only: wp

        implicit none

        type(tracers_info_type), intent(in) :: tracers_info
        integer, intent(in) :: myDim_nod2D, eDim_nod2D
        integer, intent(in) :: ulevels_nod2D(:), nlevels_nod2D(:)
        real(kind=wp), intent(in) :: rad, geo_coord_nod2D(:, :), Z_3d_n(:, :)

        integer :: row, k, nzmin, nzmax

        ! Iron tracer (ID=1019) is always BGC tracer #17 in the base block
        ! (1004..1022 shifted by 1003), so with the fixed T,S | BGC | [age] |
        ! [transit] layout it always sits at slot 21 (= 2 physical + 19th BGC
        ! slot: 1019 - 1000 = 19 -> n_base_physical + 19 = 2 + 19 = 21).
        integer, parameter :: iron_slot = 21

        ! Mask hydrothermal vent in Eastern Equatorial Pacific GO
        do row = 1, myDim_nod2D + eDim_nod2D
            !if (ulevels_nod2D(row)>1) cycle
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row) - 1
            do k = nzmin, nzmax
                ! do not take regions shallower than 2000 m into account
                if (((geo_coord_nod2D(2, row) > -12.5 * rad) .and. (geo_coord_nod2D(2, row) < 9.5 &
                        &* rad))&
                        .and. ((geo_coord_nod2D(1, row) > -106.0 * rad) .and. (geo_coord_nod2D(1, &
                        row) < -72.0 * rad))) then
                    if (abs(Z_3d_n(k, row)) < 2000.0_WP) cycle
                    tracers_info%data_pointers(iron_slot)%tracer_data(k, row) = &
                            min(0.3, tracers_info%data_pointers(iron_slot)%tracer_data(k, row)) ! OG todo: try 0.6
                end if
            end do
        end do

        ! Mask negative values
        tracers_info%data_pointers(iron_slot)%tracer_data(:, :) = &
                max(tiny, tracers_info%data_pointers(iron_slot)%tracer_data(:, :))
    end subroutine mask_hydrothermal_vents

    subroutine initialization_diagnostics(tracers_info, myDim_nod2D, ulevels_nod2D, nlevels_nod2D, &
            MPI_COMM_FESOM, mype)
        use REcoM_glovar, only: tracers_info_type
        use recom_config, only: enable_3zoo2det, enable_coccos
        use recom_declarations, only: WP, is_3zoo2det, is_coccos

        use mpi

        implicit none

        type(tracers_info_type), intent(in) :: tracers_info
        integer, intent(in) :: myDim_nod2D
        integer, intent(in) :: ulevels_nod2D(:), nlevels_nod2D(:)
        integer, intent(in) :: MPI_COMM_FESOM, mype

        integer :: MPIerr, n
        real(kind=WP) :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
        real(kind=WP) :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin
        real(kind=WP) :: locO2max, locO2min

        ! With the fixed T,S | BGC | [age] | [transit] layout, BGC tracers
        ! always start at slot 3 (right after T,S) regardless of use_transit.
        ! These fixed slot numbers are therefore now constant and correct:
        !   DIN (1001) -> slot 3, DIC (1002) -> slot 4, Alk (1003) -> slot 5,
        !   DSi (1018) -> slot 20, DFe (1019) -> slot 21, O2 (1022) -> slot 24
        ! (previously these shifted whenever use_transit was active, which was
        ! a bug: this fixed indexing is now always correct, independent of
        ! use_transit.)
        integer, parameter :: din_slot = 3
        integer, parameter :: dic_slot = 4
        integer, parameter :: alk_slot = 5
        integer, parameter :: dsi_slot = 20
        integer, parameter :: dfe_slot = 21
        integer, parameter :: o2_slot  = 24

        if (mype == 0) write(*, *) 'Tracers have been initialized as spinup from WOA/glodap' // &
                ' netcdf files'
        locDINmax = -66666
        locDINmin = 66666
        locDICmax = locDINmax
        locDICmin = locDINmin
        locAlkmax = locDINmax
        locAlkmin = locDINmin
        locDSimax = locDINmax
        locDSimin = locDINmin
        locDFemax = locDINmax
        locDFemin = locDINmin
        locO2max = locDINmax
        locO2min = locDINmin

        do n = 1, myDim_nod2d
            locDINmax = max(locDINmax, maxval(tracers_info%data_pointers(din_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDINmin = min(locDINmin, minval(tracers_info%data_pointers(din_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDICmax = max(locDICmax, maxval(tracers_info%data_pointers(dic_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDICmin = min(locDICmin, minval(tracers_info%data_pointers(dic_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locAlkmax = max(locAlkmax, maxval(tracers_info%data_pointers(alk_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locAlkmin = min(locAlkmin, minval(tracers_info%data_pointers(alk_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDSimax = max(locDSimax, maxval(tracers_info%data_pointers(dsi_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDSimin = min(locDSimin, minval(tracers_info%data_pointers(dsi_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDFemax = max(locDFemax, maxval(tracers_info%data_pointers(dfe_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locDFemin = min(locDFemin, minval(tracers_info%data_pointers(dfe_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locO2max = max(locO2max, maxval(tracers_info%data_pointers(o2_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
            locO2min = min(locO2min, minval(tracers_info%data_pointers(o2_slot)%tracer_data(&
                    ulevels_nod2D(n):nlevels_nod2D(n) - 1, n)))
        end do

        if (mype == 0) write(*, *) "Sanity check for REcoM variables after recom_init call"
        call MPI_AllREDUCE(locDINmax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. DIN. =', glo
        call MPI_AllREDUCE(locDINmin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal min init. DIN. =', glo

        call MPI_AllREDUCE(locDICmax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. DIC. =', glo
        call MPI_AllREDUCE(locDICmin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal min init. DIC. =', glo
        call MPI_AllREDUCE(locAlkmax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. Alk. =', glo
        call MPI_AllREDUCE(locAlkmin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal min init. Alk. =', glo
        call MPI_AllREDUCE(locDSimax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. DSi. =', glo
        call MPI_AllREDUCE(locDSimin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal min init. DSi. =', glo
        call MPI_AllREDUCE(locDFemax, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. DFe. =', glo
        call MPI_AllREDUCE(locDFemin, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, &
                MPIerr)
        if (mype == 0) write(*, *) '  `-> gobal min init. DFe. =', glo
        call MPI_AllREDUCE(locO2max, glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype == 0) write(*, *) '  |-> gobal max init. O2. =', glo
        call MPI_AllREDUCE(locO2min, glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype == 0) write(*, *) '  `-> gobal min init. O2. =', glo

        if (enable_3zoo2det) then
            is_3zoo2det = 1.0_WP
        else
            is_3zoo2det = 0.0_WP
        end if

        if (enable_coccos) then
            is_coccos = 1.0_WP
        else
            is_coccos = 0.0_WP
        end if
    end subroutine initialization_diagnostics

end module recom_init_interface
