! ------------
! 23.03.2023
! OG
!===============================================================================
! allocate & initialise arrays for REcoM
module recom_init_interface
contains
!
!
!_______________________________________________________________________________
    subroutine recom_init(nl, ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, Z_3d_n,   &
                          myDim_nod2d, eDim_nod2D, mype, MPI_COMM_FESOM, myDim_elem2D, &
                          eDim_elem2D, tracers_info, num_tracers, rad)

        use REcoM_declarations
        use REcoM_GloVar
        use REcoM_locVar
        use recom_config
        use REcoM_ciso

        implicit none

        integer,        intent(in)                  :: nl, mydim_nod2d, edim_nod2d, mype, num_tracers
        integer,        intent(in)                  :: mpi_comm_fesom, mydim_elem2d, edim_elem2d
        real(kind=WP), intent(in)                  :: rad
        integer,        intent(in), dimension(:)    :: ulevels_nod2d, nlevels_nod2d
        real(kind=wp),  intent(in), dimension(:, :) :: geo_coord_nod2d, z_3d_n
        type(tracers_info_type), intent(in) :: tracers_info

        !___________________________________________________________________________
        ! pointer on necessary derived types
        integer                                 :: n, k, row, nzmin, nzmax, i, id
        integer                                 :: elem_size, node_size

        ! After reading tracer namelist - validate actual IDs
        integer, dimension(num_tracers) :: tracer_id_array

        elem_size   = myDim_elem2D + eDim_elem2D
        node_size   = myDim_nod2D + eDim_nod2D

        call initialize_memory(node_size, nl, num_tracers)

        ! After reading parecomsetup namelist
        call initialize_tracer_indices

        ! Validation check here
        call validate_recom_tracers(num_tracers, mype)

        ! ... populate tracer_id_array from namelist ...
        !tracer_id_array = tracers%data(1:tracers%num_tracers)%ID
        tracer_id_array = tracers_info%ids(1:num_tracers)
        call validate_tracer_id_sequence(tracer_id_array, num_tracers, mype)

        call initialize_tracer_data(num_tracers, tracers_info)

    !------------------------------------------

        call mask_hydrothermal_vents(tracers_info, myDim_nod2D, eDim_nod2D, ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, Z_3d_n, rad)
    !------------------------------------------

        call initialization_diagnostics(tracers_info, myDim_nod2D, ulevels_nod2D, nlevels_nod2D, MPI_COMM_FESOM, mype)
    end subroutine recom_init

    subroutine initialize_memory(node_size, nl, num_tracers)
        use recom_declarations
        use recom_glovar
        use recom_locvar
        use recom_config

        implicit none

        integer, intent(in) :: node_size, nl, num_tracers

        !! *** Allocate and initialize ***

        !! * Fe and N deposition as surface boundary condition *
        allocate(GloFeDust             ( node_size ), source = 0.d0)
        allocate(AtmFeInput            ( node_size ), source = 0.d0)
        allocate(GloNDust              ( node_size ), source = 0.d0)
        allocate(AtmNInput             ( node_size ), source = 0.d0)

        !! * River nutrients as surface boundary condition *
        allocate(RiverDIN2D            ( node_size ), source = 0.d0)
        allocate(RiverDON2D            ( node_size ), source = 0.d0)
        allocate(RiverDOC2D            ( node_size ), source = 0.d0)
        allocate(RiverDSi2D            ( node_size ), source = 0.d0)
        allocate(RiverDIC2D            ( node_size ), source = 0.d0)
        allocate(RiverAlk2D            ( node_size ), source = 0.d0)
        allocate(RiverFe               ( node_size ), source = 0.d0)

        !! * Erosion nutrients as surface boundary condition *
        allocate(ErosionTON2D          ( node_size ), source = 0.d0)
        allocate(ErosionTOC2D          ( node_size ), source = 0.d0)
        allocate(ErosionTSi2D          ( node_size ), source = 0.d0)

        !! * Alkalinity restoring to climatology *
        allocate(relax_alk             ( node_size ), source = 0.d0)
        allocate(virtual_alk           ( node_size ), source = 0.d0)

        allocate(cosAI                 ( node_size ), source = 0.d0)
        allocate(GloPCO2surf           ( node_size ), source = 0.d0)
        allocate(GloCO2flux            ( node_size ), source = 0.d0)
        allocate(GloO2flux             ( node_size ), source = 0.d0)
        allocate(GloCO2flux_seaicemask ( node_size ), source = 0.d0)
        allocate(GloO2flux_seaicemask  ( node_size ), source = 0.d0)
        allocate(GlodPCO2surf          ( node_size ), source = 0.d0)
        allocate(DenitBen              ( node_size ), source = 0.d0)
        allocate(PistonVelocity        ( node_size ), source = 0.d0)
        allocate(alphaCO2              ( node_size ), source = 0.d0)
        allocate(GlodecayBenthos       ( node_size, benthos_num ), source = 0.d0)
        allocate(Benthos               ( node_size, benthos_num ), source = 0.d0)
        allocate(Benthos_tr            ( node_size, benthos_num, num_tracers ), source = 0.d0) ! kh 25.03.22 buffer per tracer index
        allocate(GloHplus              ( node_size ), source = exp(-8.d0 * log(10.d0)))

        allocate(LocBenthos            ( benthos_num ), source = 0.d0)
        allocate(decayBenthos          ( benthos_num ), source = 0.d0)     ! [1/day] Decay rate of detritus in the benthic layer
        allocate(PAR3D                 ( nl-1, node_size ), source = 0.d0)

        if (Diags) then
            !! *** Allocate 2D diagnostics ***
            allocate(NPPn    ( node_size ), source = 0.d0)
            allocate(NPPd    ( node_size ), source = 0.d0)
            allocate(NPPc    ( node_size ), source = 0.d0)
            allocate(NPPp    ( node_size ), source = 0.d0)
            allocate(GPPn    ( node_size ), source = 0.d0)
            allocate(GPPd    ( node_size ), source = 0.d0)
            allocate(GPPc    ( node_size ), source = 0.d0)
            allocate(GPPp    ( node_size ), source = 0.d0)
            allocate(NNAn    ( node_size ), source = 0.d0)
            allocate(NNAd    ( node_size ), source = 0.d0)
            allocate(NNAc    ( node_size ), source = 0.d0)
            allocate(NNAp    ( node_size ), source = 0.d0)
            allocate(Chldegn ( node_size ), source = 0.d0)
            allocate(Chldegd ( node_size ), source = 0.d0)
            allocate(Chldegc ( node_size ), source = 0.d0)
            allocate(Chldegp ( node_size ), source = 0.d0)

            allocate(grazmeso_tot(node_size), source = 0.d0)
            allocate(grazmeso_n(node_size), source = 0.d0)
            allocate(grazmeso_d(node_size), source = 0.d0)
            allocate(grazmeso_c(node_size), source = 0.d0)
            allocate(grazmeso_p(node_size), source = 0.d0)
            allocate(grazmeso_det(node_size), source = 0.d0)
            allocate(grazmeso_mic(node_size), source = 0.d0)
            allocate(grazmeso_det2(node_size), source = 0.d0)

            allocate(grazmacro_tot(node_size), source = 0.d0)
            allocate(grazmacro_n(node_size), source = 0.d0)
            allocate(grazmacro_d(node_size), source = 0.d0)
            allocate(grazmacro_c(node_size), source = 0.d0)
            allocate(grazmacro_p(node_size), source = 0.d0)
            allocate(grazmacro_mes(node_size), source = 0.d0)
            allocate(grazmacro_det(node_size), source = 0.d0)
            allocate(grazmacro_mic(node_size), source = 0.d0)
            allocate(grazmacro_det2(node_size), source = 0.d0)

            allocate(grazmicro_tot(node_size), source = 0.d0)
            allocate(grazmicro_n(node_size), source = 0.d0)
            allocate(grazmicro_d(node_size), source = 0.d0)
            allocate(grazmicro_c(node_size), source = 0.d0)
            allocate(grazmicro_p(node_size), source = 0.d0)

            !! *** Allocate 3D diagnostics ***
            allocate(respmeso     ( nl-1, node_size ), source = 0.d0)
            allocate(respmacro    ( nl-1, node_size ), source = 0.d0)
            allocate(respmicro    ( nl-1, node_size ), source = 0.d0)
            allocate(calcdiss     ( nl-1, node_size ), source = 0.d0)
            allocate(calcif       ( nl-1, node_size ), source = 0.d0)
            allocate(aggn         ( nl-1, node_size ), source = 0.d0)
            allocate(aggd         ( nl-1, node_size ), source = 0.d0)
            allocate(aggc         ( nl-1, node_size ), source = 0.d0)
            allocate(aggp         ( nl-1, node_size ), source = 0.d0)
            allocate(docexn       ( nl-1, node_size ), source = 0.d0)
            allocate(docexd       ( nl-1, node_size ), source = 0.d0)
            allocate(docexc       ( nl-1, node_size ), source = 0.d0)
            allocate(docexp       ( nl-1, node_size ), source = 0.d0)
            allocate(respn        ( nl-1, node_size ), source = 0.d0)
            allocate(respd        ( nl-1, node_size ), source = 0.d0)
            allocate(respc        ( nl-1, node_size ), source = 0.d0)
            allocate(respp        ( nl-1, node_size ), source = 0.d0)
            allocate(NPPn3D       ( nl-1, node_size ), source = 0.d0)
            allocate(NPPd3D       ( nl-1, node_size ), source = 0.d0)
            allocate(NPPc3D       ( nl-1, node_size ), source = 0.d0)
            allocate(NPPp3D       ( nl-1, node_size ), source = 0.d0)

            !! From Hannahs new temperature function (not sure if needed as diagnostic):
            allocate(TTemp_diatoms  (nl-1,node_size), source = 0.d0)
            allocate(TTemp_phyto    (nl-1,node_size), source = 0.d0)
            allocate(TTemp_cocco    (nl-1,node_size), source = 0.d0)
            allocate(TTemp_phaeo    (nl-1,node_size), source = 0.d0)

            allocate(TPhyCO2        (nl-1,node_size), source = 0.d0)
            allocate(TDiaCO2        (nl-1,node_size), source = 0.d0)
            allocate(TCoccoCO2      (nl-1,node_size), source = 0.d0)
            allocate(TPhaeoCO2      (nl-1,node_size), source = 0.d0)

            allocate(TqlimitFac_phyto     (nl-1,node_size), source = 0.d0)
            allocate(TqlimitFac_diatoms   (nl-1,node_size), source = 0.d0)
            allocate(TqlimitFac_cocco     (nl-1,node_size), source = 0.d0)
            allocate(TqlimitFac_phaeo     (nl-1,node_size), source = 0.d0)

            allocate(TCphotLigLim_diatoms    (nl-1,node_size), source = 0.d0)
            allocate(TCphotLigLim_phyto      (nl-1,node_size), source = 0.d0)
            allocate(TCphotLigLim_cocco      (nl-1,node_size), source = 0.d0)
            allocate(TCphotLigLim_phaeo      (nl-1,node_size), source = 0.d0)

            allocate(TCphot_diatoms       (nl-1,node_size), source = 0.d0)
            allocate(TCphot_phyto         (nl-1,node_size), source = 0.d0)
            allocate(TCphot_cocco         (nl-1,node_size), source = 0.d0)
            allocate(TCphot_phaeo         (nl-1,node_size), source = 0.d0)

            allocate(TSi_assimDia         (nl-1,node_size), source = 0.d0)
        end if

        !! *** Allocate 3D mocsy ***
        allocate(CO23D        ( nl-1, node_size ), source = 0.d0)
        allocate(pH3D         ( nl-1, node_size ), source = 0.d0)
        allocate(pCO23D       ( nl-1, node_size ), source = 0.d0)
        allocate(HCO33D       ( nl-1, node_size ), source = 0.d0)
        allocate(CO33D        ( nl-1, node_size ), source = 0.d0)
        allocate(OmegaC3D     ( nl-1, node_size ), source = 0.d0)
        allocate(kspc3D       ( nl-1, node_size ), source = 0.d0)
        allocate(rhoSW3D      ( nl-1, node_size ), source = 0.d0)

        !! *** Allocate ballasting ***
        allocate(rho_particle1       ( nl-1, node_size ), source = 0.d0)
        allocate(rho_particle2       ( nl-1, node_size ), source = 0.d0)
        allocate(scaling_density1_3D ( nl,   node_size ), source = 0.d0)
        allocate(scaling_density2_3D ( nl,   node_size ), source = 0.d0)
        allocate(scaling_visc_3D     ( nl,   node_size ), source = 0.d0)
        allocate(seawater_visc_3D    ( nl-1, node_size ), source = 0.d0)

        allocate(Sinkingvel1(nl,node_size), source = 0.d0)
        allocate(Sinkingvel2(nl,node_size), source = 0.d0)

        allocate(Sinkvel1_tr(nl,node_size,num_tracers), source = 0.d0)  ! OG 16.03.23
        allocate(Sinkvel2_tr(nl,node_size,num_tracers), source = 0.d0)  ! OG 16.03.23

        if (use_MEDUSA) then
            allocate(GloSed(node_size,sedflx_num), source = 0.d0)
            allocate(SinkFlx(node_size,bottflx_num), source = 0.d0)
            allocate(SinkFlx_tr(node_size,bottflx_num,num_tracers), source = 0.d0) ! kh 25.03.22 buffer sums per tracer index
            allocate(lb_flux(node_size,9), source = 0.d0)
        end if
    end subroutine initialize_memory

    subroutine initialize_tracer_data(num_tracers, tracers_info)
        use REcoM_glovar
        use REcoM_config

        implicit none

        integer, intent(in) :: num_tracers
        type(tracers_info_type), intent(in) :: tracers_info

        integer :: i, id
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

        DO i=num_tracers-bgc_num+1, num_tracers
            id=tracers_info%ids(i)

            SELECT CASE (id)

                !---------------------------------------------------------------------------
                ! Base Model: 2 Phytoplankton + 1 Zooplankton + 1 Detritus
                !---------------------------------------------------------------------------
                ! Skip: DIN, DIC, Alk, DSi and O2 are read from files
                ! Fe [mol/L] => [umol/m3] Check the units again!

                ! --- Small Phytoplankton
                CASE (1004)  ! PhyN - Phytoplankton Nitrogen
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl/chl2N_max

                CASE (1005)  ! PhyC - Phytoplankton Carbon
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl/chl2N_max/NCmax

                CASE (1006)  ! PhyChl - Phytoplankton Chlorophyll
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl

                ! --- Detritus (Non-living organic matter) ---
                CASE (1007)  ! DetN - Detrital Nitrogen
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                CASE (1008)  ! DetC - Detrital Carbon
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                ! --- Mesozooplankton (Heterotrophs) ---
                CASE (1009)  ! HetN - Heterotroph Nitrogen
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                CASE (1010)  ! HetC - Heterotroph Carbon (using Redfield ratio)
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny * Redfield

                ! --- Dissolved Organic Matter ---
                CASE (1011)  ! DON - Dissolved Organic Nitrogen
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                CASE (1012)  ! DOC - Dissolved Organic Carbon
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                ! --- Diatoms ---
                CASE (1013)  ! DiaN - Diatom Nitrogen
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl/chl2N_max

                CASE (1014)  ! DiaC - Diatom Carbon
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl/chl2N_max/NCmax

                CASE (1015)  ! DiaChl - Diatom Chlorophyll
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl

                CASE (1016)  ! DiaSi - Diatom Silica
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl/chl2N_max_d/NCmax_d/SiCmax

                CASE (1017)  ! DetSi - Detrital Silica
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                ! --- Iron (micronutrient) ---
                CASE (1019)  ! Fe - Iron (unit conversion: mol/L => umol/m3)
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tracers_info%data_pointers(i)%tracer_data(:,:)* 1.e9

                ! --- Calcium Carbonate (Calcite) ---
                CASE (1020)  ! PhyCalc - Phytoplankton Calcite
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny * Redfield

                CASE (1021)  ! DetCalc - Detrital Calcite
                    tracers_info%data_pointers(i)%tracer_data(:,:) = tiny

                !---------------------------------------------------------------------------
                ! Extended Model: Additional Zooplankton and Detritus (enable_3zoo2det)
                !---------------------------------------------------------------------------

                CASE (1023)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! Zoo2N - Macrozooplankton Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! CoccoN - Coccolithophore Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max
                    END IF

                CASE (1024)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! Zoo2C - Macrozooplankton Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny * Redfield
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! CoccoC - Coccolithophore Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max / NCmax
                    END IF

                CASE (1025)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! DetZ2N - Macrozooplankton Detrital Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! CoccoChl - Coccolithophore Chlorophyll
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl
                    END IF

                CASE (1026)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! DetZ2C - Macrozooplankton Detrital Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! PhaeoN - Phaeocystis Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max
                    END IF

                CASE (1027)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! DetZ2Si - Zooplankton 2 Detrital Silica
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! PhaeoC - Phaeocystis Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max / NCmax
                    END IF

                CASE (1028)
                    IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! DetZ2Calc - Macrozooplankton Detrital Calcite
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    ELSE IF (enable_coccos .AND. .NOT. enable_3zoo2det) THEN
                        ! PhaeoChl - Phaeocystis Chlorophyll
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl
                    END IF

                !---------------------------------------------------------------------------
                ! Extended Model: Coccolithophores with 3Zoo2Det
                !---------------------------------------------------------------------------

                CASE (1029)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! CoccoN - Coccolithophore Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max
                    ELSE IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! Zoo3N - Microzooplankton Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    END IF

                CASE (1030)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! CoccoC - Coccolithophore Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max / NCmax
                    ELSE IF (enable_3zoo2det .AND. .NOT. enable_coccos) THEN
                        ! Zoo3C - Microzooplankton Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny * Redfield
                    END IF

                CASE (1031)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! CoccoChl - Coccolithophore Chlorophyll
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl
                    END IF

                CASE (1032)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! PhaeoN - Phaeocystis Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max
                    END IF

                CASE (1033)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! PhaeoC - Phaeocystis Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl / chl2N_max / NCmax
                    END IF

                CASE (1034)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! PhaeoChl - Phaeocystis Chlorophyll
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny_chl
                    END IF

                CASE (1035)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! Zoo3N - Zooplankton 3 Nitrogen
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny
                    END IF

                CASE (1036)
                    IF (enable_coccos .AND. enable_3zoo2det) THEN
                        ! Zoo3C - Zooplankton 3 Carbon
                        tracers_info%data_pointers(i)%tracer_data(:,:) = tiny * Redfield
                    END IF

            END SELECT
        END DO
    end subroutine initialize_tracer_data

    subroutine mask_hydrothermal_vents(tracers_info, myDim_nod2D, eDim_nod2D, ulevels_nod2D, nlevels_nod2D, geo_coord_nod2D, Z_3d_n, rad)

        use REcoM_glovar
        use REcoM_config
        use recom_declarations

        implicit none

        type(tracers_info_type), intent(in) :: tracers_info
        integer, intent(in) :: myDim_nod2D, eDim_nod2D
        integer, intent(in) :: ulevels_nod2D(:), nlevels_nod2D(:)
        real(wp), intent(in) :: rad, geo_coord_nod2D(:, :), Z_3d_n(:,:)

        integer :: row, k, nzmin, nzmax

        !< Mask hydrothermal vent in Eastern Equatorial Pacific GO
        do row=1, myDim_nod2D+eDim_nod2D
            !if (ulevels_nod2D(row)>1) cycle
            nzmin = ulevels_nod2D(row)
            nzmax = nlevels_nod2D(row)-1
            do k=nzmin, nzmax
                ! do not take regions shallower than 2000 m into account
                if (((geo_coord_nod2D(2,row) > -12.5*rad) .and. (geo_coord_nod2D(2,row) < 9.5*rad))&
                    .and.((geo_coord_nod2D(1,row)> -106.0*rad) .and. (geo_coord_nod2D(1,row) < -72.0*rad))) then
                    if (abs(Z_3d_n(k,row))<2000.0_WP) cycle
                    tracers_info%data_pointers(21)%tracer_data(k,row) = min(0.3, tracers_info%data_pointers(21)%tracer_data(k,row)) ! OG todo: try 0.6
                end if
            end do
        end do

        !< Mask negative values
        tracers_info%data_pointers(21)%tracer_data(:,:) = max(tiny, tracers_info%data_pointers(21)%tracer_data(:,:))
    end subroutine mask_hydrothermal_vents

    subroutine initialization_diagnostics(tracers_info, myDim_nod2D, ulevels_nod2D, nlevels_nod2D, MPI_COMM_FESOM, mype)
        use REcoM_glovar
        use recom_declarations
        use recom_config

        use mpi

        implicit none

        type(tracers_info_type), intent(in) :: tracers_info
        integer, intent(in) :: myDim_nod2D
        integer, intent(in) :: ulevels_nod2D(:), nlevels_nod2D(:)
        integer, intent(in) :: MPI_COMM_FESOM, mype

        integer       :: MPIerr, n
        real(kind=WP) :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax, glo
        real(kind=WP) :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin
        real(kind=WP) :: locO2max, locO2min


        if(mype==0) write(*,*) 'Tracers have been initialized as spinup from WOA/glodap netcdf files'
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
        locO2max  = locDINmax
        locO2min  = locDINmin

        do n=1, myDim_nod2d
            locDINmax = max(locDINmax,maxval(tracers_info%data_pointers(3)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDINmin = min(locDINmin,minval(tracers_info%data_pointers(3)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmax = max(locDICmax,maxval(tracers_info%data_pointers(4)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDICmin = min(locDICmin,minval(tracers_info%data_pointers(4)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmax = max(locAlkmax,maxval(tracers_info%data_pointers(5)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locAlkmin = min(locAlkmin,minval(tracers_info%data_pointers(5)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimax = max(locDSimax,maxval(tracers_info%data_pointers(20)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDSimin = min(locDSimin,minval(tracers_info%data_pointers(20)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemax = max(locDFemax,maxval(tracers_info%data_pointers(21)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locDFemin = min(locDFemin,minval(tracers_info%data_pointers(21)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2max  = max(locO2max,maxval( tracers_info%data_pointers(24)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
            locO2min  = min(locO2min,minval( tracers_info%data_pointers(24)%tracer_data(ulevels_nod2D(n):nlevels_nod2D(n)-1,n)) )
        end do

        if (mype==0) write(*,*) "Sanity check for REcoM variables after recom_init call"
        call MPI_AllREDUCE(locDINmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIN. =', glo
        call MPI_AllREDUCE(locDINmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIN. =', glo

        call MPI_AllREDUCE(locDICmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DIC. =', glo
        call MPI_AllREDUCE(locDICmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DIC. =', glo
        call MPI_AllREDUCE(locAlkmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. Alk. =', glo
        call MPI_AllREDUCE(locAlkmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. Alk. =', glo
        call MPI_AllREDUCE(locDSimax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DSi. =', glo
        call MPI_AllREDUCE(locDSimin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal min init. DSi. =', glo
        call MPI_AllREDUCE(locDFemax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. DFe. =', glo
        call MPI_AllREDUCE(locDFemin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. DFe. =', glo
        call MPI_AllREDUCE(locO2max , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  |-> gobal max init. O2. =', glo
        call MPI_AllREDUCE(locO2min , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
        if (mype==0) write(*,*) '  `-> gobal min init. O2. =', glo

        if (enable_3zoo2det) then
            is_3zoo2det=1.0_WP
        else
            is_3zoo2det=0.0_WP
        endif

        if (enable_coccos) then
            is_coccos=1.0_WP
        else
            is_coccos=0.0_WP
        endif
    end subroutine initialization_diagnostics

end module
