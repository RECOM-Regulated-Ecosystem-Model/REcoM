!>
!! @par Copyright
!! This code is subject to the FESOM-REcoM - License - Agreement in it's most recent form.
!! Please see URL xxx
!!
!! @brief Module for defining variables used in REcoM, ex constant sinking velocity and
!! local time step dt
!!
!! @remarks This module contains namelist for recom
!! @author xxx, FESOM-REcoM, Bremerhaven (2019)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   2014: original implementation (V. Schourup-Kristensen)
!

module REcoM_declarations
    implicit none
    public

    save

    integer, parameter :: WP = 8 ! Working precision
    real(kind=WP), parameter :: pi = 3.14159265358979

    integer :: save_count_recom
    real(kind=wp) :: tiny_N ! Min PhyN
    real(kind=wp) :: tiny_N_d ! Min DiaN
    real(kind=wp) :: tiny_N_c ! Min CocN                 ! NEW
    real(kind=wp) :: tiny_N_p ! Min PhaN                 ! Phaeocystis
    real(kind=wp) :: tiny_C ! Min PhyC
    real(kind=wp) :: tiny_C_d ! Min DiaC
    real(kind=wp) :: tiny_C_c ! Min CocC                 ! NEW
    real(kind=wp) :: tiny_C_p ! Min PhaC                 ! Phaeocystis
    real(kind=wp) :: tiny_Si ! Min DiaSi
    !!------------------------------------------------------------------------------
    !! *** Temperature dependence of rates ***
    real(kind=wp) :: rTref ! [1/K] Reciproque value of reference temp for Arrhenius function
    real(kind=wp) :: rTloc ! [1/K] Reciproque of local ocean temp
    real(kind=wp) :: arrFunc ! []    Temp dependence of rates (also for Phaeocystis)
    real(kind=wp) :: CoccoTFunc ! []    Temp dependence of coccolithophores
    real(kind=wp) :: Temp_diatoms ! []    Temp dependence of diatoms
    real(kind=wp) :: Temp_phyto ! []    Temp dependence of small phyto
    real(kind=wp) :: Temp_cocco ! []    Temp dependence of coccolithophores
    real(kind=wp) :: Temp_phaeo ! []    Temp dependence of phaeocystis
    real(kind=wp) :: arrFuncZoo2 ! []    Temperature function for krill
    real(kind=wp) :: q10_mic ! 3Zoo
    real(kind=wp) :: q10_mic_res ! 3Zoo
    real(kind=wp) :: q10_mes ! 3Zoo
    real(kind=wp) :: q10_mes_res ! 3Zoo
    real(kind=wp) :: reminSiT
    real(kind=wp) :: O2Func ! O2remin
    !!------------------------------------------------------------------------------
    !! *** CO2 dependence of rates ! NEW CO2 ***
    real(kind=wp) :: h_depth(1) ! pH from mocsy is converted to proton concentration
    !Real(kind=wp)  :: d_CT_CL_phy            ! NEW inter For the interaction term between CO2 and
    !both temperature and light
    !Real(kind=wp)  :: d_CT_CL_dia
    !Real(kind=wp)  :: d_CT_CL_coc
    real(kind=wp) :: CoccoCO2
    real(kind=wp) :: DiaCO2
    real(kind=wp) :: PhyCO2
    real(kind=wp) :: PhaeoCO2

    !!------------------------------------------------------------------------------
    !! *** Quotas ***
    ! [mmol N/mmol C]  Quota between phytoplankton N and C (NEW changed term)
    real(kind=wp) :: quota, quota_dia, quota_cocco, quota_phaeo
    ! [mmol C/mmol N]  Reciproque of 'quota' (NEW changed term)
    real(kind=wp) :: recipQuota, recipQuota_dia, recipQuota_cocco, recipQuota_phaeo
    ! [mg ChlA/mmol C] Quota between phytoplankton ChlA and C (NEW changed term)
    real(kind=wp) :: Chl2C, Chl2C_dia, Chl2C_cocco, Chl2C_phaeo
    ! [mg ChlA/mmol C] needed for photodamage (NEW changed term)
    real(kind=wp) :: Chl2C_plast, Chl2C_plast_dia, CHL2C_plast_cocco, CHL2C_plast_phaeo
    ! [mg ChlA/mmol N] Quota between phytoplankton ChlA and N (NEW changed term)
    real(kind=wp) :: Chl2N, Chl2N_dia, Chl2N_cocco, Chl2N_phaeo
    real(kind=wp) :: qSiC
    real(kind=wp) :: qSiN
    real(kind=wp) :: recipQZoo ! [mmol C/mmol N]  Quota between heterotrophic C and N
    real(kind=wp) :: recipQZoo2 ! [mmol C/mmol N]  Quota between second zoo  C and N
    real(kind=wp) :: recipQZoo3 ! Zoo3 [mmol C/mmol N] Quota between third zoo C and N
    !!! Grazing detritus Quotas for converting
    real(kind=wp) :: recipDet ! [mmol C/mmol N]  Quota between second zoo  C and N
    real(kind=wp) :: recipDet2 ! [mmol C/mmol N]  Quota between second zoo  C and N

    !!------------------------------------------------------------------------------
    !! *** For limiter function ***
    real(kind=wp) :: qlimitFac, qlimitFacTmp ! Factor that regulates photosynthesis
    real(kind=wp), external :: recom_limiter ! Function calculating qlimitFac
    real(kind=wp) :: FeLimitFac ! [Mumol/m3] Half sat constant for iron
    ! [1/day]    Maximum rate of C-specific photosynthesis
    real(kind=wp) :: pMax, pMax_dia, pMax_cocco, pMax_phaeo
    !!------------------------------------------------------------------------------
    !! *** Light ***
    real(kind=wp) :: kappar ! [1/m]  Light attenuation coefficient modified by chla
    real(kind=wp) :: kappastar ! []
    ! []     light attenuation * deltaZ at lower and upper control volume border
    real(kind=wp) :: kdzUpper, kdzLower
    ! [mg/m3]     chl  at lower and upper control volume border
    real(kind=wp) :: chl_upper, chl_lower
    real(kind=wp) :: Chlave ! [mg/m3]     vertical average chl between two nodes
    ! [?]    light at upper and lower border of control volume
    real(kind=wp) :: Upperlight, Lowerlight
    real(kind=wp) :: PARave ! [?]    Average light in the control volumes
    !!------------------------------------------------------------------------------
    !! *** Photosynthesis ***
    ! [1/day] C-specific rate of photosynthesis
    real(kind=wp) :: Cphot, Cphot_dia, Cphot_cocco, Cphot_phaeo
    !!------------------------------------------------------------------------------
    !! *** Assimilation ***
    real(kind=wp) :: V_cm ! scaling factor for temperature dependent maximum of C-specific N-uptake
    ! Factor that regulates N-assimilation. Calc from function recom_limiter
    real(kind=wp) :: limitFacN, limitFacN_dia, limitFacN_cocco, limitFacN_phaeo
    real(kind=wp) :: limitFacSi
    ! [mmol N/(mmol C * day)] C specific N utilization rate
    real(kind=wp) :: N_assim, N_assim_dia, N_assim_Cocco, N_assim_phaeo
    real(kind=wp) :: Si_assim
    !!------------------------------------------------------------------------------
    !! *** Chlorophyll ***
    ! [mg CHL/ mmol N] CHL a synthesis regulation term
    real(kind=wp) :: ChlSynth, ChlSynth_dia, ChlSynth_cocco, ChlSynth_phaeo
    ! [1/day] Phytoplankton respiration rate
    real(kind=wp) :: phyRespRate, phyRespRate_dia, phyRespRate_cocco, phyRespRate_phaeo
    ! coefficient for damage to the photosynthetic apparatus
    real(kind=wp) :: KOchl, KOchl_dia, KOchl_cocco, KOchl_phaeo
    !!------------------------------------------------------------------------------
    !! *** Vertical only Decomposition of phytoplankton growth components ***
    real(kind=wp), allocatable, dimension(:) :: VTTemp_diatoms, VTTemp_phyto, VTTemp_cocco, &
            VTTemp_phaeo ! Vertical 1D  temperature effect on phytoplankton photosynthesis
    ! CO2 effect
    real(kind=wp), allocatable, dimension(:) :: VTPhyCO2, VTDiaCO2, VTCoccoCO2, VTPhaeoCO2
    real(kind=wp), allocatable, dimension(:) :: VTqlimitFac_phyto, VTqlimitFac_diatoms, &
            VTqlimitFac_cocco, VTqlimitFac_phaeo ! nutrient effect
    real(kind=wp), allocatable, dimension(:) :: VTCphotLigLim_phyto, VTCphotLigLim_diatoms, &
            VTCphotLigLim_cocco, VTCphotLigLim_phaeo ! light limitation
    real(kind=wp), allocatable, dimension(:) :: VTCphot_phyto, VTCphot_diatoms, VTCphot_cocco, &
            VTCphot_phaeo
    real(kind=wp), allocatable, dimension(:) :: VTSi_assimDia

    !!------------------------------------------------------------------------------
    !! *** Iron chemistry ***
    real(kind=wp), external :: iron_chemistry, iron_chemistry_2ligands
    real(kind=wp) :: logK1, logK2, Klig1, Klig2
    !!------------------------------------------------------------------------------
    !! *** Zooplankton ***
    real(kind=wp) :: DiaNsq
    real(kind=wp) :: varpzdia, fDiaN ! Part of Diatoms available for food
    real(kind=wp) :: PhyNsq
    real(kind=wp) :: varpzPhy, fPhyN ! Part of Nano available for food
    real(kind=wp) :: CoccoNsq
    real(kind=wp) :: varpzCocco, fCoccoN
    real(kind=wp) :: PhaeoNsq
    real(kind=wp) :: varpzPhaeo, fPhaeoN
    real(kind=wp) :: MicZooNsq ! NEW 3Zoo
    real(kind=wp) :: varpzMicZoo, fMicZooN ! NEW 3Zoo Part of microzooplankton available for food
    real(kind=wp) :: food, foodsq ! [(mmol N)2/m6]
    ! [mmol N / (m3 * day)] (NEW changed term)
    real(kind=wp) :: grazingFlux_phy, grazingFlux_Dia, grazingFlux_Cocco, grazingFlux_Phaeo
    real(kind=wp) :: grazingFlux_miczoo ! NEW 3Zoo
    real(kind=wp) :: grazingFlux
    real(kind=wp) :: grazEff ! NEW 3Zoo
    real(kind=wp) :: HetRespFlux ! Zooplankton respiration
    real(kind=wp) :: HetLossFlux ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
    !!------------------------------------------------------------------------------
    !! *** Second Zooplankton  ***
    real(kind=wp) :: DiaNsq2, PhyNsq2, CoccoNsq2, PhaeoNsq2, HetNsq ! NEW (changed term)
    real(kind=wp) :: varpzDia2, fDiaN2, varpzPhy2, fPhyN2, varpzCocco2, fCoccoN2, varpzPhaeo2, &
            fPhaeoN2, varpzHet, fHetN ! Part of Diatoms available for food
    real(kind=wp) :: MicZooNsq2 ! NEW Zoo3
    real(kind=wp) :: varpzMicZoo2, fMicZooN2 ! NEW Zoo3
    real(kind=wp) :: food2, foodsq2 ! [(mmol N)2/m6]
    real(kind=wp) :: grazingFlux_phy2, grazingFlux_Dia2, grazingFlux_Cocco2, grazingFlux_Phaeo2, &
            grazingFlux_het2 ! [mmol N / (m3 * day)  (NEW changed term)
    real(kind=wp) :: grazingFlux_miczoo2 ! NEW Zoo3
    real(kind=wp) :: grazingFlux2
    real(kind=wp) :: Zoo2RespFlux ! Zooplankton respiration
    real(kind=wp) :: Zoo2LossFlux ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
    real(kind=wp) :: Zoo2fecalloss_n ! [(mmol N)/(m3*day)] Second zoo fecal pellet
    real(kind=wp) :: Zoo2fecalloss_c ! [(mmol N)/(m3*day)] Second zoo fecal pellet
    real(kind=wp) :: Mesfecalloss_n ! NEW Zoo3
    real(kind=wp) :: Mesfecalloss_c ! NEW Zoo3
    real(kind=wp) :: recip_res_zoo22
    !!------------------------------------------------------------------------------
    !! *** Grazing Detritus  ***
    real(kind=wp) :: DetNsq, DetZ2Nsq, DetNsq2, DetZ2Nsq2
    ! Part of Diatoms available for food
    real(kind=wp) :: varpzDet, varpzDetZ2, varpzDet2, varpzDetZ22
    real(kind=wp) :: fDetN, fDetZ2N, fDetN2, fDetZ2N2
    real(kind=wp) :: grazingFlux_Det, grazingFlux_DetZ2 ! [mmol N / (m3 * day)]
    real(kind=wp) :: grazingFlux_Det2, grazingFlux_DetZ22 ! [mmol N / (m3 * day)]
    !!------------------------------------------------------------------------------
    !! *** Third zooplankton  ***       ! NEW 3Zoo
    real(kind=wp) :: DiaNsq3
    real(kind=wp) :: varpzDia3, fDiaN3 ! Part of diatoms available for food
    real(kind=wp) :: loss_hetfd
    real(kind=wp) :: PhyNsq3
    real(kind=wp) :: varpzPhy3, fPhyN3 ! Part of small phytoplankton available for food
    real(kind=wp) :: CoccoNsq3
    real(kind=wp) :: varpzCocco3, fCoccoN3 ! Part of coccolithophores available for food
    real(kind=wp) :: PhaeoNsq3
    real(kind=wp) :: varpzPhaeo3, fPhaeoN3 ! Part of phaeocystis available for food
    real(kind=wp) :: food3, foodsq3 ! [(mmol N)2/m6]
    ! [mmol N / (m3 * day)]
    real(kind=wp) :: grazingFlux_phy3, grazingFlux_Dia3, grazingFlux_Cocco3, grazingFlux_Phaeo3
    real(kind=wp) :: grazingFlux3
    real(kind=wp) :: MicZooRespFlux ! Zooplankton respiration
    real(kind=wp) :: MicZooLossFlux ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
    !!------------------------------------------------------------------------------
    !! *** Aggregation  ***
    real(kind=wp) :: AggregationRate ! [1/day] AggregationRate (of nitrogen)
    !!------------------------------------------------------------------------------
    !! *** Calcification  ***
    ! NEW (before it was defined as a fixed value, but now dependent on cocco and T)
    real(kind=wp) :: calc_prod_ratio_cocco
    real(kind=wp) :: calcification
    real(kind=wp) :: calc_loss_agg
    real(kind=wp) :: calc_loss_gra
    real(kind=wp) :: calc_diss
    real(kind=wp) :: calc_diss_ben ! NEW DISS
    real(kind=wp) :: calc_loss_gra2 ! zoo2 detritus
    real(kind=wp) :: calc_diss2 ! zoo2 detritus
    real(kind=wp) :: calc_loss_gra3 ! NEW Zoo3 detritus
    real(kind=wp) :: Ca ! NEW DISS (calcium ion concentration)
    real(kind=wp) :: CO3_sat ! NEW DISS (saturated CO3 concentration, calculated from kspc and Ca)
    !!------------------------------------------------------------------------------
    !! *** Diagnostics  ***
    real(kind=wp) :: recipbiostep ! 1/number of steps per recom cycle
    real(kind=wp), allocatable, dimension(:, :) :: Diags3Dloc

    ! ==================================================================
    ! SMALL PHYTOPLANKTON (n suffix)
    ! ==================================================================
    real(kind=wp) :: locNPPn, locGPPn, locNNAn, locChldegn
    real(kind=wp), allocatable, dimension(:) :: vertNPPn, vertGPPn, vertNNAn, vertChldegn
    real(kind=wp), allocatable, dimension(:) :: vertrespn
    real(kind=wp), allocatable, dimension(:) :: vertdocexn
    real(kind=wp), allocatable, dimension(:) :: vertaggn

    ! ==================================================================
    ! DIATOMS (d suffix)
    ! ==================================================================
    real(kind=wp) :: locNPPd, locGPPd, locNNAd, locChldegd
    real(kind=wp), allocatable, dimension(:) :: vertNPPd, vertGPPd, vertNNAd, vertChldegd
    real(kind=wp), allocatable, dimension(:) :: vertrespd
    real(kind=wp), allocatable, dimension(:) :: vertdocexd
    real(kind=wp), allocatable, dimension(:) :: vertaggd

    ! ==================================================================
    ! COCCOLITHOPHORES (c suffix)
    ! ==================================================================
    real(kind=wp) :: locNPPc, locGPPc, locNNAc, locChldegc
    real(kind=wp), allocatable, dimension(:) :: vertNPPc, vertGPPc, vertNNAc, vertChldegc
    real(kind=wp), allocatable, dimension(:) :: vertrespc
    real(kind=wp), allocatable, dimension(:) :: vertdocexc
    real(kind=wp), allocatable, dimension(:) :: vertaggc
    real(kind=wp), allocatable, dimension(:) :: vertcalcdiss, vertcalcif

    ! ==================================================================
    ! PHAEOCYSTIS (p suffix)
    ! ==================================================================
    real(kind=wp) :: locNPPp, locGPPp, locNNAp, locChldegp
    real(kind=wp), allocatable, dimension(:) :: vertNPPp, vertGPPp, vertNNAp, vertChldegp
    real(kind=wp), allocatable, dimension(:) :: vertrespp
    real(kind=wp), allocatable, dimension(:) :: vertdocexp
    real(kind=wp), allocatable, dimension(:) :: vertaggp

    ! ==================================================================
    ! MICROZOOPLANKTON
    ! ==================================================================
    real(kind=wp) :: locgrazmicro_tot, locgrazmicro_n, locgrazmicro_d, locgrazmicro_c, &
            locgrazmicro_p
    real(kind=wp), allocatable, dimension(:) :: vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d&
            ,&
            & vertgrazmicro_c, vertgrazmicro_p
    real(kind=wp), allocatable, dimension(:) :: vertrespmicro

    ! ==================================================================
    ! MESOZOOPLANKTON
    ! ==================================================================
    real(kind=wp) :: locgrazmeso_tot, locgrazmeso_n, locgrazmeso_d, locgrazmeso_c, locgrazmeso_p
    real(kind=wp) :: locgrazmeso_det, locgrazmeso_mic, locgrazmeso_det2
    real(kind=wp), allocatable, dimension(:) :: vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d, &
            vertgrazmeso_c, vertgrazmeso_p
    real(kind=wp), allocatable, dimension(:) :: vertgrazmeso_det, vertgrazmeso_mic, &
            vertgrazmeso_det2
    real(kind=wp), allocatable, dimension(:) :: vertrespmeso

    ! ==================================================================
    ! MACROZOOPLANKTON
    ! ==================================================================
    real(kind=wp) :: locgrazmacro_tot, locgrazmacro_n, locgrazmacro_d, locgrazmacro_c, &
            locgrazmacro_p
    real(kind=wp) :: locgrazmacro_mes, locgrazmacro_det, locgrazmacro_mic, locgrazmacro_det2
    real(kind=wp), allocatable, dimension(:) :: vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d&
            ,&
            & vertgrazmacro_c, vertgrazmacro_p
    real(kind=wp), allocatable, dimension(:) :: vertgrazmacro_mes, vertgrazmacro_det, &
            vertgrazmacro_mic, vertgrazmacro_det2
    real(kind=wp), allocatable, dimension(:) :: vertrespmacro

    !!------------------------------------------------------------------------------
    !! *** Benthos  ***
    ! [1/day] Decay rate of detritus in the benthic layer
    real(kind=wp), allocatable, dimension(:) :: decayBenthos

    ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
    real(kind=wp), allocatable, dimension(:) :: wFluxDet

    ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
    real(kind=wp), allocatable, dimension(:) :: wFluxPhy

    ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms
    real(kind=wp), allocatable, dimension(:) :: wFluxDia

    ! NEW [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of coccos
    real(kind=wp), allocatable, dimension(:) :: wFluxCocco

    ! NEW [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of Phaeocystis
    real(kind=wp), allocatable, dimension(:) :: wFluxPhaeo
    real(kind=wp) :: Vben_det ! [m/day] speed of sinking into benthos from water column
    real(kind=wp) :: Vben_det_seczoo !second zooplankton sinking benthos
    real(kind=wp) :: Vben_phy
    real(kind=wp) :: Vben_dia
    real(kind=wp) :: Vben_coc
    real(kind=wp) :: Vben_pha ! Phaeocystis
    real(kind=wp) :: Ironflux ! [umol Fe/(m2*day)] Flux of Fe from sediment to water
    !_______________________________________________________________________________
    ! Arrays added for RECOM implementation:
    !!---- PAR
    !real(kind=wp),allocatable,dimension(:)     :: PAR

    ! --> multiplication factor for surface boundary condition in
    !     bc_surface for river and erosion
    !     river on/off -->=1.0/0.0
    !     erosion on/off -->=1.0/0.0

    real(kind=wp) :: is_riverinput
    real(kind=wp) :: is_erosioninput

    real(kind=wp) :: is_3zoo2det
    real(kind=wp) :: is_coccos

    type :: recom_tracer_ids

        ! --- Small Phytoplankton
        integer :: phytoplankton_nitrogen
        integer :: phytoplankton_carbon
        integer :: phytoplankton_chlorophyll

        ! --- Detritus (Non-living organic matter) ---
        integer :: detrital_nitrogen
        integer :: detrital_carbon

        ! --- Mesozooplankton (Heterotrophs) ---
        integer :: heterotroph_nitrogen
        integer :: heterotroph_carbon

        ! --- Dissolved Organic Matter ---
        integer :: dissolved_organic_nitrogen
        integer :: dissolved_organic_carbon

        ! --- Diatoms ---
        integer :: diatom_nitrogen
        integer :: diatom_carbon
        integer :: diatom_chlorophyll
        integer :: diatom_silica

        ! --- Detrital Silica ---
        integer :: detrital_silica

        ! --- --- Iron (micronutrient) ---
        integer :: iron

        ! --- Calcium Carbonate (Calcite) ---
        integer :: phytoplankton_calcite
        integer :: detrital_calcite

        ! --- 3zoo2det extra tracers ---
        integer :: macrozooplankton_nitrogen
        integer :: macrozooplankton_carbon
        integer :: macrozooplankton_detrital_nitrogen
        integer :: macrozooplankton_detrital_carbon
        integer :: macrozooplankton_detrital_silica
        integer :: macrozooplankton_detrital_calcite
        integer :: microzooplankton_nitrogen
        integer :: microzooplankton_carbon

        ! --- coccos extra tracers ---
        integer :: coccolithophore_nitrogen
        integer :: coccolithophore_carbon
        integer :: coccolithophore_chlorophyll
        integer :: phaeocystis_nitrogen
        integer :: phaeocystis_carbon
        integer :: phaeocystis_chlorophyll

    end type recom_tracer_ids

    type(recom_tracer_ids) :: tracer_ids

end module REcoM_declarations

!
!===============================================================================
!

module recom_config
    use recom_declarations, only: wp
    implicit none
    public

    private :: wp

    save

    !! *** General constants ***

    ! *******************
    ! CASE 2phy 1zoo 1det
    ! *******************
    integer :: idin = 1, idic = 2, ialk = 3, iphyn = 4, iphyc = 5, &
            ipchl = 6, idetn = 7, idetc = 8, ihetn = 9, &
            ihetc = 10, idon = 11, idoc = 12, idian = 13, &
            idiac = 14, idchl = 15, idiasi = 16, idetsi = 17, &
            isi = 18, ife = 19, iphycal = 20, idetcal = 21, &
            ioxy = 22

    integer :: izoo2n = 23, izoo2c = 24, idetz2n = 25, &
            idetz2c = 26, idetz2si = 27, idetz2calc = 28

    ! Microzooplankton (third zooplankton group)
    integer :: imiczoon = 0 ! Microzooplankton Nitrogen (set below)
    integer :: imiczooc = 0 ! Microzooplankton Carbon (set below)

    ! ---------------------------------------------------------------------------
    ! PHYTOPLANKTON GROUPS (coccos configuration)
    ! ---------------------------------------------------------------------------
    ! Coccolithophores and Phaeocystis when enable_coccos = .true.

    integer :: icocn = 0 ! Coccolithophore Nitrogen (set below)
    integer :: icocc = 0 ! Coccolithophore Carbon (set below)
    integer :: icchl = 0 ! Coccolithophore Chlorophyll (set below)

    integer :: iphan = 0 ! Phaeocystis Nitrogen (set below)
    integer :: iphac = 0 ! Phaeocystis Carbon (set below)
    integer :: iphachl = 0 ! Phaeocystis Chlorophyll (set below)

    !=============================================================================

    integer :: ivphy = 1, ivdia = 2, ivdet = 3, ivdetsc = 4, ivcoc = 5, ivpha = 6

    !=============================================================================

    integer, dimension(8) :: recom_remin_tracer_id = [1001, 1002, 1003, 1018, 1019, 1022, 1302, &
            1402]

    ! OG
    ! Todo:  Make recom_sinking_tracer_id case sensitive
    integer, dimension(32) :: recom_sinking_tracer_id = [1007, 1008, 1017, 1021, 1004, 1005, 1020, &
            1006, 1013, 1014, 1016, 1015, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, &
            1034, 1308, 1321, 1305, 1320, 1314, 1408, 1421, 1405, 1420, 1414]

    integer, dimension(8) :: recom_det_tracer_id = [1007, 1008, 1017, 1021, 1308, 1321, 1408, 1421]
    integer, dimension(8) :: recom_phy_tracer_id = [1004, 1005, 1020, 1305, 1320, 1405, 1420, 1006]
    integer, dimension(6) :: recom_dia_tracer_id = [1013, 1014, 1314, 1414, 1016, 1015]

    ! Configuration-dependent tracer arrays (allocated during initialization)
    integer, dimension(3) :: recom_cocco_tracer_id
    integer, dimension(3) :: recom_phaeo_tracer_id
    integer, dimension(4) :: recom_det2_tracer_id

    !=============================================================================

    real(kind=wp) :: zero = 0.d0
    integer :: one = 1
    real(kind=wp) :: tiny = 2.23D-16
    real(kind=wp) :: tiny_chl = 0.00001
    real(kind=wp) :: SecondsPerDay = 86400.d0 ! [s/day]
    real(kind=wp) :: Pa2atm = 101325.d0 ! [Pa/atm]
    real(kind=wp) :: redO2C = 1.453 ! O2:C ratio Anderson and Sarmiento, 1994

    !! *** REcoM setup ***
    logical :: enable_3zoo2det = .false. ! Control extended zooplankton variables
    logical :: enable_coccos = .false. ! Control coccolithophore variables
    namelist /parecomsetup/ enable_3zoo2det, enable_coccos

    !! *** General configuration ***

    logical :: use_REcoM = .true.
    logical :: REcoM_restart = .false.

    ! NEW increased the number from 28 to 34 (added coccos and respiration)
    ! NEW 3Zoo changed from 31 to 33 ! added phaeocystis: changed from 33 to 36
    integer :: bgc_num = 36
    integer :: bgc_base_num = 22 ! standard tracers
    integer :: diags3d_num = 31 ! Number of diagnostic 3d tracers to be saved
    ! Sinking velocity, constant through the water column and positive downwards
    real(kind=wp) :: VDet = 20.d0
    real(kind=wp) :: VDet_zoo2 = 200.d0 ! Sinking velocity, constant through the water column
    !!! If the number of sinking velocities are different from 3, code needs to be changed !!!
    real(kind=wp) :: VPhy = 0.d0
    real(kind=wp) :: VDia = 0.d0
    real(kind=wp) :: VCocco = 0.d0 ! NEW
    real(kind=wp) :: VPhaeo = 0.d0 ! Phaeocystis
    logical :: allow_var_sinking = .true.
    integer :: biostep = 1 ! Number of times biology should be stepped forward for each time step

    ! Decides what routine should be used to calculate limiters in sms
    logical :: REcoM_Geider_limiter = .false.

    ! Decides if grazing should have preference for phyN or DiaN
    logical :: REcoM_Grazing_Variable_Preference = .true.
    logical :: Grazing_detritus = .false. ! Decides grazing on detritus
    logical :: het_resp_noredfield = .true. ! Decides respiratation of copepods
    logical :: diatom_mucus = .true. ! Effect of nutrient limitation on the aggregation

    ! NEW O2remin Add option for O2 dependency of organic matter remineralization
    logical :: O2dep_remin = .true.
    logical :: use_ballasting = .true. ! NEW BALL
    logical :: use_density_scaling = .true. ! NEW BALL
    logical :: use_viscosity_scaling = .true. ! NEW BALL

    ! NEW DISS Use mocsy calcite omega to compute calcite dissolution
    logical :: OmegaC_diss = .true.
    logical :: CO2lim = .true. ! NEW Use CO2 dependence of growth and calcification
    !Logical                :: inter_CT_CL           = .true.    ! NEW inter use interaction between
    !CO2 and both, temperature and light
    logical :: Diags = .true. !!!!!!!!!!!!!!!!!!!!!!Change in recom.F90 Diagnostics -> Diags
    logical :: constant_CO2 = .true.
    logical :: UseFeDust = .true. ! Turns dust input of iron off when set to.false.
    logical :: UseDustClim = .true.
    ! Use Albani dustclim field (If it is false Mahowald will be used)
    logical :: UseDustClimAlbani = .true.
    logical :: use_photodamage = .false. ! use Alvarez et al (2018) for chlorophyll degradation
    !MB More stable computation of zooplankton respiration fluxes adding a small number to HetN
    logical :: HetRespFlux_plus = .true.
    character(100) :: REcoMDataPath = &
            '/albedo/work/projects/MarESys/ogurses/input/mesh_CORE2_finaltopo_mean/'
    logical :: restore_alkalinity = .true.
    logical :: useRivers = .false.
    logical :: useRivFe = .false. ! river input of Fe
    logical :: useErosion = .false.

    ! This one only activates rivers! And in principle denitrification, but denitrification is
    ! commented out. When set to true, external sources and sinks of
    ! nitrogen are activated (Riverine, aeolian and denitrification)
    logical :: NitrogenSS = .false.
    logical :: useAeolianN = .false. ! When set to true, aeolian nitrogen deposition is activated
    ! The first year of the actual physical forcing (e.g. JRA-55) used
    integer :: firstyearoffesomcycle = 1958
    integer :: lastyearoffesomcycle = 2023 ! Last year of the actual physical forcing used
    integer :: numofCO2cycles = 1 ! Number of cycles of the forcing planned
    integer :: currentCO2cycle = 1 ! Which CO2 cycle we are currently running
    logical :: DIC_PI = .true.
    integer :: Nmocsy = 1 ! Length of the vector that is passed to mocsy (always one for recom)
    logical :: recom_debug = .false.
    logical :: ciso = .false. !MB main switch to enable/disable carbon isotopes (13|14C)
    integer :: benthos_num = 4 !MB number of sediment tracers = 8 if ciso = .true.
    logical :: use_MEDUSA = .false. ! main switch for sediment model
    integer :: sedflx_num = 0 ! number of sedimentary fluxs from MEDUSA, = 7 if ciso
    logical :: add_loopback = .false.
    real(kind=wp) :: lb_tscale = 1.d0 ! time scale to balance the burial loss
    ! number of stored sinking fluxes from the bottom layer, = 6 if C13 and = 8 if C14

    integer :: bottflx_num = 4
    logical :: use_atbox = .false. ! switch for atmospheric box model for CO2

    logical :: fe_2ligands = .false. ! consider Fe-ligand binding with two ligands
    ! use Fe-ligand parameterisation dependent on DOC and pH (Ye2020)
    logical :: fe_compl_nica = .false.

    namelist /pavariables/ use_REcoM, REcoM_restart, &
            bgc_num, diags3d_num, bgc_base_num, &
            VDet, VDet_zoo2, &
            VPhy, VDia, VCocco, &
            VPhaeo, &
            allow_var_sinking, biostep, REcoM_Geider_limiter, &
            REcoM_Grazing_Variable_Preference, &
            Grazing_detritus, &
            het_resp_noredfield, &
            diatom_mucus, &
            O2dep_remin, use_ballasting, use_density_scaling, & ! O2remin, NEW BALL
    ! BALL, DISS added OmegaC_diss, added CO2lim
            use_viscosity_scaling, OmegaC_diss, CO2lim, &
            Diags, constant_CO2, &
            UseFeDust, UseDustClim, UseDustClimAlbani, &
            use_photodamage, HetRespFlux_plus, REcoMDataPath, &
            restore_alkalinity, useRivers, useRivFe, &
            useErosion, NitrogenSS, useAeolianN, &
            firstyearoffesomcycle, lastyearoffesomcycle, numofCO2cycles, &
            currentCO2cycle, DIC_PI, Nmocsy, &
            recom_debug, ciso, benthos_num, &
            use_MEDUSA, sedflx_num, bottflx_num, &
            add_loopback, lb_tscale, use_atbox, &
            fe_2ligands, fe_compl_nica

    !!------------------------------------------------------------------------------
    !! *** Sinking ***
    real(kind=wp) :: Vdet_a = 0.0288 ! [1/day]
    real(kind=wp) :: Vcalc = 0.0216 ! [1/day] depth dependence of calc_diss

    namelist /pasinking/ Vdet_a, Vcalc
    !!------------------------------------------------------------------------------
    !! *** Initialization ***
    real(kind=wp) :: cPhyN = 0.2d0
    real(kind=wp) :: cHetN = 0.2d0
    real(kind=wp) :: cZoo2N = 0.2d0

    namelist /painitialization_N/ cPhyN, cHetN, cZoo2N
    !!------------------------------------------------------------------------------
    !! *** Temperature and Arrhenius functions ***
    real(kind=wp) :: recom_Tref = 288.15d0 ! [K]
    real(kind=wp) :: C2K = 273.15d0 !     Conversion from degrees C to K
    real(kind=wp) :: Ae = 4500.d0 ! [K] Slope of the linear part of the Arrhenius function

    !! *** Temperature variables for Blanchard function ***
    real(kind=wp) :: Tmax_phaeo = 16d0 ! [degC] For Blanchard temp fxn: maximum temperature
    real(kind=wp) :: Topt_phaeo = 7.5272d0 ! [degC] For Blanchard temp fxn: optimum temperature
    real(kind=wp) :: uopt_phaeo = 0.7328d0 ! [1/day] For Blanchard function: optimum growth date
    real(kind=wp) :: beta_phaeo = 0.7829d0 ! [unitless] For Blanchard function

    ! NEW MODIFIED parameters
    real(kind=wp) :: ord_d = -0.2216d0 ! parameters for diatom temperature function
    real(kind=wp) :: expon_d = 0.0406d0 ! diatom exponent
    real(kind=wp) :: ord_phy = -1.2154d0 ! small phyto ordonnee
    real(kind=wp) :: expon_phy = 0.0599d0 ! small phyto exponent
    real(kind=wp) :: ord_cocco = -0.2310d0 ! coccolith ordonnee
    real(kind=wp) :: expon_cocco = 0.0327d0 ! small phyto ordonnee
    real(kind=wp) :: ord_phaeo = -0.2310d0 ! phaeocystis ordonnee
    real(kind=wp) :: expon_phaeo = 0.0327d0 ! phaeocystis ordonnee

    real(kind=wp) :: reminSi = 0.02d0
    ! NEW O2remin mmol m-3; Table 1 in Cram 2018 cites DeVries & Weber 2017
    ! for a range of 0-30 mmol m-3
    real(kind=wp) :: k_o2_remin = 15.d0
    namelist /paArrhenius/ recom_Tref, C2K, Ae, Tmax_phaeo, Topt_phaeo, uopt_phaeo, beta_phaeo, &
            ord_d, expon_d, ord_phy, expon_phy, ord_cocco, expon_cocco, ord_phaeo, expon_phaeo, &
            reminSi, k_o2_remin

    !!------------------------------------------------------------------------------
    !! *** For limiter function ***
    real(kind=wp) :: NMinSlope = 50.d0
    real(kind=wp) :: SiMinSlope = 1000.d0
    real(kind=wp) :: NCmin = 0.04d0
    real(kind=wp) :: NCmin_d = 0.04d0
    real(kind=wp) :: NCmin_c = 0.04d0 ! NEW
    real(kind=wp) :: NCmin_p = 0.04d0 ! Phaeocystis
    real(kind=wp) :: SiCmin = 0.04d0
    real(kind=wp) :: k_Fe = 0.04d0
    real(kind=wp) :: k_Fe_d = 0.12d0
    real(kind=wp) :: k_Fe_c = 0.04 ! NEW
    real(kind=wp) :: k_Fe_p = 0.09 ! Phaeocystis (to be tuned)
    real(kind=wp) :: k_si = 4.d0
    real(kind=wp) :: P_cm = 3.0d0 ! [1/day]   the rate of C-specific photosynthesis
    real(kind=wp) :: P_cm_d = 3.5d0
    real(kind=wp) :: P_cm_c = 3.3d0 ! NEW
    real(kind=wp) :: P_cm_p = 3.4d0 ! NEW for Phaeocystis ( to be tuned)
    namelist /palimiter_function/ NMinSlope, SiMinSlope, NCmin, NCmin_d, NCmin_c, NCmin_p, SiCmin, &
            k_Fe, k_Fe_d, k_Fe_c, k_Fe_p, k_si, P_cm, P_cm_d, P_cm_c, P_cm_p

    !!------------------------------------------------------------------------------
    !! *** For light calculations ***
    real(kind=wp) :: k_w = 0.04d0 ! [1/m]              Light attenuation coefficient
    ! [1/m * 1/(mg Chl)] Chlorophyll specific attenuation coefficients
    real(kind=wp) :: a_chl = 0.03d0
    namelist /palight_calculations/ k_w, a_chl
    !!------------------------------------------------------------------------------
    !! *** Photosynthesis ***
    real(kind=wp) :: alfa = 0.14d0 ! [(mmol C*m2)/(mg Chl*W*day)]
    real(kind=wp) :: alfa_d = 0.19d0 ! An initial slope of the P-I curve
    real(kind=wp) :: alfa_c = 0.10d0 ! NEW
    real(kind=wp) :: alfa_p = 0.10d0 ! Phaeocystis (to be tuned)
    real(kind=wp) :: parFrac = 0.43d0
    namelist /paphotosynthesis/ alfa, alfa_d, alfa_c, alfa_p, parFrac
    !!------------------------------------------------------------------------------
    !! *** Assimilation ***
    ! scaling factor for temperature dependent maximum of C-specific N-uptake
    real(kind=wp) :: V_cm_fact = 0.7d0
    real(kind=wp) :: V_cm_fact_d = 0.7d0
    real(kind=wp) :: V_cm_fact_c = 0.7d0 ! NEW
    real(kind=wp) :: V_cm_fact_p = 0.7d0 ! Phaeocystis
    real(kind=wp) :: NMaxSlope = 1000.d0 ! Max slope for limiting function
    real(kind=wp) :: SiMaxSlope = 1000.d0
    real(kind=wp) :: NCmax = 0.2d0 ! [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    real(kind=wp) :: NCmax_d = 0.2d0
    real(kind=wp) :: NCmax_c = 0.15d0 ! NEW
    real(kind=wp) :: NCmax_p = 0.1d0 ! Phaeocystis (to be tuned)
    real(kind=wp) :: SiCmax = 0.8d0
    real(kind=wp) :: NCuptakeRatio = 0.2d0 ! [mmol N/mmol C] Maximum uptake ratio of N:C
    real(kind=wp) :: NCUptakeRatio_d = 0.2d0
    real(kind=wp) :: NCUptakeRatio_c = 0.2d0 ! NEW
    real(kind=wp) :: NCUptakeRatio_p = 0.2d0 ! Phaeocystis
    real(kind=wp) :: SiCUptakeRatio = 0.2d0
    real(kind=wp) :: k_din = 0.55d0 ! [mmol N/m3] Half-saturation constant for nitrate uptake
    real(kind=wp) :: k_din_d = 1.0d0
    real(kind=wp) :: k_din_c = 0.55d0 ! NEW
    real(kind=wp) :: k_din_p = 0.55d0 ! Phaeocystis (to be tuned)
    real(kind=wp) :: Chl2N_max = 3.15d0 ! [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    real(kind=wp) :: Chl2N_max_d = 4.2d0
    real(kind=wp) :: Chl2N_max_c = 3.5d0 ! NEW
    real(kind=wp) :: Chl2N_max_p = 3.5d0 ! Phaeocystis (to be tuned (?))
    real(kind=wp) :: res_phy = 0.01d0 ! [1/day] Maintenance respiration rate constant
    real(kind=wp) :: res_phy_d = 0.01d0
    real(kind=wp) :: res_phy_c = 0.0075d0 ! NEW
    real(kind=wp) :: res_phy_p = 0.008d0 ! Phaeocystis (to be tuned (?))
    real(kind=wp) :: biosynth = 2.33d0 ! [mmol C/mmol N] Cost of biosynthesis
    real(kind=wp) :: biosynthSi = 0.d0
    namelist /paassimilation/ V_cm_fact, V_cm_fact_d, V_cm_fact_c, V_cm_fact_p, NMaxSlope, &
            SiMaxSlope, NCmax, NCmax_d, NCmax_c, NCmax_p, SiCmax, &
            NCuptakeRatio, NCUptakeRatio_d, NCUptakeRatio_c, NCUptakeRatio_p, SiCUptakeRatio, &
            k_din &
            , k_din_d, k_din_c, k_din_p, &
            Chl2N_max, Chl2N_max_d, Chl2N_max_c, Chl2N_max_p, res_phy, res_phy_d, res_phy_c, &
            res_phy_p, biosynth, biosynthSi
    !!------------------------------------------------------------------------------
    !! *** Iron chemistry ***
    real(kind=wp) :: totalligand = 1.d0 ! [mumol/m3] order 1. Total free ligand
    ! [m3/mumol] order 100. Ligand-free iron stability constant
    real(kind=wp) :: ligandStabConst = 100.d0
    namelist /pairon_chem/ totalligand, ligandStabConst
    !!------------------------------------------------------------------------------
    !! *** Zooplankton ***
    real(kind=wp) :: graz_max = 2.4d0 ! [mmol N/(m3 * day)] Maximum grazing loss parameter
    real(kind=wp) :: epsilonr = 0.35d0 ! [(mmol N)2 /m6] Half saturation constant for grazing loss
    ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)
    real(kind=wp) :: res_het = 0.01d0
    real(kind=wp) :: Redfield = 6.625 ! [mmol C/mmol N] Redfield ratio of C:N = 106:16
    ! [1/day] Temperature dependent N degradation of extracellular organic N (EON)
    real(kind=wp) :: loss_het = 0.05d0
    real(kind=wp) :: pzDia = 0.5d0 ! Maximum diatom preference
    real(kind=wp) :: sDiaNsq = 0.d0
    real(kind=wp) :: pzPhy = 1.0d0 ! Maximum nano-phytoplankton preference
    real(kind=wp) :: sPhyNsq = 0.d0
    real(kind=wp) :: pzCocco = 0.5d0 ! NEW (value is just a guess)
    real(kind=wp) :: sCoccoNsq = 0.d0 ! NEW
    real(kind=wp) :: pzPhaeo = 1.0d0 ! Phaeocystis (to be tuned)
    real(kind=wp) :: sPhaeoNsq = 0.d0 ! Phaeocystis
    real(kind=wp) :: pzMicZoo = 1.0d0 ! NEW 3Zoo Maximum nano-phytoplankton preference
    real(kind=wp) :: sMicZooNsq = 0.d0 ! NEW 3Zoo

    ! for more stable computation of HetRespFlux (_plus).
    ! Value can be > tiny because HetRespFlux ~ hetC**2.
    real(kind=wp) :: tiny_het = 1.d-5
    namelist /pazooplankton/ graz_max, epsilonr, res_het, Redfield, loss_het, pzDia, sDiaNsq, &
            pzPhy, sPhyNsq, pzCocco, sCoccoNsq, pzPhaeo, sPhaeoNsq, pzMicZoo, sMicZooNsq, tiny_het
    !!-------------------------------------------------------------------------------
    !! *** SecondZooplankton (Macrozooplankton) ***
    real(kind=wp) :: graz_max2 = 0.1d0 ! [mmol N/(m3 * day)] Maximum grazing loss parameter
    real(kind=wp) :: epsilon2 = 0.0144d0 ! [(mmol N)2 /m6] Half saturation constant for grazing loss
    ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)
    real(kind=wp) :: res_zoo2 = 0.0107d0
    ! [1/day] Temperature dependent N degradation of extracellular organic N
    real(kind=wp) :: loss_zoo2 = 0.003d0
    real(kind=wp) :: fecal_rate_n = 0.13d0
    real(kind=wp) :: fecal_rate_c = 0.295d0
    real(kind=wp) :: fecal_rate_n_mes = 0.25d0 ! NEW 3Zoo
    real(kind=wp) :: fecal_rate_c_mes = 0.32d0 ! NEW 3Zoo
    real(kind=wp) :: pzDia2 = 1.d0 ! Maximum diatom preference
    real(kind=wp) :: sDiaNsq2 = 0.d0
    real(kind=wp) :: pzPhy2 = 0.5d0 ! Maximum diatom preference
    real(kind=wp) :: sPhyNsq2 = 0.d0
    real(kind=wp) :: pzCocco2 = 0.2d0 ! NEW (value is just a guess)
    real(kind=wp) :: sCoccoNsq2 = 0.d0 ! NEW
    real(kind=wp) :: pzPhaeo2 = 0.5d0 ! Phaeocystis (to be tuned)
    real(kind=wp) :: sPhaeoNsq2 = 0.d0 ! Phaeocystis
    real(kind=wp) :: pzHet = 0.8d0 ! Maximum diatom preference
    real(kind=wp) :: sHetNsq = 0.d0
    real(kind=wp) :: pzMicZoo2 = 0.8d0 ! NEW Zoo3 Maximum microzooplankton preference
    real(kind=wp) :: sMicZooNsq2 = 0.d0 ! NEW Zoo3
    real(kind=wp) :: t1_zoo2 = 28145.d0 ! Krill temp. function constant1
    real(kind=wp) :: t2_zoo2 = 272.5d0 ! Krill temp. function constant2
    real(kind=wp) :: t3_zoo2 = 105234.d0 ! Krill temp. function constant3
    real(kind=wp) :: t4_zoo2 = 274.15d0 ! Krill temp. function constant3
    namelist /pasecondzooplankton/ graz_max2, epsilon2, res_zoo2, &
            loss_zoo2, fecal_rate_n, fecal_rate_c, fecal_rate_n_mes, fecal_rate_c_mes, & ! NEW 3Zoo
            pzDia2, sDiaNsq2, pzPhy2, sPhyNsq2, pzCocco2, sCoccoNsq2, pzPhaeo2, sPhaeoNsq2, pzHet, &
            sHetNsq, pzMicZoo2, sMicZooNsq2, t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2
    !-------------------------------------------------------------------------------
    !! *** Third Zooplankton (Microzooplankton) ***
    ! NEW 3Zoo [mmol N/(m3 * day)] Maximum grazing loss parameter
    real(kind=wp) :: graz_max3 = 0.46d0
    ! NEW 3Zoo [(mmol N)2 /m6] Half saturation constant for grazing loss
    real(kind=wp) :: epsilon3 = 0.64d0
    ! NEW 3Zoo [1/day] Temperature dependent N degradation of extracellular organic N (EON)
    real(kind=wp) :: loss_miczoo = 0.01d0
    ! NEW 3Zoo [1/day] Respiration by heterotrophs and mortality (loss to detritus)
    real(kind=wp) :: res_miczoo = 0.01d0
    real(kind=wp) :: pzDia3 = 0.5d0 ! NEW 3Zoo Maximum diatom preference
    real(kind=wp) :: sDiaNsq3 = 0.d0 ! NEW 3Zoo
    real(kind=wp) :: pzPhy3 = 1.0d0 ! NEW 3Zoo Maximum nano-phytoplankton preference
    real(kind=wp) :: sPhyNsq3 = 0.d0 ! NEW 3Zoo
    ! NEW 3Zoo Maximum coccolithophore preference
    ! ATTENTION: This value needs to be tuned; I start with zero preference!
    real(kind=wp) :: pzCocco3 = 0.d0
    real(kind=wp) :: sCoccoNsq3 = 0.d0 ! NEW 3Zoo
    ! Phaeocystis 3Zoo Maximum phaeocystis preference (to be tuned (?))
    real(kind=wp) :: pzPhaeo3 = 1.0d0
    real(kind=wp) :: sPhaeoNsq3 = 0.d0 ! Phaeocystis 3Zoo
    namelist /pathirdzooplankton/ graz_max3, epsilon3, loss_miczoo, res_miczoo, pzDia3, sDiaNsq3, &
            pzPhy3, sPhyNsq3, pzCocco3, sCoccoNsq3, pzPhaeo3, sPhaeoNsq3

    !-------------------------------------------------------------------------------
    !! *** Detritus Grazing Params ***
    real(kind=wp) :: pzDet = 1.d0 ! Maximum small detritus prefence by first zooplankton
    real(kind=wp) :: sDetNsq = 0.d0
    real(kind=wp) :: pzDetZ2 = 1.d0 ! Maximum large detritus preference by first zooplankton
    real(kind=wp) :: sDetZ2Nsq = 0.d0
    real(kind=wp) :: pzDet2 = 1.d0 ! Maximum small detritus prefence by second zooplankton
    real(kind=wp) :: sDetNsq2 = 0.d0
    real(kind=wp) :: pzDetZ22 = 1.d0 ! Maximum large detritus preference by second zooplankton
    real(kind=wp) :: sDetZ2Nsq2 = 0.d0
    namelist /pagrazingdetritus/ pzDet, sDetNsq, pzDetZ2, sDetZ2Nsq, &
            pzDet2, sDetNsq2, pzDetZ22, sDetZ2Nsq2
    !!------------------------------------------------------------------------------
    !! *** Aggregation ***
    ! [m3/(mmol N * day)] Maximum aggregation loss parameter for DetN
    real(kind=wp) :: agg_PD = 0.165d0
    ! [m3/(mmol N * day)] Maximum aggregation loss parameter for PhyN and DiaN (plankton)
    real(kind=wp) :: agg_PP = 0.015d0
    namelist /paaggregation/ agg_PD, agg_PP
    !!------------------------------------------------------------------------------
    !! *** DIN ***
    ! [1/day] Temperature dependent N degradation of extracellular
    ! organic N (EON) (Remineralization of DON)
    real(kind=wp) :: rho_N = 0.11d0
    namelist /padin_rho_N/ rho_N
    !!------------------------------------------------------------------------------
    !! *** DIC ***
    ! [1/day] Temperature dependent C degradation of extracellular organic C (EOC)
    real(kind=wp) :: rho_C1 = 0.1d0
    namelist /padic_rho_C1/ rho_C1
    !!------------------------------------------------------------------------------
    !! *** Phytoplankton N ***
    real(kind=wp) :: lossN = 0.05d0 ! [1/day] Phytoplankton loss of organic N compounds
    real(kind=wp) :: lossN_d = 0.05d0
    real(kind=wp) :: lossN_c = 0.05d0
    real(kind=wp) :: lossN_p = 0.05d0 ! Phaeocystis
    namelist /paphytoplankton_N/ lossN, lossN_d, lossN_c, lossN_p
    !!------------------------------------------------------------------------------
    !! *** Phytoplankton C ***
    real(kind=wp) :: lossC = 0.10d0 ! [1/day] Phytoplankton loss of carbon
    real(kind=wp) :: lossC_d = 0.10d0
    real(kind=wp) :: lossC_c = 0.10d0
    real(kind=wp) :: lossC_p = 0.10d0 ! Phaeocystis
    namelist /paphytoplankton_C/ lossC, lossC_d, lossC_c, lossC_p
    !!------------------------------------------------------------------------------
    !! *** Phytoplankton ChlA ***
    real(kind=wp) :: deg_Chl = 0.25d0 ! [1/day]
    real(kind=wp) :: deg_Chl_d = 0.25d0
    real(kind=wp) :: deg_Chl_c = 0.20d0 ! (value is just a guess)
    real(kind=wp) :: deg_Chl_p = 0.25d0 ! Phaeocystis
    namelist /paphytoplankton_ChlA/ deg_Chl, deg_Chl_d, deg_Chl_c, deg_Chl_p
    !!------------------------------------------------------------------------------
    !! *** Detritus N ***
    ! 3Zoo [] Grazing efficiency (fraction of grazing flux into zooplankton pool)
    real(kind=wp) :: gfin = 0.3d0
    ! 3Zoo [] Grazing efficiency (fraction of grazing flux into second zooplankton pool)
    real(kind=wp) :: grazEff2 = 0.8d0
    ! 3Zoo [] Grazing efficiency (fraction of grazing flux into microzooplankton pool)
    real(kind=wp) :: grazEff3 = 0.8d0
    ! 3Zoo [1/day] Temperature dependent remineralisation rate of detritus
    real(kind=wp) :: reminN = 0.165d0
    namelist /padetritus_N/ gfin, grazEff2, grazEff3, reminN
    !!------------------------------------------------------------------------------
    !! *** Detritus C ***
    ! [1/day] Temperature dependent remineralisation rate of detritus
    real(kind=wp) :: reminC = 0.15d0
    ! [1/day] Temperature dependent C degradation of TEP-C
    real(kind=wp) :: rho_c2 = 0.1d0
    namelist /padetritus_C/ reminC, rho_c2
    !!------------------------------------------------------------------------------
    !! *** Heterotrophs ***
    real(kind=wp) :: lossN_z = 0.15d0
    real(kind=wp) :: lossC_z = 0.15d0
    namelist /paheterotrophs/ lossN_z, lossC_z
    !!------------------------------------------------------------------------------
    !! *** Second Zooplankton ***
    real(kind=wp) :: lossN_z2 = 0.02d0
    real(kind=wp) :: lossC_z2 = 0.02d0
    namelist /paseczooloss/ lossN_z2, lossC_z2
    !!-----------------------------------------------------------------------------
    !! *** Third Zooplankton ***
    real(kind=wp) :: lossN_z3 = 0.05d0 ! 3Zoo
    real(kind=wp) :: lossC_z3 = 0.05d0 ! 3Zoo
    namelist /pathirdzooloss/ lossN_z3, lossC_z3
    !!------------------------------------------------------------------------------
    !! *** Parameters for CO2 limitation ***

    ! Conversion factor between [mol/m3] (model) and [umol/kg] (function): (1000 * 1000) / 1024
    real(kind=wp) :: Cunits = 976.5625
    real(kind=wp) :: a_co2_phy = 1.162e+00 ! [unitless]
    real(kind=wp) :: a_co2_dia = 1.040e+00 ! [unitless]
    real(kind=wp) :: a_co2_cocco = 1.109e+00 ! [unitless]
    real(kind=wp) :: a_co2_phaeo = 1.162e+00 ! [unitless]
    real(kind=wp) :: a_co2_calc = 1.102e+00 ! [unitless]
    real(kind=wp) :: b_co2_phy = 4.888e+01 ! [mol/kg]
    real(kind=wp) :: b_co2_dia = 2.890e+01 ! [mol/kg]
    real(kind=wp) :: b_co2_cocco = 3.767e+01 ! [mol/kg]
    real(kind=wp) :: b_co2_phaeo = 4.888e+01 ! [mol/kg]
    real(kind=wp) :: b_co2_calc = 4.238e+01 ! [mol/kg]
    real(kind=wp) :: c_co2_phy = 2.255e-01 ! [kg/mol]
    real(kind=wp) :: c_co2_dia = 8.778e-01 ! [kg/mol]
    real(kind=wp) :: c_co2_cocco = 3.912e-01 ! [kg/mol]
    real(kind=wp) :: c_co2_phaeo = 2.255e-01 ! [kg/mol]
    real(kind=wp) :: c_co2_calc = 7.079e-01 ! [kg/mol]
    real(kind=wp) :: d_co2_phy = 1.023e+07 ! [kg/mol]
    real(kind=wp) :: d_co2_dia = 2.640e+06 ! [kg/mol]
    real(kind=wp) :: d_co2_cocco = 9.450e+06 ! [kg/mol]
    real(kind=wp) :: d_co2_phaeo = 1.023e+07 ! [kg/mol]
    real(kind=wp) :: d_co2_calc = 1.343e+07 ! [kg/mol]
    namelist /paco2lim/ Cunits, a_co2_phy, a_co2_dia, a_co2_cocco, a_co2_phaeo, a_co2_calc, &
            b_co2_phy, b_co2_dia, b_co2_cocco, b_co2_phaeo, b_co2_calc, &
            c_co2_phy, c_co2_dia, c_co2_cocco, c_co2_phaeo, c_co2_calc, &
            d_co2_phy, d_co2_dia, d_co2_cocco, d_co2_phaeo, d_co2_calc
    !!------------------------------------------------------------------------------
    !! *** Iron ***
    real(kind=wp) :: Fe2N = 0.033d0 ! Fe2C * 6.625 (Fe2C = 0.005d0)
    ! default was 0.14 Fe2C_benthos (=0.02125=0.68d0/32.d0) * 6.625 -
    ! will have to be tuned. [umol/m2/day]
    real(kind=wp) :: Fe2N_benthos = 0.15d0
    real(kind=wp) :: kScavFe = 0.07d0
    real(kind=wp) :: dust_sol = 0.02d0 !Dissolution of Dust for bioavaliable
    real(kind=wp) :: RiverFeConc = 100d0 ! mean DFe concentration in rivers
    namelist /pairon/ Fe2N, Fe2N_benthos, kScavFe, dust_sol, RiverFeConc
    !!------------------------------------------------------------------------------
    !! *** Calcification ***
    real(kind=wp) :: calc_prod_ratio = 0.02d0
    real(kind=wp) :: calc_diss_guts = 0.0d0
    real(kind=wp) :: calc_diss_rate = 0.005714d0 !20.d0/3500.d0
    real(kind=wp) :: calc_diss_rate2 = 0.005714d0
    ! NEW DISS value from Aumont et al. 2015, will be used with OmegaC_diss flag
    real(kind=wp) :: calc_diss_omegac = 0.197d0
    ! NEW DISS exponent in the dissolution rate of calcite, will be used with OmegaC_diss flag
    real(kind=wp) :: calc_diss_exp = 1.d0
    namelist /pacalc/ calc_prod_ratio, calc_diss_guts, calc_diss_rate, calc_diss_rate2, &
            calc_diss_omegac, calc_diss_exp ! NEW DISS added calc_diss_omegac, calc_diss_exp
    !!------------------------------------------------------------------------------
    !! *** Benthos ***
    real(kind=wp) :: decayRateBenN = 0.005d0
    real(kind=wp) :: decayRateBenC = 0.005d0
    real(kind=wp) :: decayRateBenSi = 0.005d0
    real(kind=wp) :: q_NC_Denit = 0.86d0 ! N:C quota of the denitrification process
    namelist /pabenthos_decay_rate/ decayRateBenN, decayRateBenC, decayRateBenSi, q_NC_Denit
    !!------------------------------------------------------------------------------
    !! *** CO2-flux ***
    ! 1.e-3/1024.5d0 ! Converting DIC from [mmol/m3] to [mol/kg]
    real(kind=wp) :: permil = 0.000000976
    real(kind=wp) :: permeg = 1.e-6 ! [atm/uatm] Changes units from uatm to atm
    real(kind=wp) :: Xacc = 1.e-12 ! Accuracy for ph-iteration (phacc)
    ! pressure of CO2
    real(kind=wp) :: CO2_for_spinup = 278.d0 !
    namelist /paco2_flux_param/ permil, permeg, Xacc, CO2_for_spinup
    !!------------------------------------------------------------------------------
    !! *** Alkalinity restoring ***
    real(kind=wp) :: surf_relax_Alk = 3.2e-07 !10.d0/31536000.d0
    namelist /paalkalinity_restoring/ surf_relax_Alk
    !!-----------------------------------------------------------------------------
    !! *** Ballasting ***                                              ! NEW BALL
    real(kind=wp) :: rho_POC = 1033.d0 ! kg m-3; density of POC (see Table 1 in Cram et al., 2018)
    real(kind=wp) :: rho_PON = 1033.d0 ! kg m-3; density of PON (see Table 1 in Cram et al., 2018)
    ! kg m-3; density of CaCO3 (see Table 1 in Cram et al., 2018)
    real(kind=wp) :: rho_CaCO3 = 2830.d0
    real(kind=wp) :: rho_opal = 2090.d0 ! kg m-3; density of Opal (see Table 1 in Cram et al., 2018)
    ! kg m-3; reference particle density (see Cram et al., 2018)
    real(kind=wp) :: rho_ref_part = 1230.d0
    ! kg m-3; reference seawater density (see Cram et al., 2018)
    real(kind=wp) :: rho_ref_water = 1027.d0
    ! kg m-1 s-1; reference seawater viscosity, at Temp=4 degC (see Cram et al., 2018)
    real(kind=wp) :: visc_ref_water = 0.d00158
    real(kind=wp) :: w_ref1 = 10.d0 ! m s-1; reference sinking velocity of small detritus
    real(kind=wp) :: w_ref2 = 200.d0 ! m s-1; reference sinking velocity of large detritus
    ! s-1; factor to increase sinking speed of det1 with depth, set to 0 if not wanted
    real(kind=wp) :: depth_scaling1 = 0.d015
    ! s-1; factor to increase sinking speed of det2 with depth, set to 0 if not wanted
    real(kind=wp) :: depth_scaling2 = 0.d0
    ! d-1; for numerical stability, set a maximum possible
    ! sinking velocity here (applies to both detritus classes)
    real(kind=wp) :: max_sinking_velocity = 250.d0
    namelist /paballasting/ rho_POC, rho_PON, rho_CaCO3, rho_opal, rho_ref_part, &
            rho_ref_water, visc_ref_water, w_ref1, w_ref2, depth_scaling1, &
            depth_scaling2, max_sinking_velocity

contains

    ! ---------------------------------------------------------------------------
    ! SUBROUTINE: initialize_tracer_indices
    ! ---------------------------------------------------------------------------
    ! Purpose: Set up tracer indices based on model configuration
    ! ---------------------------------------------------------------------------
    subroutine initialize_tracer_indices()
        implicit none

        if (enable_3zoo2det .and. enable_coccos) then
            ! =======================================================================
            ! CASE: 4 phytoplankton + 3 zooplankton + 2 detritus
            ! =======================================================================
            ! Phytoplankton: small phyto, diatoms, coccolithophores, phaeocystis
            ! Zooplankton: mesozoo, macrozoo, microzoo
            ! Detritus: det1, det2

            icocn = 29
            icocc = 30
            icchl = 31
            iphan = 32
            iphac = 33
            iphachl = 34
            imiczoon = 35
            imiczooc = 36

            !        allocate(recom_cocco_tracer_id(3))
            recom_cocco_tracer_id = [1029, 1030, 1031]

            !        allocate(recom_phaeo_tracer_id(3))
            recom_phaeo_tracer_id = [1032, 1033, 1034]

            !        allocate(recom_det2_tracer_id(4))
            recom_det2_tracer_id = [1025, 1026, 1027, 1028]

        else if (enable_coccos .and. .not.enable_3zoo2det) then
            ! =======================================================================
            ! CASE: 4 phytoplankton + 1 zooplankton + 1 detritus
            ! =======================================================================
            ! Phytoplankton: small phyto, diatoms, coccolithophores, phaeocystis
            ! Zooplankton: mesozoo only
            ! Detritus: det1 only

            icocn = 23
            icocc = 24
            icchl = 25
            iphan = 26
            iphac = 27
            iphachl = 28

            !        allocate(recom_cocco_tracer_id(3))
            recom_cocco_tracer_id = [1023, 1024, 1025]

            !        allocate(recom_phaeo_tracer_id(3))
            recom_phaeo_tracer_id = [1026, 1027, 1028]

        else if (enable_3zoo2det .and. .not.enable_coccos) then
            ! =======================================================================
            ! CASE: 2 phytoplankton + 3 zooplankton + 2 detritus
            ! =======================================================================
            ! Phytoplankton: small phyto, diatoms only
            ! Zooplankton: mesozoo, macrozoo, microzoo
            ! Detritus: det1, det2

            imiczoon = 29
            imiczooc = 30

            !        allocate(recom_det2_tracer_id(4))
            recom_det2_tracer_id = [1025, 1026, 1027, 1028]
        else
            ! =======================================================================
            ! CASE: 2 phytoplankton + 1 zooplankton + 1 detritus (BASE CONFIGURATION)
            ! =======================================================================
            ! Phytoplankton: small phyto, diatoms only
            ! Zooplankton: mesozoo only
            ! Detritus: det1 only
            ! (All indices already set to default values)
        end if
    end subroutine initialize_tracer_indices

    ! ==============================================================================
    ! SUBROUTINE: validate_recom_tracers
    ! ==============================================================================
    ! Purpose: Validate consistency between namelist tracer configuration and
    !          biogeochemical model setup (enable_3zoo2det, enable_coccos)
    ! ==============================================================================
    subroutine validate_recom_tracers(num_tracers, mype)
        use mpi, only: MPI_Abort, MPI_COMM_WORLD

        implicit none

        ! Arguments
        integer, intent(in) :: num_tracers ! Total number of tracers from namelist
        integer, intent(in) :: mype ! MPI rank

        ! Local variables
        integer :: expected_bgc_num
        integer :: actual_bgc_num
        integer :: expected_total_tracers
        integer :: num_physical_tracers
        integer :: MPIErr
        logical :: config_error

        ! For tracer ID validation
        integer :: i
        integer, dimension(:), allocatable :: expected_tracer_ids
        logical, dimension(:), allocatable :: tracer_found
        integer :: num_expected_tracers
        logical :: id_error

        ! Physical tracers (temperature, salinity, etc.) - typically first 2
        num_physical_tracers = 2

        ! Calculate actual BGC tracer count from namelist
        actual_bgc_num = num_tracers - num_physical_tracers

        ! ===========================================================================
        ! Determine expected BGC tracer count based on configuration
        ! ===========================================================================
        config_error = .false.

        if (enable_3zoo2det .and. enable_coccos) then
            ! ---------------------------------------------------------------------------
            ! Configuration 4: Full model (4 phyto + 3 zoo + 2 detritus)
            ! ---------------------------------------------------------------------------
            ! Base: 22 tracers (1001-1022)
            ! Additional 3zoo2det: 4 tracers for det2 (1025-1028)
            ! Additional coccos: 6 tracers for coccos (1029-1031)
            ! Additional phaeocystis: 3 tracers (1032-1034)
            ! Additional microzoo: 2 tracers (1035-1036)
            ! Total: 22 + 4 + 6 + 3 + 2 = 36 (actually 22 + 14 = 36)
            expected_bgc_num = 36

        else if (enable_coccos .and. .not.enable_3zoo2det) then
            ! ---------------------------------------------------------------------------
            ! Configuration 3: Coccos only (4 phyto + 1 zoo + 1 detritus)
            ! ---------------------------------------------------------------------------
            ! Base: 22 tracers (1001-1022)
            ! Additional coccos: 3 tracers (1023-1025)
            ! Additional phaeocystis: 3 tracers (1026-1028)
            ! Total: 22 + 6 = 28
            expected_bgc_num = 28

        else if (enable_3zoo2det .and. .not.enable_coccos) then
            ! ---------------------------------------------------------------------------
            ! Configuration 2: 3Zoo2Det only (2 phyto + 3 zoo + 2 detritus)
            ! ---------------------------------------------------------------------------
            ! Base: 22 tracers (1001-1022)
            ! Additional zoo2: 2 tracers (1023-1024)
            ! Additional det2: 4 tracers (1025-1028)
            ! Additional microzoo: 2 tracers (1029-1030)
            ! Total: 22 + 8 = 30
            expected_bgc_num = 30

        else
            ! ---------------------------------------------------------------------------
            ! Configuration 1: Base model (2 phyto + 1 zoo + 1 detritus)
            ! ---------------------------------------------------------------------------
            ! Base: 22 tracers (1001-1022)
            expected_bgc_num = 22

        end if

        expected_total_tracers = num_physical_tracers + expected_bgc_num

        ! ===========================================================================
        ! Build expected tracer ID list for current configuration
        ! ===========================================================================

        ! Determine total expected tracers
        num_expected_tracers = expected_total_tracers
        allocate(expected_tracer_ids(num_expected_tracers))
        allocate(tracer_found(num_expected_tracers))
        tracer_found = .false.

        ! Physical tracers (always present)
        expected_tracer_ids(1) = 1 ! Temperature
        expected_tracer_ids(2) = 2 ! Salinity

        ! Base BGC tracers (always present for all configurations)
        do i = 1, 22
            expected_tracer_ids(num_physical_tracers + i) = 1000 + i
        end do

        ! Configuration-specific tracers
        if (enable_3zoo2det .and. enable_coccos) then
            ! Full model: 1001-1022 (base) + 1023-1024 (zoo2) + 1025-1028 (det2) + 1029-1036
            ! (coccos+phaeo+zoo3)
            expected_tracer_ids(25) = 1023 ! Zoo2N
            expected_tracer_ids(26) = 1024 ! Zoo2C
            expected_tracer_ids(27) = 1025 ! DetZ2N
            expected_tracer_ids(28) = 1026 ! DetZ2C
            expected_tracer_ids(29) = 1027 ! DetZ2Si
            expected_tracer_ids(30) = 1028 ! DetZ2Calc
            expected_tracer_ids(31) = 1029 ! CoccoN
            expected_tracer_ids(32) = 1030 ! CoccoC
            expected_tracer_ids(33) = 1031 ! CoccoChl
            expected_tracer_ids(34) = 1032 ! PhaeoN
            expected_tracer_ids(35) = 1033 ! PhaeoC
            expected_tracer_ids(36) = 1034 ! PhaeoChl
            expected_tracer_ids(37) = 1035 ! Zoo3N
            expected_tracer_ids(38) = 1036 ! Zoo3C

        else if (enable_coccos .and. .not.enable_3zoo2det) then
            ! Coccos only: 1001-1022 (base) + 1023-1028 (coccos+phaeo)
            expected_tracer_ids(25) = 1023 ! CoccoN
            expected_tracer_ids(26) = 1024 ! CoccoC
            expected_tracer_ids(27) = 1025 ! CoccoChl
            expected_tracer_ids(28) = 1026 ! PhaeoN
            expected_tracer_ids(29) = 1027 ! PhaeoC
            expected_tracer_ids(30) = 1028 ! PhaeoChl

        else if (enable_3zoo2det .and. .not.enable_coccos) then
            ! 3Zoo2Det only: 1001-1022 (base) + 1023-1030 (zoo2+det2+zoo3)
            expected_tracer_ids(25) = 1023 ! Zoo2N
            expected_tracer_ids(26) = 1024 ! Zoo2C
            expected_tracer_ids(27) = 1025 ! DetZ2N
            expected_tracer_ids(28) = 1026 ! DetZ2C
            expected_tracer_ids(29) = 1027 ! DetZ2Si
            expected_tracer_ids(30) = 1028 ! DetZ2Calc
            expected_tracer_ids(31) = 1029 ! Zoo3N
            expected_tracer_ids(32) = 1030 ! Zoo3C
        end if
        ! else: base configuration only needs tracers 1, 2, 1001-1022

        ! ===========================================================================
        ! Perform validation checks
        ! ===========================================================================

        if (mype == 0) then
            write(*, *) ''
            write(*, *) '=========================================================================='
            write(*, *) 'REcoM TRACER CONFIGURATION VALIDATION'
            write(*, *) '=========================================================================='
            write(*, *) 'Model configuration:'
            write(*, *) '  enable_3zoo2det = ', enable_3zoo2det
            write(*, *) '  enable_coccos   = ', enable_coccos
            write(*, *) ''
            write(*, *) 'Tracer counts:'
            write(*, *) '  Physical tracers (T, S, ...)      = ', num_physical_tracers
            write(*, *) '  Expected BGC tracers              = ', expected_bgc_num
            write(*, *) '  Expected TOTAL tracers            = ', expected_total_tracers
            write(*, *) '  Actual tracers from namelist      = ', num_tracers
            write(*, *) '  Actual BGC tracers from namelist  = ', actual_bgc_num
            write(*, *) ''
        end if

        ! Check for inconsistencies
        if (actual_bgc_num /= expected_bgc_num) then
            config_error = .true.
            if (mype == 0) then
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'ERROR: TRACER COUNT MISMATCH!'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'The number of BGC tracers in the namelist does not match'
                write(*, *) 'the expected count for the current configuration.'
                write(*, *) ''
                write(*, *) '  Expected BGC tracers: ', expected_bgc_num
                write(*, *) '  Actual BGC tracers:   ', actual_bgc_num
                write(*, *) '  Difference:           ', actual_bgc_num - expected_bgc_num
                write(*, *) ''
                write(*, *) 'Required tracer IDs for current configuration:'
                write(*, *) '  Base tracers (always):  1001-1022 (22 tracers)'

                if (enable_3zoo2det .and. .not.enable_coccos) then
                    write(*, *) '  3Zoo2Det extension:     1023-1030 (8 tracers)'
                    write(*, *) '    - Zoo2N, Zoo2C:       1023-1024'
                    write(*, *) '    - DetZ2 pool:         1025-1028'
                    write(*, *) '    - MicZooN, MicZooC:   1029-1030'
                else if (enable_coccos .and. .not.enable_3zoo2det) then
                    write(*, *) '  Coccos extension:       1023-1028 (6 tracers)'
                    write(*, *) '    - CoccoN, C, Chl:     1023-1025'
                    write(*, *) '    - PhaeoN, C, Chl:     1026-1028'
                else if (enable_3zoo2det .and. enable_coccos) then
                    write(*, *) '    - Zoo2N, Zoo2C:       1023-1024'
                    write(*, *) '  3Zoo2Det extension:     1025-1028 (4 tracers for det2)'
                    write(*, *) '  Coccos extension:       1029-1034 (6 tracers)'
                    write(*, *) '    - CoccoN, C, Chl:     1029-1031'
                    write(*, *) '    - PhaeoN, C, Chl:     1032-1034'
                    write(*, *) '  MicroZoo extension:     1035-1036 (2 tracers)'
                end if

                write(*, *) ''
                write(*, *) 'ACTION REQUIRED:'
                write(*, *) '  1. Check your namelist.config tracer_list section'
                write(*, *) '  2. Ensure enable_3zoo2det and enable_coccos match your setup'
                write(*, *) '  3. Add/remove tracers to match the expected configuration'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
            end if
        else
            ! Validation passed
            if (mype == 0) then
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'VALIDATION PASSED: Tracer configuration is consistent!'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
            end if
        end if

        ! ===========================================================================
        ! Additional sanity check: verify bgc_num variable matches
        ! ===========================================================================
        if (bgc_num /= expected_bgc_num) then
            if (mype == 0) then
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'WARNING: bgc_num variable inconsistency!'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'The bgc_num parameter does not match the expected value.'
                write(*, *) '  Current bgc_num value: ', bgc_num
                write(*, *) '  Expected value:        ', expected_bgc_num
                write(*, *) ''
                write(*, *) 'This may indicate that bgc_num was not updated after changing'
                write(*, *) 'enable_3zoo2det or enable_coccos flags.'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
            end if
            config_error = .true.
        end if

        ! ===========================================================================
        ! Validate tracer IDs: Check for correct IDs and detect clashes
        ! ===========================================================================
        id_error = .false.

        ! This check requires access to the actual tracer IDs from the namelist
        ! We'll validate against the expected list
        if (mype == 0) then
            write(*, *) '=========================================================================='
            write(*, *) 'VALIDATING TRACER IDs'
            write(*, *) '=========================================================================='
            write(*, *) 'Expected tracer ID sequence:'
            write(*, *) ''

            ! Display expected IDs in a readable format
            write(*, *) 'Physical tracers:'
            write(*, *) '  ', expected_tracer_ids(1:num_physical_tracers)
            write(*, *) ''
            write(*, *) 'Base BGC tracers (1001-1022):'
            write(*, *) '  ', expected_tracer_ids(3:24)
            write(*, *) ''

            if (expected_bgc_num > 22) then
                write(*, *) 'Extended configuration tracers:'
                write(*, *) '  ', expected_tracer_ids(25:num_expected_tracers)
                write(*, *) ''
            end if

            write(*, *) 'CRITICAL: The tracer IDs in your namelist MUST match this sequence'
            write(*, *) '          exactly, in the same order!'
            write(*, *) ''
            write(*, *) 'Common errors to avoid:'
            write(*, *) '  - Using wrong tracer ID numbers (e.g., 1023 instead of 1025)'
            write(*, *) '  - Tracer ID clashes between configurations'
            write(*, *) '  - Incorrect order of tracer IDs in namelist'
            write(*, *) '  - Missing or duplicate tracer IDs'
            write(*, *) ''

            ! Configuration-specific warnings
            if (enable_3zoo2det .and. enable_coccos) then
                write(*, *) 'IMPORTANT for FULL MODEL (3zoo2det + coccos):'
                !  write(*,*) '  - Tracers 1023-1024 are NOT used (reserved for other configs)'
                write(*, *) '  - Zoo2 uses:         1023-1024'
                write(*, *) '  - Det2 pool uses:    1025-1028'
                write(*, *) '  - Coccos uses:       1029-1031'
                write(*, *) '  - Phaeocystis uses:  1032-1034'
                write(*, *) '  - Microzooplankton:  1035-1036'
                write(*, *) ''
            else if (enable_coccos .and. .not.enable_3zoo2det) then
                write(*, *) 'IMPORTANT for COCCOS-ONLY configuration:'
                write(*, *) '  - Coccos uses:       1023-1025 (NOT 1029-1031)'
                write(*, *) '  - Phaeocystis uses:  1026-1028 (NOT 1032-1034)'
                write(*, *) '  - Tracers 1029+ are NOT used in this configuration'
                write(*, *) ''
            else if (enable_3zoo2det .and. .not.enable_coccos) then
                write(*, *) 'IMPORTANT for 3ZOO2DET-ONLY configuration:'
                write(*, *) '  - Zoo2 uses:         1023-1024'
                write(*, *) '  - Det2 pool uses:    1025-1028'
                write(*, *) '  - Microzoo uses:     1029-1030 (NOT 1035-1036)'
                write(*, *) '  - Tracers 1031+ are NOT used in this configuration'
                write(*, *) ''
            else
                write(*, *) 'IMPORTANT for BASE configuration:'
                write(*, *) '  - Only tracers 1-2, 1001-1022 should be present'
                write(*, *) '  - Tracers 1023+ are NOT used in base configuration'
                write(*, *) ''
            end if

            write(*, *) '=========================================================================='
            write(*, *) ''
        end if

        ! ===========================================================================
        ! Check for tracer ID clashes based on configuration
        ! ===========================================================================
        if (mype == 0) then
            write(*, *) '=========================================================================='
            write(*, *) 'CHECKING FOR TRACER ID CONFLICTS'
            write(*, *) '=========================================================================='

            ! Warn about potential clashes between configurations
            if (enable_3zoo2det .and. enable_coccos) then
                write(*, *) 'Full model configuration active.'
                !  write(*,*) 'Ensure you are NOT using tracer IDs 1023-1024 in your namelist!'
                !  write(*,*) 'These are reserved for configurations WITHOUT full model.'
            else if (enable_coccos) then
                write(*, *) 'Coccos-only configuration active.'
                write(*, *) 'Coccos MUST use IDs 1023-1025 (NOT 1029-1031).'
                write(*, *) 'Phaeocystis MUST use IDs 1026-1028 (NOT 1032-1034).'
            else if (enable_3zoo2det) then
                write(*, *) '3Zoo2Det-only configuration active.'
                write(*, *) 'Microzoo MUST use IDs 1029-1030 (NOT 1035-1036).'
            end if

            write(*, *) ''
            ! write(*,*) 'No automated clash detection available without tracer array access.'
            write(*, *) 'Please manually verify your namelist tracer_list against the'
            write(*, *) 'expected sequence shown above.'
            write(*, *) '=========================================================================='
            write(*, *) ''
        end if

        ! ===========================================================================
        ! Stop execution if configuration error detected
        ! ===========================================================================
        if (config_error) then
            if (mype == 0) then
                write(*, *) ''
                write(*, *) '******************************************************************'
                write(*, *) '***  FATAL ERROR: MODEL CONFIGURATION INCONSISTENCY DETECTED   ***'
                write(*, *) '***  MODEL EXECUTION STOPPED                                   ***'
                write(*, *) '******************************************************************'
                write(*, *) ''
            end if
            deallocate(expected_tracer_ids, tracer_found)
            ! Stop execution (use appropriate stop routine for your model)
            call MPI_ABORT(MPI_COMM_WORLD, 1, MPIErr)
            stop
        end if

        ! Clean up
        deallocate(expected_tracer_ids, tracer_found)

    end subroutine validate_recom_tracers

    ! ==============================================================================
    ! SUBROUTINE: validate_tracer_id_sequence
    ! ==============================================================================
    ! Purpose: Validate that actual tracer IDs from namelist match expected sequence
    !          Call this after reading the tracer namelist
    ! ==============================================================================
    subroutine validate_tracer_id_sequence(tracer_ids, num_tracers, mype)
        use mpi, only: MPI_Abort, MPI_COMM_WORLD

        implicit none

        ! Arguments
        integer, dimension(:), intent(in) :: tracer_ids ! Actual IDs from namelist
        integer, intent(in) :: num_tracers ! Number of tracers
        integer, intent(in) :: mype ! MPI rank

        ! Local variables
        integer :: i, j
        integer, dimension(:), allocatable :: expected_ids
        integer :: MPIErr
        logical :: error_found
        logical :: duplicate_found
        integer :: num_physical_tracers

        error_found = .false.
        duplicate_found = .false.
        num_physical_tracers = 2

        ! Allocate expected IDs array
        allocate(expected_ids(num_tracers))

        ! Build expected ID sequence
        expected_ids(1) = 1
        expected_ids(2) = 2

        do i = 1, 22
            expected_ids(num_physical_tracers + i) = 1000 + i
        end do

        if (enable_3zoo2det .and. enable_coccos) then
            ! Full model configuration
            expected_ids(25:30) = [1023, 1024, 1025, 1026, 1027, 1028]
            expected_ids(31:36) = [1029, 1030, 1031, 1032, 1033, 1034]
            expected_ids(37:38) = [1035, 1036]

        else if (enable_coccos .and. .not.enable_3zoo2det) then
            expected_ids(25:30) = [1023, 1024, 1025, 1026, 1027, 1028]

        else if (enable_3zoo2det .and. .not.enable_coccos) then
            expected_ids(25:32) = [1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030]
        end if

        ! ===========================================================================
        ! Check 1: Compare actual vs expected tracer IDs
        ! ===========================================================================
        if (mype == 0) then
            write(*, *) ''
            write(*, *) '=========================================================================='
            write(*, *) 'VALIDATING TRACER ID SEQUENCE FROM NAMELIST'
            write(*, *) '=========================================================================='
        end if

        do i = 1, num_tracers
            if (tracer_ids(i) /= expected_ids(i)) then
                error_found = .true.
                if (mype == 0) then
                    write(*, *) 'ERROR at position ', i, ':'
                    write(*, *) '  Expected tracer ID: ', expected_ids(i)
                    write(*, *) '  Found tracer ID:    ', tracer_ids(i)
                    write(*, *) ''
                end if
            end if
        end do

        ! ===========================================================================
        ! Check 2: Detect duplicate tracer IDs
        ! ===========================================================================
        do i = 1, num_tracers - 1
            do j = i + 1, num_tracers
                if (tracer_ids(i) == tracer_ids(j)) then
                    duplicate_found = .true.
                    if (mype == 0) then
                        write(*, *) 'ERROR: Duplicate tracer ID detected!'
                        write(*, *) '  Tracer ID ', tracer_ids(i), ' appears at positions ', i, &
                                ' and ', j
                        write(*, *) ''
                    end if
                end if
            end do
        end do

        ! ===========================================================================
        ! Check 3: Detect forbidden tracer IDs for current configuration
        ! ===========================================================================
        !if (enable_3zoo2det .and. enable_coccos) then
        ! Check for forbidden IDs 1023-1024 in full model
        !do i = 1, num_tracers
        !if (tracer_ids(i) == 1023 .or. tracer_ids(i) == 1024) then
        !error_found = .true.
        !if (mype == 0) then
        !write(*,*) 'ERROR: Forbidden tracer ID in full model configuration!'
        !write(*,*) '  Tracer ID ', tracer_ids(i), ' at position ', i
        !write(*,*) '  IDs 1023-1024 are NOT used when both flags are enabled'
        !write(*,*) ''
        !end if
        !end if
        !end do
        !end if

        ! ===========================================================================
        ! Report results
        ! ===========================================================================
        if (error_found .or. duplicate_found) then
            if (mype == 0) then
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'TRACER ID VALIDATION FAILED!'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
                write(*, *) 'Expected tracer ID sequence for current configuration:'
                write(*, *) expected_ids
                write(*, *) ''
                write(*, *) 'Actual tracer ID sequence from namelist:'
                write(*, *) tracer_ids
                write(*, *) ''
                write(*, *) 'ACTION REQUIRED:'
                write(*, *) '  Correct the tracer IDs in your namelist.config file'
                write(*, *) '  Ensure the sequence matches exactly as expected'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
                write(*, *) '******************************************************************'
                write(*, *) '***  FATAL ERROR: INVALID TRACER ID SEQUENCE                   ***'
                write(*, *) '***  MODEL EXECUTION STOPPED                                   ***'
                write(*, *) '******************************************************************'
                write(*, *) ''
            end if
            deallocate(expected_ids)
            ! Stop execution (use appropriate stop routine for your model)
            call MPI_ABORT(MPI_COMM_WORLD, 1, MPIErr)
            stop
        else
            if (mype == 0) then
                write(*, *) '======================================================================&
                        &===='
                write(*, *) 'TRACER ID VALIDATION PASSED!'
                write(*, *) 'All tracer IDs match expected sequence - no clashes detected.'
                write(*, *) '======================================================================&
                        &===='
                write(*, *) ''
            end if
        end if

        deallocate(expected_ids)

    end subroutine validate_tracer_id_sequence

end module recom_config
!===============================================================================
! For arrays needed for the whole 2D or 3D domain, but only needed in REcoM
!-------------------------------------------------------------------------------
module REcoM_GloVar
    use recom_declarations, only: wp

    implicit none
    public

    private :: wp

    save

    ! 4 types of benthos-tracers with size [4 n2d]
    real(kind=wp), allocatable, dimension(:, :) :: Benthos
    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global
    ! sums when running the tracer loop in parallel
    real(kind=wp), allocatable, dimension(:, :, :) :: Benthos_tr

    ! [umol/m2/s] Monthly 2D field of iron soluted in surface water from dust
    real(kind=wp), allocatable, dimension(:) :: GloFeDust
    ! [mmol/m2/s] 10-year mean 2D fields of nitrogen soluted in surface water from dust
    real(kind=wp), allocatable, dimension(:) :: GloNDust
    ! [uatm] Atmospheric CO2 partial pressure. One value for the whole planet for each month
    real(kind=wp), dimension(12) :: AtmCO2

    ! [umol/m2/s] Includes ice, but is, other than that identlical to GloFeDust
    real(kind=wp), allocatable, dimension(:) :: AtmFeInput
    ! [umol/m2/s] Includes ice, but is, other than that identlical to GloNDust
    real(kind=wp), allocatable, dimension(:) :: AtmNInput
    ! [uatm] Surface ocean CO2 partial pressure
    real(kind=wp), allocatable, dimension(:) :: GloPCO2surf
    ! [mmol/m2/day] Positive downwards
    real(kind=wp), allocatable, dimension(:) :: GloCO2flux
    ! [mmol/m2/day] Positive downwards
    real(kind=wp), allocatable, dimension(:) :: GloO2flux
    ! [mmol/m2/day] Positive downwards
    real(kind=wp), allocatable, dimension(:) :: GloCO2flux_seaicemask
    ! [mmol/m2/day] Positive downwards
    real(kind=wp), allocatable, dimension(:) :: GloO2flux_seaicemask
    ! [mol/kg] Concentrations of H-plus ions in the surface ocean
    real(kind=wp), allocatable, dimension(:) :: GloHplus

    ! MOCSY: [mol/m3] Aqueous CO2 concentration for all depths
    real(kind=wp), allocatable, dimension(:, :) :: CO23D
    ! MOCSY: total scale
    real(kind=wp), allocatable, dimension(:, :) :: pH3D
    ! MOCSY: [uatm] CO2 partial pressure
    real(kind=wp), allocatable, dimension(:, :) :: pCO23D
    ! MOCSY: [mol/m3] Bicarbonate ion concentration
    real(kind=wp), allocatable, dimension(:, :) :: HCO33D
    ! DISS: [mol/m3] Carbonate ion concentration
    real(kind=wp), allocatable, dimension(:, :) :: CO33D
    ! DISS: calcite saturation state
    real(kind=wp), allocatable, dimension(:, :) :: OmegaC3D
    ! DISS: [mol^2/kg^2] stoichiometric solubility product of calcite
    real(kind=wp), allocatable, dimension(:, :) :: kspc3D
    ! DISS: [mol/m3] in-situ density of seawater
    real(kind=wp), allocatable, dimension(:, :) :: rhoSW3D

    ! BALL: density of particle class 1
    real(kind=wp), allocatable, dimension(:, :) :: rho_particle1
    ! BALL: density of particle class 2
    real(kind=wp), allocatable, dimension(:, :) :: rho_particle2
    ! BALL: scaling factor
    real(kind=wp), allocatable, dimension(:, :) :: scaling_density1_3D
    ! BALL: scaling factor
    real(kind=wp), allocatable, dimension(:, :) :: scaling_density2_3D
    ! BALL: scaling factor
    real(kind=wp), allocatable, dimension(:, :) :: scaling_visc_3D
    ! BALL: scaling factor
    real(kind=wp), allocatable, dimension(:, :) :: seawater_visc_3D

    ! [mmol/m2/day] ocean-atmosphere
    real(kind=wp), allocatable, dimension(:) :: GlodPCO2surf
    ! [1/day] Decay rate of detritus in the benthic layer saved for oce_ale_tracer.F90
    real(kind=wp), allocatable, dimension(:, :) :: GlodecayBenthos
    ! [m s-1]
    real(kind=wp), allocatable, dimension(:) :: PistonVelocity
    ! [mol L-1 atm-1]
    real(kind=wp), allocatable, dimension(:) :: alphaCO2

    real(kind=wp), allocatable, dimension(:, :) :: GlowFluxDet
    real(kind=wp), allocatable, dimension(:, :) :: GlowFluxPhy
    real(kind=wp), allocatable, dimension(:, :) :: GlowFluxDia
    real(kind=wp), allocatable, dimension(:, :) :: GlowFluxCocco
    real(kind=wp), allocatable, dimension(:, :) :: GlowFluxPhaeo

    real(kind=wp), allocatable, dimension(:, :) :: diags2D
    real(kind=wp), allocatable, dimension(:) :: NPPn
    real(kind=wp), allocatable, dimension(:) :: NPPd
    real(kind=wp), allocatable, dimension(:) :: GPPn
    real(kind=wp), allocatable, dimension(:) :: GPPd
    real(kind=wp), allocatable, dimension(:) :: NNAn
    real(kind=wp), allocatable, dimension(:) :: NNAd
    real(kind=wp), allocatable, dimension(:) :: Chldegn
    real(kind=wp), allocatable, dimension(:) :: Chldegd
    real(kind=wp), allocatable, dimension(:) :: NPPc
    real(kind=wp), allocatable, dimension(:) :: GPPc
    real(kind=wp), allocatable, dimension(:) :: NNAc
    real(kind=wp), allocatable, dimension(:) :: Chldegc
    ! Phaeocystis
    real(kind=wp), allocatable, dimension(:) :: NPPp
    real(kind=wp), allocatable, dimension(:) :: GPPp
    real(kind=wp), allocatable, dimension(:) :: NNAp
    real(kind=wp), allocatable, dimension(:) :: Chldegp
    real(kind=wp), allocatable, dimension(:) :: grazmeso_tot
    real(kind=wp), allocatable, dimension(:) :: grazmeso_n
    real(kind=wp), allocatable, dimension(:) :: grazmeso_d
    real(kind=wp), allocatable, dimension(:) :: grazmeso_c
    real(kind=wp), allocatable, dimension(:) :: grazmeso_p
    real(kind=wp), allocatable, dimension(:) :: grazmeso_det
    real(kind=wp), allocatable, dimension(:) :: grazmeso_mic
    real(kind=wp), allocatable, dimension(:) :: grazmeso_det2
    real(kind=wp), allocatable, dimension(:) :: grazmacro_tot
    real(kind=wp), allocatable, dimension(:) :: grazmacro_n
    real(kind=wp), allocatable, dimension(:) :: grazmacro_d
    real(kind=wp), allocatable, dimension(:) :: grazmacro_c
    real(kind=wp), allocatable, dimension(:) :: grazmacro_p
    real(kind=wp), allocatable, dimension(:) :: grazmacro_mes
    real(kind=wp), allocatable, dimension(:) :: grazmacro_det
    real(kind=wp), allocatable, dimension(:) :: grazmacro_mic
    real(kind=wp), allocatable, dimension(:) :: grazmacro_det2
    real(kind=wp), allocatable, dimension(:) :: grazmicro_tot
    real(kind=wp), allocatable, dimension(:) :: grazmicro_n
    real(kind=wp), allocatable, dimension(:) :: grazmicro_d
    real(kind=wp), allocatable, dimension(:) :: grazmicro_c
    real(kind=wp), allocatable, dimension(:) :: grazmicro_p
    real(kind=wp), allocatable, dimension(:, :) :: respmeso
    real(kind=wp), allocatable, dimension(:, :) :: respmacro
    real(kind=wp), allocatable, dimension(:, :) :: respmicro
    real(kind=wp), allocatable, dimension(:, :) :: calcdiss
    real(kind=wp), allocatable, dimension(:, :) :: calcif
    real(kind=wp), allocatable, dimension(:, :) :: aggn
    real(kind=wp), allocatable, dimension(:, :) :: aggd
    real(kind=wp), allocatable, dimension(:, :) :: aggc

    ! Phaeocystis
    real(kind=wp), allocatable, dimension(:, :) :: aggp

    real(kind=wp), allocatable, dimension(:, :) :: docexn
    real(kind=wp), allocatable, dimension(:, :) :: docexd

    real(kind=wp), allocatable, dimension(:, :) :: docexc
    ! Phaeocystis

    real(kind=wp), allocatable, dimension(:, :) :: docexp
    real(kind=wp), allocatable, dimension(:, :) :: respn
    real(kind=wp), allocatable, dimension(:, :) :: respd
    real(kind=wp), allocatable, dimension(:, :) :: respc

    ! Phaeocystis
    real(kind=wp), allocatable, dimension(:, :) :: respp

    real(kind=wp), allocatable, dimension(:, :) :: NPPn3D
    real(kind=wp), allocatable, dimension(:, :) :: NPPd3D
    real(kind=wp), allocatable, dimension(:, :) :: NPPc3D

    ! Phaeocystis
    real(kind=wp), allocatable, dimension(:, :) :: NPPp3D

    ! my new variables to track
    real(kind=wp), allocatable, dimension(:, :) :: TTemp_diatoms

    ! new Temperature effect
    real(kind=wp), allocatable, dimension(:, :) :: TTemp_phyto

    ! new
    real(kind=wp), allocatable, dimension(:, :) :: TTemp_cocco

    ! new
    real(kind=wp), allocatable, dimension(:, :) :: TTemp_phaeo

    ! new CO2 effect
    real(kind=wp), allocatable, dimension(:, :) :: TPhyCO2
    real(kind=wp), allocatable, dimension(:, :) :: TDiaCO2
    real(kind=wp), allocatable, dimension(:, :) :: TCoccoCO2
    real(kind=wp), allocatable, dimension(:, :) :: TPhaeoCO2
    ! new nutrient limitation
    real(kind=wp), allocatable, dimension(:, :) :: TqlimitFac_phyto

    real(kind=wp), allocatable, dimension(:, :) :: TqlimitFac_diatoms
    real(kind=wp), allocatable, dimension(:, :) :: TqlimitFac_cocco
    real(kind=wp), allocatable, dimension(:, :) :: TqlimitFac_phaeo

    ! new light limitation
    real(kind=wp), allocatable, dimension(:, :) :: TCphotLigLim_phyto

    ! new
    real(kind=wp), allocatable, dimension(:, :) :: TCphot_phyto

    ! new light limitation
    real(kind=wp), allocatable, dimension(:, :) :: TCphotLigLim_diatoms

    real(kind=wp), allocatable, dimension(:, :) :: TCphot_diatoms

    ! new light limitation
    real(kind=wp), allocatable, dimension(:, :) :: TCphotLigLim_cocco

    real(kind=wp), allocatable, dimension(:, :) :: TCphot_cocco

    ! new light limitation
    real(kind=wp), allocatable, dimension(:, :) :: TCphotLigLim_phaeo

    real(kind=wp), allocatable, dimension(:, :) :: TCphot_phaeo

    ! tracking the assimilation of Si by Diatoms
    real(kind=wp), allocatable, dimension(:, :) :: TSi_assimDia

    ! Benthic denitrification Field in 2D [n2d 1]
    real(kind=wp), allocatable, dimension(:) :: DenitBen

    !  for using MEDUSA

    ! Diagnostics in 2D [4 n2d] or [6 n2d] with ciso
    real(kind=wp), allocatable, dimension(:, :) :: SinkFlx

    ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global
    ! sums when running the tracer loop in parallel
    real(kind=wp), allocatable, dimension(:, :, :) :: SinkFlx_tr

    ! Diagnostics for vertical sinking
    real(kind=wp), allocatable, dimension(:, :) :: Sinkingvel1

    ! Diagnostics for vertical sinking
    real(kind=wp), allocatable, dimension(:, :) :: Sinkingvel2

    ! Sinking speed of particle class 1 OG 16.03.23
    real(kind=wp), allocatable, dimension(:, :, :) :: Sinkvel1_tr

    ! Sinking speed of particle class 2 OG 16.03.23
    real(kind=wp), allocatable, dimension(:, :, :) :: Sinkvel2_tr

    ! Yearly input into bottom water from sediments [n2d 5] or [n2d 7] with ciso
    real(kind=wp), allocatable, dimension(:, :) :: GloSed
    ! Yearly burial from medusa: [n2d 5] or [n2d 9] with ciso_14
    real(kind=wp), allocatable, dimension(:, :) :: lb_flux

    ! atmospheric box model:
    ! atmospheric CO2 mixing ratio (mole fraction)
    real(kind=wp), allocatable, dimension(:) :: x_co2atm

    ! Surface alkalinity field used for restoring
    real(kind=wp), allocatable, dimension(:) :: Alk_surf
    real(kind=wp), allocatable, dimension(:) :: relax_alk
    real(kind=wp), allocatable, dimension(:) :: virtual_alk

    ! Light in the water column [nl-1 n2d]
    real(kind=wp), allocatable, dimension(:, :) :: PAR3D
    real(kind=wp), allocatable, dimension(:) :: RiverineLonOrig, RiverineLatOrig, RiverineDINOrig, &
    ! Variables to save original values for riverine nutrients
            RiverineDONOrig, RiverineDOCOrig, RiverineDSiOrig
    real(kind=wp), allocatable, dimension(:) :: RiverDIN2D, RiverDON2D, RiverDOC2D, RiverDSi2D, &
            RiverAlk2D, RiverDIC2D, RiverFe
    real(kind=wp), allocatable, dimension(:) :: ErosionTSi2D, ErosionTON2D, ErosionTOC2D
    !! Cobeta, Cos(Angle of incidence)
    real(kind=wp), allocatable, dimension(:) :: cosAI

    type :: tracer_data_pointer
        real(kind=wp), dimension(:, :), pointer :: tracer_data => null()
    end type tracer_data_pointer

    type :: tracers_info_type
        integer, dimension(:), allocatable :: ids
        logical, dimension(:), allocatable :: ltra_diag
        type(tracer_data_pointer), dimension(:), allocatable :: data_pointers
    end type tracers_info_type

end module REcoM_GloVar

!===============================================================================
! For variables saved locally for each column and then used in REcoM
!-------------------------------------------------------------------------------
module REcoM_locVar
    use recom_declarations, only: wp

    implicit none
    public

    private :: wp

    ! Storing the values for benthos in current watercolumn: N,C,Si and Calc
    real(kind=wp), allocatable, dimension(:) :: LocBenthos

    ! [mol/kg] Concentrations of H-plus ions in the surface node
    real(kind=wp) :: Hplus

    ! [uatm] Partial pressure of CO2 in surface layer at current 2D node
    real(kind=wp) :: pCO2surf(1)

    ! [mmol/m2/day] Flux of CO2 into the ocean
    real(kind=wp) :: dflux(1)

    ! [mmol/m2/day] Flux of O2 into the ocean
    real(kind=wp) :: oflux(1)

    ! [mmol/m2/s] Flux of O2 into the ocean
    real(kind=wp) :: o2ex(1)

    ! Wind strength above current 2D node, change array size if used with mocsy input vector longer
    ! than one
    real(kind=wp) :: ULoc(1)

    ! [uatm] difference of oceanic pCO2 minus atmospheric pCO2
    real(kind=wp) :: dpCO2surf(1)

    ! mocsy output
    ! --------------------------------------------------------------------------------------------
    ! air-to-sea flux of CO2 [mol/(m^2 * s)]
    real(kind=wp) :: co2flux(1)

    ! time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
    real(kind=wp) :: co2ex(1)

    ! difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
    real(kind=wp) :: dpco2(1)

    ! pH on total scale
    real(kind=wp) :: ph(1)

    ! oceanic partial pressure of CO2 (uatm)
    real(kind=wp) :: pco2(1)

    ! oceanic fugacity of CO2 (uatm)
    real(kind=wp) :: fco2(1)

    ! aqueous CO2 concentration [mol/m^3]
    real(kind=wp) :: co2(1)

    ! bicarbonate (HCO3-) concentration [mol/m^3]
    real(kind=wp) :: hco3(1)

    ! carbonate (CO3--) concentration [mol/m^3]
    real(kind=wp) :: co3(1)

    ! Omega for aragonite, i.e., the aragonite saturation state
    real(kind=wp) :: OmegaA(1)

    ! Omega for calcite, i.e., the   calcite saturation state
    real(kind=wp) :: OmegaC(1)

    ! BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
    real(kind=wp) :: BetaD(1)

    ! rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
    real(kind=wp) :: rhoSW(1)

    ! pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
    real(kind=wp) :: p(1)

    ! in-situ temperature [degrees C]
    real(kind=wp) :: tempis(1)

    ! depth converted to positive values, needed in the mocsy routine
    real(kind=wp) :: dpos(1)

    ! gas transfer velocity (piston velocity) for CO2 [m/s]
    real(kind=wp) :: kw660(1)

    ! CO2 solubility
    real(kind=wp) :: K0(1)

    ! air-to-sea flux of CO2 [mmol/m2/s]
    real(kind=wp) :: co2flux_seaicemask(1)

    ! air-to-sea flux of CO2 [mmol/m2/s]
    real(kind=wp) :: o2flux_seaicemask(1)

    ! mocsy output entire depth range
    ! --------------------------------------------------------------------------------------------
    ! ! NEW MOCSY
    ! NEW MOCSY pH on total scale
    real(kind=wp) :: ph_depth(1)

    ! NEW MOCSY oceanic partial pressure of CO2 (uatm)
    real(kind=wp) :: pco2_depth(1)

    ! NEW MOCSY oceanic fugacity of CO2 (uatm)
    real(kind=wp) :: fco2_depth(1)

    ! NEW MOCSY aqueous CO2 concentration [mol/m^3]
    real(kind=wp) :: co2_depth(1)

    ! NEW MOCSY bicarbonate (HCO3-) concentration [mol/m^3]
    real(kind=wp) :: hco3_depth(1)

    ! NEW MOCSY carbonate (CO3--) concentration [mol/m^3]
    real(kind=wp) :: co3_depth(1)

    ! NEW MOCSY Omega for aragonite, i.e., the aragonite saturation state
    real(kind=wp) :: OmegaA_depth(1)

    ! NEW MOCSY Omega for calcite, i.e., the   calcite saturation state
    real(kind=wp) :: OmegaC_depth(1)

    ! NEW MOCSY BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
    real(kind=wp) :: BetaD_depth(1)

    ! NEW DISS  stoichiometric solubility product of calcite (mol^2/kg^2)
    real(kind=wp) :: kspc_depth(1)

    ! NEW MOCSY rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
    real(kind=wp) :: rhoSW_depth(1)

    ! NEW MOCSY pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth
    ! if [db]
    real(kind=wp) :: p_depth(1)

    ! NEW MOCSY in-situ temperature [degrees C]
    real(kind=wp) :: tempis_depth(1)

    ! NEW MOCSY helper value to calculate the timesteps for the carbonate system (every 7th day)
    integer :: logfile_outfreq_7

    ! NEW MOCSY helper value to calculate the timesteps for the carbonate system (every 30th day)
    integer :: logfile_outfreq_30

    !-------------------------------------------------------------------------------

    ! Common block: Species
    real(kind=wp) :: bt, dic_molal, talk_molal

    ! Common block: Equilibrium_constants
    real(kind=wp) :: k1, k2, kw, kb, ff

    ! [umol/m2/s]
    real(kind=wp) :: FeDust

    ! [mmol/m2/s]
    real(kind=wp) :: NDust

    ! Used to calculate flux of DIC in REcoM 0 -> 1
    real(kind=wp) :: Loc_ice_conc(1)

    ! [uatm]
    real(kind=wp) :: LocAtmCO2(1)

    ! (changed it from 8 to 12)
    real(kind=wp) :: LocDiags2D(12)

    real(kind=wp) :: LocRiverDIN, LocRiverDON, LocRiverDOC, LocRiverDSi, LocRiverDIC, LocRiverAlk

    real(kind=wp) :: res_zoo2_a, res_zoo2_f
    ! grazingfluxcarbon
    real(kind=wp) :: grazingFluxcarbonzoo2

    ! Zoo3
    real(kind=wp) :: grazingFluxcarbon_mes

    ! (added to make the calcification dependent on the temperature, after Krumhardt et al.
    ! 2017/2019)
    real(kind=wp) :: PICPOCtemp

    ! (to make calcification dependent on CO2)
    real(kind=wp) :: PICPOCCO2

    ! (to make calcification dependent on N-limitation)
    real(kind=wp) :: PICPOCN

    ! (added to make the calcification dependent on nutrients (N, Fe), after Krumhardt et al.
    ! 2017/2019)
    real(kind=wp) :: calc_prod_final

    integer :: currentCO2year

end module REcoM_LocVar
!===============================================================================
! Specific declarations related to carbon isotope simulations
!-------------------------------------------------------------------------------
module REcoM_ciso
    use recom_declarations, only: wp

    implicit none
    public

    private :: wp

    save

    ! Options for carbon isotope simulations (see namelist.recom)
    ! Initial fractionation of bulk organic matter
    logical :: ciso_init = .false.

    ! Include radiocarbon (-> 31 or 38 tracers)
    logical :: ciso_14 = .false.

    ! Include organic radiocarbon (-> 38 tracers)
    logical :: ciso_organic_14 = .false.

    real(kind=wp) :: delta_co2_13 = -6.61
    real(kind=wp) :: big_delta_co2_14(3) = [0., 0., 0.]

    ! Decay constant of carbon-14
    real(kind=wp) :: lambda_14 = 3.8561e-12

    ! for revised atbox 14CO2 implementation
    logical :: atbox_spinup = .true.
    real(kind=wp) :: cosmic_14_init = 2.0 ! Initial 14C production flux (atoms / s / cm**2)

    namelist /paciso/ ciso_init, ciso_14, ciso_organic_14, &
            lambda_14, delta_co2_13, big_delta_co2_14, &
            atbox_spinup, cosmic_14_init

    ! Extensions of other modules or subroutines
    ! Module REcoM_constants: ciso tracer indices
    integer :: idic_13, iphyc_13, idetc_13, ihetc_13, idoc_13, idiac_13, iphycal_13, idetcal_13, &
            idic_14, iphyc_14, idetc_14, ihetc_14, idoc_14, idiac_14, iphycal_14, idetcal_14

    ! Module REcoM_declarations:
    ! quotas
    real(kind=wp) :: quota_13, quota_14, quota_dia_13, quota_dia_14, &
    ! reciprocal quotas
            recipQuota_13, recipQuota_14, recipQuota_dia_13, recipQuota_dia_14, &
            recipQZoo_13, recipQZoo_14

    ! zooplankton respiration fluxes
    real(kind=wp) :: HetRespFlux_13, HetRespFlux_14

    ! calcification
    real(kind=wp) :: calcification_13, calcification_14, &
            calc_loss_agg_13, calc_loss_agg_14, &
            calc_loss_gra_13, calc_loss_gra_14, &
            calc_diss_13, calc_diss_14

    ! Module REcoM_GloVar:
    ! [uatm] Atmospheric 13CO2 partial pressure. One value for the whole planet for each month
    real(kind=wp), dimension(12) :: AtmCO2_13

    ! [uatm] Atmospheric 14CO2 partial pressure. Three latitude zones for each month
    real(kind=wp), dimension(3, 12) :: AtmCO2_14

    ! [uatm] Surface ocean 13|14CO2 partial pressure
    real(kind=wp), allocatable, dimension(:) :: GloPCO2surf_13, GloPCO2surf_14

    ! [mmol/m2/day] Positive downwards
    real(kind=wp), allocatable, dimension(:) :: GloCO2flux_13, GloCO2flux_14
    real(kind=wp), allocatable, dimension(:) :: GloCO2flux_seaicemask_13, GloCO2flux_seaicemask_14
    real(kind=wp), allocatable, dimension(:) :: RiverineDOCOrig_13, RiverineDOCOrig_14, &
            RiverDOC2D_13, RiverDOC2D_14

    ! Module REcoM_LocVar:
    ! [uatm] Partial pressure of 13|14CO2 in surface layer at current 2D node
    real(kind=wp) :: pCO2surf_13(1), pCO2surf_14(1), &

    ! mocsy output: air-to-sea flux of 13|14CO2 [mol/(m^2 * s)]
            co2flux_13(1), co2flux_14(1), &

    ! air-to-sea flux of CO2 [mmol/m2/s]
            co2flux_seaicemask_13(1), co2flux_seaicemask_14(1)

    ! [uatm]
    real(kind=wp) :: LocAtmCO2_13(1), LocAtmCO2_14(1)

    real(kind=wp) :: LocRiverDOC_13, LocRiverDOC_14

    ! Subroutines REcoM_main & REcoM_extra:
    real(kind=wp) :: lat_val ! nodal latitude (of atmospheric input)

    ! Subroutine REcoM_extra:
    real(kind=wp) :: delta_co2_14 ! atmospheric Delta14CO2

    ! Subroutine REcoM_forcing:
    real(kind=wp) :: co2sat, & ! dissolved CO2 at saturation (CO2*air) [mol / m**3]
            kwco2 ! piston velocity of CO2

    ! Subroutine REcoM_sms:
    ! [mmol/m3] Dissolved Inorganic 13|14Carbon
    real(kind=wp) :: DIC_13, DIC_14, &

    ! [mmol/m3] Intracellular conc of 13|14Carbon in small phytoplankton
            PhyC_13, PhyC_14, &

    ! [mmol/m3] Conc of 13|14C in Detritus
            DetC_13, DetC_14, &

    ! [mmol/m3] Conc of 13|14C in heterotrophs
            HetC_13, HetC_14, &

    ! [mmol/m3] Extracellular Organic 13|14C conc
            EOC_13, EOC_14, &

    ! [mmol/m3] Intracellular conc of 13|14Carbon in diatoms
            DiaC_13, DiaC_14, &

    ! [mmol/m3] Conc of 13|14C in calcite of phytoplankton
            PhyCalc_13, PhyCalc_14, &

    ! [mmol/m3] Conc of 13|14C in calcite of detritus
            DetCalc_13, DetCalc_14

    ! Vertical profiles of photosynthesis rates, fesom1: 46 -> 47 in fesom2
    real(kind=wp), allocatable, dimension(:) :: Cphot_z, Cphot_dia_z

    ! Subroutine REcoM_init:
    ! auxiliary initial
    real(kind=wp), allocatable, dimension(:, :) :: delta_dic_13_init, &

    ! d|Delta13|14C
            delta_dic_14_init, &

    ! fields
            big_delta_dic_14_init

    ! Atmospheric box model (global variables):
    ! atmospheric CO2 mixing ratio (mole fraction)
    real(kind=wp), allocatable, dimension(:) :: x_co2atm_13, x_co2atm_14, &

    ! cosmogenic 14 production (mol / s)
            cosmic_14

    ! conversion factor
    real(kind=wp) :: production_rate_to_flux_14, &

    ! 13|14CO2 / 12CO2 spinup ratios
            r_atm_spinup_13, r_atm_spinup_14

    ! Specific factors related the carbon-isotopic composition
    ! Isotopic ratios
    ! atmospheric CO2
    real(kind=wp) :: r_atm_13, r_atm_14, &

    ! dissolved CO2
            r_co2s_13, r_co2s_14, &

    ! DIC in seawater
            r_dic_13, r_dic_14, &

    ! nanophytoplankton
            r_phyc_13, r_phyc_14, &

    ! diatoms
            r_diac_13, r_diac_14, &

    ! initial ratios of organic matter
            r_iorg_13 = 0.975, &
            r_iorg_14 = 0.950

    ! Fractionation factors
    ! gas transfer (kinetic fractionation,
    real(kind=wp) :: alpha_k_13 = 0.99912, &

    ! mean values for 5-21C by Zhang et al., 1995)
            alpha_k_14 = 0.99824, &

    ! dissolution of CO2 in sewater (equilibrium fractionation)
            alpha_aq_13, alpha_aq_14, &

    ! hydrolysis / dissociation of CO2 <-> DIC (equilibrium fract.)
            alpha_dic_13, alpha_dic_14, &

    ! photosynthesis of nanophytoplankton
            alpha_p_13, alpha_p_14, &

    ! photosynthesis of diatoms
            alpha_p_dia_13, alpha_p_dia_14, &

    ! calcification (Romanek et al., 1992: 1.001, 1.002)
            alpha_calc_13 = 1.000, &
            alpha_calc_14 = 1.000, &
    ! dissolution of calcite (Romanek et al., 1992: 0.999, 0.998)
            alpha_dcal_13 = 1.000, &
            alpha_dcal_14 = 1.000
    ! Radioactive decay constant of carbon-14
    ! t1/2 = 5700 years (Bé et al., 2013; recommended by Orr et al., 2017, for OMIP-BGC)
    ! if 1 year := 365.25 days:  lambda_14 = 3.8534e-12 / second
    ! if 1 year := 365.00 days:  lambda_14 = 3.8561e-12 / second
    ! if 1 year := 360    days:  lambda_14 = 3.9096e-12 / second

    ! Tracer IDs to be considered in decay calculations (oce_ale_tracer.F90)
    integer, dimension(8) :: c14_tracer_id = [1402, 1405, 1408, 1410, 1412, 1414, 1420, 1421]

contains

    subroutine recom_ciso_airsea(tempc, co3, dic)
        !   ----------------------------------------------------------------------------------
        !     Subroutine to calculate carbon-isotopic fractionation during air-sea exchange
        !   ----------------------------------------------------------------------------------
        !
        !     Input variables:
        !     tempc              lokal temperature in C
        !     co3                carbonate ion concentration
        !     dic                total carbon concentration
        !
        !     Output variables, defined in module REcoM_ciso:
        !     alpha_k_13,14      kinetic fract. factors for gas transfer
        !     alpha_aq_13,14     equilib. fract. factors for dissolution
        !     alpha_dic_13,14    equilib. fract. factors for DIC <-> CO2
        !
        !     Internal variables:
        !     epsilon_aq_13,14   equilib. fractionation for dissolution
        !     epsilon_dic_13,14  equilib. fractionation for DIC <-> CO2
        !     fco3               total carbon fraction
        !
        !     mbutzin, 2016 - 2019.

        !     Declarations
        implicit none

        real(kind=wp), intent(in) :: tempc, co3, dic
        real(kind=wp) :: epsilon_aq_13, epsilon_dic_13, fco3

        !     Calculation of carbon-isotopic fractionation factors, where
        !
        !     alpha_xy   = Rx / Ry               = fractionation factor
        !     epsilon_xy = (alpha_xy - 1) * 1000 = fractionation (in per mill)
        !     epsilon_14 = 2 * epsilon_13 => alpha_14 = 2 * alpha_13 - 1.

        !     We use parametrisations and numerical values determined for carbon-13
        !     by Zhang et al. (1995).

        !     Kinetic fractionation during gas transfer, mean values between 5 and 21C
        !     (values are defined in module REcoM_ciso)
        !     epsilon_k_13 = -0.86 => alpha_k_13 =  0.99914, alpha_k_14 =  0.99828

        !     Equilibrium fractionation during gas dissolution
        !
        epsilon_aq_13 = 0.0049 * tempc - 1.31
        alpha_aq_13 = 1. + 0.001 * epsilon_aq_13

        !     Equilibrium fractionation between DIC and CO2
        !
        !     The equilibrium fractionation between DIC and CO2 cannot be simply
        !     calculated from the fractionation factors for HCO3, CO3 and CO2star.
        !     Here, we employ an empirical function involving fCO3 = [CO3] / DIC
        !     assuming that fCO3 is the same for all carbon isotopes
        fco3 = co3 / dic
        epsilon_dic_13 = (0.014 * fco3 - 0.107) * tempc + 10.53
        alpha_dic_13 = 1. + 0.001 * epsilon_dic_13

        !     Fractionation of radiocarbon
        if (ciso_organic_14) then
            alpha_aq_14 = 2. * alpha_aq_13 - 1.
            alpha_dic_14 = 2. * alpha_dic_13 - 1.
        else
            !       no fractionation in the inorganic approximation
            alpha_aq_14 = 1.
            alpha_dic_14 = 1.
        end if

        return
    end subroutine recom_ciso_airsea

    !   ----------------------------------------------------------------------------------

    subroutine recom_ciso_photo(co2st)
        !   ----------------------------------------------------------------------------------
        !        Subroutine calculating carbon-isotopic fractionation during photosynthesis
        !   ----------------------------------------------------------------------------------
        !     Input:
        !     dissolved CO2 (co2st) in mol / m**3
        !
        !     Output:
        !     isotopic fractionation factors for phytoplankton and diatoms due to
        !     photosynthesis (alpha_p_13|14, declared at the head of the module)
        !
        !     Note that we are interested in effective values (implictly including the
        !     fractionation of dissolved CO2) which are actually derived in field studies
        !     or lab experiments. Young et al. 2013, eq. (5) with values from paragraph [35]
        !
        !     Here, we follow Young et al. 2013, eq. (5) with values from paragraph [35]
        !     eps_p = eps_pm * (1. - rho / co2aq) = 17.6 * (1 - 2.02 / co2aq)
        !     where co2aq is in umol / L
        !
        !     mbutzin, 2017 - 2021.

        implicit none
        real(kind=wp), intent(in) :: co2st
        real(kind=wp) :: co2aq

        !     Convert dissolved CO2 from mol / m**3 to umol / L and prevent from division by zero
        co2aq = max(1.d-8, co2st * 1000.)

        !     Fractionation wrt carbon-13
        alpha_p_13 = max(1., 1. + 0.001 * (17.6 * (1 - 2.02 / co2aq)))
        alpha_p_dia_13 = alpha_p_13

        !     Fractionation wrt carbon-14
        alpha_p_14 = 2. * alpha_p_13 - 1.
        alpha_p_dia_14 = 2. * alpha_p_dia_13 - 1.

        return
    end subroutine recom_ciso_photo

    !   ----------------------------------------------------------------------------------

    function lat_zone(lat_n)
        !   ----------------------------------------------------------------------------------
        !   Assign latitude zones from nodal latitude values
        !   ----------------------------------------------------------------------------------

        implicit none
        integer :: lat_zone

        !     Input: Latitude value corresponding to node n
        real(kind=wp), intent(in) :: lat_n

        !     Binning of latitudes to three zones
        if (lat_n > 30.) then ! Northern Hemisphere polewards of 30°N
            lat_zone = 1
        else if (lat_n < -30.) then ! Southern Hemisphere polewards of 30°S
            lat_zone = 3
        else ! (Sub-) Tropical zone
            lat_zone = 2
        end if

        return
    end function lat_zone

    function wind_10(windstr_x, windstr_y)
        !   ----------------------------------------------------------------------------------
        !    computes wind speed at 10 m height "wind10" from wind stress fields tau_x, tau_y
        !    as long as wind10 is not properly passed from ECHAM in coupled simulations.
        !    We follow Peixoto & Oort (1992, Eq. (10.28), (10,29)) and Charnock (1955);
        !    also see MPI report 349 (2003), Eq. (5.7).
        !   ----------------------------------------------------------------------------------
        implicit none

        real(kind=wp) :: wind_10

        !     Input
        real(kind=wp), intent(in) :: windstr_x, windstr_y

        !     Internal variables and parameters
        !     Zonal and meridional velocities at 10 m height
        real(kind=wp) :: u_10, v_10
        !     Zonal and meridional friction velocities
        real(kind=wp) :: u_fric, v_fric
        !     Zonal and meridional roughness lengths
        real(kind=wp) :: l_rough_x, l_rough_y
        !     Inverse von-Karman constant (0.4), Charnock constant (0.018) divided by g, inverse
        ! density of air (1.3), log(10)
        real(kind=wp), parameter :: inv_karm = 2.5, charn_g = 0.00173, inv_dens_air = 0.76923, &
                log_10 = 2.30258

        !     Calculate friction velocities (Peixoto & Oort, 1992, Eq. (10.28))
        u_fric = sqrt(abs(windstr_x) * inv_dens_air)
        v_fric = sqrt(abs(windstr_y) * inv_dens_air)

        !     Calculate roughness lengths (MPI report 349, 2003, Eq. (5.7), quoting Charnock, 1955)
        l_rough_x = max((charn_g * u_fric ** 2), 1.5e-5)
        l_rough_y = max((charn_g * v_fric ** 2), 1.5e-5)

        !     Calculate wind speed at 10 m (Peixoto & Oort, 1992, Eq. (10.29))
        u_10 = inv_karm * u_fric * (log_10 - log(l_rough_x))
        v_10 = inv_karm * v_fric * (log_10 - log(l_rough_y))

        wind_10 = sqrt(u_10 ** 2 + v_10 ** 2)

        return
    end function wind_10
    !   ----------------------------------------------------------------------------------

end module REcoM_ciso

module recom_diags_management
    use recom_config, only: enable_3zoo2det, enable_coccos, grazing_detritus

    implicit none

    private

    public :: allocate_and_init_diags
    public :: update_2d_diags
    public :: update_3d_diags
    public :: deallocate_diags

contains

    ! ==============================================================================
    ! SUBROUTINE: allocate_and_init_diags
    ! Purpose: Allocate and initialize all diagnostic arrays for a water column
    ! ==============================================================================
    subroutine allocate_and_init_diags(nl)

        use REcoM_declarations, only: vertaggc, vertaggd, vertaggn, vertaggp, vertcalcdiss, &
                vertcalcif, vertchldegc, vertchldegd, vertchldegn, vertchldegp, vertdocexc, &
                vertdocexd, vertdocexn, vertdocexp, vertgppc, vertgppd, vertgppn, vertgppp, &
                vertgrazmacro_c, vertgrazmacro_d, vertgrazmacro_det, vertgrazmacro_det2, &
                vertgrazmacro_mes, vertgrazmacro_mic, vertgrazmacro_n, vertgrazmacro_p, &
                vertgrazmacro_tot, vertgrazmeso_c, vertgrazmeso_d, vertgrazmeso_det, &
                vertgrazmeso_det2, vertgrazmeso_mic, vertgrazmeso_n, vertgrazmeso_p, &
                vertgrazmicro_c, vertgrazmicro_d, vertgrazmicro_n, vertgrazmicro_p, &
                vertnnac, vertnnad, vertnnan, vertnnap, vertnppc, vertnppd, vertnppn, vertnppp, &
                vertrespc, vertrespd, vertrespmacro, vertrespmeso, vertrespmicro, vertrespn, &
                vertrespp, vtcoccoco2, vtcphot_cocco, vtcphot_diatoms, vtcphot_phaeo, &
                vtcphotliglim_cocco, vtcphotliglim_diatoms, vtcphotliglim_phaeo, &
                vtdiaco2, vtphaeoco2, vtphyco2, vtqlimitfac_cocco, vtqlimitfac_diatoms, &
                vtqlimitfac_phaeo, vtqlimitfac_phyto, vtsi_assimdia, vttemp_cocco, vttemp_diatoms, &
                vertgrazmeso_tot, vertgrazmicro_tot, vttemp_phaeo, vttemp_phyto, &
                vtcphot_phyto, vtcphotliglim_phyto

        implicit none

        integer, intent(in) :: nl ! Number of vertical levels

        ! --------------------------------------------------------------------------
        ! Small Phytoplankton
        ! --------------------------------------------------------------------------
        allocate(vertNPPn(nl - 1), vertGPPn(nl - 1), vertNNAn(nl - 1), vertChldegn(nl - 1))
        allocate(vertrespn(nl - 1), vertdocexn(nl - 1), vertaggn(nl - 1))

        vertNPPn = 0.d0
        vertGPPn = 0.d0
        vertNNAn = 0.d0
        vertChldegn = 0.d0
        vertrespn = 0.d0
        vertdocexn = 0.d0
        vertaggn = 0.d0

        ! --------------------------------------------------------------------------
        ! Diatoms
        ! --------------------------------------------------------------------------
        allocate(vertNPPd(nl - 1), vertGPPd(nl - 1), vertNNAd(nl - 1), vertChldegd(nl - 1))
        allocate(vertrespd(nl - 1), vertdocexd(nl - 1), vertaggd(nl - 1))

        vertNPPd = 0.d0
        vertGPPd = 0.d0
        vertNNAd = 0.d0
        vertChldegd = 0.d0
        vertrespd = 0.d0
        vertdocexd = 0.d0
        vertaggd = 0.d0

        ! --------------------------------------------------------------------------
        ! Coccolithophores (if enabled)
        ! --------------------------------------------------------------------------
        if (enable_coccos) then
            allocate(vertNPPc(nl - 1), vertGPPc(nl - 1), vertNNAc(nl - 1), vertChldegc(nl - 1))
            allocate(vertrespc(nl - 1), vertdocexc(nl - 1), vertaggc(nl - 1))
            allocate(vertcalcdiss(nl - 1), vertcalcif(nl - 1))

            vertNPPc = 0.d0
            vertGPPc = 0.d0
            vertNNAc = 0.d0
            vertChldegc = 0.d0
            vertrespc = 0.d0
            vertdocexc = 0.d0
            vertaggc = 0.d0
            vertcalcdiss = 0.d0
            vertcalcif = 0.d0

            ! ----------------------------------------------------------------------
            ! Phaeocystis
            ! ----------------------------------------------------------------------
            allocate(vertNPPp(nl - 1), vertGPPp(nl - 1), vertNNAp(nl - 1), vertChldegp(nl - 1))
            allocate(vertrespp(nl - 1), vertdocexp(nl - 1), vertaggp(nl - 1))

            vertNPPp = 0.d0
            vertGPPp = 0.d0
            vertNNAp = 0.d0
            vertChldegp = 0.d0
            vertrespp = 0.d0
            vertdocexp = 0.d0
            vertaggp = 0.d0
        else
            ! Allocate calcification arrays even if coccos are disabled
            allocate(vertcalcdiss(nl - 1), vertcalcif(nl - 1))
            vertcalcdiss = 0.d0
            vertcalcif = 0.d0
        end if

        ! --------------------------------------------------------------------------
        ! Zooplankton Grazing (if enabled)
        ! --------------------------------------------------------------------------
        if (Grazing_detritus) then

            if (enable_3zoo2det) then
                ! Microzooplankton
                allocate(vertgrazmicro_tot(nl - 1), &
                        vertgrazmicro_n(nl - 1), &
                        vertgrazmicro_d(nl - 1))
                allocate(vertrespmicro(nl - 1))

                vertgrazmicro_tot = 0.d0
                vertgrazmicro_n = 0.d0
                vertgrazmicro_d = 0.d0
                vertrespmicro = 0.d0

                ! Mesozooplankton
                allocate(vertgrazmeso_tot(nl - 1), vertgrazmeso_n(nl - 1), vertgrazmeso_d(nl - 1))
                allocate(vertgrazmeso_det(nl - 1), vertgrazmeso_mic(nl - 1), vertgrazmeso_det2(nl &
                        &- 1))
                allocate(vertrespmeso(nl - 1))

                vertgrazmeso_tot = 0.d0
                vertgrazmeso_n = 0.d0
                vertgrazmeso_d = 0.d0
                vertgrazmeso_det = 0.d0
                vertgrazmeso_mic = 0.d0
                vertgrazmeso_det2 = 0.d0
                vertrespmeso = 0.d0

                if (enable_coccos) then
                    allocate(vertgrazmicro_c(nl - 1), vertgrazmicro_p(nl - 1))
                    allocate(vertgrazmeso_c(nl - 1), vertgrazmeso_p(nl - 1))

                    vertgrazmicro_c = 0.d0
                    vertgrazmicro_p = 0.d0
                    vertgrazmeso_c = 0.d0
                    vertgrazmeso_p = 0.d0
                end if
            end if

            ! Macrozooplankton
            allocate(vertgrazmacro_tot(nl - 1), vertgrazmacro_n(nl - 1), vertgrazmacro_d(nl - 1))
            allocate(vertgrazmacro_mes(nl - 1), vertgrazmacro_det(nl - 1))
            allocate(vertgrazmacro_mic(nl - 1), vertgrazmacro_det2(nl - 1))
            allocate(vertrespmacro(nl - 1))

            vertgrazmacro_tot = 0.d0
            vertgrazmacro_n = 0.d0
            vertgrazmacro_d = 0.d0
            vertgrazmacro_mes = 0.d0
            vertgrazmacro_det = 0.d0
            vertgrazmacro_mic = 0.d0
            vertgrazmacro_det2 = 0.d0
            vertrespmacro = 0.d0

            if (enable_coccos) then
                allocate(vertgrazmacro_c(nl - 1), vertgrazmacro_p(nl - 1))
                vertgrazmacro_c = 0.d0
                vertgrazmacro_p = 0.d0
            end if
        end if

        ! --------------------------------------------------------------------------
        ! Temperature and Photosynthesis Tracking Variables
        ! --------------------------------------------------------------------------

        allocate(VTPhyCO2(nl - 1), VTDiaCO2(nl - 1))
        VTPhyCO2 = 0.d0
        VTDiaCO2 = 0.d0

        allocate(VTCphotLigLim_phyto(nl - 1), VTCphotLigLim_diatoms(nl - 1))
        VTCphotLigLim_phyto = 0.d0
        VTCphotLigLim_diatoms = 0.d0

        allocate(VTCphot_phyto(nl - 1), VTCphot_diatoms(nl - 1))
        VTCphot_phyto = 0.d0
        VTCphot_diatoms = 0.d0

        if (enable_coccos) then

            allocate(VTTemp_diatoms(nl - 1), VTTemp_phyto(nl - 1))
            VTTemp_diatoms = 0.d0
            VTTemp_phyto = 0.d0

            allocate(VTqlimitFac_phyto(nl - 1), VTqlimitFac_diatoms(nl - 1))
            VTqlimitFac_phyto = 0.d0
            VTqlimitFac_diatoms = 0.d0

            allocate(VTSi_assimDia(nl - 1))
            VTSi_assimDia = 0.d0

            ! --------------------------------------------------------------------------
            ! Coccolithophores and Phaeocystis Tracking (if enabled)
            ! --------------------------------------------------------------------------

            allocate(VTTemp_cocco(nl - 1), VTTemp_phaeo(nl - 1))
            VTTemp_cocco = 0.d0
            VTTemp_phaeo = 0.d0

            allocate(VTCoccoCO2(nl - 1), VTPhaeoCO2(nl - 1))
            VTCoccoCO2 = 0.d0
            VTPhaeoCO2 = 0.d0

            allocate(VTqlimitFac_cocco(nl - 1), VTqlimitFac_phaeo(nl - 1))
            VTqlimitFac_cocco = 0.d0
            VTqlimitFac_phaeo = 0.d0

            allocate(VTCphotLigLim_cocco(nl - 1), VTCphotLigLim_phaeo(nl - 1))
            VTCphotLigLim_cocco = 0.d0
            VTCphotLigLim_phaeo = 0.d0

            allocate(VTCphot_cocco(nl - 1), VTCphot_phaeo(nl - 1))
            VTCphot_cocco = 0.d0
            VTCphot_phaeo = 0.d0

        end if

    end subroutine allocate_and_init_diags

    ! ==============================================================================
    ! SUBROUTINE: update_2d_diags
    ! Purpose: Transfer local diagnostic values to 2D global arrays
    ! ==============================================================================
    subroutine update_2d_diags(n)

        use recom_glovar, only: nppc, gppc, nnac, chldegc, nppp, gppp, nnap, chldegp, &
                grazmeso_tot, grazmeso_n, grazmeso_d, grazmeso_det, grazmeso_c, grazmeso_p, &
                grazmeso_mic, grazmeso_det2, grazmacro_tot, grazmacro_n, grazmacro_d, &
                grazmacro_det, grazmacro_mic, grazmacro_det2, grazmacro_c, grazmacro_p, &
                grazmicro_n, grazmicro_d, grazmicro_c, grazmicro_p, chldegd, chldegn, gppd, gppn, &
                grazmicro_tot, grazmacro_mes, nnad, nnan, nppd, nppn

        use REcoM_declarations, only: locchldegc, locchldegd, locchldegn, locchldegp, locgppc, &
                locgppd, locgppn, locgppp, locgrazmacro_c, locgrazmacro_d, locgrazmacro_det, &
                locgrazmacro_det2, locgrazmacro_mes, locgrazmacro_mic, locgrazmacro_n, &
                locgrazmacro_tot, locgrazmeso_c, locgrazmeso_d, locgrazmeso_det, locgrazmeso_det2, &
                locgrazmeso_mic, locgrazmeso_n, locgrazmeso_p, locgrazmeso_tot, locgrazmicro_c, &
                locgrazmicro_d, locgrazmicro_n, locgrazmicro_p, locgrazmicro_tot, locnnac, &
                locgrazmacro_p, locnnad, locnnan, locnnap, locnppc, locnppd, locnppn, locnppp

        implicit none

        integer, intent(in) :: n ! Node index

        ! --------------------------------------------------------------------------
        ! Small Phytoplankton
        ! --------------------------------------------------------------------------
        NPPn(n) = locNPPn
        GPPn(n) = locGPPn
        NNAn(n) = locNNAn
        Chldegn(n) = locChldegn

        ! --------------------------------------------------------------------------
        ! Diatoms
        ! --------------------------------------------------------------------------
        NPPd(n) = locNPPd
        GPPd(n) = locGPPd
        NNAd(n) = locNNAd
        Chldegd(n) = locChldegd

        ! --------------------------------------------------------------------------
        ! Coccolithophores and Phaeocystis (if enabled)
        ! --------------------------------------------------------------------------
        if (enable_coccos) then
            NPPc(n) = locNPPc
            GPPc(n) = locGPPc
            NNAc(n) = locNNAc
            Chldegc(n) = locChldegc

            NPPp(n) = locNPPp
            GPPp(n) = locGPPp
            NNAp(n) = locNNAp
            Chldegp(n) = locChldegp
        end if

        ! --------------------------------------------------------------------------
        ! Zooplankton Grazing (if enabled)
        ! --------------------------------------------------------------------------
        if (Grazing_detritus) then
            ! Mesozooplankton
            grazmeso_tot(n) = locgrazmeso_tot
            grazmeso_n(n) = locgrazmeso_n
            grazmeso_d(n) = locgrazmeso_d
            grazmeso_det(n) = locgrazmeso_det

            if (enable_coccos) then
                grazmeso_c(n) = locgrazmeso_c
                grazmeso_p(n) = locgrazmeso_p
            end if

            if (enable_3zoo2det) then
                grazmeso_mic(n) = locgrazmeso_mic
                grazmeso_det2(n) = locgrazmeso_det2

                ! Macrozooplankton
                grazmacro_tot(n) = locgrazmacro_tot
                grazmacro_n(n) = locgrazmacro_n
                grazmacro_d(n) = locgrazmacro_d
                grazmacro_mes(n) = locgrazmacro_mes
                grazmacro_det(n) = locgrazmacro_det
                grazmacro_mic(n) = locgrazmacro_mic
                grazmacro_det2(n) = locgrazmacro_det2

                if (enable_coccos) then
                    grazmacro_c(n) = locgrazmacro_c
                    grazmacro_p(n) = locgrazmacro_p
                end if

                ! Microzooplankton
                grazmicro_tot(n) = locgrazmicro_tot
                grazmicro_n(n) = locgrazmicro_n
                grazmicro_d(n) = locgrazmicro_d

                if (enable_coccos) then
                    grazmicro_c(n) = locgrazmicro_c
                    grazmicro_p(n) = locgrazmicro_p
                end if
            end if
        end if

    end subroutine update_2d_diags

    ! ==============================================================================
    ! SUBROUTINE: update_3d_diags
    ! Purpose: Transfer vertical profile diagnostic values to 3D global arrays
    ! ==============================================================================
    subroutine update_3d_diags(n, nzmax)

        use recom_glovar, only: aggn, aggd, aggc, aggp, calcdiss, calcif, &
                nppc3D, nppd3D, nppn3D, nppp3D, respn, respd, respc, respp, respmeso, &
                respmicro, tphyco2, tdiaco2, TCphotLigLim_phyto, TCphot_phyto, TCphot_diatoms, &
                TCphotLigLim_diatoms, TTemp_phyto, TqlimitFac_phyto, TTemp_diatoms, &
                TSi_assimDia, TCoccoCO2, TphaeoCO2, TTemp_cocco, TTemp_phaeo, TqlimitFac_cocco, &
                TCphotLigLim_cocco, TCphot_cocco, TqlimitFac_phaeo, TCphotLigLim_phaeo, &
                docexd, docexn, docexp, TCphot_phaeo, docexc, respmacro, TqlimitFac_diatoms

        use REcoM_declarations, only: vertaggc, vertaggd, vertaggn, vertaggp, vertcalcdiss, &
                vertrespc, vertrespd, vertrespmacro, vertrespmeso, vertrespmicro, vertrespn, &
                vertrespp, vtcoccoco2, vtcphot_cocco, vtcphot_diatoms, vtcphot_phaeo, &
                vtcphotliglim_cocco, vtcphotliglim_diatoms, vtcphotliglim_phaeo, &
                vtdiaco2, vtphaeoco2, vtphyco2, vtqlimitfac_cocco, vtqlimitfac_diatoms, &
                vtqlimitfac_phaeo, vtqlimitfac_phyto, vtsi_assimdia, vttemp_cocco, vttemp_diatoms, &
                vtcphot_phyto, vertcalcif, vertdocexc, vttemp_phaeo, vttemp_phyto, &
                vtcphotliglim_phyto, vertdocexd, vertdocexn, vertdocexp, vertnppc, vertnppd, &
                vertnppn, vertnppp

        implicit none

        integer, intent(in) :: n ! Node index
        integer, intent(in) :: nzmax ! Maximum vertical level for this node

        ! --------------------------------------------------------------------------
        ! Small Phytoplankton
        ! --------------------------------------------------------------------------
        aggn(1:nzmax, n) = vertaggn(1:nzmax)
        docexn(1:nzmax, n) = vertdocexn(1:nzmax)
        respn(1:nzmax, n) = vertrespn(1:nzmax)
        NPPn3D(1:nzmax, n) = vertNPPn(1:nzmax)

        ! --------------------------------------------------------------------------
        ! Diatoms
        ! --------------------------------------------------------------------------
        aggd(1:nzmax, n) = vertaggd(1:nzmax)
        docexd(1:nzmax, n) = vertdocexd(1:nzmax)
        respd(1:nzmax, n) = vertrespd(1:nzmax)
        NPPd3D(1:nzmax, n) = vertNPPd(1:nzmax)

        ! --------------------------------------------------------------------------
        ! Coccolithophores and Phaeocystis (if enabled)
        ! --------------------------------------------------------------------------
        if (enable_coccos) then
            aggc(1:nzmax, n) = vertaggc(1:nzmax)
            docexc(1:nzmax, n) = vertdocexc(1:nzmax)
            respc(1:nzmax, n) = vertrespc(1:nzmax)
            NPPc3D(1:nzmax, n) = vertNPPc(1:nzmax)

            aggp(1:nzmax, n) = vertaggp(1:nzmax)
            docexp(1:nzmax, n) = vertdocexp(1:nzmax)
            respp(1:nzmax, n) = vertrespp(1:nzmax)
            NPPp3D(1:nzmax, n) = vertNPPp(1:nzmax)
        end if

        ! --------------------------------------------------------------------------
        ! Calcification
        ! --------------------------------------------------------------------------
        calcdiss(1:nzmax, n) = vertcalcdiss(1:nzmax)
        calcif(1:nzmax, n) = vertcalcif(1:nzmax)

        ! --------------------------------------------------------------------------
        ! Zooplankton Respiration
        ! --------------------------------------------------------------------------
        respmeso(1:nzmax, n) = vertrespmeso(1:nzmax)

        if (enable_3zoo2det) then
            respmacro(1:nzmax, n) = vertrespmacro(1:nzmax)
            respmicro(1:nzmax, n) = vertrespmicro(1:nzmax)
        end if

        TPhyCO2(1:nzmax, n) = VTPhyCO2(1:nzmax)
        TDiaCO2(1:nzmax, n) = VTDiaCO2(1:nzmax)
        TCphotLigLim_phyto(1:nzmax, n) = VTCphotLigLim_phyto(1:nzmax)
        TCphot_phyto(1:nzmax, n) = VTCphot_phyto(1:nzmax)
        TCphotLigLim_diatoms(1:nzmax, n) = VTCphotLigLim_diatoms(1:nzmax)
        TCphot_diatoms(1:nzmax, n) = VTCphot_diatoms(1:nzmax)

        if (enable_coccos) then
            ! --------------------------------------------------------------------------
            ! Temperature and Photosynthesis Tracking - Phytoplankton
            ! --------------------------------------------------------------------------
            TTemp_phyto(1:nzmax, n) = VTTemp_phyto(1:nzmax)
            TqlimitFac_phyto(1:nzmax, n) = VTqlimitFac_phyto(1:nzmax)

            ! --------------------------------------------------------------------------
            ! Temperature and Photosynthesis Tracking - Diatoms
            ! --------------------------------------------------------------------------
            TTemp_diatoms(1:nzmax, n) = VTTemp_diatoms(1:nzmax)
            TqlimitFac_diatoms(1:nzmax, n) = VTqlimitFac_diatoms(1:nzmax)
            TSi_assimDia(1:nzmax, n) = VTSi_assimDia(1:nzmax)

            ! --------------------------------------------------------------------------
            ! Temperature and Photosynthesis Tracking - Coccos/Phaeo (if enabled)
            ! --------------------------------------------------------------------------

            TTemp_cocco(1:nzmax, n) = VTTemp_cocco(1:nzmax)
            TCoccoCO2(1:nzmax, n) = VTCoccoCO2(1:nzmax)
            TqlimitFac_cocco(1:nzmax, n) = VTqlimitFac_cocco(1:nzmax)
            TCphotLigLim_cocco(1:nzmax, n) = VTCphotLigLim_cocco(1:nzmax)
            TCphot_cocco(1:nzmax, n) = VTCphot_cocco(1:nzmax)

            TTemp_phaeo(1:nzmax, n) = VTTemp_phaeo(1:nzmax)
            TPhaeoCO2(1:nzmax, n) = VTPhaeoCO2(1:nzmax)
            TqlimitFac_phaeo(1:nzmax, n) = VTqlimitFac_phaeo(1:nzmax)
            TCphotLigLim_phaeo(1:nzmax, n) = VTCphotLigLim_phaeo(1:nzmax)
            TCphot_phaeo(1:nzmax, n) = VTCphot_phaeo(1:nzmax)
        end if

    end subroutine update_3d_diags

    ! ==============================================================================
    ! SUBROUTINE: deallocate_diags
    ! Purpose: Deallocate all diagnostic arrays
    ! ==============================================================================
    subroutine deallocate_diags()

        use REcoM_declarations, only: vertaggc, vertaggd, vertaggn, vertaggp, vertcalcdiss, &
                vertcalcif, vertchldegc, vertchldegd, vertchldegn, vertchldegp, vertdocexc, &
                vertdocexd, vertdocexn, vertdocexp, vertgppc, vertgppd, vertgppn, vertgppp, &
                vertgrazmacro_c, vertgrazmacro_d, vertgrazmacro_det, vertgrazmacro_det2, &
                vertgrazmacro_mes, vertgrazmacro_mic, vertgrazmacro_n, vertgrazmacro_p, &
                vertgrazmacro_tot, vertgrazmeso_c, vertgrazmeso_d, vertgrazmeso_det, &
                vertgrazmeso_det2, vertgrazmeso_mic, vertgrazmeso_n, vertgrazmeso_p, &
                vertgrazmicro_c, vertgrazmicro_d, vertgrazmicro_n, vertgrazmicro_p, &
                vertnnac, vertnnad, vertnnan, vertnnap, vertnppc, vertnppd, vertnppn, vertnppp, &
                vertrespc, vertrespd, vertrespmacro, vertrespmeso, vertrespmicro, vertrespn, &
                vertrespp, vtcoccoco2, vtcphot_cocco, vtcphot_diatoms, vtcphot_phaeo, &
                vtcphotliglim_cocco, vtcphotliglim_diatoms, vtcphotliglim_phaeo, &
                vtdiaco2, vtphaeoco2, vtphyco2, vtqlimitfac_cocco, vtqlimitfac_diatoms, &
                vtqlimitfac_phaeo, vtqlimitfac_phyto, vtsi_assimdia, vttemp_cocco, vttemp_diatoms, &
                vttemp_phaeo, vttemp_phyto, vtcphot_phyto, vertgrazmeso_tot, vertgrazmicro_tot, &
                vtcphotliglim_phyto

        implicit none

        ! --------------------------------------------------------------------------
        ! Small Phytoplankton
        ! --------------------------------------------------------------------------
        deallocate(vertNPPn, vertGPPn, vertNNAn, vertChldegn)
        deallocate(vertaggn, vertdocexn, vertrespn)
        deallocate(VTPhyCO2, VTCphotLigLim_phyto, VTCphot_phyto)

        ! --------------------------------------------------------------------------
        ! Diatoms
        ! --------------------------------------------------------------------------
        deallocate(vertNPPd, vertGPPd, vertNNAd, vertChldegd)
        deallocate(vertaggd, vertdocexd, vertrespd)
        deallocate(VTDiaCO2, VTCphotLigLim_diatoms, VTCphot_diatoms)

        if (enable_coccos) then
            deallocate(VTTemp_phyto, VTqlimitFac_phyto)
            deallocate(VTTemp_diatoms, VTqlimitFac_diatoms)
            deallocate(VTSi_assimDia)

            ! --------------------------------------------------------------------------
            ! Coccolithophores and Phaeocystis (if enabled)
            ! --------------------------------------------------------------------------
            deallocate(vertNPPc, vertGPPc, vertNNAc, vertChldegc)
            deallocate(vertaggc, vertdocexc, vertrespc)
            deallocate(vertcalcdiss, vertcalcif)
            deallocate(VTTemp_cocco, VTCoccoCO2, VTqlimitFac_cocco)
            deallocate(VTCphotLigLim_cocco, VTCphot_cocco)

            deallocate(vertNPPp, vertGPPp, vertNNAp, vertChldegp)
            deallocate(vertaggp, vertdocexp, vertrespp)
            deallocate(VTTemp_phaeo, VTPhaeoCO2, VTqlimitFac_phaeo)
            deallocate(VTCphotLigLim_phaeo, VTCphot_phaeo)
        else
            deallocate(vertcalcdiss, vertcalcif)
        end if

        ! --------------------------------------------------------------------------
        ! Zooplankton Grazing (if enabled)
        ! --------------------------------------------------------------------------
        if (Grazing_detritus) then
            deallocate(vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d)
            deallocate(vertgrazmeso_det, vertrespmeso)

            if (enable_coccos) then
                deallocate(vertgrazmeso_c, vertgrazmeso_p)
            end if

            if (enable_3zoo2det) then
                deallocate(vertgrazmeso_mic, vertgrazmeso_det2)

                deallocate(vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d)
                deallocate(vertgrazmacro_mes, vertgrazmacro_det)
                deallocate(vertgrazmacro_mic, vertgrazmacro_det2)
                deallocate(vertrespmacro)

                if (enable_coccos) then
                    deallocate(vertgrazmacro_c, vertgrazmacro_p)
                end if

                deallocate(vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d)
                deallocate(vertrespmicro)

                if (enable_coccos) then
                    deallocate(vertgrazmicro_c, vertgrazmicro_p)
                end if
            end if
        end if

    end subroutine deallocate_diags

end module recom_diags_management
