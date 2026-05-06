module recom_sms_module
    implicit none
    private

    public :: REcoM_sms
    public :: iron_chemistry

contains

    subroutine REcoM_sms(n, Nn, state, thick, SurfSR, sms, Temp, Sali_depth, CO2_watercolumn, &
            pH_watercolumn, pCO2_watercolumn, HCO3_watercolumn, CO3_watercolumn, &
            OmegaC_watercolumn, kspc_watercolumn, rhoSW_watercolumn, Loc_slp, &
            zF, PAR, Latd, daynew, dt, kappa, mstep, MPI_COMM_FESOM, mype, &
            myDim_nod2D, eDim_nod2D, nl, geo_coord_nod2D)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        use recom_declarations, only: wp, decaybenthos, vertrespp, vertrespc, vertrespd, &
                vertrespn, vertcalcif, vertdocexp, vertdocexc, vertdocexd, vertdocexn, vertaggp, &
                vertaggc, vertaggd, vertaggn, vertcalcdiss, vertrespmicro, vertrespmacro, &
                vertrespmeso, vertgrazmicro_p, vertgrazmicro_c, vertgrazmicro_d, vertgrazmicro_n, &
                vertgrazmicro_tot, vertgrazmacro_det, vertgrazmacro_mic, vertgrazmacro_det2, &
                vertgrazmacro_mes, vertgrazmacro_p, vertgrazmacro_c, vertgrazmacro_d, &
                vertgrazmacro_n, vertgrazmeso_det, vertgrazmeso_mic, vertgrazmeso_det2, &
                vertgrazmeso_p, vertgrazmeso_c, vertgrazmeso_d, vertgrazmeso_n, vertgrazmeso_tot, &
                vertgrazmacro_tot, aggregationrate, arrfunc, arrfunczoo2, ca, calc_diss, &
                calc_diss2, calc_diss_ben, calc_loss_agg, calc_loss_gra, calc_loss_gra2, &
                calc_loss_gra3, calcification, chl2c, chl2c_cocco, chl2c_dia, chl2c_phaeo, &
                chl2c_plast, chl2c_plast_cocco, chl2c_plast_dia, chl2c_plast_phaeo, chl2n, &
                chl2n_cocco, chl2n_dia, chl2n_phaeo, chl_lower, chl_upper, chlave, chlsynth, &
                chlsynth_cocco, chlsynth_dia, chlsynth_phaeo, co3_sat, coccoco2, &
                cphot, cphot_cocco, cphot_dia, cphot_phaeo, diaco2, fcoccon, fcoccon2, fcoccon3, &
                fdetn, fdetn2, fdetz2n, fdetz2n2, fdian, fdian2, fdian3, felimitfac, fhetn, &
                fmiczoon, fmiczoon2, food, food2, food3, foodsq, foodsq2, foodsq3, fphaeon, &
                fphaeon2, fphaeon3, fphyn, fphyn2, fphyn3, grazeff, grazingflux, grazingflux2, &
                grazingflux_cocco, grazingflux_cocco2, grazingflux_cocco3, grazingflux_det, &
                grazingflux_det2, grazingflux_detz2, grazingflux_detz22, grazingflux_dia, &
                grazingflux_dia2, grazingflux_dia3, grazingflux_het2, grazingflux_miczoo, &
                grazingflux_miczoo2, grazingflux_phaeo, grazingflux_phaeo2, grazingflux_phaeo3, &
                grazingflux_phy, grazingflux_phy2, grazingflux_phy3, hetlossflux, hetrespflux, &
                is_3zoo2det, is_coccos, kappastar, kdzlower, kdzupper, kochl, kochl_cocco, &
                kochl_phaeo, limitfacn, limitfacn_cocco, limitfacn_dia, limitfacn_phaeo, &
                lowerlight, mesfecalloss_c, mesfecalloss_n, miczoolossflux, miczoorespflux, &
                n_assim_cocco, n_assim_dia, n_assim_phaeo, o2func, parave, phaeoco2, phyco2, &
                phyresprate, phyresprate_cocco, phyresprate_dia, phyresprate_phaeo, pmax, &
                pmax_dia, pmax_phaeo, q10_mes, q10_mes_res, q10_mic, q10_mic_res, qlimitfac, &
                qlimitfactmp, qsic, qsin, quota, quota_cocco, quota_dia, quota_phaeo, &
                kochl_dia, grazingflux3, limitfacsi, n_assim, pmax_cocco, recip_res_zoo22, &
                recipbiostep, recipdet, recipdet2, recipquota, recipquota_cocco, recipquota_dia, &
                recipquota_phaeo, recipqzoo, recipqzoo2, recipqzoo3, reminsit, rtloc, rtref, &
                temp_cocco, temp_diatoms, temp_phaeo, temp_phyto, tiny_c, tiny_c_c, tiny_c_d, &
                tiny_c_p, tiny_n, tiny_n_c, tiny_n_d, tiny_n_p, tiny_si, v_cm, varpzcocco, &
                varpzcocco2, varpzcocco3, varpzdet, varpzdet2, varpzdetz2, varpzdetz22, varpzdia, &
                varpzdia2, varpzdia3, varpzhet, varpzmiczoo, varpzmiczoo2, varpzphaeo, &
                varpzphaeo3, varpzphy, varpzphy2, varpzphy3, zoo2fecalloss_c, zoo2fecalloss_n, &
                zoo2lossflux, zoo2respflux, vttemp_phyto, vttemp_diatoms, vttemp_cocco, &
                h_depth, vtphyco2, vtdiaco2, vtcoccoco2, vtphaeoco2, vtcphotliglim_phyto, &
                vtcphot_phyto, vtcphotliglim_diatoms, vtcphot_diatoms, vtcphotliglim_cocco, &
                vtcphot_cocco, vtcphotliglim_phaeo, vtcphot_phaeo, vtsi_assimdia, &
                vertnppn, vertnppd, vertnppc, vertnppp, vertgppn, vertgppd, vertgppc, vertgppc, &
                vertgppp, vertgppp, vertnnan, vertnnan, vertnnad, vertnnad, vertnnac, vertnnac, &
                vertnnap, vertnnap, vertchldegn, vertchldegn, vertchldegd, vertchldegd, &
                si_assim, varpzphaeo2, vttemp_phaeo, vertchldegc, vertchldegc, vertchldegp, &
                vertchldegp

        use recom_locvar, only: betad_depth, co2_depth, co3_depth, dpos, fco2_depth, &
                grazingfluxcarbon_mes, grazingfluxcarbonzoo2, hco3_depth, kspc_depth, locbenthos, &
                locriverdoc, logfile_outfreq_30, logfile_outfreq_7, omegaa_depth, omegac_depth, &
                pco2_depth, p_depth, ph_depth, picpocco2, picpocn, picpoctemp, res_zoo2_a, &
                res_zoo2_f, rhosw_depth, tempis_depth

        use recom_config, only: bgc_num, a_chl, a_co2_calc, a_co2_cocco, a_co2_dia, a_co2_phaeo, &
                a_co2_phy, ae, agg_pd, agg_pp, alfa, alfa_c, alfa_d, alfa_p, b_co2_calc, &
                b_co2_cocco, b_co2_dia, b_co2_phaeo, b_co2_phy, beta_phaeo, biostep, biosynth, &
                biosynthsi, c2k, c_co2_calc, c_co2_cocco, c_co2_dia, c_co2_phaeo, c_co2_phy, &
                d_co2_cocco, d_co2_dia, d_co2_phaeo, d_co2_phy, decayratebenc, decayratebenn, &
                decayratebensi, deg_chl, deg_chl_c, deg_chl_d, deg_chl_p, diags, diatom_mucus, &
                expon_phy, fe2n, fecal_rate_c, fecal_rate_c_mes, fecal_rate_n, fecal_rate_n_mes, &
                graz_max, graz_max2, graz_max3, grazeff2, grazeff3, grazing_detritus, &
                het_resp_noredfield, ialk, icchl, icocc, icocn, idchl, idetc, idetcal, idetn, &
                idetz2c, idetz2calc, idetz2n, idetz2si, idiac, idian, idiasi, idic, idin, idoc, &
                ife, ihetc, ihetn, imiczooc, imiczoon, ioxy, ipchl, iphac, iphachl, iphan, iphyc, &
                k_fe_d, k_fe_p, k_o2_remin, k_si, k_w, kscavfe, ligandstabconst, loss_het, &
                loss_miczoo, loss_zoo2, lossc, lossc_c, lossc_d, lossc_p, lossc_z, lossc_z2, &
                lossn, lossn_c, lossn_d, lossn_p, lossn_z, lossn_z2, lossn_z3, ncmax, ncmax_c, &
                ncuptakeratio_d, ncuptakeratio_p, nmaxslope, nminslope, nmocsy, o2dep_remin, &
                omegac_diss, one, ord_cocco, ord_d, ord_phy, p_cm, p_cm_d, pa2atm, pzcocco, &
                pzmiczoo2, pzphaeo, pzphaeo2, pzphaeo3, pzphy, pzphy2, pzphy3, &
                recom_grazing_variable_preference, recom_tref, redfield, redo2c, reminc, reminn, &
                reminsi, res_het, res_miczoo, res_phy, res_phy_c, res_phy_d, res_phy_p, res_zoo2, &
                rho_c1, rho_n, secondsperday, sedflx_num, sicmax, sicmin, sicuptakeratio, &
                calc_diss_exp, gfin, idetsi, idon, k_fe_c, lossc_z3, pzcocco2, simaxslope, &
                siminslope, t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2, tiny, tiny_chl, tiny_het, &
                topt_phaeo, totalligand, uopt_phaeo, use_medusa, use_photodamage, v_cm_fact, &
                v_cm_fact_c, v_cm_fact_d, v_cm_fact_p, vdet, vdet_a, zero, tmax_phaeo, &
                calc_diss_guts, calc_diss_omegac, calc_diss_rate, calc_diss_rate2, &
                calc_prod_ratio, chl2n_max, chl2n_max_c, chl2n_max_d, chl2n_max_p, ciso, co2lim, &
                cunits, d_co2_calc, enable_3zoo2det, enable_coccos, epsilon2, epsilon3, epsilonr, &
                expon_cocco, expon_d, iphycal, iphyn, isi, izoo2c, izoo2n, k_din, k_din_c, &
                k_din_d, k_din_p, k_fe, ncmax_d, ncmax_p, ncmin, ncmin_c, ncmin_d, ncmin_p, &
                ncuptakeratio, ncuptakeratio_c, pzcocco3, pzdet, pzdet2, pzdetz2, pzdetz22, &
                pzdia, pzdia2, pzdia3, pzhet, pzmiczoo

        use recoM_ciso, only: alpha_calc_13, alpha_calc_14, alpha_dcal_13, alpha_dcal_14, &
                calc_diss_13, calc_diss_14, calc_loss_agg_13, calc_loss_agg_14, calc_loss_gra_13, &
                calc_loss_gra_14, calcification_13, calcification_14, ciso_14, ciso_organic_14, &
                eoc_14, hetc_13, hetc_14, hetrespflux_13, hetrespflux_14, idetc_13, idetc_14, &
                idetcal_13, idetcal_14, idiac_13, idiac_14, idic_13, idic_14, idoc_13, idoc_14, &
                ihetc_13, ihetc_14, iphyc_13, iphyc_14, iphycal_13, iphycal_14, phyc_13, phyc_14, &
                phycalc_13, phycalc_14, quota_13, quota_14, quota_dia_13, quota_dia_14, r_iorg_13, &
                r_iorg_14, recipquota_13, recipquota_14, recipquota_dia_13, recipquota_dia_14, &
                recipqzoo_13, recipqzoo_14, detc_13, detc_14, detcalc_13, detcalc_14, diac_13, &
                diac_14, dic_13, dic_14, eoc_13

        use mvars, only: vars_sprac
        use recom_glovar, only: cosAI
        use recom_extra, only: krill_resp
        use recom_update_sms, only: update_sms

        implicit none

        integer, intent(in) :: n, daynew, mype, myDim_nod2D, eDim_nod2D, nl, mstep
        integer, intent(in) :: MPI_COMM_FESOM
        integer, intent(in) :: Nn !< Total number of nodes in the vertical

        real(kind=wp), intent(in) :: SurfSR, dt !< [W/m2] ShortWave radiation at surface
        real(kind=wp), intent(in) :: Loc_slp ![Pa] sea-level pressure
        real(kind=wp), intent(in) :: Latd(1) ! latitude in degree

        !< [m] Vertical distance between two nodes = Thickness
        real(kind=wp), intent(in), dimension(nl - 1) :: thick
        real(kind=wp), intent(in), dimension(nl - 1) :: Temp !< [degrees C] Ocean temperature

        !< NEW MOCSY Salinity for the whole water column
        real(kind=wp), intent(in), dimension(nl - 1) :: Sali_depth
        real(kind=wp), intent(in), dimension(nl) :: zF !< [m] Depth of fluxes
        real(kind=WP), intent(in), dimension(:, :) :: geo_coord_nod2D

        real(kind=wp), intent(inout) :: kappa
        real(kind=wp), intent(inout), dimension(nl - 1) :: CO2_watercolumn !< [mol/m3]
        real(kind=wp), intent(inout), dimension(nl - 1) :: pH_watercolumn !< on total scale
        real(kind=wp), intent(inout), dimension(nl - 1) :: pCO2_watercolumn !< [uatm]
        real(kind=wp), intent(inout), dimension(nl - 1) :: HCO3_watercolumn !< [mol/m3]
        real(kind=wp), intent(inout), dimension(nl - 1) :: CO3_watercolumn !< [mol/m3]

        !< calcite saturation state
        real(kind=wp), intent(inout), dimension(nl - 1) :: OmegaC_watercolumn

        !< stoichiometric solubility product [mol^2/kg^2]
        real(kind=wp), intent(inout), dimension(nl - 1) :: kspc_watercolumn

        !< in-situ density of seawater [kg/m3]
        real(kind=wp), intent(inout), dimension(nl - 1) :: rhoSW_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: PAR

        !< ChlA conc in phytoplankton [mg/m3]
        real(kind=wp), intent(inout), dimension(nl - 1, bgc_num) :: state
        real(kind=wp), intent(inout), dimension(nl - 1, bgc_num) :: sms !< Source-Minus-Sinks term

        real(kind=wp) :: dt_d !< Size of time steps [day]
        real(kind=wp) :: dt_b !< Size of time steps [day]

        real(kind=wp) :: recip_hetN_plus !< MB's addition to heterotrophic respiration

        !< [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)
        real(kind=wp) :: recip_res_het
        real(kind=wp) :: Sink_Vel
        real(kind=wp) :: aux
        integer :: k, step

        real(kind=wp) :: Patm_depth(1)

        ! MOCSY temperature for the whole water column for mocsy minimum defined as -2
        real(kind=wp) :: REcoM_T_depth(1)
        real(kind=wp) :: REcoM_S_depth(1)
        real(kind=wp) :: REcoM_DIC_depth(1)
        real(kind=wp) :: REcoM_Alk_depth(1)
        real(kind=wp) :: REcoM_Si_depth(1)
        real(kind=wp) :: REcoM_Phos_depth(1)
        real(kind=wp) :: mocsy_step_per_day

        ! --- Biogeochemical state variables ---
        real(kind=wp) :: &
                DIN, & ! [mmol/m3] Dissolved inorganic nitrogen
                DIC, & ! [mmol/m3] Dissolved inorganic carbon
                Alk, & ! [mmol/m3] Total alkalinity
                PhyN, & ! [mmol/m3] Phytoplankton nitrogen (small)
                PhyC, & ! [mmol/m3] Phytoplankton carbon (small)
                PhyChl, & ! [mg/m3] Phytoplankton chlorophyll
                DetN, & ! [mmol/m3] Detrital nitrogen
                DetC, & ! [mmol/m3] Detrital carbon
                HetN, & ! [mmol/m3] Heterotroph nitrogen
                HetC, & ! [mmol/m3] Heterotroph carbon
                DON, & ! [mmol/m3] Dissolved organic nitrogen
                EOC, & ! [mmol/m3] Extracellular organic carbon
                DiaN, & ! [mmol/m3] Diatom nitrogen
                DiaC, & ! [mmol/m3] Diatom carbon
                DiaChl, & ! [mg/m3] Diatom chlorophyll
                DiaSi, & ! [mmol/m3] Diatom silicate
                DetSi, & ! [mmol/m3] Detrital silicate
                Si, & ! [mmol/m3] Dissolved silicate
                Fe, & ! [mmol/m3] Dissolved iron
                PhyCalc, & ! [mmol/m3] Phytoplankton calcite
                DetCalc, & ! [mmol/m3] Detrital calcite
                FreeFe, & ! [mmol/m3] Free iron
                O2 ! [mmol/m3] Dissolved oxygen

        ! Coccolithophore variables (conditionally used based on namelist)
        real(kind=wp) :: &
                CoccoN, & ! [mmol/m3] Coccolithophore nitrogen
                CoccoC, & ! [mmol/m3] Coccolithophore carbon
                CoccoChl, & ! [mg/m3] Coccolithophore chlorophyll
                PhaeoN, & ! [mmol/m3] Phaeocystis nitrogen
                PhaeoC, & ! [mmol/m3] Phaeocystis carbon
                PhaeoChl ! [mg/m3] Phaeocystis chlorophyll

        ! Extended zooplankton variables (conditionally used based on namelist)
        real(kind=wp) :: &
                Zoo2N, & ! [mmol/m3] Zooplankton type 2 nitrogen
                Zoo2C, & ! [mmol/m3] Zooplankton type 2 carbon
                DetZ2N, & ! [mmol/m3] Zooplankton detritus nitrogen
                DetZ2C, & ! [mmol/m3] Zooplankton detritus carbon
                DetZ2Si, & ! [mmol/m3] Zooplankton detritus silicate
                DetZ2Calc, & ! [mmol/m3] Zooplankton detritus calcite
                MicZooN, & ! [mmol/m3] Microzooplankton nitrogen
                MicZooC ! [mmol/m3] Microzooplankton carbon

        ! ===========================================================================
        ! VARIABLE DECLARATIONS AND INITIALIZATION
        ! ===========================================================================

        !===============================================================================
        ! INITIALIZATION AND SETUP
        !===============================================================================
        ! Initializes arrays, calculates minimum thresholds, and sets up time stepping.
        !
        ! Variables:
        !   sms             : Source-minus-sink array for tracer updates [mmol m-3 day-1]
        !   zero            : Double precision zero (0.0d0)
        !   tiny            : Very small positive number (numerical stability)
        !   tiny_chl        : Minimum chlorophyll concentration [mgChl m-3]
        !
        ! Source-Minus-Sink (SMS) Array:
        !   - Accumulates all biogeochemical fluxes
        !   - Dimension: (vertical_levels, number_of_tracers)
        !   - Updated throughout biological calculations
        !   - Applied to state variables at end of time step
        !   - Units: concentration per day [mmol m-3 day-1]
        !-------------------------------------------------------------------------------

        ! Initialize source-minus-sink array
        ! Sets all biogeochemical fluxes to zero at start of time step
        sms = zero

        !===============================================================================
        ! MINIMUM THRESHOLD VALUES
        !===============================================================================
        ! Calculates minimum allowable concentrations for all biological state variables.
        ! Prevents division by zero and ensures numerical stability.
        !
        ! Threshold Calculation Strategy:
        !   - Based on physiological maximum ratios
        !   - Works backward from minimum chlorophyll
        !   - Ensures stoichiometric consistency
        !   - Species-specific values
        !
        ! Rationale for Minimum Thresholds:
        !   - Division by zero prevention in quota calculations
        !   - Numerical stability in resource limitation terms
        !   - Prevents spurious negative values
        !   - Represents detection limits or "ghost populations"
        !
        ! Typical Minimum Values:
        !   - Chlorophyll: ~0.001-0.01 mgChl m-3
        !   - Nitrogen: ~0.001-0.01 mmolN m-3
        !   - Carbon: ~0.01-0.1 mmolC m-3
        !-------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------
        ! Small Phytoplankton Thresholds
        !-------------------------------------------------------------------------------
        ! Variables:
        !   tiny_N      : Minimum small phyto nitrogen [mmolN m-3]
        !   tiny_C      : Minimum small phyto carbon [mmolC m-3]
        !   tiny_chl    : Minimum chlorophyll (set externally) [mgChl m-3]
        !   chl2N_max   : Maximum Chlorophyll:Nitrogen ratio [mgChl mmolN-1]
        !   NCmax       : Maximum Nitrogen:Carbon quota [mmolN mmolC-1]
        !
        ! Calculation Logic:
        !   1. Start with minimum observable Chl (tiny_chl)
        !   2. Calculate minimum N using maximum Chl:N ratio
        !   3. Calculate minimum C using maximum N:C quota
        !
        ! Typical Values:
        !   chl2N_max = 3.15 mgChl/mmolN (high Chl per N, low light acclimation)
        !   NCmax = 0.2 mmolN/mmolC (luxury N consumption maximum)

        ! Minimum nitrogen based on minimum chlorophyll and maximum Chl:N ratio
        tiny_N = tiny_chl / chl2N_max

        ! Minimum carbon based on minimum nitrogen and maximum N:C quota
        tiny_C = tiny_N / NCmax

        !-------------------------------------------------------------------------------
        ! Diatom Thresholds
        !-------------------------------------------------------------------------------
        ! Variables:
        !   tiny_N_d    : Minimum diatom nitrogen [mmolN m-3]
        !   tiny_C_d    : Minimum diatom carbon [mmolC m-3]
        !   tiny_Si     : Minimum diatom silicate [mmolSi m-3]
        !   chl2N_max_d : Maximum diatom Chl:N ratio [mgChl mmolN-1]
        !   NCmax_d     : Maximum diatom N:C quota [mmolN mmolC-1]
        !   SiCmax      : Maximum diatom Si:C quota [mmolSi mmolC-1]
        !
        ! Typical Values:
        !   chl2N_max_d = 4.2 mgChl/mmolN (diatoms can have higher Chl:N)
        !   NCmax_d = 0.2 mmolN/mmolC
        !   SiCmax = 0.8 mmolSi/mmolC (heavily silicified frustules)
        !
        ! Silicon Requirement:
        !   - Unique to diatoms (frustule formation)
        !   - Calculated from minimum carbon and maximum Si:C ratio

        ! Minimum diatom nitrogen
        tiny_N_d = tiny_chl / chl2N_max_d

        ! Minimum diatom carbon
        tiny_C_d = tiny_N_d / NCmax_d

        ! Minimum silicate (based on diatom carbon and maximum Si:C quota)
        tiny_Si = tiny_C_d / SiCmax

        !-------------------------------------------------------------------------------
        ! Coccolithophore and Phaeocystis Thresholds (Optional)
        !-------------------------------------------------------------------------------
        ! Only calculated when 4-plankton functional type model is enabled

        if (enable_coccos) then

            ! Coccolithophore thresholds
            ! Variables:
            !   tiny_N_c    : Minimum cocco nitrogen [mmolN m-3]
            !   tiny_C_c    : Minimum cocco carbon [mmolC m-3]
            !   chl2N_max_c : Maximum cocco Chl:N ratio [mgChl mmolN-1]
            !   NCmax_c     : Maximum cocco N:C quota [mmolN mmolC-1]

            tiny_N_c = tiny_chl / chl2N_max_c
            tiny_C_c = tiny_N_c / NCmax_c

            ! Phaeocystis thresholds
            ! Variables:
            !   tiny_N_p    : Minimum Phaeo nitrogen [mmolN m-3]
            !   tiny_C_p    : Minimum Phaeo carbon [mmolC m-3]
            !   chl2N_max_p : Maximum Phaeo Chl:N ratio [mgChl mmolN-1]
            !   NCmax_p     : Maximum Phaeo N:C quota [mmolN mmolC-1]

            tiny_N_p = tiny_chl / chl2N_max_p
            tiny_C_p = tiny_N_p / NCmax_p

        end if

        ! Reciprocal of heterotroph respiration rate
        ! Used in Redfield-based respiration calculations
        ! Variable:
        !   recip_res_het : Reciprocal respiration parameter [day]
        !   res_het       : Heterotroph respiration rate [day-1]
        ! Typical value: res_het = 0.01 day-1 (combined respiration + mortality)
        recip_res_het = 1.d0 / res_het

        ! Atmospheric pressure conversion for carbonate chemistry
        ! MOCSY requires pressure in atmospheres
        ! Variables:
        !   Patm_depth : Atmospheric pressure [atm]
        !   Loc_slp    : Local sea level pressure [Pa]
        !   Pa2atm     : Pascal to atmosphere conversion factor [Pa atm-1]
        Patm_depth = Loc_slp / Pa2atm

        !===============================================================================
        ! TIME STEPPING SETUP
        !===============================================================================
        ! Converts time steps and sets up sub-cycling for biogeochemistry.
        !
        ! Time Step Hierarchy:
        !   1. FESOM physics time step (dt) [seconds]
        !   2. Daily time step (dt_d) [days]
        !   3. REcoM biogeochemistry sub-step (dt_b) [days]
        !
        ! Sub-Cycling Rationale:
        !   - Physics: Slow (hours to days)
        !   - Fast biology: Minutes to hours (phytoplankton growth, grazing)
        !   - Allows finer temporal resolution for biological processes
        !   - Improves numerical stability for stiff biological systems
        !
        ! Variables:
        !   rTref          : Reciprocal reference temperature [K-1]
        !   recom_Tref     : Reference temperature (typically 288.15 K = 15degC) [K]
        !   dt             : FESOM physics time step [seconds]
        !   dt_d           : Physics time step in days [days]
        !   dt_b           : REcoM biogeochemistry sub-time step [days]
        !   SecondsPerDay  : Conversion factor (86400) [s day-1]
        !   biostep        : Number of biogeochemistry steps per physics step [-]
        !
        ! Typical Configuration:
        !   dt = 3600 s (1 hour physics step)
        !   biostep = 4 (15-minute biogeochemistry steps)
        !   dt_b = 0.25/24 = 0.0104 days
        !-------------------------------------------------------------------------------

        ! Convert FESOM physics time step to days
        dt_d = dt / SecondsPerDay

        ! Calculate biogeochemistry sub-time step
        ! Divides physics step into smaller biogeochemical steps
        dt_b = dt_d / real(biostep)

        ! Reciprocal reference temperature for Arrhenius calculations
        ! Reference: 288.15 K (15degC) used in temperature dependence functions
        rTref = real(one) / recom_Tref

        !===============================================================================
        ! MAIN TIME INTEGRATION LOOP
        !===============================================================================
        ! Iterates through biogeochemical sub-time steps within each physics time step.
        ! Allows fine temporal resolution for fast biological processes.
        !
        ! Loop Structure:
        !   - Outer loop: Biogeochemistry sub-steps (biostep iterations)
        !   - Inner loop: Vertical layers (surface to bottom)
        !
        ! Variables:
        !   step     : Current biogeochemistry sub-step counter [-]
        !   biostep  : Total number of sub-steps per physics step [-]
        !   kdzUpper : Cumulative light attenuation from surface [dimensionless]
        !   sms      : Source-minus-sink flux array [mmol m-3 day-1]
        !   tiny     : Threshold for negligible fluxes [mmol m-3 day-1]
        !
        ! Numerical Cleanup:
        !   - Removes negligible SMS values before applying
        !   - Prevents accumulation of numerical noise
        !   - Improves computational efficiency
        !-------------------------------------------------------------------------------

        do step = one, biostep

            ! Reset upper light attenuation (top of the cell) at start of each sub-step
            ! Light attenuation integrates downward through water column
            kdzUpper = 0.d0

            ! Clean up negligible SMS values to prevent numerical issues
            ! Sets very small fluxes to exactly zero
            if (any(abs(sms(:, :)) <= tiny)) sms(:, :) = zero

            !===========================================================================
            ! VERTICAL LOOP THROUGH WATER COLUMN
            !===========================================================================
            ! Processes each vertical layer from surface to bottom.
            ! Updates state variables and enforces constraints.
            !
            ! Variables:
            !   k  : Current vertical level index [-]
            !   nn : Number of vertical levels in water column [-]
            !
            ! Note: Alternative loop configurations commented out:
            !   - nzmin, nzmax: Min/max active levels (for dry cells)
            !   - myDim_nod2D: Horizontal dimension (for 3D unstructured grids)
            !---------------------------------------------------------------------------

            do k = one, nn
                ! Alternative loop structures (commented out):
                ! do n = 1, myDim_nod2D
                !     Nn = nlevels_nod2D(n) - 1
                !     nzmin = ulevels_nod2D(row)
                !     nzmax = nlevels_nod2D(row)

                !-----------------------------------------------------------------------
                ! DISSOLVED INORGANIC NUTRIENTS
                !-----------------------------------------------------------------------
                ! Update nutrient concentrations: state(previous) + SMS(fluxes)
                ! Enforce minimum values for numerical stability
                !
                ! Variables:
                !   DIN : Dissolved inorganic nitrogen (NO3- + NH4+) [mmolN m-3]
                !   Si  : Dissolved silicate (Si(OH)4) [mmolSi m-3]
                !   Fe  : Dissolved iron (bioavailable Fe) [mmolFe m-3]
                !
                ! max() function ensures non-negative concentrations
                !-----------------------------------------------------------------------

                DIN = max(tiny, state(k, idin) + sms(k, idin))
                Si = max(tiny, state(k, isi) + sms(k, isi))
                Fe = max(tiny, state(k, ife) + sms(k, ife))

                !-----------------------------------------------------------------------
                ! CARBON SYSTEM VARIABLES
                !-----------------------------------------------------------------------
                ! Variables:
                !   DIC : Dissolved inorganic carbon (CO2 + HCO3- + CO3--) [mmolC m-3]
                !   ALK : Total alkalinity [meq m-3]
                !   O2  : Dissolved oxygen [mmolO2 m-3]
                !-----------------------------------------------------------------------

                DIC = max(tiny, state(k, idic) + sms(k, idic))
                ALK = max(tiny, state(k, ialk) + sms(k, ialk))
                O2 = max(tiny, state(k, ioxy) + sms(k, ioxy))

                !-----------------------------------------------------------------------
                ! DISSOLVED ORGANIC MATTER
                !-----------------------------------------------------------------------
                ! Variables:
                !   DON : Dissolved organic nitrogen (labile + semi-labile) [mmolN m-3]
                !   EOC : Dissolved organic carbon (labile + semi-labile) [mmolC m-3]
                !
                ! Note: EOC naming convention (Enhanced Organic Carbon) is historical
                !-----------------------------------------------------------------------

                DON = max(tiny, state(k, idon) + sms(k, idon))
                EOC = max(tiny, state(k, idoc) + sms(k, idoc))

                !-----------------------------------------------------------------------
                ! SMALL PHYTOPLANKTON
                !-----------------------------------------------------------------------
                ! General phytoplankton functional type
                ! Variables:
                !   PhyN    : Small phyto nitrogen [mmolN m-3]
                !   PhyC    : Small phyto carbon [mmolC m-3]
                !   PhyChl  : Small phyto chlorophyll [mgChl m-3]
                !   PhyCalc : Small phyto calcite (if calcifying) [mmolC m-3]
                !-----------------------------------------------------------------------

                PhyN = max(tiny_N, state(k, iphyn) + sms(k, iphyn))
                PhyC = max(tiny_C, state(k, iphyc) + sms(k, iphyc))
                PhyChl = max(tiny_chl, state(k, ipchl) + sms(k, ipchl))
                PhyCalc = max(tiny, state(k, iphycal) + sms(k, iphycal))

                !-----------------------------------------------------------------------
                ! DIATOMS (SILICIFYING PHYTOPLANKTON)
                !-----------------------------------------------------------------------
                ! Large phytoplankton with silica frustules
                ! Variables:
                !   DiaN   : Diatom nitrogen [mmolN m-3]
                !   DiaC   : Diatom carbon [mmolC m-3]
                !   DiaChl : Diatom chlorophyll [mgChl m-3]
                !   DiaSi  : Diatom silicate (frustule) [mmolSi m-3]
                !-----------------------------------------------------------------------

                DiaN = max(tiny_N_d, state(k, idian) + sms(k, idian))
                DiaC = max(tiny_C_d, state(k, idiac) + sms(k, idiac))
                DiaChl = max(tiny_chl, state(k, idchl) + sms(k, idchl))
                DiaSi = max(tiny_si, state(k, idiasi) + sms(k, idiasi))

                if (enable_coccos) then

                    !-------------------------------------------------------------------
                    ! COCCOLITHOPHORES (CALCIFYING PHYTOPLANKTON)
                    !-------------------------------------------------------------------
                    ! Variables:
                    !   CoccoN   : Cocco nitrogen [mmolN m-3]
                    !   CoccoC   : Cocco carbon [mmolC m-3]
                    !   CoccoChl : Cocco chlorophyll [mgChl m-3]
                    !-------------------------------------------------------------------

                    CoccoN = max(tiny_N_c, state(k, icocn) + sms(k, icocn))
                    CoccoC = max(tiny_C_c, state(k, icocc) + sms(k, icocc))
                    CoccoChl = max(tiny_chl, state(k, icchl) + sms(k, icchl))

                    !-------------------------------------------------------------------
                    ! PHAEOCYSTIS (COLONIAL PHYTOPLANKTON)
                    !-------------------------------------------------------------------
                    ! Variables:
                    !   PhaeoN   : Phaeo nitrogen [mmolN m-3]
                    !   PhaeoC   : Phaeo carbon [mmolC m-3]
                    !   PhaeoChl : Phaeo chlorophyll [mgChl m-3]
                    !-------------------------------------------------------------------

                    PhaeoN = max(tiny_N_p, state(k, iphan) + sms(k, iphan))
                    PhaeoC = max(tiny_C_p, state(k, iphac) + sms(k, iphac))
                    PhaeoChl = max(tiny_chl, state(k, iphachl) + sms(k, iphachl))

                end if

                !-----------------------------------------------------------------------
                ! HETEROTROPHS (ZOOPLANKTON)
                !-----------------------------------------------------------------------
                ! Primary grazers (mesozooplankton)
                ! Variables:
                !   HetN : Mesozooplankton nitrogen [mmolN m-3]
                !   HetC : Mesozooplankton carbon [mmolC m-3]
                !-----------------------------------------------------------------------

                HetN = max(tiny, state(k, ihetn) + sms(k, ihetn))
                HetC = max(tiny, state(k, ihetc) + sms(k, ihetc))

                if (enable_3zoo2det) then

                    !-------------------------------------------------------------------
                    ! ADDITIONAL ZOOPLANKTON (3-ZOOPLANKTON MODEL)
                    !-------------------------------------------------------------------
                    ! Macrozooplankton (e.g., krill)
                    ! Variables:
                    !   Zoo2N : Macrozooplankton nitrogen [mmolN m-3]
                    !   Zoo2C : Macrozooplankton carbon [mmolC m-3]
                    Zoo2N = max(tiny, state(k, izoo2n) + sms(k, izoo2n))
                    Zoo2C = max(tiny, state(k, izoo2c) + sms(k, izoo2c))

                    ! Microzooplankton (e.g., ciliates, heterotrophic dinoflagellates)
                    ! Variables:
                    !   MicZooN : Microzooplankton nitrogen [mmolN m-3]
                    !   MicZooC : Microzooplankton carbon [mmolC m-3]
                    MicZooN = max(tiny, state(k, imiczoon) + sms(k, imiczoon))
                    MicZooC = max(tiny, state(k, imiczooc) + sms(k, imiczooc))

                end if

                !-----------------------------------------------------------------------
                ! DETRITUS (DEAD ORGANIC MATTER)
                !-----------------------------------------------------------------------
                ! Slow-sinking detritus pools
                ! Variables:
                !   DetN    : Detrital nitrogen [mmolN m-3]
                !   DetC    : Detrital carbon [mmolC m-3]
                !   DetSi   : Detrital silicate [mmolSi m-3]
                !   DetCalc : Detrital calcite [mmolC m-3]
                !-----------------------------------------------------------------------

                DetN = max(tiny, state(k, idetn) + sms(k, idetn))
                DetC = max(tiny, state(k, idetc) + sms(k, idetc))
                DetSi = max(tiny, state(k, idetsi) + sms(k, idetsi))
                DetCalc = max(tiny, state(k, idetcal) + sms(k, idetcal))

                if (enable_3zoo2det) then

                    !-------------------------------------------------------------------
                    ! FAST-SINKING DETRITUS (FECAL PELLETS)
                    !-------------------------------------------------------------------
                    ! Large, rapidly sinking particles
                    ! Variables:
                    !   DetZ2N    : Fast detritus nitrogen [mmolN m-3]
                    !   DetZ2C    : Fast detritus carbon [mmolC m-3]
                    !   DetZ2Si   : Fast detritus silicate [mmolSi m-3]
                    !   DetZ2Calc : Fast detritus calcite [mmolC m-3]
                    !-------------------------------------------------------------------

                    DetZ2N = max(tiny, state(k, idetz2n) + sms(k, idetz2n))
                    DetZ2C = max(tiny, state(k, idetz2c) + sms(k, idetz2c))
                    DetZ2Si = max(tiny, state(k, idetz2si) + sms(k, idetz2si))
                    DetZ2Calc = max(tiny, state(k, idetz2calc) + sms(k, idetz2calc))

                end if

                !-----------------------------------------------------------------------
                ! FREE IRON INITIALIZATION
                !-----------------------------------------------------------------------
                ! Free iron will be calculated later from total iron budget
                ! Accounts for scavenging, complexation, and biological uptake
                !-----------------------------------------------------------------------

                FreeFe = zero

                !=======================================================================
                ! PHYSICAL ENVIRONMENT CONSTRAINTS FOR CARBONATE CHEMISTRY
                !=======================================================================
                ! Validates and constrains physical parameters for MOCSY carbonate
                ! system calculations. Ensures inputs are within valid ranges.
                !
                ! MOCSY Valid Ranges (Lueker K1/K2 formulation):
                !   Temperature: 2-35degC
                !   Salinity: 19-43 psu
                !
                ! Rationale for Constraints:
                !   - Equilibrium constants are empirical fits
                !   - Extrapolation outside valid range introduces errors
                !   - Numerical instability at extreme values
                !   - Ice formation creates low-salinity issues
                !-----------------------------------------------------------------------

                !-----------------------------------------------------------------------
                ! Temperature Constraints
                !-----------------------------------------------------------------------
                ! Variables:
                !   REcoM_T_depth : Constrained temperature for MOCSY [degC]
                !   Temp(k)       : Actual temperature at depth k [degC]
                !
                ! Constraints:
                !   Minimum: 2degC (prevents extrapolation below calibration range)
                !   Maximum: 40degC (safety limit, ocean rarely exceeds 35degC)
                !
                ! Note: K1/K2 Lueker formulation valid for 2-35degC
                !-----------------------------------------------------------------------

                REcoM_T_depth = max(2.d0, Temp(k)) ! Apply minimum
                REcoM_T_depth = min(REcoM_T_depth, 40.d0) ! Apply maximum

                !-----------------------------------------------------------------------
                ! Salinity Constraints
                !-----------------------------------------------------------------------
                ! Variables:
                !   REcoM_S_depth : Constrained salinity for MOCSY [psu]
                !   Sali_depth(k) : Actual salinity at depth k [psu]
                !
                ! Constraints:
                !   Minimum: 21 psu (increased from 19 to avoid numerical issues)
                !   Maximum: 43 psu (upper limit of calibration range)
                !
                ! Problematic Conditions:
                !   - Salinity 19-21 psu with ice concentration > 97%
                !   - Causes numerical instability in MOCSY
                !   - Conservative constraint (21 psu minimum) prevents issues
                !
                ! Note: Brackish water and ice-covered regions require special care
                !-----------------------------------------------------------------------

                REcoM_S_depth = max(21.d0, Sali_depth(k)) ! Apply minimum
                REcoM_S_depth = min(REcoM_S_depth, 43.d0) ! Apply maximum

                !-----------------------------------------------------------------------
                ! Unit Conversions for MOCSY
                !-----------------------------------------------------------------------
                ! MOCSY requires concentrations in mol/m3 (not mmol/m3)
                ! Conversion factor: 1e-3 (mmol -> mol)
                !
                ! Variables (output):
                !   REcoM_DIC_depth  : DIC for MOCSY [mol m-3]
                !   REcoM_Alk_depth  : Alkalinity for MOCSY [mol m-3]
                !   REcoM_Si_depth   : Silicate for MOCSY [mol m-3]
                !   REcoM_Phos_depth : Phosphate for MOCSY [mol m-3]
                !
                ! Sources:
                !   state(k, idic) + sms(k, idic) : DIC [mmol m-3]
                !   state(k, ialk) + sms(k, ialk) : Alkalinity [mmol m-3]
                !   state(k, isi)  + sms(k, isi)  : Silicate [mmol m-3]
                !   state(k, idin) + sms(k, idin) : DIN -> Phosphate via Redfield
                !-----------------------------------------------------------------------

                ! Dissolved inorganic carbon [mol m-3]
                REcoM_DIC_depth = max(tiny * 1e-3, state(k, idic) * 1e-3 + sms(k, idic) * 1e-3)

                ! Total alkalinity [mol m-3]
                REcoM_Alk_depth = max(tiny * 1e-3, state(k, ialk) * 1e-3 + sms(k, ialk) * 1e-3)

                ! Silicate [mol m-3]
                REcoM_Si_depth = max(tiny * 1e-3, state(k, isi) * 1e-3 + sms(k, isi) * 1e-3)

                ! Phosphate [mol m-3]
                ! Calculated from nitrogen using Redfield ratio (N:P = 16:1)
                ! Model tracks nitrogen but MOCSY needs phosphate
                REcoM_Phos_depth = max(tiny * 1e-3, state(k, idin) * 1e-3 + sms(k, idin) * 1e-3) &
                        / 16.d0

                ! ===================================================================
                ! CELLULAR QUOTAS AND RATIOS CALCULATIONS
                ! ===================================================================

                !===============================================================================
                ! Small Phytoplankton Quotas
                !===============================================================================
                ! Calculates stoichiometric ratios for the small phytoplankton functional type.
                ! Represents diverse group of small flagellates and picophytoplankton.
                !
                ! Variables:
                !   quota           : Nitrogen:Carbon quota [mmolN mmolC-1]
                !   recipquota      : Carbon:Nitrogen ratio [mmolC mmolN-1]
                !   Chl2C           : Chlorophyll:Carbon ratio [mgChl mmolC-1]
                !   Chl2N           : Chlorophyll:Nitrogen ratio [mgChl mmolN-1]
                !   CHL2C_plast     : Plastidic Chlorophyll:Carbon ratio [mgChl mmolC-1]
                !   PhyN            : Small phytoplankton nitrogen [mmolN m-3]
                !   PhyC            : Small phytoplankton carbon [mmolC m-3]
                !   PhyChl          : Small phytoplankton chlorophyll [mgChl m-3]
                !   NCmin           : Minimum N:C quota (subsistence quota) [mmolN mmolC-1]
                !
                ! Quota Interpretation:
                !   - High quota (N:C > 0.15): Nutrient replete, luxury consumption
                !   - Medium quota (N:C ≈ 0.10): Balanced growth
                !   - Low quota (N:C < 0.06): Severely N-limited, near subsistence
                !   - Minimum quota (NCmin ≈ 0.04): Zero growth threshold
                !
                ! Plastidic Chlorophyll Concept:
                !   - Total Chl includes storage and structural chlorophyll
                !   - Plastidic Chl represents functional photosynthetic apparatus
                !   - Correction factor: quota/(quota - NCmin)
                !   - Higher correction when quota approaches minimum (more Chl in chloroplasts)
                !-------------------------------------------------------------------------------

                ! Nitrogen:Carbon quota (cellular N:C ratio)
                ! Controls growth rate via Droop limitation
                quota = PhyN / PhyC

                ! Carbon:Nitrogen ratio (reciprocal)
                ! Used for converting N-based fluxes to carbon
                recipquota = real(one) / quota

                ! Chlorophyll:Carbon ratio
                ! Reflects photoacclimation state (higher in low light)
                Chl2C = PhyChl / PhyC

                ! Chlorophyll:Nitrogen ratio
                ! Links photosynthetic machinery to nitrogen investment
                Chl2N = PhyChl / PhyN

                ! Plastidic Chlorophyll:Carbon ratio
                ! Estimates chlorophyll in active photosynthetic apparatus
                ! Correction accounts for non-photosynthetic N (structural proteins, storage)
                ! Formula: Chl2C x (quota / (quota - NCmin))
                ! As quota -> NCmin, more N is in photosynthetic machinery
                CHL2C_plast = Chl2C * (quota / (quota - NCmin))

                !===============================================================================
                ! Diatom Quotas
                !===============================================================================
                ! Calculates stoichiometric ratios for diatoms (large phytoplankton with
                ! silica frustules). Includes unique silicon quotas.
                !
                ! Variables:
                !   quota_dia       : Diatom N:C quota [mmolN mmolC-1]
                !   recipQuota_dia  : Diatom C:N ratio [mmolC mmolN-1]
                !   Chl2C_dia       : Diatom Chl:C ratio [mgChl mmolC-1]
                !   Chl2N_dia       : Diatom Chl:N ratio [mgChl mmolN-1]
                !   CHL2C_plast_dia : Diatom plastidic Chl:C ratio [mgChl mmolC-1]
                !   qSiC            : Diatom Si:C quota [mmolSi mmolC-1]
                !   qSiN            : Diatom Si:N quota [mmolSi mmolN-1]
                !   DiaN            : Diatom nitrogen [mmolN m-3]
                !   DiaC            : Diatom carbon [mmolC m-3]
                !   DiaChl          : Diatom chlorophyll [mgChl m-3]
                !   DiaSi           : Diatom silicon [mmolSi m-3]
                !   NCmin_d         : Minimum diatom N:C quota [mmolN mmolC-1]
                !
                ! Silicon Quota Significance:
                !   - Required for frustule (shell) formation
                !   - Typical Si:C ≈ 0.13 (Brzezinski 1985)
                !   - Low Si:C -> thin frustules, increased sinking mortality
                !   - High Si:C -> thick frustules, enhanced sinking
                !   - Si limitation can occur even when N is abundant
                !
                ! Diatom-Specific Features:
                !   - Generally lower Chl:C than small phytoplankton (package effect)
                !   - Higher maximum growth rates when nutrient replete
                !   - Bloom-forming under high-nutrient conditions
                !-------------------------------------------------------------------------------

                ! Nitrogen:Carbon quota
                quota_dia = DiaN / DiaC

                ! Carbon:Nitrogen ratio (reciprocal)
                recipQuota_dia = real(one) / quota_dia

                ! Chlorophyll:Carbon ratio
                ! Generally lower than small phytoplankton due to large cell size (package effect)
                Chl2C_dia = DiaChl / DiaC

                ! Chlorophyll:Nitrogen ratio
                Chl2N_dia = DiaChl / DiaN

                ! Plastidic Chlorophyll:Carbon ratio
                ! Corrected for non-photosynthetic nitrogen allocation
                CHL2C_plast_dia = Chl2C_dia * (quota_dia / (quota_dia - NCmin_d))

                ! Silicon:Carbon quota
                ! Critical for frustule formation and diatom physiology
                ! Low Si:C indicates silicon limitation
                qSiC = DiaSi / DiaC

                ! Silicon:Nitrogen quota
                ! Alternative measure of silicon status relative to cellular nitrogen
                qSiN = DiaSi / DiaN

                !===============================================================================
                ! Additional Phytoplankton Quotas (OPTIONAL)
                !===============================================================================
                ! Calculates quotas for coccolithophores and Phaeocystis when enabled.
                ! These groups have distinct biogeochemical roles.
                !
                ! Coccolithophores:
                !   - Calcifying phytoplankton (produce CaCO3 plates)
                !   - Warm-water adapted
                !   - Important for carbonate counter-pump
                !
                ! Phaeocystis:
                !   - Colonial phytoplankton (can form large blooms)
                !   - Produces mucilaginous matrix
                !   - Cold-water species (polar and temperate)
                !   - High aggregation potential
                !-------------------------------------------------------------------------------

                if (enable_coccos) then

                    !===========================================================================
                    ! Coccolithophore Quotas
                    !===========================================================================
                    ! Calcifying phytoplankton with calcium carbonate plates (coccoliths)
                    !
                    ! Variables:
                    !   quota_cocco       : Cocco N:C quota [mmolN mmolC-1]
                    !   recipQuota_cocco  : Cocco C:N ratio [mmolC mmolN-1]
                    !   Chl2C_cocco       : Cocco Chl:C ratio [mgChl mmolC-1]
                    !   Chl2N_cocco       : Cocco Chl:N ratio [mgChl mmolN-1]
                    !   CHL2C_plast_cocco : Cocco plastidic Chl:C ratio [mgChl mmolC-1]
                    !   CoccoN            : Coccolithophore nitrogen [mmolN m-3]
                    !   CoccoC            : Coccolithophore carbon [mmolC m-3]
                    !   CoccoChl          : Coccolithophore chlorophyll [mgChl m-3]
                    !   NCmin_c           : Minimum cocco N:C quota [mmolN mmolC-1]
                    !
                    ! Note: Additional calcite quotas (CaCO3:C) calculated in calcification module
                    !---------------------------------------------------------------------------

                    ! Nitrogen:Carbon quota
                    quota_cocco = CoccoN / CoccoC

                    ! Carbon:Nitrogen ratio
                    recipQuota_cocco = real(one) / quota_cocco

                    ! Chlorophyll:Carbon ratio
                    Chl2C_cocco = CoccoChl / CoccoC

                    ! Chlorophyll:Nitrogen ratio
                    Chl2N_cocco = CoccoChl / CoccoN

                    ! Plastidic Chlorophyll:Carbon ratio
                    CHL2C_plast_cocco = Chl2C_cocco * (quota_cocco / (quota_cocco - NCmin_c))

                    !===========================================================================
                    ! Phaeocystis Quotas
                    !===========================================================================
                    ! Colonial phytoplankton that forms large blooms in cold waters
                    !
                    ! Variables:
                    !   quota_phaeo       : Phaeo N:C quota [mmolN mmolC-1]
                    !   recipQuota_phaeo  : Phaeo C:N ratio [mmolC mmolN-1]
                    !   Chl2C_phaeo       : Phaeo Chl:C ratio [mgChl mmolC-1]
                    !   Chl2N_phaeo       : Phaeo Chl:N ratio [mgChl mmolN-1]
                    !   CHL2C_plast_phaeo : Phaeo plastidic Chl:C ratio [mgChl mmolC-1]
                    !   PhaeoN            : Phaeocystis nitrogen [mmolN m-3]
                    !   PhaeoC            : Phaeocystis carbon [mmolC m-3]
                    !   PhaeoChl          : Phaeocystis chlorophyll [mgChl m-3]
                    !   NCmin_p           : Minimum Phaeo N:C quota [mmolN mmolC-1]
                    !
                    ! Ecological Notes:
                    !   - Forms colonial mucilaginous matrix (contributes to DOM)
                    !   - Can dominate Arctic/Antarctic spring blooms
                    !   - Enhanced aggregation and export potential
                    !---------------------------------------------------------------------------

                    ! Nitrogen:Carbon quota
                    quota_phaeo = PhaeoN / PhaeoC

                    ! Carbon:Nitrogen ratio
                    recipQuota_phaeo = real(one) / quota_phaeo

                    ! Chlorophyll:Carbon ratio
                    Chl2C_phaeo = PhaeoChl / PhaeoC

                    ! Chlorophyll:Nitrogen ratio
                    Chl2N_phaeo = PhaeoChl / PhaeoN

                    ! Plastidic Chlorophyll:Carbon ratio
                    CHL2C_plast_phaeo = Chl2C_phaeo * (quota_phaeo / (quota_phaeo - NCmin_p))

                end if

                !===============================================================================
                ! Zooplankton and Detritus Quotas
                !===============================================================================
                ! Calculates carbon:nitrogen ratios for consumers and detrital pools.
                ! These ratios are more constrained than phytoplankton (less variable).
                !
                ! Zooplankton C:N ratios:
                !   - Typically near Redfield ratio (C:N ≈ 6.6)
                !   - Less variable than phytoplankton (homeostatic regulation)
                !   - Important for grazer nutrition and trophic transfer efficiency
                !
                ! Detritus C:N ratios:
                !   - Reflects source material composition
                !   - Can increase with depth (preferential N remineralization)
                !   - Affects remineralization stoichiometry
                !
                ! Variables:
                !   recipQZoo       : Mesozooplankton C:N ratio [mmolC mmolN-1]
                !   recipQZoo2      : Macrozooplankton C:N ratio [mmolC mmolN-1]
                !   recipQZoo3      : Microzooplankton C:N ratio [mmolC mmolN-1]
                !   recipDet        : Slow-sinking detritus C:N ratio [mmolC mmolN-1]
                !   recipDet2       : Fast-sinking detritus C:N ratio [mmolC mmolN-1]
                !   recip_hetN_plus : Stable divisor for respiration calculations [mmolN-1 m3]
                !   HetC, HetN      : Mesozooplankton carbon and nitrogen [mmol m-3]
                !   Zoo2C, Zoo2N    : Macrozooplankton carbon and nitrogen [mmol m-3]
                !   MicZooC, MicZooN: Microzooplankton carbon and nitrogen [mmol m-3]
                !   DetC, DetN      : Detritus carbon and nitrogen [mmol m-3]
                !   DetZ2C, DetZ2N  : Fast-sinking detritus carbon and nitrogen [mmol m-3]
                !   tiny_het        : Small number to prevent division by zero [mmolN m-3]
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Mesozooplankton (Primary Heterotroph) Quotas
                !-------------------------------------------------------------------------------
                ! Primary grazer consuming phytoplankton and microzooplankton

                ! Carbon:Nitrogen ratio
                ! Used for converting nitrogen-based grazing to carbon fluxes
                recipQZoo = HetC / HetN

                ! Stable reciprocal for respiration calculations
                ! Prevents division by zero when zooplankton biomass is very low
                ! Used in Redfield-based respiration formulation
                recip_hetN_plus = 1.d0 / (HetN + tiny_het)

                !-------------------------------------------------------------------------------
                ! Detritus Quotas
                !-------------------------------------------------------------------------------
                ! Organic matter pools with variable C:N ratios

                if (Grazing_detritus) then
                    ! Slow-sinking detritus C:N ratio
                    ! Source: Unassimilated food, mortality, aggregation
                    ! Generally close to Redfield but can be elevated
                    recipDet = DetC / DetN
                end if

                if (enable_3zoo2det) then

                    !---------------------------------------------------------------------------
                    ! Additional Zooplankton Quotas (3-Zooplankton Model)
                    !---------------------------------------------------------------------------

                    ! Macrozooplankton (e.g., krill) C:N ratio
                    ! Larger organisms with potentially different stoichiometry
                    recipQZoo2 = Zoo2C / Zoo2N

                    ! Microzooplankton (e.g., ciliates, heterotrophic dinoflagellates) C:N ratio
                    ! Smallest heterotrophs, link to microbial loop
                    recipQZoo3 = MicZooC / MicZooN

                    if (Grazing_detritus) then
                        !-----------------------------------------------------------------------
                        ! Fast-Sinking Detritus (Fecal Pellets)
                        !-----------------------------------------------------------------------
                        ! Large, rapidly sinking particles
                        ! Source: Zooplankton fecal pellets, large aggregates
                        ! Important for biological pump and carbon export

                        recipDet2 = DetZ2C / DetZ2N
                    end if

                end if

                ! ===================================================================
                ! CARBON ISOTOPE TRACERS (if enabled)
                ! ===================================================================

                if (ciso) then
                    ! 13C isotope tracers
                    DIC_13 = max(tiny, state(k, idic_13) + sms(k, idic_13))
                    PhyC_13 = max(tiny_C, state(k, iphyc_13) + sms(k, iphyc_13))
                    DetC_13 = max(tiny, state(k, idetc_13) + sms(k, idetc_13))
                    HetC_13 = max(tiny, state(k, ihetc_13) + sms(k, ihetc_13))
                    EOC_13 = max(tiny, state(k, idoc_13) + sms(k, idoc_13))
                    DiaC_13 = max(tiny_C, state(k, idiac_13) + sms(k, idiac_13))
                    PhyCalc_13 = max(tiny, state(k, iphycal_13) + sms(k, iphycal_13))
                    DetCalc_13 = max(tiny, state(k, idetcal_13) + sms(k, idetcal_13))

                    ! 13C calcite dissolution with fractionation
                    calc_diss_13 = alpha_dcal_13 * calc_diss

                    ! 13C quotas
                    quota_13 = PhyN / PhyC_13
                    recipQuota_13 = real(one) / quota_13
                    quota_dia_13 = DiaN / DiaC_13
                    recipQuota_dia_13 = real(one) / quota_dia_13
                    recipQZoo_13 = HetC_13 / HetN

                    ! 14C radiocarbon tracers (if enabled)
                    if (ciso_14) then
                        DIC_14 = max(tiny, state(k, idic_14) + sms(k, idic_14))

                        if (ciso_organic_14) then
                            PhyC_14 = max(tiny_C, state(k, iphyc_14) + sms(k, iphyc_14))
                            DetC_14 = max(tiny, state(k, idetc_14) + sms(k, idetc_14))
                            HetC_14 = max(tiny, state(k, ihetc_14) + sms(k, ihetc_14))
                            EOC_14 = max(tiny, state(k, idoc_14) + sms(k, idoc_14))
                            DiaC_14 = max(tiny_C, state(k, idiac_14) + sms(k, idiac_14))
                            PhyCalc_14 = max(tiny, state(k, iphycal_14) + sms(k, iphycal_14))
                            DetCalc_14 = max(tiny, state(k, idetcal_14) + sms(k, idetcal_14))

                            calc_diss_14 = alpha_dcal_14 * calc_diss

                            quota_14 = PhyN / PhyC_14
                            recipQuota_14 = real(one) / quota_14
                            quota_dia_14 = DiaN / DiaC_14
                            recipQuota_dia_14 = real(one) / quota_dia_14
                            recipQZoo_14 = HetC_14 / HetN
                        end if ! ciso_organic_14
                    end if ! ciso_14
                end if ! ciso

                !===============================================================================
                ! TEMPERATURE DEPENDENCE OF METABOLIC RATES
                !===============================================================================
                ! Calculates how temperature affects biological rates using multiple
                ! formulations appropriate for different organism groups.
                !
                ! General Principle:
                !   - Metabolic rates increase exponentially with temperature (Q10 rule)
                !   - Different organisms have different thermal optima and tolerances
                !   - Cold-adapted vs warm-adapted species
                !
                ! Variables (Arrhenius):
                !   rTloc       : Inverse of local absolute temperature [K-1]
                !   arrFunc     : Arrhenius temperature function [-]
                !   Temp(k)     : Temperature at depth k [degC]
                !   C2K         : Celsius to Kelvin conversion (273.15) [K]
                !   Ae          : Activation energy parameter [K]
                !   rTref       : Inverse reference temperature (1/288.15 K at 15degC) [K-1]
                !
                ! Alternative Formulations (commented in code):
                !   - Eppley (1972): Log-linear for phytoplankton
                !   - Li (1980): Parabolic curve
                !   - Ahlgren (1987): Optimum curve
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Standard Arrhenius Function
                !-------------------------------------------------------------------------------
                ! General metabolic rate temperature dependence (Schourup 2013, Eq. A54)
                ! Exponential increase with temperature (no upper thermal limit)
                !
                ! Equation:
                !   f(T) = exp(-Ae × (1/T - 1/Tref))
                !
                ! Where:
                !   - Ae: Slope of Arrhenius plot (activation energy/gas constant)
                !   - Tref: Reference temperature (typically 15degC = 288.15 K)
                !   - T: Absolute temperature [K]
                !
                ! Used for: Most metabolic processes without strong thermal limits
                !   (e.g., remineralization, basal metabolism)

                ! Calculate inverse absolute temperature
                rTloc = real(one) / (Temp(k) + C2K)

                ! Calculate Arrhenius function
                arrFunc = exp(-Ae * (rTloc - rTref))

                !-------------------------------------------------------------------------------
                ! Phytoplankton-Specific Temperature Functions
                !-------------------------------------------------------------------------------
                ! Species-specific exponential and optimum-curve temperature responses
                ! Tuned for 4-plankton functional type version (small phyto, diatoms, coccos,
                ! Phaeocystis)
                !
                ! Note: New functions require calibration if adapted to 2-plankton version
                !-------------------------------------------------------------------------------

                ! Old coccolithophore temperature function (commented out):
                ! Power law form from Fielding (2013) based on observational growth rates
                ! if (enable_coccos) then
                !     CoccoTFunc = max(0.1419d0 * Temp(k)**0.8151d0, tiny)
                ! endif

                if (enable_coccos) then

                    !---------------------------------------------------------------------------
                    ! Small Phytoplankton
                    !---------------------------------------------------------------------------
                    ! Exponential temperature response: f(T) = exp(a + b×T)
                    ! Represents diverse group with broad thermal tolerance
                    ! Monotonic increase with temperature (no upper limit in this formulation)

                    Temp_phyto = exp(ord_phy + expon_phy * Temp(k))
                    VTTemp_phyto(k) = Temp_phyto ! Store for diagnostics

                    !---------------------------------------------------------------------------
                    ! Diatoms
                    !---------------------------------------------------------------------------
                    ! Large phytoplankton with silica frustules
                    ! Exponential form with different parameters than small phyto
                    ! Generally favored by cooler, nutrient-rich conditions

                    Temp_diatoms = exp(ord_d + expon_d * Temp(k))
                    VTTemp_diatoms(k) = Temp_diatoms ! Store for diagnostics

                    !---------------------------------------------------------------------------
                    ! Coccolithophores
                    !---------------------------------------------------------------------------
                    ! Calcifying phytoplankton with minimum temperature threshold
                    ! Cold-intolerant: minimal growth below 5degC
                    ! Exponential increase above threshold temperature
                    !
                    ! Ecological rationale: Coccos typically dominate in warm, stratified waters

                    if (Temp(k) < 5.0) then
                        ! Below threshold: minimal metabolic activity
                        Temp_cocco = tiny
                    else
                        ! Above threshold: exponential response
                        Temp_cocco = exp(ord_cocco + expon_cocco * Temp(k))
                        Temp_cocco = max(Temp_cocco, tiny) ! Ensure positive values
                    end if
                    VTTemp_cocco(k) = Temp_cocco ! Store for diagnostics

                    !---------------------------------------------------------------------------
                    ! Phaeocystis
                    !---------------------------------------------------------------------------
                    ! Colonial phytoplankton with bell-shaped temperature response
                    ! Blanchard function from Grimaud et al. (2017)
                    ! Has optimal temperature with decline at high and low temperatures
                    !
                    ! Equation: f(T) = uopt × ((Tmax-T)/(Tmax-Topt))^β ×
                    ! exp(-β×(Topt-T)/(Tmax-Topt))
                    !
                    ! Where:
                    !   - uopt: Maximum growth rate at optimal temperature [day-1]
                    !   - Topt: Optimal temperature [degC]
                    !   - Tmax: Maximum temperature (growth = 0) [degC]
                    !   - β: Shape parameter (steepness of curve) [-]
                    !
                    ! Ecological rationale: Phaeocystis blooms occur at specific temperature ranges
                    ! (typically cold-temperate waters, 0-10degC)

                    Temp_phaeo = uopt_phaeo &
                            * ((Tmax_phaeo - Temp(k)) / (Tmax_phaeo - Topt_phaeo)) ** beta_phaeo &
                            * exp(-beta_phaeo * (Topt_phaeo - Temp(k)) / (Tmax_phaeo - Topt_phaeo))
                    Temp_phaeo = max(Temp_phaeo, tiny) ! Ensure positive values
                    VTTemp_phaeo(k) = Temp_phaeo ! Store for diagnostics

                end if

                !-------------------------------------------------------------------------------
                ! Zooplankton Temperature Dependencies
                !-------------------------------------------------------------------------------
                ! Temperature functions for different zooplankton types with thermal limits
                ! Q10 formulations: exponential increase with ~doubling per 10degC

                if (enable_3zoo2det) then

                    !---------------------------------------------------------------------------
                    ! Macrozooplankton (Krill) Temperature Function
                    !---------------------------------------------------------------------------
                    ! Sigmoid function with upper thermal limit
                    ! Accounts for thermal stress at high temperatures
                    !
                    ! Equation: f(T) = exp(t1/t2 - t1/T) / (1 + exp(t3/t4 - t3/T))
                    !
                    ! Numerator: Exponential increase with temperature
                    ! Denominator: Sigmoid decline at high temperatures (thermal stress)
                    !
                    ! Ecological rationale: Macrozooplankton have defined thermal niches
                    ! (e.g., Antarctic krill prefer cold water, decline above ~4degC)

                    arrFuncZoo2 = exp(t1_zoo2 / t2_zoo2 - t1_zoo2 * rTloc) / &
                            (1.0 + exp(t3_zoo2 / t4_zoo2 - t3_zoo2 * rTloc))

                    !---------------------------------------------------------------------------
                    ! Q10 Temperature Coefficients
                    !---------------------------------------------------------------------------
                    ! Q10 formulation: rate = Q10^(T/10)
                    ! Simple exponential increase with temperature
                    !
                    ! Q10 values:
                    !   - ~1.02-1.04: Moderate temperature sensitivity (typical for metabolism)
                    !   - ~1.09: Higher sensitivity (respiration processes)
                    !
                    ! Ecological interpretation:
                    !   - Smaller organisms (microzooplankton) often have higher Q10
                    !   - Respiration Q10 > growth Q10 (maintenance costs increase faster)

                    q10_mes = 1.0242 ** (Temp(k)) ! Mesozooplankton metabolism
                    q10_mic = 1.04 ** (Temp(k)) ! Microzooplankton metabolism
                    q10_mes_res = 1.0887 ** (Temp(k)) ! Mesozooplankton respiration
                    q10_mic_res = 1.0897 ** (Temp(k)) ! Microzooplankton respiration

                end if

                !-------------------------------------------------------------------------------
                ! Silicate Dissolution Temperature Dependence
                !-------------------------------------------------------------------------------
                ! Temperature effect on biogenic silica (diatom frustule) dissolution
                ! Higher temperatures accelerate chemical dissolution kinetics
                !
                ! Variables:
                !   reminSiT : Temperature-dependent Si dissolution rate [day-1]
                !   reminSi  : Minimum dissolution rate [day-1]
                !
                ! Exponential formulation: 2.6× increase per 10degC
                ! Reference temperature: 10degC
                !

                reminSiT = max(0.023d0 * 2.6d0 ** ((Temp(k) - 10.0) / 10.0), reminSi)

                ! Alternative Kamatani (1982) function  (commented out):
                ! reminSiT = min(1.32e16 * exp(-11200.d0 * rTloc), reminSi)

                !===============================================================================
                ! 2. OXYGEN DEPENDENCE OF REMINERALIZATION
                !===============================================================================
                ! Calculates oxygen limitation effects on aerobic organic matter decomposition.
                ! Important for oxygen minimum zones (OMZs) and suboxic/anoxic environments.
                !
                ! Variables:
                !   O2Func      : Oxygen limitation factor [0-1, 0=anoxic, 1=oxic]
                !   O2          : Dissolved oxygen concentration [mmolO2 m-3]
                !   k_o2_remin  : Half-saturation for O2-limited remineralization [mmolO2 m-3]
                !   O2dep_remin : Flag to enable O2-dependent remineralization [logical]
                !
                ! Michaelis-Menten Formulation:
                !   f(O2) = O2 / (k_O2 + O2)
                !
                ! Parameter Value:
                !   k_o2_remin = 15 mmolO2 m-3 (half-saturation constant)
                !   Range: 0-30 mmolO2 m-3 based on DeVries & Weber (2017), cited in Cram (2018)
                !
                ! Ecological/Biogeochemical Significance:
                !   - Aerobic respiration dominates in oxic waters (O2 > 30 mmol m-3)
                !   - Suboxic/anoxic metabolism (denitrification, sulfate reduction) in OMZs
                !   - Reduced remineralization efficiency in low-O2 environments
                !   - Important for nutrient cycling and carbon export in OMZs
                !
                ! Note: When O2 < 0.1 mmol m-3, consider switching to anaerobic pathways
                !       (not implemented in this version)
                !-------------------------------------------------------------------------------

                ! Default: no oxygen limitation (fully oxic conditions)
                O2Func = 1.d0

                if (O2dep_remin) then
                    ! Enable oxygen-dependent remineralization
                    ! Michaelis-Menten type limitation
                    ! Becomes significant when O2 < ~30 mmol m-3
                    O2Func = O2 / (k_o2_remin + O2)
                end if

                !===============================================================================
                ! LIGHT AVAILABILITY CALCULATION
                !===============================================================================
                ! Calculates photosynthetically available radiation (PAR) through the water
                ! column using Beer-Lambert law with chlorophyll-based attenuation.
                !
                ! Light Attenuation Components:
                !   1. Water attenuation (k_w): Clear water absorption and scattering
                !   2. Chlorophyll attenuation (a_chl): Phytoplankton self-shading
                !
                ! Variables:
                !   PARave          : Average PAR at depth k [W m-2]
                !   PAR(k)          : Stored PAR for layer k [W m-2]
                !   SurfSR          : Surface solar radiation [W m-2]
                !   chl_upper       : Chlorophyll at upper layer boundary [mgChl m-3]
                !   chl_lower       : Chlorophyll at lower layer boundary [mgChl m-3]
                !   Chlave          : Average chlorophyll in layer [mgChl m-3]
                !   kappa           : Total attenuation coefficient [m-1]
                !   kappastar       : Angle-corrected attenuation coefficient [m-1]
                !   k_w             : Water attenuation coefficient [m-1]
                !   a_chl           : Chlorophyll-specific attenuation [(m2 mgChl-1)]
                !   cosAI(n)        : Cosine of solar zenith angle [-]
                !   thick(k)        : Layer thickness [m]
                !   kdzLower        : Cumulative optical depth [dimensionless]
                !   kdzUpper        : Cumulative optical depth at upper boundary [dimensionless]
                !
                ! Beer-Lambert Law:
                !   I(z) = I0 × exp(-κ×z)
                !   where κ = k_w + a_chl×[Chl]
                !
                ! Self-Shading Effect:
                !   - High chlorophyll reduces light penetration
                !   - Limits bloom depth and total biomass
                !   - Creates trade-off between cell density and light availability
                !
                ! Ecological Significance:
                !   - Defines euphotic zone depth (1% surface light)
                !   - Controls vertical distribution of primary production
                !   - Self-shading is key negative feedback on bloom magnitude
                !-------------------------------------------------------------------------------

                if (k == 1) then

                    !===========================================================================
                    ! SURFACE LAYER INITIALIZATION
                    !===========================================================================
                    ! Surface layer receives full incident solar radiation
                    ! Initialize chlorophyll and optical depth for subsurface calculations
                    !---------------------------------------------------------------------------

                    ! Surface PAR equals incident solar radiation
                    PARave = max(tiny, SurfSR)
                    PAR(k) = PARave

                    ! Initialize surface chlorophyll for attenuation calculation
                    ! Sum all phytoplankton functional types
                    chl_upper = (PhyChl + DiaChl) ! Base groups (always present)

                    if (enable_coccos) then
                        ! Add coccolithophores and Phaeocystis if enabled
                        chl_upper = chl_upper + CoccoChl + PhaeoChl
                    end if
                else

                    !===========================================================================
                    ! SUBSURFACE LIGHT ATTENUATION
                    !===========================================================================
                    ! Calculate light penetration through water column using Beer-Lambert law
                    ! with chlorophyll-based self-shading
                    !---------------------------------------------------------------------------

                    !---------------------------------------------------------------------------
                    ! Calculate Current Layer Chlorophyll
                    !---------------------------------------------------------------------------

                    chl_lower = PhyChl + DiaChl

                    if (enable_coccos) then
                        chl_lower = chl_lower + CoccoChl + PhaeoChl
                    end if

                    ! Average chlorophyll between layer boundaries
                    ! Assumes linear interpolation within layer
                    Chlave = (chl_upper + chl_lower) * 0.5

                    !---------------------------------------------------------------------------
                    ! Calculate Attenuation Coefficient
                    !---------------------------------------------------------------------------

                    ! Total attenuation coefficient
                    ! k_w: Clear water absorption (~0.04 m-1 in ocean)
                    ! a_chl: Chlorophyll-specific attenuation (~0.03-0.05 m2 mgChl-1)
                    kappa = k_w + a_chl * Chlave

                    ! Correct for solar zenith angle (path length through water)
                    ! Lower sun angle -> longer path -> more attenuation
                    kappastar = kappa / cosAI(n)

                    ! Cumulative optical depth (dimensionless)
                    ! Integrates attenuation over depth
                    kdzLower = kdzUpper + kappastar * thick(k - 1)

                    !---------------------------------------------------------------------------
                    ! Calculate Light at Layer
                    !---------------------------------------------------------------------------

                    ! Beer-Lambert law: exponential decay with optical depth
                    Lowerlight = SurfSR * exp(-kdzLower)
                    Lowerlight = max(tiny, Lowerlight) ! Ensure positive value

                    ! Store PAR for this layer
                    PARave = Lowerlight
                    PAR(k) = PARave

                    ! Update variables for next layer
                    chl_upper = chl_lower ! Current lower becomes next upper
                    kdzUpper = kdzLower ! Current cumulative depth for next layer

                end if

                !===============================================================================
                ! MARINE CARBONATE SYSTEM CALCULATIONS (MOCSY)
                !===============================================================================
                ! Calculates complete marine carbonate chemistry using the MOCSY package
                ! (Marine Ocean Carbon System Solver).
                !
                ! This module calculates:
                !   1. Carbonate system speciation (CO2, HCO3-, CO3--)
                !   2. pH and partial pressure of CO2 (pCO2)
                !   3. Carbonate saturation states (Omega for calcite and aragonite)
                !   4. Solubility products and seawater properties
                !
                ! Key Features:
                !   - Adaptive update frequency (depth-dependent)
                !   - Euphotic zone: Weekly updates (high biological activity)
                !   - Deep waters: Monthly updates (slower changes)
                !   - Complete thermodynamic consistency
                !   - Pressure correction for depth
                !
                ! MOCSY Package:
                !   - Developed for OCMIP5 project
                !   - Solves carbonate system from two known parameters
                !   - Accounts for temperature, salinity, pressure effects
                !   - Multiple equilibrium constant formulations available
                !
                ! Update Strategy:
                !   - Initialize on first time step (mstep = 1)
                !   - Euphotic zone (PAR > 1% surface): 7-day updates
                !   - Deep waters (PAR < 1% surface): 30-day updates
                !   - Rationale: Biological activity drives rapid changes near surface
                !
                ! Input Parameters:
                !   - Temperature, Salinity (from physical model)
                !   - DIC, Alkalinity (from biogeochemical tracers)
                !   - Silicate, Phosphate (affects equilibrium constants)
                !   - Atmospheric pressure, Latitude (for gas exchange)
                !
                ! Output Variables:
                !   - pH, pCO2, fCO2 (CO2 partial and fugacity)
                !   - CO2, HCO3-, CO3-- (carbonate species concentrations)
                !   - OmegaC, OmegaA (calcite and aragonite saturation)
                !   - Solubility products, seawater density
                !
                ! Ecological/Biogeochemical Significance:
                !   - Controls CO2 uptake/release (air-sea exchange)
                !   - Regulates calcification and dissolution rates
                !   - Affects phytoplankton carbon acquisition
                !   - Critical for ocean acidification studies
                !
                ! References:
                !   - MOCSY: http://ocmip5.ipsl.jussieu.fr/mocsy/
                !   - Orr & Epitalon (2015) - MOCSY 2.0 user guide
                !===============================================================================

                !===============================================================================
                ! PREPARATION AND INITIALIZATION
                !===============================================================================
                ! Prepares input data and initializes carbonate system on first time step.
                !
                ! Variables:
                !   dpos(1)                 : Depth for pressure calculations [m, positive down]
                !   zF(k)                   : Model depth coordinate [m, negative down]
                !   mstep                   : Model time step counter [-]
                !   k                       : Vertical layer index [-]
                !
                ! Input Arrays for MOCSY:
                !   REcoM_T_depth           : Temperature [degC, potential temperature]
                !   REcoM_S_depth           : Salinity [psu, practical salinity]
                !   REcoM_Alk_depth         : Total alkalinity [mol m-3]
                !   REcoM_DIC_depth         : Dissolved inorganic carbon [mol m-3]
                !   REcoM_Si_depth          : Silicate concentration [mol m-3]
                !   REcoM_Phos_depth        : Phosphate concentration [mol m-3]
                !   Patm_depth              : Atmospheric pressure [atm]
                !   Latd                    : Latitude [degrees]
                !   Nmocsy                  : Number of points (1 for single depth)
                !
                ! Note: Depth convention conversion required (model uses negative depths)
                !-------------------------------------------------------------------------------

                ! Convert model depth coordinate to positive depth for MOCSY
                ! Model convention: zF(k) is negative (e.g., -100 m)
                ! MOCSY convention: depth is positive (e.g., 100 m)
                dpos(1) = -zF(k)

                ! Calculate update frequencies based on model time step
                mocsy_step_per_day = 1.0 / dt_b
                logfile_outfreq_7 = int(mocsy_step_per_day * 7.0) ! Steps in 7 days
                logfile_outfreq_30 = int(mocsy_step_per_day * 30.0) ! Steps in 30 days

                !===============================================================================
                ! INITIAL CARBONATE SYSTEM CALCULATION
                !===============================================================================
                ! Calculates complete carbonate system on first model time step.
                ! Provides initial conditions for all carbonate chemistry variables.
                !
                ! Output Variables (from MOCSY):
                !   ph_depth(1)             : pH on total scale [-]
                !   pco2_depth(1)           : Partial pressure of CO2 [μatm]
                !   fco2_depth(1)           : Fugacity of CO2 [μatm]
                !   co2_depth(1)            : Dissolved CO2 concentration [mol m-3]
                !   hco3_depth(1)           : Bicarbonate concentration [mol m-3]
                !   co3_depth(1)            : Carbonate ion concentration [mol m-3]
                !   OmegaA_depth(1)         : Aragonite saturation state [-]
                !   OmegaC_depth(1)         : Calcite saturation state [-]
                !   kspc_depth(1)           : Calcite solubility product [mol2 kg-2]
                !   BetaD_depth(1)          : Revelle factor (buffer capacity) [-]
                !   rhoSW_depth(1)          : Seawater density [kg m-3]
                !   p_depth(1)              : Pressure [bar]
                !   tempis_depth(1)         : In situ temperature [degC]
                !
                ! Water Column Storage:
                !   CO2_watercolumn(k)      : Stored CO2 for biological calculations [mol m-3]
                !   pH_watercolumn(k)       : Stored pH [-]
                !   pCO2_watercolumn(k)     : Stored pCO2 [μatm]
                !   HCO3_watercolumn(k)     : Stored bicarbonate [mol m-3]
                !   CO3_watercolumn(k)      : Stored carbonate [mol m-3]
                !   OmegaC_watercolumn(k)   : Stored calcite saturation [-]
                !   kspc_watercolumn(k)     : Stored solubility product [mol2 kg-2]
                !   rhoSW_watercolumn(k)    : Stored seawater density [kg m-3]
                !
                ! MOCSY Options:
                !   optCON='mol/m3'  : Concentration units (mol/m3)
                !   optT='Tpot   '   : Temperature is potential temperature
                !   optP='m '        : Pressure given as depth in meters
                !   optB='u74'       : Boron:Salinity ratio from Uppström (1974)
                !   optK1K2='l  '    : Carbonic acid constants from Lueker et al. (2000)
                !   optKf='dg'       : HF constant from Dickson & Goyet (1994)
                !   optGAS='Pinsitu' : Pressure is in situ (accounts for depth)
                !   optS='Sprc'      : Salinity on practical scale
                !-------------------------------------------------------------------------------

                !===============================================================================
                ! ADAPTIVE CARBONATE SYSTEM UPDATE FREQUENCY
                !===============================================================================
                ! Determines how often to recalculate carbonate system based on depth.
                ! More frequent updates where biological activity drives rapid changes.
                !
                ! Variables:
                !   mocsy_step_per_day  : Number of model time steps per day [-]
                !   dt_b                : Model time step [days]
                !   logfile_outfreq_7   : Number of steps in 7 days [-]
                !   logfile_outfreq_30  : Number of steps in 30 days [-]
                !   PARave              : Average photosynthetically active radiation [W m-2]
                !   SurfSR              : Surface solar radiation [W m-2]
                !
                ! Update Strategy:
                !   - Euphotic zone (PAR > 1% surface): 7-day updates
                !     * High photosynthesis rates alter DIC and pH rapidly
                !     * Important for accurate phytoplankton CO2 responses
                !   - Deep waters (PAR < 1% surface): 30-day updates
                !     * Slower changes dominated by remineralization and mixing
                !     * Reduces computational cost while maintaining accuracy
                !
                ! Computational Cost Considerations:
                !   - Carbonate system solving is computationally expensive
                !   - Adaptive frequency balances accuracy and performance
                !   - Typical speedup: 4× faster than daily updates everywhere
                !-------------------------------------------------------------------------------

                !===============================================================================
                ! EUPHOTIC ZONE UPDATES (WEEKLY)
                !===============================================================================
                ! Frequent updates in sunlit surface waters where biological activity is high.
                ! Euphotic zone defined as PAR > 1% of surface irradiance.
                !
                ! Biological Drivers:
                !   - Photosynthesis removes DIC, increases pH
                !   - Respiration/remineralization adds DIC, decreases pH
                !   - Calcification removes alkalinity
                !   - Rapid daily and seasonal cycles
                !
                ! Why 7-day updates?
                !   - Captures weekly-scale biological dynamics
                !   - Adequate for phytoplankton bloom progression
                !   - Reasonable computational cost
                !-------------------------------------------------------------------------------

                !===============================================================================
                ! DEEP WATER UPDATES (MONTHLY)
                !===============================================================================
                ! Less frequent updates in dark deep waters where changes are slower.
                ! Deep waters defined as PAR < 1% of surface irradiance.
                !
                ! Physical/Chemical Drivers:
                !   - Slow remineralization of sinking organic matter
                !   - Calcite dissolution (below saturation horizon)
                !   - Mixing and advection
                !   - No photosynthesis to drive rapid changes
                !
                ! Why 30-day updates?
                !   - Changes occur on monthly to seasonal timescales
                !   - Dominated by physical transport and slow remineralization
                !   - Significant computational savings with minimal accuracy loss
                !
                ! Note: Below permanent pycnocline, even longer update intervals
                !       could be justified (e.g., seasonal)
                !-------------------------------------------------------------------------------

                if (mstep == 1 .or. &
                        PARave > 0.01 * SurfSR .and. mod(mstep, logfile_outfreq_7) == 0 .or. &
                        PARave < 0.01 * SurfSR .and. mod(mstep, logfile_outfreq_30) == 0) then

                    ! Monthly updates in deep waters (low biological activity)
                    call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, &
                            co3_depth, &
                            OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, &
                            rhoSW_depth, p_depth, tempis_depth, &
                            REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, &
                            REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, dpos, Latd, Nmocsy, &
                            optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', &
                            optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')

                    ! Update water column arrays with new carbonate chemistry
                    CO2_watercolumn(k) = co2_depth(1)
                    pH_watercolumn(k) = ph_depth(1)
                    pCO2_watercolumn(k) = pco2_depth(1)
                    HCO3_watercolumn(k) = hco3_depth(1)
                    CO3_watercolumn(k) = co3_depth(1)
                    OmegaC_watercolumn(k) = OmegaC_depth(1)
                    kspc_watercolumn(k) = kspc_depth(1)
                    rhoSW_watercolumn(k) = rhoSW_depth(1)
                end if

                !===============================================================================
                ! CO2 EFFECTS AND CALCITE DISSOLUTION
                !===============================================================================
                ! Simulates ocean acidification impacts on phytoplankton growth and calcium
                ! carbonate dissolution in the water column.
                !
                ! This module calculates:
                !   1. CO2/pH effects on phytoplankton photosynthesis (species-specific)
                !   2. Calcite dissolution rates (saturation-state or depth-dependent)
                !
                ! Key Features:
                !   - Multi-component CO2 response (HCO3- benefit, CO2 toxicity, pH stress)
                !   - Species-specific sensitivities to ocean acidification
                !   - Two dissolution mechanisms (thermodynamic vs empirical)
                !   - Velocity-dependent dissolution for sinking particles
                !
                ! CO2 Response Function:
                !   f(CO2) = a×[HCO3-]/(b+[HCO3-]) - exp(-c×[CO2]) - d×[H+]
                !
                !   Three competing effects:
                !     1. HCO3- availability (substrate for carbon fixation) - POSITIVE
                !     2. High CO2 concentration (toxicity/stress) - NEGATIVE
                !     3. Low pH (proton stress on enzymes) - NEGATIVE
                !
                ! Phytoplankton Groups:
                !   - Small phytoplankton: Moderate sensitivity
                !   - Diatoms: Generally tolerant
                !   - Coccolithophores: High sensitivity (calcifying organisms)
                !   - Phaeocystis: Variable sensitivity
                !
                ! Dissolution Mechanisms:
                !   A) Saturation-state dependent (OmegaC_diss = TRUE)
                !      - Thermodynamically based on carbonate saturation
                !      - Aumont et al. (2015) parameterization
                !   B) Depth-dependent (OmegaC_diss = FALSE)
                !      - Empirical dissolution rates scaled by sinking velocity
                !
                ! Ecological/Biogeochemical Significance:
                !   - Ocean acidification reduces calcification and alters competitive balance
                !   - Dissolution releases CO2 and alkalinity back to seawater
                !   - Critical for understanding climate change impacts on marine ecosystems
                !
                ! References:
                !   - Aumont et al. (2015) - Saturation-state dependent dissolution
                !   - Schourup-Kristensen et al. (2013) - REcoM model description
                !===============================================================================

                !===============================================================================
                ! CO2 EFFECTS ON PHYTOPLANKTON GROWTH
                !===============================================================================
                ! Calculates how ocean acidification (elevated CO2, reduced pH) affects
                ! phytoplankton photosynthesis rates through a complex response function.
                !
                ! Response Function Components:
                !   1. Michaelis-Menten HCO3- uptake (carbonate benefit)
                !   2. Exponential CO2 inhibition (high CO2 toxicity)
                !   3. Linear H+ inhibition (pH stress on cellular processes)
                !
                ! Variables (Small Phytoplankton):
                !   PhyCO2              : CO2 effect modifier for small phyto [0-3]
                !   a_co2_phy           : HCO3- uptake parameter [-]
                !   b_co2_phy           : HCO3- half-saturation [mmolC m-3]
                !   c_co2_phy           : CO2 inhibition coefficient [m3 mmolC-1]
                !   d_co2_phy           : H+ stress coefficient [L mol-1]
                !   HCO3_watercolumn(k) : Bicarbonate concentration [mmolC m-3]
                !   CO2_watercolumn(k)  : Dissolved CO2 concentration [mmolC m-3]
                !   pH_watercolumn(k)   : Water column pH [-]
                !   Cunits              : Concentration units conversion factor [-]
                !   h_depth(1)          : Proton concentration at surface [mol L-1]
                !   VTPhyCO2(k)         : Diagnostic output for CO2 effect [-]
                !
                ! Constraints:
                !   - Upper limit: 3× enhancement (April 2022 modification)
                !   - Lower limit: 0 (no negative values, July 2022 modification)
                !
                ! Note: Similar calculations for diatoms, coccolithophores, and Phaeocystis
                !       with species-specific parameters reflecting different sensitivities
                !-------------------------------------------------------------------------------

                ! Convert pH to proton concentration for calculations
                ! pH = -log10[H+], therefore [H+] = 10^(-pH)
                h_depth(1) = 10.d0 ** (-ph_depth(1))

                ! Note: Cunits conversion not needed for [H+] because pH is already in mol/L

                !-------------------------------------------------------------------------------
                ! Small Phytoplankton CO2 Response
                !-------------------------------------------------------------------------------
                ! Moderate sensitivity to ocean acidification
                ! Represents diverse group with varied carbon acquisition strategies

                PhyCO2 = a_co2_phy * HCO3_watercolumn(k) * Cunits &
                        / (b_co2_phy + HCO3_watercolumn(k) * Cunits) &
                        - exp(-c_co2_phy * CO2_watercolumn(k) * Cunits) &
                        - d_co2_phy * 10.d0 ** (-pH_watercolumn(k))

                ! Apply empirical constraints based on observations
                PhyCO2 = min(PhyCO2, 3.d0) ! Upper limit: maximum 3x enhancement
                PhyCO2 = max(0.d0, PhyCO2) ! Lower limit: prevent negative growth response

                ! Store for diagnostics and output
                VTPhyCO2(k) = PhyCO2

                !-------------------------------------------------------------------------------
                ! Diatoms CO2 Response
                !-------------------------------------------------------------------------------
                ! Generally tolerant to elevated CO2
                ! Efficient carbon concentrating mechanisms (CCMs)

                DiaCO2 = a_co2_dia * HCO3_watercolumn(k) * Cunits &
                        / (b_co2_dia + HCO3_watercolumn(k) * Cunits) &
                        - exp(-c_co2_dia * CO2_watercolumn(k) * Cunits) &
                        - d_co2_dia * 10.d0 ** (-pH_watercolumn(k))

                ! Apply constraints
                DiaCO2 = min(DiaCO2, 3.d0) ! Upper limit: 3x enhancement
                DiaCO2 = max(0.d0, DiaCO2) ! Lower limit: no negative effect

                ! Store for diagnostics
                VTDiaCO2(k) = DiaCO2

                if (enable_coccos) then

                    !---------------------------------------------------------------------------
                    ! Coccolithophores CO2 Response
                    !---------------------------------------------------------------------------
                    ! Calcifying phytoplankton - highly sensitive to ocean acidification
                    ! Both photosynthesis and calcification affected by carbonate chemistry
                    ! May be disadvantaged under future high-CO2 conditions

                    CoccoCO2 = a_co2_cocco * HCO3_watercolumn(k) * Cunits &
                            / (b_co2_cocco + HCO3_watercolumn(k) * Cunits) &
                            - exp(-c_co2_cocco * CO2_watercolumn(k) * Cunits) &
                            - d_co2_cocco * 10.d0 ** (-pH_watercolumn(k))

                    ! Apply constraints
                    CoccoCO2 = min(CoccoCO2, 3.d0) ! Upper limit: 3x enhancement
                    CoccoCO2 = max(0.d0, CoccoCO2) ! Lower limit: no negative effect

                    ! Store for diagnostics
                    VTCoccoCO2(k) = CoccoCO2

                    !---------------------------------------------------------------------------
                    ! Phaeocystis CO2 Response
                    !---------------------------------------------------------------------------
                    ! Colonial phytoplankton with variable CO2 sensitivity
                    ! Response may depend on bloom stage and environmental conditions

                    PhaeoCO2 = a_co2_phaeo * HCO3_watercolumn(k) * Cunits &
                            / (b_co2_phaeo + HCO3_watercolumn(k) * Cunits) &
                            - exp(-c_co2_phaeo * CO2_watercolumn(k) * Cunits) &
                            - d_co2_phaeo * 10.d0 ** (-pH_watercolumn(k))

                    ! Apply constraints
                    PhaeoCO2 = min(PhaeoCO2, 3.d0) ! Upper limit: 3× enhancement
                    PhaeoCO2 = max(0.d0, PhaeoCO2) ! Lower limit: no negative effect

                    ! Store for diagnostics
                    VTPhaeoCO2(k) = PhaeoCO2

                end if

                !===============================================================================
                ! CALCITE DISSOLUTION
                !===============================================================================
                ! Calculates dissolution rates of calcium carbonate (CaCO3) in seawater.
                ! Dissolution depends on carbonate saturation state or depth (two modes).
                !
                ! Carbonate Saturation State (Omega):
                !   Ω = [CO3²⁻] / [CO3²⁻]sat
                !   - Ω > 1: Supersaturated (favors precipitation, slow dissolution)
                !   - Ω < 1: Undersaturated (favors dissolution)
                !
                ! Variables (Saturation-dependent mode):
                !   calc_diss           : Primary detritus dissolution rate [day-1]
                !   calc_diss2          : Secondary detritus dissolution rate [day-1]
                !   calc_diss_ben       : Benthic dissolution rate [day-1]
                !   calc_diss_omegac    : Base dissolution rate coefficient [day-1]
                !   calc_diss_exp       : Dissolution order exponent [-]
                !   Ca                  : Calcium ion concentration [mol kg-1]
                !   CO3_sat             : Saturated carbonate concentration [mol m-3]
                !   CO3_watercolumn(k)  : Actual carbonate concentration [mol m-3]
                !   kspc_watercolumn(k) : Solubility product for calcite [mol2 kg-2]
                !   rhoSW_watercolumn(k): Seawater density [kg m-3]
                !   Sali_depth(k)       : Salinity [psu]
                !
                ! Variables (Depth-dependent mode):
                !   calc_diss_rate      : Primary dissolution rate constant [day-1]
                !   calc_diss_rate2     : Secondary dissolution rate constant [day-1]
                !   Sink_Vel            : Particle sinking velocity [m day-1]
                !   Vdet_a              : Depth-dependent sinking coefficient [day-1]
                !   Vdet                : Base sinking velocity [m day-1]
                !   zF(k)               : Depth at layer k [m]
                !
                ! Note: Dissolution releases CO2 and alkalinity, affecting carbonate chemistry
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Calculate Sinking Velocity
                !-------------------------------------------------------------------------------
                ! Sinking velocity increases with depth (particle compaction, reduced drag)

                Sink_Vel = Vdet_a * abs(zF(k)) + Vdet

                if (OmegaC_diss) then

                    !===========================================================================
                    ! SATURATION-STATE DEPENDENT DISSOLUTION
                    !===========================================================================
                    ! Thermodynamically-based dissolution using carbonate saturation state
                    ! Based on Aumont et al. (2015) parameterization
                    !---------------------------------------------------------------------------

                    !---------------------------------------------------------------------------
                    ! Calculate Calcium Ion Concentration
                    !---------------------------------------------------------------------------
                    ! Conservative element: scales linearly with salinity
                    ! Reference: 0.02128 mol/kg at salinity 35 psu

                    Ca = (0.02128d0 / 40.078d0) * Sali_depth(k) / 1.80655d0

                    !---------------------------------------------------------------------------
                    ! Calculate Saturated Carbonate Ion Concentration
                    !---------------------------------------------------------------------------
                    ! [CO3²⁻]sat = Ksp / [Ca²+]
                    ! Where Ksp is the solubility product for calcite
                    ! Convert from mol/kg to mol/m3 using seawater density

                    CO3_sat = (kspc_watercolumn(k) / Ca) * rhoSW_watercolumn(k)

                    !---------------------------------------------------------------------------
                    ! Calculate Dissolution Rate
                    !---------------------------------------------------------------------------
                    ! Dissolution increases with undersaturation (Ω < 1)
                    ! Power law relationship: rate ∝ (1 - Ω)^n
                    ! Exponent (n) typically 1-4 depending on calcite form

                    calc_diss = calc_diss_omegac &
                            * max(zero, (1.0 - (CO3_watercolumn(k) / CO3_sat))) ** (calc_diss_exp)

                    ! Apply same dissolution rate to all detritus types
                    ! Assumes similar calcite characteristics regardless of source
                    if (enable_3zoo2det) then
                        calc_diss2 = calc_diss ! Fast-sinking detritus
                    end if
                    calc_diss_ben = calc_diss ! Benthic detritus

                else

                    !===========================================================================
                    ! DEPTH-DEPENDENT DISSOLUTION
                    !===========================================================================
                    ! Empirical dissolution rates scaled by sinking velocity
                    ! Simpler approach when carbonate chemistry is not fully resolved
                    !---------------------------------------------------------------------------

                    !---------------------------------------------------------------------------
                    ! Primary Detritus Dissolution
                    !---------------------------------------------------------------------------
                    ! Dissolution rate scales with sinking velocity
                    ! Faster sinking -> less time for dissolution per unit depth
                    ! Reference velocity: 20 m/day

                    calc_diss = calc_diss_rate * Sink_Vel / 20.d0

                    if (enable_3zoo2det) then
                        !-----------------------------------------------------------------------
                        ! Secondary Detritus Dissolution (Fast-sinking)
                        !-----------------------------------------------------------------------
                        ! May have different dissolution characteristics than primary detritus
                        ! (e.g., more compact, different organic coating)

                        calc_diss2 = calc_diss_rate2 * Sink_Vel / 20.d0
                    end if

                    !---------------------------------------------------------------------------
                    ! Benthic Dissolution
                    !---------------------------------------------------------------------------
                    ! Calcite reaching the sediment interface
                    ! Uses same rate as primary detritus

                    calc_diss_ben = calc_diss_rate * Sink_Vel / 20.d0

                end if

                !===============================================================================
                ! PHOTOSYNTHESIS LIMITATION FACTORS
                !===============================================================================
                ! Calculates how nutrients, iron, and temperature limit maximum photosynthesis
                ! rates using intracellular quota-based regulation (Droop model).
                !
                ! Limitation Approach:
                !   - Intracellular regulation through N:C, Si:C quotas
                !   - Michaelis-Menten iron limitation
                !   - Liebig's law: most limiting factor controls growth
                !
                ! Variables (Small Phytoplankton):
                !   qlimitFac    : Combined nutrient limitation factor [0-1, 0=limited, 1=replete]
                !   feLimitFac   : Iron limitation factor [0-1]
                !   quota        : Phytoplankton N:C ratio [mmolN mmolC-1]
                !   NCmin        : Minimum N:C ratio (subsistence quota) [mmolN mmolC-1]
                !   NMinSlope    : Steepness of limitation curve [-]
                !   Fe           : Dissolved iron concentration [mmolFe m-3]
                !   k_Fe         : Half-saturation for Fe uptake [mmolFe m-3]
                !   pMax         : Maximum photosynthesis rate [day-1]
                !   P_cm         : Maximum rate constant [day-1]
                !   Temp_phyto   : Temperature function for small phyto [-]
                !   arrFunc      : Arrhenius temperature function [-]
                !
                ! recom_limiter Function:
                !   Returns limitation factor based on quota:
                !   - Returns ~0 when quota near minimum (severely limited)
                !   - Returns ~1 when quota is high (nutrient replete)
                !   - Smooth transition controlled by slope parameter
                !
                ! Note: Similar calculations for diatoms, coccolithophores, and Phaeocystis
                !       with species-specific parameters and additional Si limitation for diatoms
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Small Phytoplankton
                !-------------------------------------------------------------------------------
                ! Limited by nitrogen and iron
                ! Represents small flagellates, small cyanobacteria, etc.

                ! Nitrogen limitation based on intracellular N:C quota
                qlimitFac = recom_limiter(NMinSlope, NCmin, quota)

                ! Iron limitation using Michaelis-Menten kinetics
                ! Iron is often the limiting micronutrient in HNLC regions
                feLimitFac = Fe / (k_Fe + Fe)

                ! Apply Liebig's law: most limiting nutrient controls photosynthesis
                qlimitFac = min(qlimitFac, feLimitFac)

                ! Calculate maximum photosynthesis rate with temperature correction
                if (enable_coccos) then
                    ! Use species-specific temperature function (when cocco module active)
                    pMax = qlimitFac * Temp_phyto
                else
                    ! Use standard Arrhenius temperature function
                    pMax = P_cm * qlimitFac * arrFunc
                end if

                !-------------------------------------------------------------------------------
                ! Diatoms
                !-------------------------------------------------------------------------------
                ! Limited by nitrogen, silicon, and iron
                ! Large phytoplankton with silica frustules (shells)

                ! Nitrogen limitation
                qlimitFac = recom_limiter(NMinSlope, NCmin_d, quota_dia)

                ! Silicon limitation (unique to diatoms)
                ! Required for frustule formation - critical for diatom growth
                qlimitFacTmp = recom_limiter(SiMinSlope, SiCmin, qSiC)
                qlimitFac = min(qlimitFac, qlimitFacTmp)

                ! Iron limitation
                feLimitFac = Fe / (k_Fe_d + Fe)
                qlimitFac = min(qlimitFac, feLimitFac)

                ! Calculate maximum photosynthesis rate
                if (enable_coccos) then
                    pMax_dia = qlimitFac * Temp_diatoms
                else
                    pMax_dia = P_cm_d * qlimitFac * arrFunc
                end if

                !-------------------------------------------------------------------------------
                ! Coccolithophores (Optional)
                !-------------------------------------------------------------------------------
                ! Limited by nitrogen and iron
                ! Calcifying phytoplankton that produce calcite plates

                if (enable_coccos) then

                    ! Nitrogen limitation
                    qlimitFac = recom_limiter(NMinSlope, NCmin_c, quota_cocco)

                    ! Iron limitation
                    feLimitFac = Fe / (k_Fe_c + Fe)
                    qlimitFac = min(qlimitFac, feLimitFac)

                    ! Calculate maximum photosynthesis rate
                    pMax_cocco = qlimitFac * Temp_cocco

                    !---------------------------------------------------------------------------
                    ! Phaeocystis (Optional)
                    !---------------------------------------------------------------------------
                    ! Limited by nitrogen and iron
                    ! Colonial phytoplankton that can form large blooms

                    ! Nitrogen limitation
                    qlimitFac = recom_limiter(NMinSlope, NCmin_p, quota_phaeo)

                    ! Iron limitation
                    feLimitFac = Fe / (k_Fe_p + Fe)
                    qlimitFac = min(qlimitFac, feLimitFac)

                    ! Calculate maximum photosynthesis rate
                    pMax_phaeo = qlimitFac * Temp_phaeo

                end if

                !===============================================================================
                ! LIGHT-DEPENDENT PHOTOSYNTHESIS RATE CALCULATIONS
                !===============================================================================
                ! Calculates actual photosynthesis rates using photosynthesis-irradiance (P-I)
                ! curves. Uses exponential saturation model (no photoinhibition).
                !
                ! P-I Curve Model:
                !   P = Pmax × (1 - exp(-α × Chl:C × PAR / Pmax))
                !
                ! Where:
                !   P     : Actual photosynthesis rate [day-1]
                !   Pmax  : Maximum rate (nutrient and temperature limited) [day-1]
                !   α     : Initial slope of P-I curve (photosynthetic efficiency) [-]
                !   Chl:C : Chlorophyll to carbon ratio [mgChl mmolC-1]
                !   PAR   : Photosynthetically active radiation [W m-2]
                !
                ! Variables (Small Phytoplankton):
                !   Cphot               : Carbon-specific photosynthesis rate [day-1]
                !   PARave              : Average PAR in mixed layer [W m-2]
                !   alfa                : Initial slope parameter [-]
                !   Chl2C               : Chlorophyll:Carbon ratio [mgChl mmolC-1]
                !   PhyCO2              : CO2 limitation factor [0-1]
                !   VTCphotLigLim_phyto : Light limitation factor for diagnostics [0-1]
                !   VTCphot_phyto       : Final photosynthesis rate for diagnostics [day-1]
                !
                ! Safety Checks:
                !   - Check for darkness (pMax < tiny)
                !   - Check for NaN values (PARave /= PARave)
                !   - Check for valid Chl:C ratios
                !
                ! Note: Similar calculations for all phytoplankton types
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Small Phytoplankton and Diatom Photosynthesis
                !-------------------------------------------------------------------------------
                call calculate_photosynthesis_rate(pMax, PARave, Chl2C, PhyCO2, alfa, Cphot, &
                        VTCphotLigLim_phyto(k))
                call calculate_photosynthesis_rate(pMax_dia, PARave, Chl2C_dia, DiaCO2, &
                        alfa_d, Cphot_dia, VTCphotLigLim_diatoms(k))
                !! Store photosynthesis rate for diagnostics and output
                VTCphot_phyto(k) = Cphot
                VTCphot_diatoms(k) = Cphot_dia

                if (enable_coccos) then
                    !-------------------------------------------------------------------------------
                    ! Coccolithophore and Phaeocystis Photosynthesis (Optional)
                    !-------------------------------------------------------------------------------
                    call calculate_photosynthesis_rate(pMax_cocco, PARave, Chl2C_cocco, CoccoCO2, &
                            alfa_c, Cphot_cocco, VTCphotLigLim_cocco(k))
                    call calculate_photosynthesis_rate(pMax_phaeo, PARave, Chl2C_phaeo, PhaeoCO2, &
                            alfa_p, Cphot_phaeo, VTCphotLigLim_phaeo(k))

                    !! Store photosynthesis rate for diagnostics and output
                    VTCphot_cocco(k) = Cphot_cocco
                    VTCphot_phaeo(k) = Cphot_phaeo
                end if

                !===============================================================================
                ! CHLOROPHYLL DEGRADATION
                !===============================================================================
                ! Calculates chlorophyll degradation rates with optional photodamage effects.
                ! Chlorophyll degrades due to senescence and light-induced damage.
                !
                ! Two Modes:
                !   A) Base degradation: Constant rate (use_photodamage = FALSE)
                !   B) Photodamage-dependent: Light-dependent rate (use_photodamage = TRUE)
                !
                ! Variables:
                !   KOchl           : Chlorophyll degradation rate for small phyto [day-1]
                !   KOchl_dia       : Chlorophyll degradation rate for diatoms [day-1]
                !   KOchl_cocco     : Chlorophyll degradation rate for coccos [day-1]
                !   KOchl_phaeo     : Chlorophyll degradation rate for Phaeocystis [day-1]
                !   deg_Chl         : Base degradation rate constant [day-1]
                !   CHL2C_plast     : Chlorophyll:Carbon ratio in plastids [mgChl mmolC-1]
                !   alfa            : P-I curve initial slope [-]
                !   PARave          : Average PAR [W m-2]
                !   pMax            : Maximum photosynthesis rate [day-1]
                !
                ! Photodamage Model:
                !   - High light damages photosystem II and degrades chlorophyll
                !   - Uses same exponential form as P-I curve (saturating damage)
                !   - Minimum degradation rate (10% of base) in darkness
                !   - No upper cap (commented out safety constraint)
                !
                ! Note: Photodamage increases chlorophyll turnover at high light
                !       This helps prevent photoinhibition and photooxidative stress
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Set Base Chlorophyll Degradation Rates
                !-------------------------------------------------------------------------------
                ! Constant degradation rates independent of light (senescence)

                KOchl = deg_Chl ! Small phytoplankton
                KOchl_dia = deg_Chl_d ! Diatoms

                if (enable_coccos) then
                    KOchl_cocco = deg_Chl_c ! Coccolithophores
                    KOchl_phaeo = deg_Chl_p ! Phaeocystis
                end if

                if (use_photodamage) then

                    !===========================================================================
                    ! PHOTODAMAGE-DEPENDENT DEGRADATION
                    !===========================================================================
                    ! Light-dependent chlorophyll degradation using saturation model
                    ! Higher light intensity increases chlorophyll turnover
                    !---------------------------------------------------------------------------

                    !---------------------------------------------------------------------------
                    ! Small Phytoplankton Chlorophyll Loss
                    !---------------------------------------------------------------------------

                    call calculate_chlorophyll_degradation(pMax, PARave, CHL2C_plast, &
                            deg_Chl, alfa, KOchl)

                    !---------------------------------------------------------------------------
                    ! Diatom Chlorophyll Loss
                    !---------------------------------------------------------------------------

                    call calculate_chlorophyll_degradation(pMax_dia, PARave, CHL2C_plast_dia, &
                            deg_Chl_d, alfa_d, KOchl_dia)

                    if (enable_coccos) then

                        !-----------------------------------------------------------------------
                        ! Coccolithophore Chlorophyll Loss
                        !-----------------------------------------------------------------------

                        call calculate_chlorophyll_degradation(pMax_cocco, PARave, &
                                CHL2C_plast_cocco, deg_Chl_c, &
                                alfa_c, KOchl_cocco)

                        !-----------------------------------------------------------------------
                        ! Phaeocystis Chlorophyll Loss
                        !-----------------------------------------------------------------------

                        call calculate_chlorophyll_degradation(pMax_phaeo, PARave, &
                                CHL2C_plast_phaeo, deg_Chl_p, &
                                alfa_p, KOchl_phaeo)

                    end if ! enable_coccos

                end if ! use_photodamage

                !---------------------------------------------------------------------------
                ! ERROR CHECKING AND DEBUGGING
                ! Check for NaN values in chlorophyll degradation rates
                !---------------------------------------------------------------------------

                ! Small phytoplankton
                if (ieee_is_nan(KOchl)) then
                    print*, 'ERROR: KOchl is NaN'
                    print*, '  deg_Chl =', deg_Chl
                    print*, '  alfa =', alfa
                    print*, '  CHL2C_plast =', CHL2C_plast
                    print*, '  PARave =', PARave
                    print*, '  pMax =', pMax
                    stop
                end if

                ! Diatoms
                if (ieee_is_nan(KOchl_dia)) then
                    print*, 'ERROR: KOchl_dia is NaN'
                    print*, '  deg_Chl_d =', deg_Chl_d
                    print*, '  alfa_d =', alfa_d
                    print*, '  CHL2C_plast_dia =', CHL2C_plast_dia
                    print*, '  PARave =', PARave
                    print*, '  pMax_dia =', pMax_dia
                    stop
                end if

                ! Additional species (if enabled)
                if (enable_coccos) then
                    if (ieee_is_nan(KOchl_cocco)) then
                        print*, 'ERROR: KOchl_cocco is NaN'
                        print*, '  deg_Chl_c =', deg_Chl_c
                        print*, '  alfa_c =', alfa_c
                        print*, '  CHL2C_plast_cocco =', CHL2C_plast_cocco
                        print*, '  PARave =', PARave
                        print*, '  pMax_cocco =', pMax_cocco
                        stop
                    end if

                    if (ieee_is_nan(KOchl_phaeo)) then
                        print*, 'ERROR: KOchl_phaeo is NaN'
                        print*, '  deg_Chl_p =', deg_Chl_p
                        print*, '  alfa_p =', alfa_p
                        print*, '  CHL2C_plast_phaeo =', CHL2C_plast_phaeo
                        print*, '  PARave =', PARave
                        print*, '  pMax_phaeo =', pMax_phaeo
                        stop
                    end if
                end if

                !===============================================================================
                ! PHYTOPLANKTON ASSIMILATION SECTION
                !===============================================================================
                ! Computes nutrient uptake, chlorophyll synthesis, and respiration rates
                ! for multiple phytoplankton functional types following Geider et al. 1998
                !===============================================================================

                !===============================================================================
                ! NITROGEN ASSIMILATION
                !===============================================================================
                ! Calculates nitrogen uptake rates for all phytoplankton groups based on
                ! nutrient availability, cell quota status, and maximum uptake capacity
                !
                ! Key Parameters:
                !   V_cm_fact       : Scaling factor for C-specific N uptake [-]
                !   NCmax           : Maximum cell quota of nitrogen (N:C) [mmolN mmolC-1]
                !   NMaxSlope       : Maximum slope for quota limiting function [-]
                !   NCuptakeRatio   : Maximum uptake ratio N:C [mmolN mmolC-1]
                !   k_din           : Half-saturation constant for DIN uptake [mmolN m-3]
                !   pMax            : Maximum photosynthesis rate [day-1]
                !
                ! Note: Cell quota limiting function prevents luxury N uptake when internal
                !       N:C ratio approaches maximum capacity
                !-------------------------------------------------------------------------------

                ! --- Small phytoplankton Nitrogen Uptake ---
                V_cm = V_cm_fact
                limitFacN = recom_limiter(NMaxSlope, quota, NCmax)
                N_assim = V_cm * pMax * NCuptakeRatio * limitFacN * (DIN / (DIN + k_din))

                ! --- Diatom Nitrogen Uptake ---
                V_cm = V_cm_fact_d
                limitFacN_dia = recom_limiter(NMaxSlope, quota_dia, NCmax_d)
                N_assim_dia = V_cm * pMax_dia * NCUptakeRatio_d * limitFacN_dia * DIN &
                        / (DIN + k_din_d)

                ! --- Optional Coccolithophore and Phaeocystis Groups ---
                if (enable_coccos) then
                    ! Coccolithophore nitrogen uptake
                    V_cm = V_cm_fact_c
                    limitFacN_cocco = recom_limiter(NMaxSlope, quota_cocco, NCmax_c)
                    N_assim_cocco = V_cm * pMax_cocco * NCUptakeRatio_c * limitFacN_cocco * &
                            DIN / (DIN + k_din_c)

                    ! Phaeocystis nitrogen uptake
                    V_cm = V_cm_fact_p
                    limitFacN_phaeo = recom_limiter(NMaxSlope, quota_phaeo, NCmax_p)
                    N_assim_phaeo = V_cm * pMax_phaeo * NCUptakeRatio_p * limitFacN_phaeo * &
                            DIN / (DIN + k_din_p)
                end if

                !===============================================================================
                ! SILICON ASSIMILATION (DIATOMS ONLY)
                !===============================================================================
                ! Calculates silicate uptake for diatom frustule formation
                !
                ! Key Parameters:
                !   SiCUptakeRatio  : Maximum uptake ratio Si:C [mmolSi mmolC-1]
                !   SiCmax          : Maximum cell quota of silicon (Si:C) [mmolSi mmolC-1]
                !   SiMaxSlope      : Maximum slope for Si quota limiting function [-]
                !   k_si            : Half-saturation constant for Si uptake [mmolSi m-3]
                !   P_cm_d          : Maximum photosynthesis rate for diatoms [day-1]
                !
                ! Note: Silicon uptake is coupled to nitrogen status - high N:C ratios indicate
                !       low intracellular energy reserves, limiting energy-intensive Si uptake
                !-------------------------------------------------------------------------------

                limitFacSi = recom_limiter(SiMaxSlope, qSiC, SiCmax) * limitFacN_dia

                if (.not.enable_coccos) then
                    ! Standard silicon assimilation formulation
                    Si_assim = V_cm_fact_d * P_cm_d * arrFunc * SiCUptakeRatio * limitFacSi * &
                            Si / (Si + k_si)
                else
                    ! Alternative formulation with temperature dependence
                    Si_assim = V_cm_fact_d * Temp_diatoms * SiCUptakeRatio * limitFacSi * &
                            Si / (Si + k_si)
                    VTSi_assimDia(k) = Si_assim
                end if

                !===============================================================================
                ! 3. IRON CHEMISTRY
                !===============================================================================
                ! Computes free (bioavailable) iron concentration from total dissolved iron
                ! accounting for complexation with organic ligands
                !
                ! Variables:
                !   Fe                  : Total dissolved iron [µmol m-3]
                !   totalligand         : Total organic ligand concentration [µmol m-3]
                !   ligandStabConst     : Conditional stability constant [M-1]
                !   freeFe              : Free (inorganic) iron concentration [µmol m-3]
                !-------------------------------------------------------------------------------

                freeFe = iron_chemistry(Fe, totalligand, ligandStabConst)

                !===============================================================================
                ! 4. CHLOROPHYLL SYNTHESIS
                !===============================================================================
                ! Calculates chlorophyll production coupled to nitrogen assimilation
                ! following photoacclimation theory (Geider et al. 1998)
                !
                ! Key Parameters:
                !   Chl2N_max       : Maximum Chl:N ratio [mg Chl mmolN-1]
                !   Chl2C           : Chlorophyll to carbon ratio [mg Chl mmolC-1]
                !   alfa            : Initial slope of P-I curve [mmolC (mg Chl)-1 m2 µmol-1
                ! day-1]
                !   Cphot           : Carbon-specific photosynthesis rate [day-1]
                !   PARave          : Depth-averaged photosynthetically available radiation [µmol
                ! m-2 s-1]
                !
                ! Note: Chlorophyll synthesis is down-regulated when light is sufficient,
                !       preventing over-investment in light-harvesting machinery
                !-------------------------------------------------------------------------------

                ! --- Small phytoplankton Chlorophyll Synthesis ---
                chlSynth = zero
                if (PARave >= tiny .and. .not.ieee_is_nan(PARave)) then
                    chlSynth = N_assim * Chl2N_max * &
                            min(real(one), Cphot / (alfa * Chl2C * PARave))
                end if

                ! --- Diatom Chlorophyll Synthesis ---
                ChlSynth_dia = zero
                if (PARave >= tiny .and. .not.ieee_is_nan(PARave)) then
                    ChlSynth_dia = N_assim_dia * Chl2N_max_d * &
                            min(real(one), Cphot_dia / (alfa_d * Chl2C_dia * PARave))
                end if

                ! --- Optional Coccolithophore and Phaeocystis Chlorophyll Synthesis ---
                if (enable_coccos) then
                    ! Coccolithophore chlorophyll synthesis
                    ChlSynth_cocco = zero
                    if (PARave >= tiny .and. .not.ieee_is_nan(PARave)) then
                        ChlSynth_cocco = N_assim_cocco * Chl2N_max_c * &
                                min(real(one), Cphot_cocco / (alfa_c * Chl2C_cocco * PARave))
                    end if

                    ! Phaeocystis chlorophyll synthesis
                    ChlSynth_phaeo = zero
                    if (PARave >= tiny .and. .not.ieee_is_nan(PARave)) then
                        ChlSynth_phaeo = N_assim_phaeo * Chl2N_max_p * &
                                min(real(one), Cphot_phaeo / (alfa_p * Chl2C_phaeo * PARave))
                    end if
                end if

                !===============================================================================
                ! 5. PHYTOPLANKTON RESPIRATION RATES
                !===============================================================================
                ! Computes carbon loss through maintenance respiration and biosynthetic costs
                !
                ! Key Parameters:
                !   res_phy         : Maintenance respiration rate constant [day-1]
                !   biosynth        : Biosynthetic cost of nitrogen assimilation [mmolC mmolN-1]
                !   biosynthSi      : Biosynthetic cost of silicon assimilation [mmolC mmolSi-1]
                !
                ! Components:
                !   1. Maintenance respiration: Quota-dependent baseline metabolic cost
                !   2. Biosynthetic costs: Additional respiration for nutrient assimilation
                !
                ! Equation: R = res_phy * limitFacN + biosynth * N_assim + biosynthSi * Si_assim
                !-------------------------------------------------------------------------------

                ! --- Small phytoplankton Respiration ---
                phyRespRate = res_phy * limitFacN + biosynth * N_assim

                ! --- Diatom Respiration (includes silicon biosynthesis cost) ---
                phyRespRate_dia = res_phy_d * limitFacN_dia + biosynth * N_assim_dia + &
                        biosynthSi * Si_assim

                ! --- Optional Coccolithophore and Phaeocystis Respiration ---
                if (enable_coccos) then
                    phyRespRate_cocco = res_phy_c * limitFacN_cocco + biosynth * N_assim_cocco
                    phyRespRate_phaeo = res_phy_p * limitFacN_phaeo + biosynth * N_assim_phaeo
                end if

                !===============================================================================
                ! MESOZOOPLANKTON GRAZING
                !===============================================================================
                ! Simulates mesozooplankton grazing on multiple prey types using a Holling
                ! Type III functional response with food-dependent preferences and efficiency.
                !
                ! This module calculates:
                !   1. Food availability and grazing preferences (fixed or variable)
                !   2. Total grazing flux with Holling Type III response
                !   3. Distribution of grazing among prey types
                !   4. Food-dependent grazing efficiency
                !   5. Carbon flux from grazed prey to mesozooplankton
                !
                ! Key Features:
                !   - Variable or fixed prey preferences
                !   - Optional prey types (coccolithophores, microzooplankton, detritus)
                !   - Temperature-dependent grazing rate (Q10 or Arrhenius)
                !   - Food-dependent assimilation efficiency
                !
                ! References:
                !   - Schourup-Kristensen et al. (2013) - REcoM model description
                !===============================================================================

                !===============================================================================
                ! 1. FOOD AVAILABILITY AND GRAZING PREFERENCES
                !===============================================================================
                ! Calculates which prey types are available and their relative preferences.
                ! Two modes:
                !   - Variable preferences: Adjust based on relative prey abundance
                !   - Fixed preferences: Use constant maximum preference values
                !
                ! Variables:
                !   pzPhy, pzDia    : Max preference for small phyto and diatoms [-]
                !   pzCocco, pzPhaeo: Max preference for coccoliths and Phaeocystis [-]
                !   pzDet, pzDetZ2  : Max preference for slow/fast sinking detritus [-]
                !   pzMicZoo        : Max preference for microzooplankton [-]
                !   PhyN, DiaN      : Small phytoplankton and diatom nitrogen [mmolN m-3]
                !   CoccoN, PhaeoN  : Coccolithophore and Phaeocystis nitrogen [mmolN m-3]
                !   DetN, DetZ2N    : Slow and fast sinking detritus nitrogen [mmolN m-3]
                !   MicZooN         : Microzooplankton nitrogen [mmolN m-3]
                !   varpz*          : Variable preferences (calculated from availability) [-]
                !   f*N             : Available food pools (preference × concentration) [mmolN
                ! m-3]
                !   aux             : Total weighted food availability [mmolN m-3]
                !
                ! Logic:
                !   IF variable preferences: varpz_i = (pz_i × prey_i) / Σ(pz_j × prey_j)
                !   ELSE: Use fixed maximum preferences (pz_i)
                !   Food pools: f_i = preference_i × prey_i
                !-------------------------------------------------------------------------------

                if (REcoM_Grazing_Variable_Preference) then
                    !---------------------------------------------------------------------------
                    ! VARIABLE PREFERENCE MODE
                    ! Preferences scale with relative abundance of each prey type
                    !---------------------------------------------------------------------------

                    ! Calculate total weighted food availability (denominator)
                    aux = pzPhy * PhyN + pzDia * DiaN

                    if (Grazing_detritus) then
                        aux = aux + pzDet * DetN
                    end if

                    if (enable_3zoo2det) then
                        if (Grazing_detritus) aux = aux + pzDetZ2 * DetZ2N ! Fast-sinking detritus
                        aux = aux + pzMicZoo * MicZooN ! Microzooplankton
                    end if

                    if (enable_coccos) then
                        aux = aux + pzCocco * CoccoN + pzPhaeo * PhaeoN
                    end if

                    ! Calculate variable preferences (normalized by total availability)
                    varpzPhy = (pzPhy * PhyN) / aux
                    varpzDia = (pzDia * DiaN) / aux

                    if (Grazing_detritus) then
                        varpzDet = (pzDet * DetN) / aux
                    end if

                    if (enable_3zoo2det) then
                        if (Grazing_detritus) varpzDetZ2 = (pzDetZ2 * DetZ2N) / aux
                        varpzMicZoo = (pzMicZoo * MicZooN) / aux
                    end if

                    if (enable_coccos) then
                        varpzCocco = (pzCocco * CoccoN) / aux
                        varpzPhaeo = (pzPhaeo * PhaeoN) / aux
                    end if

                    ! Calculate available food pools (preference × concentration)
                    fPhyN = varpzPhy * PhyN
                    fDiaN = varpzDia * DiaN

                    if (Grazing_detritus) then
                        fDetN = varpzDet * DetN
                    end if

                    if (enable_3zoo2det) then
                        if (Grazing_detritus) fDetZ2N = varpzDetZ2 * DetZ2N
                        fMicZooN = varpzMicZoo * MicZooN
                    end if

                    if (enable_coccos) then
                        fCoccoN = varpzCocco * CoccoN
                        fPhaeoN = varpzPhaeo * PhaeoN
                    end if

                else
                    !---------------------------------------------------------------------------
                    ! FIXED PREFERENCE MODE
                    ! Use constant maximum preference values
                    !---------------------------------------------------------------------------

                    fPhyN = pzPhy * PhyN
                    fDiaN = pzDia * DiaN

                    if (Grazing_detritus) then
                        fDetN = pzDet * DetN
                    end if

                    if (enable_3zoo2det) then
                        if (Grazing_detritus) fDetZ2N = pzDetZ2 * DetZ2N
                        fMicZooN = pzMicZoo * MicZooN
                    end if

                    if (enable_coccos) then
                        fCoccoN = pzCocco * CoccoN
                        fPhaeoN = pzPhaeo * PhaeoN
                    end if

                end if ! REcoM_Grazing_Variable_Preference

                !===============================================================================
                ! 2. TOTAL GRAZING FLUX (HOLLING TYPE III)
                !===============================================================================
                ! Calculates total grazing rate using a sigmoidal (Type III) functional response.
                ! This creates a threshold effect where grazing accelerates at higher food levels.
                !
                ! Variables:
                !   food            : Total available food [mmolN m-3]
                !   foodsq          : Squared food concentration [mmolN2 m-6]
                !   grazingFlux     : Total N grazing rate [mmolN m-3 day-1]
                !   Graz_max        : Maximum specific grazing rate [day-1]
                !   epsilonr        : Half-saturation constant squared [mmolN2 m-6]
                !   HetN            : Mesozooplankton nitrogen concentration [mmolN m-3]
                !   q10_mes         : Q10 temperature function for mesozooplankton [-]
                !   arrFunc         : Arrhenius temperature function [-]
                !
                ! Equation: Holling Type III
                !   grazingFlux = (Graz_max × food²) / (epsilonr + food²) × HetN × T_func
                !
                ! Note: Uses Q10 when 3-zoo/2-detritus enabled, otherwise uses Arrhenius
                !-------------------------------------------------------------------------------

                ! Sum all available food pools
                food = fPhyN + fDiaN

                if (Grazing_detritus) then
                    food = food + fDetN
                end if

                if (enable_3zoo2det) then
                    if (Grazing_detritus) food = food + fDetZ2N
                    food = food + fMicZooN
                end if

                if (enable_coccos) then
                    food = food + fCoccoN + fPhaeoN
                end if

                ! Calculate grazing flux with Holling Type III functional response
                foodsq = food ** 2

                if (enable_3zoo2det) then
                    grazingFlux = (Graz_max * foodsq) / (epsilonr + foodsq) * HetN * q10_mes
                else
                    grazingFlux = (Graz_max * foodsq) / (epsilonr + foodsq) * HetN * arrFunc
                end if

                !===============================================================================
                ! 3. GRAZING FLUX DISTRIBUTION
                !===============================================================================
                ! Partitions total grazing among prey types proportional to their availability.
                !
                ! Variables:
                !   grazingFlux_phy    : Grazing on small phytoplankton [mmolN m-3 day-1]
                !   grazingFlux_Dia    : Grazing on diatoms [mmolN m-3 day-1]
                !   grazingFlux_Det    : Grazing on slow-sinking detritus [mmolN m-3 day-1]
                !   grazingFlux_DetZ2  : Grazing on fast-sinking detritus [mmolN m-3 day-1]
                !   grazingFlux_miczoo : Grazing on microzooplankton [mmolN m-3 day-1]
                !   grazingFlux_Cocco  : Grazing on coccolithophores [mmolN m-3 day-1]
                !   grazingFlux_Phaeo  : Grazing on Phaeocystis [mmolN m-3 day-1]
                !
                ! Equation for each prey type i:
                !   grazingFlux_i = grazingFlux × (f_i / total_food)
                !-------------------------------------------------------------------------------

                grazingFlux_phy = grazingFlux * fPhyN / food
                grazingFlux_Dia = grazingFlux * fDiaN / food

                if (Grazing_detritus) then
                    grazingFlux_Det = grazingFlux * fDetN / food
                end if

                if (enable_3zoo2det) then
                    if (Grazing_detritus) grazingFlux_DetZ2 = grazingFlux * fDetZ2N / food
                    grazingFlux_miczoo = grazingFlux * fMicZooN / food
                end if

                if (enable_coccos) then
                    grazingFlux_Cocco = grazingFlux * fCoccoN / food
                    grazingFlux_Phaeo = grazingFlux * fPhaeoN / food
                end if

                !===============================================================================
                ! 4. GRAZING EFFICIENCY AND CARBON FLUX
                !===============================================================================
                ! Calculates food-dependent assimilation efficiency and converts grazed nitrogen
                ! to carbon flux using prey-specific C:N ratios.
                !
                ! Variables:
                !   grazEff               : Grazing/assimilation efficiency [-]
                !   gfin                  : Baseline grazing efficiency [-]
                !   grazingFluxcarbon_mes : Total carbon flux to mesozooplankton [mmolC m-3 day-1]
                !   recipQuota            : Small phytoplankton C:N ratio [mmolC mmolN-1]
                !   recipQuota_Dia        : Diatom C:N ratio [mmolC mmolN-1]
                !   recipQuota_Cocco      : Coccolithophore C:N ratio [mmolC mmolN-1]
                !   recipQuota_phaeo      : Phaeocystis C:N ratio [mmolC mmolN-1]
                !   recipDet, recipDet2   : Detritus C:N ratios [mmolC mmolN-1]
                !   recipQZoo3            : Microzooplankton C:N ratio [mmolC mmolN-1]
                !
                ! Grazing Efficiency Equation:
                !   grazEff = gfin + 1/(0.2×food + 2)
                !   Higher food -> higher efficiency (asymptotes to gfin + 0.5)
                !
                ! Carbon Flux Equation:
                !   C_flux = Σ(grazingFlux_i × C:N_ratio_i × grazEff)
                !-------------------------------------------------------------------------------

                ! Calculate food-dependent grazing efficiency
                ! Increases with food availability, representing improved assimilation at higher
                ! rations
                grazEff = gfin + 1.0 / (0.2 * food + 2.0)

                ! Convert grazed nitrogen to carbon flux using prey C:N ratios
                grazingFluxcarbon_mes = (grazingFlux_phy * recipQuota * grazEff) + &
                        (grazingFlux_Dia * recipQuota_Dia * grazEff)

                if (Grazing_detritus) then
                    grazingFluxcarbon_mes = grazingFluxcarbon_mes + &
                            (grazingFlux_Det * recipDet * grazEff)
                end if

                if (enable_3zoo2det) then
                    if (Grazing_detritus) then
                        grazingFluxcarbon_mes = grazingFluxcarbon_mes + &
                                (grazingFlux_DetZ2 * recipDet2 * grazEff)
                    end if
                    grazingFluxcarbon_mes = grazingFluxcarbon_mes + &
                            (grazingFlux_miczoo * recipQZoo3 * grazEff)
                end if

                if (enable_coccos) then
                    grazingFluxcarbon_mes = grazingFluxcarbon_mes + &
                            (grazingFlux_Cocco * recipQuota_Cocco * grazEff) + &
                            (grazingFlux_Phaeo * recipQuota_phaeo * grazEff)
                end if

                !===============================================================================
                ! MACROZOOPLANKTON GRAZING
                !===============================================================================
                ! Simulates macrozooplankton (second zooplankton) grazing on multiple prey types
                ! using a Holling Type II functional response with food-dependent preferences.
                !
                ! This module calculates:
                !   1. Food availability and grazing preferences (fixed or variable)
                !   2. Total grazing flux with Holling Type II response
                !   3. Distribution of grazing among prey types
                !   4. Carbon assimilation from grazed prey using C:N ratios
                !
                ! Key Features:
                !   - Carnivorous feeding: grazes on mesozooplankton and microzooplankton
                !   - Herbivorous feeding: grazes on phytoplankton, diatoms, coccolithophores
                !   - Detritivorous feeding: optional grazing on slow/fast-sinking detritus
                !   - Variable or fixed prey preferences
                !   - Temperature-dependent grazing rate (Arrhenius)
                !   - Constant assimilation efficiency (grazEff2)
                !
                ! Prey Types:
                !   - Primary: Mesozooplankton, microzooplankton, small phyto, diatoms
                !   - Optional: Coccolithophores, Phaeocystis, slow/fast-sinking detritus
                !
                ! References:
                !   - Schourup-Kristensen et al. (2013) - REcoM model description
                !===============================================================================

                if (enable_3zoo2det) then

                    !===========================================================================
                    ! 1. FOOD AVAILABILITY AND GRAZING PREFERENCES
                    !===========================================================================
                    ! Calculates which prey types are available and their relative preferences.
                    ! Two modes:
                    !   - Variable preferences: Adjust based on relative prey abundance
                    !   - Fixed preferences: Use constant maximum preference values
                    !
                    ! Variables:
                    !   pzPhy2, pzDia2   : Max preference for small phyto and diatoms [-]
                    !   pzCocco2, pzPhaeo2: Max preference for coccoliths and Phaeocystis [-]
                    !   pzHet            : Max preference for mesozooplankton [-]
                    !   pzMicZoo2        : Max preference for microzooplankton [-]
                    !   pzDet2, pzDetZ22 : Max preference for slow/fast sinking detritus [-]
                    !   PhyN, DiaN       : Small phytoplankton and diatom nitrogen [mmolN m-3]
                    !   CoccoN, PhaeoN   : Coccolithophore and Phaeocystis nitrogen [mmolN m-3]
                    !   HetN             : Mesozooplankton nitrogen [mmolN m-3]
                    !   MicZooN          : Microzooplankton nitrogen [mmolN m-3]
                    !   DetN, DetZ2N     : Slow and fast sinking detritus nitrogen [mmolN m-3]
                    !   varpz*2, varpzHet: Variable preferences (calculated from availability) [-]
                    !   f*N2, fHetN      : Available food pools (preference × concentration)
                    ! [mmolN
                    ! m-3]
                    !   aux              : Total weighted food availability [mmolN m-3]
                    !
                    ! Logic:
                    !   IF variable preferences: varpz_i = (pz_i × prey_i) / Σ(pz_j × prey_j)
                    !   ELSE: Use fixed maximum preferences (pz_i)
                    !   Food pools: f_i = preference_i × prey_i
                    !---------------------------------------------------------------------------

                    if (REcoM_Grazing_Variable_Preference) then
                        !-----------------------------------------------------------------------
                        ! VARIABLE PREFERENCE MODE
                        ! Preferences scale with relative abundance of each prey type
                        !-----------------------------------------------------------------------

                        ! Calculate total weighted food availability (denominator)
                        ! Core prey: phytoplankton, diatoms, meso- and microzooplankton
                        aux = pzPhy2 * PhyN + pzDia2 * DiaN + pzHet * HetN + pzMicZoo2 * MicZooN

                        ! Add detritus pools if detrital grazing is enabled
                        if (Grazing_detritus) then
                            aux = aux + pzDet2 * DetN + pzDetZ22 * DetZ2N
                        end if

                        ! Add coccolithophores and Phaeocystis if enabled
                        if (enable_coccos) then
                            aux = aux + pzCocco2 * CoccoN + pzPhaeo2 * PhaeoN
                        end if

                        ! Calculate variable preferences (normalized by total availability)
                        ! Each preference = (max_pref × prey_conc) / total_weighted_food
                        varpzPhy2 = (pzPhy2 * PhyN) / aux
                        varpzDia2 = (pzDia2 * DiaN) / aux
                        varpzMicZoo2 = (pzMicZoo2 * MicZooN) / aux
                        varpzHet = (pzHet * HetN) / aux

                        if (enable_coccos) then
                            varpzCocco2 = (pzCocco2 * CoccoN) / aux
                            varpzPhaeo2 = (pzPhaeo2 * PhaeoN) / aux
                        end if

                        if (Grazing_detritus) then
                            varpzDet2 = (pzDet2 * DetN) / aux
                            varpzDetZ22 = (pzDetZ22 * DetZ2N) / aux
                        end if

                        ! Calculate available food pools (preference × concentration)
                        fPhyN2 = varpzPhy2 * PhyN
                        fDiaN2 = varpzDia2 * DiaN
                        fMicZooN2 = varpzMicZoo2 * MicZooN
                        fHetN = varpzHet * HetN

                        if (enable_coccos) then
                            fCoccoN2 = varpzCocco2 * CoccoN
                            fPhaeoN2 = varpzPhaeo2 * PhaeoN
                        end if

                        if (Grazing_detritus) then
                            fDetN2 = varpzDet2 * DetN
                            fDetZ2N2 = varpzDetZ22 * DetZ2N
                        end if

                    else
                        !-----------------------------------------------------------------------
                        ! FIXED PREFERENCE MODE
                        ! Use constant maximum preference values
                        !-----------------------------------------------------------------------

                        fPhyN2 = pzPhy2 * PhyN
                        fDiaN2 = pzDia2 * DiaN
                        fMicZooN2 = pzMicZoo2 * MicZooN
                        fHetN = pzHet * HetN

                        if (enable_coccos) then
                            fCoccoN2 = pzCocco2 * CoccoN
                            fPhaeoN2 = pzPhaeo2 * PhaeoN
                        end if

                        if (Grazing_detritus) then
                            fDetN2 = pzDet2 * DetN
                            fDetZ2N2 = pzDetZ22 * DetZ2N
                        end if

                    end if ! REcoM_Grazing_Variable_Preference

                    !===========================================================================
                    ! 2. TOTAL GRAZING FLUX (HOLLING TYPE II)
                    !===========================================================================
                    ! Calculates total grazing rate using a hyperbolic (Type II) functional
                    ! response. This creates saturating grazing at high food concentrations.
                    !
                    ! Variables:
                    !   food2             : Total available food [mmolN m-3]
                    !   foodsq2           : Squared food concentration [mmolN2 m-6]
                    !   grazingFlux2      : Total N grazing rate [mmolN m-3 day-1]
                    !   Graz_max2         : Maximum specific grazing rate [day-1]
                    !   epsilon2          : Half-saturation constant squared [mmolN2 m-6]
                    !   Zoo2N             : Macrozooplankton nitrogen concentration [mmolN m-3]
                    !   arrFuncZoo2       : Arrhenius temperature function for macrozooplankton [-]
                    !
                    ! Equation: Holling Type II
                    !   grazingFlux2 = (Graz_max2 × food²) / (epsilon2 + food²) × Zoo2N ×
                    ! T_func
                    !
                    ! Note: Uses Arrhenius temperature dependency (arrFuncZoo2)
                    !---------------------------------------------------------------------------

                    ! Sum all available food pools
                    food2 = fPhyN2 + fDiaN2 + fHetN + fMicZooN2

                    if (Grazing_detritus) then
                        food2 = food2 + fDetN2 + fDetZ2N2
                    end if

                    if (enable_coccos) then
                        food2 = food2 + fCoccoN2 + fPhaeoN2
                    end if

                    ! Calculate grazing flux with Holling Type II functional response
                    ! Type II uses squared food (similar to Type III but with different parameters)
                    foodsq2 = food2 ** 2
                    grazingFlux2 = (Graz_max2 * foodsq2) / (epsilon2 + foodsq2) * Zoo2N &
                            * arrFuncZoo2

                    !===========================================================================
                    ! 3. GRAZING FLUX DISTRIBUTION
                    !===========================================================================
                    ! Partitions total grazing among prey types proportional to their availability.
                    !
                    ! Variables:
                    !   grazingFlux_phy2    : Grazing on small phytoplankton [mmolN m-3 day-1]
                    !   grazingFlux_Dia2    : Grazing on diatoms [mmolN m-3 day-1]
                    !   grazingFlux_miczoo2 : Grazing on microzooplankton [mmolN m-3 day-1]
                    !   grazingFlux_het2    : Grazing on mesozooplankton [mmolN m-3 day-1]
                    !   grazingFlux_Cocco2  : Grazing on coccolithophores [mmolN m-3 day-1]
                    !   grazingFlux_Phaeo2  : Grazing on Phaeocystis [mmolN m-3 day-1]
                    !   grazingFlux_Det2    : Grazing on slow-sinking detritus [mmolN m-3 day-1]
                    !   grazingFlux_DetZ22  : Grazing on fast-sinking detritus [mmolN m-3 day-1]
                    !
                    ! Equation for each prey type i:
                    !   grazingFlux_i = grazingFlux2 × (f_i / total_food)
                    !---------------------------------------------------------------------------

                    ! Distribute total grazing among prey types proportionally
                    grazingFlux_phy2 = (grazingFlux2 * fPhyN2) / food2
                    grazingFlux_Dia2 = (grazingFlux2 * fDiaN2) / food2
                    grazingFlux_miczoo2 = (grazingFlux2 * fMicZooN2) / food2
                    grazingFlux_het2 = (grazingFlux2 * fHetN) / food2

                    if (enable_coccos) then
                        grazingFlux_Cocco2 = (grazingFlux2 * fCoccoN2) / food2
                        grazingFlux_Phaeo2 = (grazingFlux2 * fPhaeoN2) / food2
                    end if

                    if (Grazing_detritus) then
                        grazingFlux_Det2 = (grazingFlux2 * fDetN2) / food2
                        grazingFlux_DetZ22 = (grazingFlux2 * fDetZ2N2) / food2
                    end if

                    !===========================================================================
                    ! 4. CARBON ASSIMILATION
                    !===========================================================================
                    ! Converts grazed nitrogen to assimilated carbon flux using prey-specific
                    ! C:N ratios and constant grazing efficiency.
                    !
                    ! Variables:
                    !   grazingFluxcarbonzoo2 : Total carbon flux to macrozooplankton [mmolC m-3
                    ! day-1]
                    !   grazEff2              : Macrozooplankton grazing efficiency (constant) [-]
                    !   recipQuota            : Small phytoplankton C:N ratio [mmolC mmolN-1]
                    !   recipQuota_Dia        : Diatom C:N ratio [mmolC mmolN-1]
                    !   recipQuota_Cocco      : Coccolithophore C:N ratio [mmolC mmolN-1]
                    !   recipQuota_Phaeo      : Phaeocystis C:N ratio [mmolC mmolN-1]
                    !   recipQZoo             : Mesozooplankton C:N ratio [mmolC mmolN-1]
                    !   recipQZoo3            : Microzooplankton C:N ratio [mmolC mmolN-1]
                    !   recipDet, recipDet2   : Detritus C:N ratios [mmolC mmolN-1]
                    !
                    ! Carbon Flux Equation:
                    !   C_flux = Σ(grazingFlux_i × C:N_ratio_i × grazEff2)
                    !
                    ! Note: Unlike mesozooplankton, macrozooplankton uses constant efficiency
                    !       (grazEff2) rather than food-dependent efficiency
                    !---------------------------------------------------------------------------

                    ! Convert grazed nitrogen to carbon flux using prey C:N ratios
                    ! Start with core prey types (always present)
                    grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) + &
                            (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) + &
                            (grazingFlux_het2 * recipQZoo * grazEff2) + &
                            (grazingFlux_miczoo2 * recipQZoo3 * grazEff2)

                    ! Add detritus contribution if detrital grazing is enabled
                    if (Grazing_detritus) then
                        grazingFluxcarbonzoo2 = grazingFluxcarbonzoo2 + &
                                (grazingFlux_Det2 * recipDet * grazEff2) + &
                                (grazingFlux_DetZ22 * recipDet2 * grazEff2)
                    end if

                    ! Add coccolithophore and Phaeocystis contribution if enabled
                    if (enable_coccos) then
                        grazingFluxcarbonzoo2 = grazingFluxcarbonzoo2 + &
                                (grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2) + &
                                (grazingFlux_Phaeo2 * recipQuota_Phaeo * grazEff2)
                    end if

                end if ! enable_3zoo2det

                !===============================================================================
                ! MICROZOOPLANKTON GRAZING
                !===============================================================================
                ! Simulates microzooplankton (third zooplankton) grazing on phytoplankton prey
                ! using a Holling Type II functional response with food-dependent preferences.
                !
                ! This module calculates:
                !   1. Food availability and grazing preferences (fixed or variable)
                !   2. Total grazing flux with Holling Type II response
                !   3. Distribution of grazing among prey types
                !   4. Carbon assimilation from grazed prey (calculated elsewhere)
                !
                ! Key Features:
                !   - Strictly herbivorous: only grazes on phytoplankton groups
                !   - No carnivory: does not graze on other zooplankton
                !   - No detritivory: does not graze on detritus
                !   - Variable or fixed prey preferences
                !   - Temperature-dependent grazing rate (Q10)
                !   - Smallest zooplankton size class (fastest response to phytoplankton)
                !
                ! Prey Types:
                !   - Primary: Small phytoplankton, diatoms
                !   - Optional: Coccolithophores, Phaeocystis
                !
                ! Ecological Role:
                !   - Links small phytoplankton to higher trophic levels
                !   - Prey for mesozooplankton and macrozooplankton
                !   - Rapid response to phytoplankton blooms
                !
                ! References:
                !   - Schourup-Kristensen et al. (2013) - REcoM model description
                !===============================================================================

                if (enable_3zoo2det) then

                    !===========================================================================
                    ! 1. FOOD AVAILABILITY AND GRAZING PREFERENCES
                    !===========================================================================
                    ! Calculates which phytoplankton prey types are available and their
                    ! relative preferences.
                    ! Two modes:
                    !   - Variable preferences: Adjust based on relative prey abundance
                    !   - Fixed preferences: Use constant maximum preference values
                    !
                    ! Variables:
                    !   pzPhy3, pzDia3      : Max preference for small phyto and diatoms [-]
                    !   pzCocco3, pzPhaeo3  : Max preference for coccoliths and Phaeocystis [-]
                    !   PhyN, DiaN          : Small phytoplankton and diatom nitrogen [mmolN m-3]
                    !   CoccoN, PhaeoN      : Coccolithophore and Phaeocystis nitrogen [mmolN m-3]
                    !   varpzPhy3, varpzDia3: Variable preferences (calculated from availability)
                    ! [-]
                    !   varpzCocco3, varpzPhaeo3: Variable preferences for coccos/phaeo [-]
                    !   fPhyN3, fDiaN3      : Available food pools (preference × concentration)
                    ! [mmolN m-3]
                    !   fCoccoN3, fPhaeoN3  : Available coccolithophore/Phaeocystis pools [mmolN
                    ! m-3]
                    !   aux                 : Total weighted food availability [mmolN m-3]
                    !
                    ! Logic:
                    !   IF variable preferences: varpz_i = (pz_i × prey_i) / Σ(pz_j × prey_j)
                    !   ELSE: Use fixed maximum preferences (pz_i)
                    !   Food pools: f_i = preference_i × prey_i
                    !
                    ! Note: Microzooplankton have simpler feeding than meso/macrozooplankton
                    !       - No carnivory (no zooplankton prey)
                    !       - No detritivory (no detritus grazing)
                    !---------------------------------------------------------------------------

                    if (REcoM_Grazing_Variable_Preference) then
                        !-----------------------------------------------------------------------
                        ! VARIABLE PREFERENCE MODE
                        ! Preferences scale with relative abundance of each phytoplankton type
                        !-----------------------------------------------------------------------

                        ! Calculate total weighted food availability (denominator)
                        ! Core phytoplankton prey: small phytoplankton and diatoms
                        aux = pzPhy3 * PhyN + pzDia3 * DiaN

                        ! Add coccolithophores and Phaeocystis if enabled
                        if (enable_coccos) then
                            aux = aux + pzCocco3 * CoccoN + pzPhaeo3 * PhaeoN
                        end if

                        ! Calculate variable preferences (normalized by total availability)
                        ! Each preference = (max_pref × prey_conc) / total_weighted_food
                        varpzPhy3 = (pzPhy3 * PhyN) / aux
                        varpzDia3 = (pzDia3 * DiaN) / aux

                        if (enable_coccos) then
                            varpzCocco3 = (pzCocco3 * CoccoN) / aux
                            varpzPhaeo3 = (pzPhaeo3 * PhaeoN) / aux
                        end if

                        ! Calculate available food pools (preference × concentration)
                        fPhyN3 = varpzPhy3 * PhyN
                        fDiaN3 = varpzDia3 * DiaN

                        if (enable_coccos) then
                            fCoccoN3 = varpzCocco3 * CoccoN
                            fPhaeoN3 = varpzPhaeo3 * PhaeoN
                        end if

                    else
                        !-----------------------------------------------------------------------
                        ! FIXED PREFERENCE MODE
                        ! Use constant maximum preference values
                        !-----------------------------------------------------------------------

                        fPhyN3 = pzPhy3 * PhyN
                        fDiaN3 = pzDia3 * DiaN

                        if (enable_coccos) then
                            fCoccoN3 = pzCocco3 * CoccoN
                            fPhaeoN3 = pzPhaeo3 * PhaeoN
                        end if

                    end if ! REcoM_Grazing_Variable_Preference

                    !===========================================================================
                    ! 2. TOTAL GRAZING FLUX (HOLLING TYPE II)
                    !===========================================================================
                    ! Calculates total grazing rate using a hyperbolic (Type II) functional
                    ! response. This creates saturating grazing at high phytoplankton
                    ! concentrations.
                    !
                    ! Variables:
                    !   food3        : Total available phytoplankton food [mmolN m-3]
                    !   foodsq3      : Squared food concentration [mmolN2 m-6]
                    !   grazingFlux3 : Total N grazing rate [mmolN m-3 day-1]
                    !   Graz_max3    : Maximum specific grazing rate [day-1]
                    !   epsilon3     : Half-saturation constant squared [mmolN2 m-6]
                    !   MicZooN      : Microzooplankton nitrogen concentration [mmolN m-3]
                    !   q10_mic      : Q10 temperature function for microzooplankton [-]
                    !
                    ! Equation: Holling Type II
                    !   grazingFlux3 = (Graz_max3 × food³) / (epsilon3 + food²) × MicZooN ×
                    ! T_func
                    !
                    ! Note: Uses Q10 temperature dependency (q10_mic) for rapid metabolic response
                    !       Smaller organisms typically have higher Q10 sensitivity
                    !---------------------------------------------------------------------------

                    ! Sum all available phytoplankton food pools
                    food3 = fPhyN3 + fDiaN3

                    if (enable_coccos) then
                        food3 = food3 + fCoccoN3 + fPhaeoN3
                    end if

                    ! Calculate grazing flux with Holling Type II functional response
                    ! Type II creates hyperbolic saturation: grazing rate increases rapidly
                    ! at low food, then saturates at high food concentrations
                    foodsq3 = food3 ** 2
                    grazingFlux3 = (Graz_max3 * foodsq3) / (epsilon3 + foodsq3) * MicZooN * q10_mic

                    !===========================================================================
                    ! 3. GRAZING FLUX DISTRIBUTION
                    !===========================================================================
                    ! Partitions total grazing among phytoplankton prey types proportional
                    ! to their availability.
                    !
                    ! Variables:
                    !   grazingFlux_phy3   : Grazing on small phytoplankton [mmolN m-3 day-1]
                    !   grazingFlux_Dia3   : Grazing on diatoms [mmolN m-3 day-1]
                    !   grazingFlux_Cocco3 : Grazing on coccolithophores [mmolN m-3 day-1]
                    !   grazingFlux_Phaeo3 : Grazing on Phaeocystis [mmolN m-3 day-1]
                    !
                    ! Equation for each prey type i:
                    !   grazingFlux_i = grazingFlux3 × (f_i / total_food)
                    !
                    ! Note: This proportional distribution ensures conservation of mass
                    !       The sum of all distributed fluxes equals the total grazing flux
                    !---------------------------------------------------------------------------

                    ! Distribute total grazing among phytoplankton prey types proportionally
                    grazingFlux_phy3 = (grazingFlux3 * fPhyN3) / food3
                    grazingFlux_Dia3 = (grazingFlux3 * fDiaN3) / food3

                    if (enable_coccos) then
                        grazingFlux_Cocco3 = (grazingFlux3 * fCoccoN3) / food3
                        grazingFlux_Phaeo3 = (grazingFlux3 * fPhaeoN3) / food3
                    end if

                    !===========================================================================
                    ! 4. CARBON ASSIMILATION
                    !===========================================================================
                    ! Carbon flux calculation for microzooplankton is handled elsewhere in the
                    ! model code using the nitrogen grazing fluxes calculated above and
                    ! prey-specific C:N ratios.
                    !
                    ! The carbon assimilation follows the general pattern:
                    !   C_flux = Σ(grazingFlux_i × C:N_ratio_i × grazEff3)
                    !
                    ! Where:
                    !   grazEff3 = Microzooplankton grazing efficiency (constant) [-]
                    !   C:N ratios for each phytoplankton prey type [mmolC mmolN-1]
                    !
                    ! This calculation is performed in the carbon balance section of the model
                    ! to maintain consistency with other carbon flux calculations.
                    !===========================================================================

                end if ! enable_3zoo2det

                !===============================================================================
                ! MESOZOOPLANKTON RESPIRATION
                !===============================================================================
                ! Calculates mesozooplankton carbon respiration using either simple temperature-
                ! dependent respiration or Redfield-based respiration that accounts for C:N
                ! stoichiometry deviations.
                !
                ! Two respiration modes:
                !   A) Simple mode: Constant specific rate with temperature correction
                !   B) Redfield mode: Respirate excess carbon above Redfield C:N ratio
                !
                ! Variables:
                !   HetRespFlux      : Mesozooplankton respiration rate [mmolC m-3 day-1]
                !   res_het          : Baseline respiration rate constant [day-1]
                !   HetC             : Mesozooplankton carbon concentration [mmolC m-3]
                !   HetN             : Mesozooplankton nitrogen concentration [mmolN m-3]
                !   q10_mes_res      : Q10 temperature function for meso respiration [-]
                !   arrFunc          : Arrhenius temperature function [-]
                !   recip_hetN_plus  : Reciprocal of (HetN + small number) [mmolN-1 m3]
                !   redfield         : Redfield C:N ratio (106:16 ≈ 6.625) [mmolC mmolN-1]
                !   recip_res_het    : Reciprocal respiration parameter [day-1]
                !
                ! Logic:
                !   IF simple mode: RespFlux = res_het × T_func × HetC
                !   ELSE Redfield mode: RespFlux = recip_res_het × T_func × (C:N - Redfield) ×
                ! HetC
                !   Redfield mode respirates excess carbon when C:N > Redfield ratio
                !
                ! Note: Redfield-based respiration maintains stoichiometric homeostasis
                !-------------------------------------------------------------------------------

                if (het_resp_noredfield) then
                    !---------------------------------------------------------------------------
                    ! SIMPLE RESPIRATION MODE
                    ! Constant specific rate with temperature dependency
                    !---------------------------------------------------------------------------

                    if (enable_3zoo2det) then
                        ! Use Q10 temperature function (exponential temperature sensitivity)
                        HetRespFlux = res_het * q10_mes_res * HetC
                    else
                        ! Use Arrhenius temperature function
                        HetRespFlux = res_het * arrFunc * HetC
                    end if

                else
                    !---------------------------------------------------------------------------
                    ! REDFIELD-BASED RESPIRATION MODE
                    ! Respirate excess carbon above Redfield C:N ratio
                    !---------------------------------------------------------------------------

                    ! Calculate respiration based on deviation from Redfield ratio
                    ! When C:N > Redfield, organism has excess carbon that must be respired
                    ! This maintains stoichiometric balance in the organism
                    HetRespFlux = recip_res_het * arrFunc * (HetC * recip_hetN_plus - redfield) &
                            * HetC

                    ! Ensure non-negative flux (no respiration when C:N ≤ Redfield)
                    HetRespFlux = max(zero, HetRespFlux)

                end if

                !===============================================================================
                ! CARBON ISOTOPE TRACKING (OPTIONAL)
                !===============================================================================
                ! Tracks respiration fluxes of carbon isotopes (13C and 14C) for paleoclimate
                ! and carbon cycle studies.
                !
                ! Variables:
                !   HetRespFlux_13  : 13C respiration flux [mmol13C m-3 day-1]
                !   HetRespFlux_14  : 14C respiration flux [mmol14C m-3 day-1]
                !   HetC_13         : Mesozooplankton 13C concentration [mmol13C m-3]
                !   HetC_14         : Mesozooplankton 14C concentration [mmol14C m-3]
                !   HetC            : Total mesozooplankton carbon [mmolC m-3]
                !   ciso            : Flag for 13C isotope tracking [logical]
                !   ciso_14         : Flag for 14C isotope tracking [logical]
                !   ciso_organic_14 : Flag for 14C in organic matter [logical]
                !
                ! Equation:
                !   IsotopeFlux = TotalFlux × (Isotope_C / Total_C)
                !   Assumes isotope ratios are conserved during respiration
                !-------------------------------------------------------------------------------

                if (ciso) then
                    ! Track 13C respiration proportional to 13C:12C ratio in biomass
                    HetRespFlux_13 = HetRespFlux * HetC_13 / HetC

                    if (ciso_14 .and. ciso_organic_14) then
                        ! Track 14C respiration (for radiocarbon dating applications)
                        ! Used to study carbon residence times and ocean circulation
                        HetRespFlux_14 = HetRespFlux * HetC_14 / HetC
                    end if
                end if

                !===============================================================================
                ! MESOZOOPLANKTON MORTALITY
                !===============================================================================
                ! Calculates density-dependent mortality using quadratic formulation.
                ! Higher densities lead to disproportionately higher mortality.
                !
                ! Variables:
                !   hetLossFlux : Mesozooplankton mortality flux [mmolN m-3 day-1]
                !   loss_het    : Mortality rate constant [day-1 (mmolN m-3)-1]
                !   HetN        : Mesozooplankton nitrogen concentration [mmolN m-3]
                !
                ! Equation: Quadratic Mortality
                !   Mortality = loss_het × HetN²
                !
                ! Ecological Rationale:
                !   - Linear term represents background mortality (disease, senescence)
                !   - Quadratic term represents density-dependent processes:
                !     * Increased predation pressure at high densities
                !     * Disease transmission (increases with encounter rate)
                !     * Intraspecific competition for resources
                !
                ! Note: Mortality products go to detritus pool (recycling pathway)
                !-------------------------------------------------------------------------------

                hetLossFlux = loss_het * HetN * HetN

                if (enable_3zoo2det) then

                    !===========================================================================
                    ! MACROZOOPLANKTON RESPIRATION (KRILL)
                    !===========================================================================
                    ! Calculates macrozooplankton respiration with feeding-dependent stress
                    ! response. Poor feeding conditions increase metabolic costs.
                    !
                    ! Variables:
                    !   Zoo2RespFlux         : Macrozooplankton respiration [mmolC m-3 day-1]
                    !   res_zoo2             : Baseline respiration rate [day-1]
                    !   res_zoo2_f           : Feeding success modifier [-]
                    !   res_zoo2_a           : Additional respiration factor [-]
                    !   recip_res_zoo22      : Combined respiration coefficient [day-1]
                    !   Zoo2C                : Macrozooplankton carbon concentration [mmolC m-3]
                    !   grazingFluxcarbonzoo2: Carbon ingestion rate [mmolC m-3 day-1]
                    !
                    ! Feeding Success Response:
                    !   - IF specific ingestion < 0.1 day-1: Stress response activated
                    !     * Increases metabolic costs (searching, maintenance)
                    !     * res_zoo2_f scales linearly with feeding rate
                    !   - ELSE: Normal feeding conditions (res_zoo2_f = 1.0)
                    !
                    ! Total Respiration:
                    !   RespFlux = res_zoo2 × (1 + res_zoo2_f + res_zoo2_a) × Zoo2C
                    !
                    ! Note: Calls external subroutine krill_resp() for detailed metabolism
                    !---------------------------------------------------------------------------

                    ! Call detailed krill respiration subroutine
                    ! Handles additional physiological processes (e.g., molting, reproduction)
                    call krill_resp(n, daynew, geo_coord_nod2D)

                    ! Calculate feeding success modifier
                    ! Low feeding rates trigger stress response with elevated respiration
                    if ((grazingFluxcarbonzoo2 / Zoo2C) <= 0.1) then
                        ! Stress response: metabolic costs increase linearly with feeding deficit
                        ! Factor of 0.1 converts to percentage basis
                        res_zoo2_f = 0.1 * (grazingFluxcarbonzoo2 / Zoo2C * 100.0)
                    else
                        ! Normal feeding: no additional metabolic stress
                        res_zoo2_f = 1.0
                    end if

                    ! Calculate combined respiration coefficient
                    ! Includes baseline + feeding stress + additional factors
                    recip_res_zoo22 = res_zoo2 * (1.0 + res_zoo2_f + res_zoo2_a)

                    ! Calculate total respiration flux
                    Zoo2RespFlux = recip_res_zoo22 * Zoo2C

                    !===========================================================================
                    ! MACROZOOPLANKTON MORTALITY
                    !===========================================================================
                    ! Quadratic density-dependent mortality for macrozooplankton.
                    !
                    ! Variables:
                    !   Zoo2LossFlux : Macrozooplankton mortality flux [mmolN m-3 day-1]
                    !   loss_zoo2    : Mortality rate constant [day-1 (mmolN m-3)-1]
                    !   zoo2N        : Macrozooplankton nitrogen concentration [mmolN m-3]
                    !
                    ! Equation: Quadratic Mortality
                    !   Mortality = loss_zoo2 × zoo2N²
                    !---------------------------------------------------------------------------

                    Zoo2LossFlux = loss_zoo2 * zoo2N * zoo2N

                    !===========================================================================
                    ! MICROZOOPLANKTON RESPIRATION AND MORTALITY
                    !===========================================================================
                    ! Calculates microzooplankton respiration and mortality.
                    ! Simpler than macro/mesozooplankton (no feeding-dependent stress).
                    !
                    ! Variables:
                    !   MicZooRespFlux : Microzooplankton respiration [mmolC m-3 day-1]
                    !   MicZooLossFlux : Microzooplankton mortality [mmolN m-3 day-1]
                    !   res_miczoo     : Respiration rate constant [day-1]
                    !   loss_miczoo    : Mortality rate constant [day-1 (mmolN m-3)-1]
                    !   q10_mic_res    : Q10 temperature function for micro respiration [-]
                    !   MicZooC        : Microzooplankton carbon concentration [mmolC m-3]
                    !   MicZooN        : Microzooplankton nitrogen concentration [mmolN m-3]
                    !
                    ! Respiration Equation:
                    !   RespFlux = res_miczoo × q10_mic_res × MicZooC
                    !
                    ! Mortality Equation:
                    !   Mortality = loss_miczoo × MicZooN²
                    !---------------------------------------------------------------------------

                    ! Temperature-dependent respiration (Q10 formulation)
                    ! Smaller organisms typically have higher temperature sensitivity
                    MicZooRespFlux = res_miczoo * q10_mic_res * MicZooC

                    ! Quadratic density-dependent mortality
                    MicZooLossFlux = loss_miczoo * MicZooN * MicZooN

                end if ! enable_3zoo2det

                !===============================================================================
                ! FECAL PELLET PRODUCTION
                !===============================================================================
                ! Calculates fecal pellet production from mesozooplankton and macrozooplankton.
                ! Fecal pellets are key to the biological carbon pump due to fast sinking.
                !
                ! Fecal pellets represent:
                !   - Undigested/unassimilated food material
                !   - Packaged waste with high sinking velocity (100-1000 m/day)
                !   - Major pathway for carbon export to deep ocean
                !   - Important food source for deep-sea organisms
                !
                ! Variables:
                !   Zoo2fecalloss_n     : Macrozooplankton fecal N [mmolN m-3 day-1]
                !   Zoo2fecalloss_c     : Macrozooplankton fecal C [mmolC m-3 day-1]
                !   mesfecalloss_n      : Mesozooplankton fecal N [mmolN m-3 day-1]
                !   mesfecalloss_c      : Mesozooplankton fecal C [mmolC m-3 day-1]
                !   fecal_rate_n        : Fecal production rate for N [fraction of ingestion]
                !   fecal_rate_c        : Fecal production rate for C [fraction of ingestion]
                !   fecal_rate_n_mes    : Mesozooplankton N fecal rate [-]
                !   fecal_rate_c_mes    : Mesozooplankton C fecal rate [-]
                !   grazingFlux2        : Macrozooplankton N ingestion [mmolN m-3 day-1]
                !   grazingFlux         : Mesozooplankton N ingestion [mmolN m-3 day-1]
                !   grazingFluxcarbonzoo2   : Macro C ingestion [mmolC m-3 day-1]
                !   grazingFluxcarbon_mes   : Meso C ingestion [mmolC m-3 day-1]
                !
                ! Equation:
                !   Fecal_flux = fecal_rate × ingestion_flux
                !   Fecal rate typically 10-30% of ingestion (remainder is assimilated)
                !
                ! Note: Fecal C:N ratio may differ from prey due to selective digestion
                !-------------------------------------------------------------------------------

                if (enable_3zoo2det) then

                    !---------------------------------------------------------------------------
                    ! Macrozooplankton Fecal Pellet Production
                    !---------------------------------------------------------------------------
                    ! Larger zooplankton produce larger, faster-sinking fecal pellets

                    Zoo2fecalloss_n = fecal_rate_n * grazingFlux2 ! Nitrogen
                    Zoo2fecalloss_c = fecal_rate_c * grazingFluxcarbonzoo2 ! Carbon

                    !---------------------------------------------------------------------------
                    ! Mesozooplankton Fecal Pellet Production
                    !---------------------------------------------------------------------------
                    ! Smaller fecal pellets than macrozooplankton but still fast-sinking

                    mesfecalloss_n = fecal_rate_n_mes * grazingFlux ! Nitrogen
                    mesfecalloss_c = fecal_rate_c_mes * grazingFluxcarbon_mes ! Carbon

                end if ! enable_3zoo2det

                !===============================================================================
                ! PHYTOPLANKTON AND DETRITUS AGGREGATION
                !===============================================================================
                ! Calculates particle aggregation that forms larger, faster-sinking particles.
                ! Critical process for biological carbon pump and export production.
                !
                ! Aggregation mechanisms:
                !   - Particle collision and sticking (Brownian motion, shear, settling)
                !   - Diatom mucus production enhances aggregation (nutrient stress response)
                !   - Transparent exopolymer particles (TEP) glue particles together
                !
                ! Variables:
                !   aggregationrate : Total particle aggregation rate [mmolN m-3 day-1]
                !   agg_PP          : Phytoplankton-phytoplankton aggregation rate [day-1]
                !   agg_PD          : Phytoplankton-detritus aggregation rate [day-1]
                !   PhyN, DiaN      : Phytoplankton nitrogen concentrations [mmolN m-3]
                !   CoccoN, PhaeoN  : Coccolithophore and Phaeocystis N [mmolN m-3]
                !   DetN, DetZ2N    : Slow and fast-sinking detritus N [mmolN m-3]
                !   qlimitFac       : Nutrient limitation factor (0=replete, 1=limited) [-]
                !   feLimitFac      : Iron limitation factor [-]
                !   quota_dia       : Diatom N:C quota [mmolN mmolC-1]
                !   qSiC            : Diatom Si:C quota [mmolSi mmolC-1]
                !
                ! Note: Aggregation products enter fast-sinking detritus pool (DetZ2N)
                !-------------------------------------------------------------------------------

                !---------------------------------------------------------------------------
                ! Diatom Mucus-Enhanced Aggregation (Optional)
                !---------------------------------------------------------------------------
                ! Nutrient limitation triggers diatom mucus production, enhancing aggregation
                ! and forming large, rapidly-sinking marine snow aggregates.

                if (diatom_mucus) then

                    !-----------------------------------------------------------------------
                    ! Calculate nutrient limitation factors
                    !-----------------------------------------------------------------------

                    ! Nitrogen limitation factor
                    ! Uses limiter function: returns 0 when quota is high (replete)
                    !                        returns 1 when quota is low (limited)
                    qlimitFac = recom_limiter(NMinSlope, NCmin_d, quota_dia)

                    ! Silicon limitation factor
                    ! Diatoms require silica for frustule (shell) formation
                    qlimitFacTmp = recom_limiter(SiMinSlope, SiCmin, qSiC)
                    qlimitFac = min(qlimitFac, qlimitFacTmp) ! Most limiting nutrient

                    ! Iron limitation factor
                    ! Iron is essential for photosynthesis and often limiting in ocean
                    feLimitFac = Fe / (k_Fe_d + Fe)
                    qlimitFac = min(qlimitFac, feLimitFac) ! Most limiting of all nutrients

                    ! Calculate mucus-enhanced aggregation rate
                    ! Aggregation increases with nutrient stress (high limitation factor)
                    ! Factor (1 - qlimitFac): 0 when replete, 1 when severely limited
                    aggregationrate = agg_PP * (1.0 - qlimitFac) * DiaN

                else

                    !-----------------------------------------------------------------------
                    ! Simple aggregation without nutrient limitation effect
                    !-----------------------------------------------------------------------
                    aggregationrate = agg_PP * DiaN

                end if

                !---------------------------------------------------------------------------
                ! Total Aggregation Rate
                !---------------------------------------------------------------------------
                ! Sum contributions from all particle types
                ! Each particle type can aggregate with itself and other particles

                ! Add small phytoplankton and detritus aggregation
                aggregationrate = aggregationrate + agg_PD * DetN + agg_PP * PhyN

                if (enable_3zoo2det) then
                    ! Add fast-sinking detritus aggregation
                    aggregationrate = aggregationrate + agg_PD * DetZ2N
                end if

                if (enable_coccos) then
                    ! Add coccolithophore and Phaeocystis aggregation
                    ! These can form large blooms with high aggregation potential
                    aggregationrate = aggregationrate + agg_PP * CoccoN + agg_PP * PhaeoN
                end if

                !===============================================================================
                ! MARINE CALCIFICATION
                !===============================================================================
                ! Simulates the formation and dissolution of calcium carbonate (CaCO3) in
                ! marine environments, primarily through phytoplankton calcification.
                !
                ! This module calculates:
                !   1. Calcification rates (organism-specific or general)
                !   2. Environmental controls on calcification (temperature, CO2, nutrients)
                !   3. Loss processes (grazing, aggregation)
                !   4. Carbon isotope fractionation (13C, 14C - optional)
                !
                ! Key Features:
                !   - Explicit coccolithophore calcification (temperature, CO2, nutrient effects)
                !   - General small phytoplankton calcification (simplified)
                !   - Ocean acidification sensitivity (CO2/carbonate chemistry effects)
                !   - PIC:POC ratio variability (environmental controls on calcite production)
                !   - Isotope fractionation during calcification
                !
                ! Calcification Scenarios:
                !   A) enable_coccos = TRUE:  Explicit coccolithophore model
                !      - Temperature-dependent PIC:POC ratios (Krumhardt et al. 2017, 2019)
                !      - CO2/carbonate chemistry effects (ocean acidification)
                !      - Nitrogen limitation effects
                !   B) enable_coccos = FALSE: General small phytoplankton
                !      - Fixed fraction of photosynthesis goes to calcification
                !
                ! Ecological/Biogeochemical Significance:
                !   - Produces ballast that increases particle sinking (biological pump)
                !   - Releases CO2 during calcification (carbonate counter-pump)
                !   - Major component of ocean alkalinity cycle
                !   - Sensitive to ocean acidification (climate change impact)
                !
                ! References:
                !   - Krumhardt et al. (2017, 2019) - Coccolithophore parameterization
                !   - Gürses et al. (2023) - REcoM2 model description
                !===============================================================================

                !===============================================================================
                ! CALCIFICATION RATE CALCULATION
                !===============================================================================
                ! Calculates the rate at which calcium carbonate (CaCO3) is formed by marine
                ! organisms, with environmental modulation of the PIC:POC ratio.
                !
                ! PIC = Particulate Inorganic Carbon (CaCO3)
                ! POC = Particulate Organic Carbon (photosynthetic biomass)
                !
                ! Variables (Coccolithophore mode):
                !   calcification   : CaCO3 formation rate [mmolC m-3 day-1]
                !   Cphot_cocco     : Coccolithophore C-specific photosynthesis [day-1]
                !   CoccoC          : Coccolithophore carbon concentration [mmolC m-3]
                !   PICPOCtemp      : Temperature modifier for PIC:POC ratio [-]
                !   PICPOCCO2       : CO2/carbonate chemistry modifier for PIC:POC [-]
                !   PICPOCN         : Nitrogen limitation modifier for PIC:POC [-]
                !   Temp(k)         : Temperature at depth k [degC]
                !   HCO3_watercolumn: Bicarbonate concentration [mmolC m-3]
                !   CO2_watercolumn : Dissolved CO2 concentration [mmolC m-3]
                !   pH_watercolumn  : Water column pH [-]
                !   DIN             : Dissolved inorganic nitrogen [mmolN m-3]
                !   k_din_c         : Half-saturation for DIN effect on calcification [mmolN m-3]
                !   a,b,c,d_co2_calc: CO2 effect parameters [-]
                !   Cunits          : Conversion factor for units [-]
                !
                ! Variables (General phytoplankton mode):
                !   calc_prod_ratio : Fixed calcification:photosynthesis ratio [-]
                !   Cphot           : Small phytoplankton photosynthesis rate [day-1]
                !   PhyC            : Small phytoplankton carbon [mmolC m-3]
                !-------------------------------------------------------------------------------

                if (enable_coccos) then

                    !===========================================================================
                    ! COCCOLITHOPHORE-SPECIFIC CALCIFICATION
                    !===========================================================================
                    ! Coccolithophores are marine phytoplankton that produce calcite plates
                    ! (coccoliths) as an external protective covering.
                    !
                    ! Environmental Controls:
                    !   - Temperature: Optimum calcification at warm temperatures (>10.6degC)
                    !   - CO2: Ocean acidification reduces calcification efficiency
                    !   - Nitrogen: High nutrient availability reduces PIC:POC ratio
                    !---------------------------------------------------------------------------

                    !---------------------------------------------------------------------------
                    ! Temperature Effect on PIC:POC Ratio
                    !---------------------------------------------------------------------------
                    ! Based on Krumhardt et al. (2017, 2019) parameterization
                    ! Warmer temperatures favor higher calcification efficiency

                    if (Temp(k) < 10.6) then
                        ! Linear relationship for temperatures below 10.6degC
                        ! Calcification efficiency increases with temperature
                        PICPOCtemp = 0.104d0 * Temp(k) - 0.108d0
                    else
                        ! Maximum efficiency at temperatures ≥ 10.6degC
                        ! Represents optimal conditions for coccolith formation
                        PICPOCtemp = 1.0d0
                    end if

                    ! Prevent negative values (can occur at very low temperatures)
                    PICPOCtemp = max(tiny, PICPOCtemp)

                    !---------------------------------------------------------------------------
                    ! CO2/Carbonate Chemistry Effect on PIC:POC Ratio
                    !---------------------------------------------------------------------------
                    ! Complex function representing ocean acidification impacts
                    ! Three components:
                    !   1. HCO3- availability (substrate for calcification)
                    !   2. High CO2 inhibition (ocean acidification stress)
                    !   3. Low pH inhibition (direct pH stress on calcification)

                    PICPOCCO2 = a_co2_calc * HCO3_watercolumn(k) * Cunits &
                            / (b_co2_calc + HCO3_watercolumn(k) * Cunits) &
                            - exp(-c_co2_calc * CO2_watercolumn(k) * Cunits) &
                            - d_co2_calc * 10.d0 ** (-pH_watercolumn(k))

                    ! Apply empirical constraints to CO2 effect
                    PICPOCCO2 = min(PICPOCCO2, 3.d0) ! Upper limit (April 2022 modification)
                    PICPOCCO2 = max(0.d0, PICPOCCO2) ! Lower limit (July 2022 modification)

                    !---------------------------------------------------------------------------
                    ! Nitrogen Limitation Effect on PIC:POC Ratio
                    !---------------------------------------------------------------------------
                    ! Higher DIN availability reduces calcification efficiency
                    ! Ecological rationale: Under nutrient replete conditions, organisms
                    ! allocate resources preferentially to organic carbon (growth) rather
                    ! than calcite plates (protection)

                    PICPOCN = -0.31 * (DIN / (DIN + k_din_c)) + 1.31
                    PICPOCN = max(tiny, PICPOCN) ! Prevent negative values

                    !---------------------------------------------------------------------------
                    ! Calculate Final Coccolithophore Calcification Rate
                    !---------------------------------------------------------------------------
                    ! Combine photosynthesis rate with all environmental modifiers
                    ! Base rate × Temperature effect × Nutrient effect

                    calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN

                    ! Apply CO2 limitation if ocean acidification sensitivity is enabled
                    if (CO2lim) then
                        calcification = calcification * PICPOCCO2
                    end if

                else

                    !===========================================================================
                    ! GENERAL SMALL PHYTOPLANKTON CALCIFICATION
                    !===========================================================================
                    ! Simplified approach: calcification as a fixed fraction of photosynthesis
                    ! Used when coccolithophores are not explicitly represented
                    ! Represents calcification by other small calcifying organisms
                    !---------------------------------------------------------------------------

                    calcification = calc_prod_ratio * Cphot * PhyC

                end if

                !===============================================================================
                ! CALCIFICATION LOSS PROCESSES
                !===============================================================================
                ! Calculates removal of calcified material through grazing and aggregation.
                ! These processes transfer CaCO3 to zooplankton or sinking particles.
                !
                ! Variables:
                !   calc_loss_agg   : CaCO3 loss through aggregation [mmolC m-3 day-1]
                !   calc_loss_gra   : CaCO3 loss to primary grazer [mmolC m-3 day-1]
                !   calc_loss_gra2  : CaCO3 loss to secondary grazer [mmolC m-3 day-1]
                !   calc_loss_gra3  : CaCO3 loss to tertiary grazer [mmolC m-3 day-1]
                !   PhyCalc         : Phytoplankton calcite content [mmolC m-3]
                !   aggregationrate : Particle aggregation rate [mmolN m-3 day-1]
                !   recipQuota      : Small phytoplankton C:N ratio [mmolC mmolN-1]
                !   recipQuota_Cocco: Coccolithophore C:N ratio [mmolC mmolN-1]
                !   grazingFlux_*   : Grazing rates on different prey [mmolN m-3 day-1]
                !   aux             : Auxiliary conversion factor [mmolC mmolN-1]
                !
                ! Note: CaCO3 in grazed material can dissolve in zooplankton guts or
                !       provide ballast for fecal pellets (enhancing sinking)
                !-------------------------------------------------------------------------------

                !-------------------------------------------------------------------------------
                ! Aggregation Loss
                !-------------------------------------------------------------------------------
                ! CaCO3 is incorporated into aggregates (marine snow) that sink rapidly
                ! Provides ballast effect: increases particle density and sinking velocity

                calc_loss_agg = aggregationrate * PhyCalc

                if (enable_coccos) then

                    !---------------------------------------------------------------------------
                    ! Coccolithophore Grazing Losses
                    !---------------------------------------------------------------------------
                    ! Calculate CaCO3 consumed by grazers feeding on coccolithophores
                    ! Requires conversion from nitrogen-based grazing to carbon flux

                    ! Calculate auxiliary conversion factor (CaCO3:N ratio in coccolithophores)
                    aux = recipQuota_Cocco / (CoccoC + tiny) * PhyCalc

                    ! Primary grazer (mesozooplankton) consumption
                    calc_loss_gra = grazingFlux_Cocco * aux

                    if (enable_3zoo2det) then
                        ! Secondary grazer (macrozooplankton) consumption
                        calc_loss_gra2 = grazingFlux_Cocco2 * aux

                        ! Tertiary grazer (microzooplankton) consumption
                        calc_loss_gra3 = grazingFlux_Cocco3 * aux
                    end if

                else

                    !---------------------------------------------------------------------------
                    ! Small Phytoplankton Grazing Losses
                    !---------------------------------------------------------------------------
                    ! Calculate CaCO3 consumed by grazers feeding on calcifying small phytoplankton

                    ! Calculate auxiliary conversion factor (CaCO3:N ratio in small phytoplankton)
                    aux = recipQuota / (PhyC + tiny) * PhyCalc

                    ! Primary grazer (mesozooplankton) consumption
                    calc_loss_gra = grazingFlux_phy * aux

                    if (enable_3zoo2det) then
                        ! Secondary grazer (macrozooplankton) consumption
                        calc_loss_gra2 = grazingFlux_phy2 * aux

                        ! Tertiary grazer (microzooplankton) consumption
                        calc_loss_gra3 = grazingFlux_phy3 * aux
                    end if

                end if

                !===============================================================================
                ! CARBON ISOTOPE FRACTIONATION
                !===============================================================================
                ! Handles carbon-13 and carbon-14 isotope fractionation during calcification
                ! and loss processes. Used for paleoclimate reconstructions and carbon cycle
                ! studies.
                !
                ! Isotope Fractionation:
                !   - Light isotopes (12C) are preferentially incorporated during calcification
                !   - Fractionation factors (alpha) quantify this discrimination
                !   - Different fractionation for 13C and 14C
                !
                ! Variables:
                !   calcification_13    : 13C calcification rate [mmol13C m-3 day-1]
                !   calcification_14    : 14C calcification rate [mmol14C m-3 day-1]
                !   alpha_calc_13       : 13C fractionation factor during calcification [-]
                !   alpha_calc_14       : 14C fractionation factor during calcification [-]
                !   calc_loss_agg_13/14 : Isotope losses through aggregation [mmol m-3 day-1]
                !   calc_loss_gra_13/14 : Isotope losses through grazing [mmol m-3 day-1]
                !   PhyCalc_13/14       : Phytoplankton calcite isotope content [mmol m-3]
                !   PhyC_13/14          : Phytoplankton organic carbon isotopes [mmol m-3]
                !   recipQuota_13/14    : Isotope-specific C:N ratios [mmol mmol-1]
                !
                ! Applications:
                !   - 13C: Paleoclimate proxies, carbon source tracing
                !   - 14C: Radiocarbon dating, carbon residence times
                !
                ! Note: Only executed if carbon isotope tracking is enabled (ciso = TRUE)
                !-------------------------------------------------------------------------------

                if (ciso) then

                    !===========================================================================
                    ! Carbon-13 Isotope Calculations
                    !===========================================================================
                    ! Track 13C through calcification and loss processes
                    !---------------------------------------------------------------------------

                    ! Calcification with isotopic fractionation
                    ! Alpha factor < 1 means light isotope (12C) is preferentially incorporated
                    calcification_13 = calcification * alpha_calc_13

                    ! Isotopic losses through aggregation
                    ! Assumes no fractionation during physical aggregation process
                    calc_loss_agg_13 = aggregationRate * PhyCalc_13

                    ! Isotopic losses through grazing
                    ! Requires isotope-specific conversion factor
                    calc_loss_gra_13 = grazingFlux_phy * recipQuota_13 / (PhyC_13 + tiny) &
                            * PhyCalc_13

                    if (ciso_14 .and. ciso_organic_14) then

                        !=======================================================================
                        ! Carbon-14 Isotope Calculations
                        !=======================================================================
                        ! Track radiocarbon (14C) through calcification and loss processes
                        ! Used for determining carbon residence times and age dating
                        !-----------------------------------------------------------------------

                        ! Calcification with 14C fractionation
                        ! 14C fractionation is approximately twice that of 13C
                        calcification_14 = calc_prod_ratio * Cphot * PhyC_14 * alpha_calc_14

                        ! 14C losses through aggregation
                        calc_loss_agg_14 = aggregationRate * PhyCalc_14

                        ! 14C losses through grazing
                        calc_loss_gra_14 = grazingFlux_phy * recipQuota_14 / (PhyC_14 + tiny) &
                                * PhyCalc_14

                    end if ! ciso_14 .and. ciso_organic_14

                end if ! ciso

                call update_sms(k, dt_b, dt, MPI_COMM_FESOM, sms, state, PhyN, PhyC, PhyChl, &
                        PhyCalc, DiaN, DiaC, DiaChl, DiaSi, CoccoN, CoccoC, CoccoChl, PhaeoN, &
                        PhaeoC, PhaeoChl, HetN, HetC, Zoo2N, Zoo2C, MicZooN, MicZooC, DetN, DetC, &
                        DetSi, DetCalc, DetZ2N, DetZ2C, DetZ2Si, DetZ2Calc, DON, EOC, FreeFe)

                !-------------------------------------------------------------------------------
                ! DIAGNOSTIC ACCUMULATION INITIALIZATION
                !-------------------------------------------------------------------------------
                ! Calculate averaging weight for accumulating diagnostics over biogeochemical
                ! sub-time steps.
                !
                ! Variables:
                !   recipbiostep : Reciprocal of number of bio sub-steps (averaging weight) [-]
                !   biostep      : Number of biogeochemistry steps per physics step [-]
                !
                ! Averaging Approach:
                !   Each sub-step contributes: (1/biostep) × rate
                !   After biostep iterations: average rate over physics time step
                !
                ! Example:
                !   If biostep = 4 (4 bio steps per physics step)
                !   recipbiostep = 0.25
                !   Each contribution weighted by 0.25
                !   Sum of 4 contributions = time-averaged rate
                !-------------------------------------------------------------------------------

                recipbiostep = 1.d0 / real(biostep)

                if (Diags) then

                    !===========================================================================
                    ! PRIMARY PRODUCTION DIAGNOSTICS
                    !===========================================================================
                    ! Net Primary Production (NPP) = Gross Production - Autotrophic Respiration
                    ! Gross Primary Production (GPP) = Photosynthetic carbon fixation only
                    !
                    ! NPP represents carbon available for:
                    !   - Growth (biomass increase)
                    !   - Exudation (DOC production)
                    !   - Grazing by herbivores
                    !
                    ! Units: mmolC m-3 day-1
                    !---------------------------------------------------------------------------

                    !===========================================================================
                    ! NET PRIMARY PRODUCTION (NPP)
                    !===========================================================================
                    ! NPP = Photosynthesis - Respiration
                    ! Carbon available after meeting metabolic maintenance costs
                    !
                    ! Variables (Small phytoplankton):
                    !   vertNPPn(k)   : Time-averaged NPP for small phyto [mmolC m-3 day-1]
                    !   Cphot         : Carbon-specific photosynthesis rate [day-1]
                    !   PhyC          : Small phytoplankton carbon [mmolC m-3]
                    !   PhyRespRate   : Respiration rate [day-1]
                    !
                    ! Ecological Interpretation:
                    !   High NPP: Nutrient replete, optimal growth conditions
                    !   Low NPP: Nutrient/light limited, high respiration costs
                    !   Negative NPP: Respiration exceeds photosynthesis (rare, stress)
                    !
                    ! Validation Targets:
                    !   - 14C uptake measurements
                    !   - Satellite-based NPP estimates (VGPM, CbPM, CAFE)
                    !   - Time series stations (HOT, BATS)
                    !   - Typical range: 0.1-10 mmolC m-3 day-1 (surface)
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton NPP
                    vertNPPn(k) = vertNPPn(k) + ( &
                            +Cphot * PhyC & ! Photosynthetic production
                            - PhyRespRate * PhyC & ! Minus respiration costs
                            ) * recipbiostep

                    ! Diatom NPP
                    vertNPPd(k) = vertNPPd(k) + ( &
                            +Cphot_dia * DiaC &
                            - PhyRespRate_dia * DiaC &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore NPP
                        vertNPPc(k) = vertNPPc(k) + ( &
                                +Cphot_cocco * CoccoC &
                                - PhyRespRate_cocco * CoccoC &
                                ) * recipbiostep

                        ! Phaeocystis NPP
                        vertNPPp(k) = vertNPPp(k) + ( &
                                +Cphot_phaeo * PhaeoC &
                                - PhyRespRate_phaeo * PhaeoC &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! GROSS PRIMARY PRODUCTION (GPP)
                    !===========================================================================
                    ! GPP = Total photosynthetic carbon fixation (before respiration)
                    ! Represents maximum potential carbon uptake
                    !
                    ! Relationship:
                    !   GPP = NPP + Respiration
                    !   Growth efficiency = NPP / GPP (typically 0.6-0.8)
                    !
                    ! Uses:
                    !   - Calculate carbon use efficiency
                    !   - Compare to oxygen evolution (photosynthetic quotient)
                    !   - Understand temperature effects (GPP and R differ)
                    !
                    ! Typical Values:
                    !   Surface: 1-20 mmolC m-3 day-1
                    !   Deep chlorophyll max: 0.1-5 mmolC m-3 day-1
                    !   Below euphotic zone: <0.01 mmolC m-3 day-1
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton GPP
                    vertGPPn(k) = vertGPPn(k) + ( &
                            +Cphot * PhyC & ! Total photosynthetic fixation
                            ) * recipbiostep

                    ! Diatom GPP
                    vertGPPd(k) = vertGPPd(k) + ( &
                            +Cphot_dia * DiaC &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore GPP
                        vertGPPc(k) = vertGPPc(k) + ( &
                                +Cphot_cocco * CoccoC &
                                ) * recipbiostep

                        ! Phaeocystis GPP
                        vertGPPp(k) = vertGPPp(k) + ( &
                                +Cphot_phaeo * PhaeoC &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! NET NITROGEN ASSIMILATION
                    !===========================================================================
                    ! Net nitrogen uptake = Assimilation - Exudation losses
                    ! Represents nitrogen incorporation into biomass
                    !
                    ! Variables (Small phytoplankton):
                    !   vertNNAn(k)   : Net N assimilation for small phyto [mmolN m-3 day-1]
                    !   N_assim       : N-specific assimilation rate [mmolN mmolC-1 day-1]
                    !   PhyC          : Small phytoplankton carbon [mmolC m-3]
                    !   lossN         : N exudation rate constant [day-1]
                    !   limitFacN     : Nutrient limitation factor [0-1]
                    !   PhyN          : Small phytoplankton nitrogen [mmolN m-3]
                    !
                    ! Ecological Interpretation:
                    !   Positive: Net nitrogen accumulation (growth)
                    !   Negative: Net nitrogen loss (stress-induced exudation)
                    !   Zero: Balanced uptake and loss (steady state)
                    !
                    ! Relationship to C:N Ratio:
                    !   Net N assimilation / NPP = change in N:C quota
                    !   High ratio -> decreasing C:N (N accumulation)
                    !   Low ratio -> increasing C:N (N limitation)
                    !
                    ! Validation:
                    !   - 15N tracer studies
                    !   - Nutrient depletion experiments
                    !   - Typical range: 0.01-2 mmolN m-3 day-1
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton net N assimilation
                    vertNNAn(k) = vertNNAn(k) + ( &
                            +N_assim * PhyC & ! Nitrogen uptake from DIN
                            - lossN * limitFacN * PhyN & ! Minus exudation (stress-induced)
                            ) * recipbiostep

                    ! Diatom net N assimilation
                    vertNNAd(k) = vertNNAd(k) + ( &
                            +N_assim_dia * DiaC &
                            - lossN * limitFacN_dia * DiaN &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore net N assimilation
                        vertNNAc(k) = vertNNAc(k) + ( &
                                +N_assim_cocco * CoccoC &
                                - lossN * limitFacN_cocco * CoccoN &
                                ) * recipbiostep

                        ! Phaeocystis net N assimilation
                        vertNNAp(k) = vertNNAp(k) + ( &
                                +N_assim_phaeo * PhaeoC &
                                - lossN * limitFacN_phaeo * PhaeoN &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! CHLOROPHYLL DEGRADATION RATES
                    !===========================================================================
                    ! Tracks chlorophyll turnover rates (photodamage + senescence)
                    ! Important for understanding chlorophyll:carbon dynamics
                    !
                    ! Variables:
                    !   vertChldegn(k) : Small phyto Chl degradation rate [day-1]
                    !   KOchl          : Chlorophyll degradation rate coefficient [day-1]
                    !
                    ! Rate Components:
                    !   - Base degradation: Constant senescence (~0.01-0.05 day-1)
                    !   - Photodamage: Light-dependent additional degradation
                    !   - Total: KOchl = base + photodamage component
                    !
                    ! Ecological Significance:
                    !   - High rates: Photodamage stress, high light
                    !   - Low rates: Low light, minimal photoinhibition
                    !   - Affects Chl:C ratio and satellite retrievals
                    !
                    ! Applications:
                    !   - Understand Chl:C variability
                    !   - Validate photoacclimation dynamics
                    !   - Interpret satellite chlorophyll trends
                    !
                    ! Note: Changed from gross N-assimilation diagnostic (previous version)
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton Chl degradation
                    vertChldegn(k) = vertChldegn(k) + ( &
                            +KOchl & ! Degradation rate [day-1]
                            ) * recipbiostep

                    ! Diatom Chl degradation
                    vertChldegd(k) = vertChldegd(k) + ( &
                            +KOchl_dia &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore Chl degradation
                        vertChldegc(k) = vertChldegc(k) + ( &
                                +KOchl_cocco &
                                ) * recipbiostep

                        ! Phaeocystis Chl degradation
                        vertChldegp(k) = vertChldegp(k) + ( &
                                +KOchl_phaeo &
                                ) * recipbiostep
                    end if

                    !--------------------------------------------------------------------------------------------------------------------------------------

                    ! GRAZING FLUXES
                    ! Only for the case with detritus grazing, not without detritus grazing,
                    ! because
                    ! this output is probably anyway not needed as a default.
                    ! diagnostics, combined from Onur and Cara, modified by Miriam

                    if (Grazing_detritus) then

                        !===========================================================================
                        ! MESOZOOPLANKTON GRAZING
                        !===========================================================================
                        ! Tracks grazing by mesozooplankton (100-2000 μm) on multiple prey types
                        ! Central trophic link connecting primary producers to higher consumers
                        !
                        ! Variables:
                        !   vertgrazmeso_tot(k)   : Total assimilated grazing [mmolC m-3 day-1]
                        !   vertgrazmeso_n(k)     : Grazing on nanophytoplankton [mmolC m-3 day-1]
                        !   vertgrazmeso_d(k)     : Grazing on diatoms [mmolC m-3 day-1]
                        !   vertgrazmeso_c(k)     : Grazing on coccolithophores [mmolC m-3 day-1]
                        !   vertgrazmeso_p(k)     : Grazing on Phaeocystis [mmolC m-3 day-1]
                        !   vertgrazmeso_det(k)   : Grazing on detritus 1 [mmolC m-3 day-1]
                        !   vertgrazmeso_det2(k)  : Grazing on detritus 2 [mmolC m-3 day-1]
                        !   vertgrazmeso_mic(k)   : Grazing on microzooplankton [mmolC m-3 day-1]
                        !
                        ! Prey Composition:
                        !   PHY    : Nanophytoplankton (2-20 μm) - primary prey
                        !   DIA    : Diatoms (>20 μm) - silicified cells
                        !   COCCO  : Coccolithophores - calcifying cells (optional)
                        !   PHAEO  : Phaeocystis colonies (optional)
                        !   DET    : Detritus type 1 - particulate organic matter
                        !   DET2   : Detritus type 2 - larger particles (optional)
                        !   MICZOO : Microzooplankton - intraguild predation (optional)
                        !
                        ! Ecological Significance:
                        !   - Grazing efficiency (grazEff): 30-70%, rest to fecal pellets
                        !   - Prey selection via Holling Type II/III functional response
                        !   - Major pathway for carbon export via fecal pellets
                        !   - Controls phytoplankton bloom magnitude/duration
                        !
                        ! Validation:
                        !   - Gut content analysis & pigment measurements
                        !   - Dilution experiments (Landry-Hassett method)
                        !   - Typical rates: 0.01-10 mmolC m-3 day-1
                        !   - Clearance rates: 1-100 mL ind-1 day-1
                        !---------------------------------------------------------------------------

                        ! Total assimilated grazing (with efficiency applied)
                        vertgrazmeso_tot(k) = vertgrazmeso_tot(k) + ( &
                                +grazingFlux_phy * recipQuota * grazEff & ! Nanophytoplankton
                                + grazingFlux_Dia * recipQuota_Dia * grazEff & ! Diatoms
                                + grazingFlux_Det * recipDet * grazEff & ! Detritus 1
                                ) * recipbiostep

                        if (enable_coccos) then
                            vertgrazmeso_tot(k) = vertgrazmeso_tot(k) + ( &
                            ! Coccolithophores
                                    +grazingFlux_Cocco * recipQuota_Cocco * grazEff &
                            ! Phaeocystis
                                    + grazingFlux_Phaeo * recipQuota_Phaeo * grazEff &
                                    ) * recipbiostep
                        end if

                        if (enable_3zoo2det) then
                            vertgrazmeso_tot(k) = vertgrazmeso_tot(k) + ( &
                                    +GrazingFlux_DetZ2 * recipDet2 * grazEff & ! Detritus 2
                                    + grazingFlux_miczoo * recipQZoo3 * grazEff & ! Microzooplankton
                                    ) * recipbiostep
                        end if

                        ! Prey-specific mortality (loss terms, no efficiency applied)
                        ! These track carbon removal from each prey population

                        ! Small phytoplankton mortality
                        vertgrazmeso_n(k) = vertgrazmeso_n(k) + ( &
                                +grazingFlux_phy * recipQuota &
                                ) * recipbiostep

                        ! Diatom mortality
                        vertgrazmeso_d(k) = vertgrazmeso_d(k) + ( &
                                +grazingFlux_dia * recipQuota_dia &
                                ) * recipbiostep

                        if (enable_coccos) then
                            ! Coccolithophore mortality
                            vertgrazmeso_c(k) = vertgrazmeso_c(k) + ( &
                                    +grazingFlux_Cocco * recipQuota_cocco &
                                    ) * recipbiostep

                            ! Phaeocystis mortality
                            vertgrazmeso_p(k) = vertgrazmeso_p(k) + ( &
                                    +grazingFlux_Phaeo * recipQuota_Phaeo &
                                    ) * recipbiostep
                        end if

                        ! Detritus 1 consumption
                        vertgrazmeso_det(k) = vertgrazmeso_det(k) + ( &
                                +grazingFlux_Det * recipDet &
                                ) * recipbiostep

                        if (enable_3zoo2det) then
                            ! Microzooplankton mortality (intraguild predation)
                            vertgrazmeso_mic(k) = vertgrazmeso_mic(k) + ( &
                                    +grazingFlux_miczoo * recipQZoo3 &
                                    ) * recipbiostep

                            ! Detritus 2 consumption
                            vertgrazmeso_det2(k) = vertgrazmeso_det2(k) + ( &
                                    +GrazingFlux_DetZ2 * recipDet2 &
                                    ) * recipbiostep
                        end if

                        !===========================================================================
                        ! MACROZOOPLANKTON GRAZING (KRILL)
                        !===========================================================================
                        ! Tracks grazing by macrozooplankton (2-20 mm), often dominated by krill
                        ! Top mesozooplankton predator with omnivorous feeding strategy
                        !
                        ! Variables:
                        !   vertgrazmacro_tot(k)  : Total assimilated grazing [mmolC m-3 day-1]
                        !   vertgrazmacro_n(k)    : Grazing on small phytoplankton [mmolC m-3
                        ! day-1]
                        !   vertgrazmacro_d(k)    : Grazing on diatoms [mmolC m-3 day-1]
                        !   vertgrazmacro_c(k)    : Grazing on coccolithophores [mmolC m-3 day-1]
                        !   vertgrazmacro_p(k)    : Grazing on Phaeocystis [mmolC m-3 day-1]
                        !   vertgrazmacro_mes(k)  : Grazing on mesozooplankton [mmolC m-3 day-1]
                        !   vertgrazmacro_mic(k)  : Grazing on microzooplankton [mmolC m-3 day-1]
                        !   vertgrazmacro_det(k)  : Grazing on detritus 1 [mmolC m-3 day-1]
                        !   vertgrazmacro_det2(k) : Grazing on detritus 2 [mmolC m-3 day-1]
                        !
                        ! Prey Composition:
                        !   PHY    : Small phytoplankton - supplementary prey
                        !   DIA    : Diatoms - preferred prey in productive waters
                        !   COCCO  : Coccolithophores (optional)
                        !   PHAEO  : Phaeocystis colonies (optional)
                        !   HET    : Mesozooplankton - carnivory/cannibalism
                        !   MICZOO : Microzooplankton - smaller zooplankton
                        !   DET    : Detritus 1 - opportunistic feeding
                        !   DET2   : Detritus 2 - larger particles
                        !
                        ! Ecological Significance:
                        !   - Grazing efficiency (grazEff2): typically 30-60%
                        !   - Major prey for fish, seabirds, marine mammals
                        !   - Produces large, fast-sinking fecal pellets
                        !   - Vertical migration enhances carbon export
                        !   - Key species: Antarctic krill (Euphausia superba)
                        !
                        ! Validation:
                        !   - Net tows & acoustic surveys for biomass
                        !   - Feeding experiments with size-fractionated prey
                        !   - Typical rates: 0.1-50 mmolC m-3 day-1
                        !   - Individual ingestion: 10-100% body C day-1
                        !---------------------------------------------------------------------------

                        if (enable_3zoo2det) then

                            ! Total assimilated grazing (with efficiency applied)
                            vertgrazmacro_tot(k) = vertgrazmacro_tot(k) + ( &
                            ! Small phytoplankton
                                    +grazingFlux_phy2 * recipQuota * grazEff2 &
                                    + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 & ! Diatoms
                                    + grazingFlux_het2 * recipQZoo * grazEff2 & ! Mesozooplankton
                            ! Microzooplankton
                                    + grazingFlux_miczoo2 * recipQZoo3 * grazEff2 &
                                    + grazingFlux_Det2 * recipDet * grazEff2 & ! Detritus 1
                                    + grazingFlux_DetZ22 * recipDet2 * grazEff2 & ! Detritus 2
                                    ) * recipbiostep

                            if (enable_coccos) then
                                vertgrazmacro_tot(k) = vertgrazmacro_tot(k) + ( &
                                ! Coccolithophores
                                        +grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 &
                                ! Phaeocystis
                                        + grazingFlux_Phaeo2 * recipQuota_Phaeo * grazEff2 &
                                        ) * recipbiostep
                            end if

                            ! Prey-specific mortality (loss terms, no efficiency applied)
                            ! These track carbon removal from each prey population

                            ! Small phytoplankton mortality
                            vertgrazmacro_n(k) = vertgrazmacro_n(k) + ( &
                                    +grazingFlux_phy2 * recipQuota &
                                    ) * recipbiostep

                            ! Diatom mortality
                            vertgrazmacro_d(k) = vertgrazmacro_d(k) + ( &
                                    +grazingFlux_Dia2 * recipQuota_Dia &
                                    ) * recipbiostep

                            if (enable_coccos) then
                                ! Coccolithophore mortality
                                vertgrazmacro_c(k) = vertgrazmacro_c(k) + ( &
                                        +grazingFlux_Cocco2 * recipQuota_cocco &
                                        ) * recipbiostep

                                ! Phaeocystis mortality
                                vertgrazmacro_p(k) = vertgrazmacro_p(k) + ( &
                                        +grazingFlux_Phaeo2 * recipQuota_Phaeo &
                                        ) * recipbiostep
                            end if

                            ! Mesozooplankton mortality (carnivory)
                            vertgrazmacro_mes(k) = vertgrazmacro_mes(k) + ( &
                                    +grazingFlux_het2 * recipQZoo &
                                    ) * recipbiostep

                            ! Detritus 1 consumption
                            vertgrazmacro_det(k) = vertgrazmacro_det(k) + ( &
                                    +grazingFlux_Det2 * recipDet &
                                    ) * recipbiostep

                            ! Microzooplankton mortality
                            vertgrazmacro_mic(k) = vertgrazmacro_mic(k) + ( &
                                    +grazingFlux_miczoo2 * recipQZoo3 &
                                    ) * recipbiostep

                            ! Detritus 2 consumption
                            vertgrazmacro_det2(k) = vertgrazmacro_det2(k) + ( &
                                    +GrazingFlux_DetZ22 * recipDet2 &
                                    ) * recipbiostep

                        end if

                        !===========================================================================
                        ! MICROZOOPLANKTON GRAZING
                        !===========================================================================
                        ! Tracks grazing by microzooplankton (20-200 μm), mainly
                        ! ciliates/dinoflagellates
                        ! Critical link between picoplankton and mesozooplankton
                        !
                        ! Variables:
                        !   vertgrazmicro_tot(k) : Total assimilated grazing [mmolC m-3 day-1]
                        !   vertgrazmicro_n(k)   : Grazing on nanophytoplankton [mmolC m-3 day-1]
                        !   vertgrazmicro_d(k)   : Grazing on diatoms [mmolC m-3 day-1]
                        !   vertgrazmicro_c(k)   : Grazing on coccolithophores [mmolC m-3 day-1]
                        !   vertgrazmicro_p(k)   : Grazing on Phaeocystis [mmolC m-3 day-1]
                        !
                        ! Prey Composition:
                        !   PHY   : Snall phytoplankton (2-20 μm) - primary prey
                        !   DIA   : Small diatoms - when available
                        !   COCCO : Coccolithophores - calcifying prey (optional)
                        !   PHAEO : Phaeocystis colonies/single cells (optional)
                        !
                        ! Ecological Significance:
                        !   - Grazing efficiency (grazEff3): typically 40-70%
                        !   - Consumes 60-100% of primary production in oligotrophic waters
                        !   - High growth rates (1-2 doublings day-1 at 20°C)
                        !   - Size-selective feeding (optimal prey 10-50% predator size)
                        !   - Regenerates nutrients in upper water column
                        !
                        ! Validation:
                        !   - Dilution experiments (most common method)
                        !   - Epifluorescence microscopy for biomass
                        !   - Typical rates: 0.01-5 mmolC m-3 day-1
                        !   - Clearance rates: 10-10,000 nL ind-1 hour-1
                        !---------------------------------------------------------------------------

                        if (enable_3zoo2det) then

                            ! Total assimilated grazing (with efficiency applied)
                            vertgrazmicro_tot(k) = vertgrazmicro_tot(k) + ( &
                            ! Small phytoplankton
                                    +grazingFlux_phy3 * recipQuota * grazEff3 &
                                    + grazingFlux_Dia3 * recipQuota_Dia * grazEff3 & ! Diatoms
                                    ) * recipbiostep

                            if (enable_coccos) then
                                vertgrazmicro_tot(k) = vertgrazmicro_tot(k) + ( &
                                ! Coccolithophores
                                        +grazingFlux_Cocco3 * recipQuota_Cocco * grazEff3 &
                                ! Phaeocystis
                                        + grazingFlux_Phaeo3 * recipQuota_Phaeo * grazEff3 &
                                        ) * recipbiostep
                            end if

                            ! Prey-specific mortality (loss terms, no efficiency applied)
                            ! These track carbon removal from each prey population

                            ! Small phytoplankton mortality
                            vertgrazmicro_n(k) = vertgrazmicro_n(k) + ( &
                                    +grazingFlux_phy3 * recipQuota &
                                    ) * recipbiostep

                            ! Diatom mortality
                            vertgrazmicro_d(k) = vertgrazmicro_d(k) + ( &
                                    +grazingFlux_Dia3 * recipQuota_Dia &
                                    ) * recipbiostep

                            if (enable_coccos) then
                                ! Coccolithophore mortality
                                vertgrazmicro_c(k) = vertgrazmicro_c(k) + ( &
                                        +grazingFlux_Cocco3 * recipQuota_cocco &
                                        ) * recipbiostep

                                ! Phaeocystis mortality
                                vertgrazmicro_p(k) = vertgrazmicro_p(k) + ( &
                                        +grazingFlux_Phaeo3 * recipQuota_Phaeo &
                                        ) * recipbiostep
                            end if

                        end if

                    end if ! Grazing_detritus

                    !===========================================================================
                    ! ZOOPLANKTON RESPIRATION
                    !===========================================================================
                    ! Tracks heterotrophic respiration (metabolic CO2 release)
                    ! Key component of carbon cycling and food web energetics
                    !
                    ! Variables:
                    !   vertrespmeso(k)  : Mesozooplankton respiration [mmolC m-3 day-1]
                    !   HetRespFlux      : Meso respiration flux [mmolC m-3 day-1]
                    !
                    ! Ecological Significance:
                    !   - Represents metabolic carbon loss
                    !   - Typically 20-40% of ingested carbon
                    !   - Temperature-dependent (Q10 ≈ 2-3)
                    !   - Increases with activity level
                    !
                    ! Validation:
                    !   - Incubation experiments (O2 consumption)
                    !   - ETS (electron transport system) measurements
                    !   - Typical range: 0.01-1 mmolC m-3 day-1
                    !---------------------------------------------------------------------------

                    ! Mesozooplankton respiration
                    vertrespmeso(k) = vertrespmeso(k) + ( &
                            +HetRespFlux & ! Carbon respiration flux
                            ) * recipbiostep

                    if (enable_3zoo2det) then
                        ! Macrozooplankton (krill) respiration
                        vertrespmacro(k) = vertrespmacro(k) + ( &
                                +Zoo2RespFlux &
                                ) * recipbiostep

                        ! Microzooplankton respiration
                        vertrespmicro(k) = vertrespmicro(k) + ( &
                                +MicZooRespFlux &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! CALCITE DISSOLUTION
                    !===========================================================================
                    ! Tracks calcium carbonate dissolution in water column
                    ! Critical for carbonate counter-pump and alkalinity cycling
                    !
                    ! Variables:
                    !   vertcalcdiss(k) : Calcite dissolution rate [mmolC m-3 day-1]
                    !   calc_diss       : Dissolution rate coefficient [day-1]
                    !   DetCalc         : Detrital calcite concentration [mmolC m-3]
                    !
                    ! Process:
                    !   CaCO3 -> Ca²⁺ + CO3²⁻
                    !   Releases DIC (+1) and alkalinity (+2)
                    !
                    ! Controls:
                    !   - Saturation state (Ω < 1 favors dissolution)
                    !   - Temperature (higher T -> faster kinetics)
                    !   - Pressure (higher P -> lower Ω)
                    !
                    ! Depth Pattern:
                    !   - Surface: Minimal (supersaturated, Ω > 1)
                    !   - Intermediate: Moderate (Ω ≈ 1)
                    !   - Deep: High (undersaturated, Ω < 1)
                    !
                    ! Applications:
                    !   - Ocean acidification impacts
                    !   - Carbon export efficiency
                    !   - Sediment CaCO3 preservation
                    !---------------------------------------------------------------------------

                    vertcalcdiss(k) = vertcalcdiss(k) + ( &
                            +calc_diss * DetCalc & ! Dissolution flux
                            ) * recipbiostep

                    !===========================================================================
                    ! PARTICLE AGGREGATION
                    !===========================================================================
                    ! Tracks formation of large particles through aggregation
                    ! Key process for biological pump and carbon export
                    !
                    ! Variables:
                    !   vertaggn(k)       : Small phyto aggregation [mmolC m-3 day-1]
                    !   aggregationrate   : Aggregation rate coefficient [day-1]
                    !   PhyC              : Small phytoplankton carbon [mmolC m-3]
                    !
                    ! Aggregation Mechanisms:
                    !   - Brownian motion (small particles)
                    !   - Differential settling (size-dependent)
                    !   - Shear aggregation (turbulence)
                    !   - Mucus/TEP gluing (especially diatoms)
                    !
                    ! Enhancement Factors:
                    !   - High phytoplankton concentration
                    !   - Nutrient limitation (mucus production)
                    !   - Turbulence (increased collision rate)
                    !   - Sticky exudates (TEP, polysaccharides)
                    !
                    ! Ecological Significance:
                    !   - Increases sinking velocity (export)
                    !   - Forms marine snow
                    !   - Provides food for deep-sea organisms
                    !   - Removes surface biomass
                    !
                    ! Typical Rates:
                    !   - Low aggregation: <0.01 day-1
                    !   - Moderate: 0.01-0.1 day-1
                    !   - High (bloom collapse): >0.1 day-1
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton aggregation
                    vertaggn(k) = vertaggn(k) + ( &
                            +aggregationrate * PhyC & ! Aggregate formation
                            ) * recipbiostep

                    ! Diatom aggregation
                    vertaggd(k) = vertaggd(k) + ( &
                            +aggregationrate * DiaC &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore aggregation
                        vertaggc(k) = vertaggc(k) + ( &
                                +aggregationrate * CoccoC &
                                ) * recipbiostep

                        ! Phaeocystis aggregation
                        vertaggp(k) = vertaggp(k) + ( &
                                +aggregationrate * PhaeoC &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! DOC EXCRETION
                    !===========================================================================
                    ! Tracks dissolved organic carbon release by phytoplankton
                    ! Important for microbial loop and carbon cycling
                    !
                    ! Variables:
                    !   vertdocexn(k) : Small phyto DOC excretion [mmolC m-3 day-1]
                    !   lossC         : Carbon exudation rate constant [day-1]
                    !   limitFacN     : Nutrient limitation factor [0-1]
                    !   phyC          : Small phytoplankton carbon [mmolC m-3]
                    !
                    ! Exudation Mechanisms:
                    !   - Passive leakage (diffusion across membranes)
                    !   - Overflow metabolism (excess photosynthate)
                    !   - Stress response (nutrient limitation)
                    !   - Viral lysis (cell breakage)
                    !
                    ! Nutrient Stress Effect:
                    !   - High limitFacN (replete): Low exudation
                    !   - Low limitFacN (limited): High exudation
                    !   - C:N imbalance drives carbon overflow
                    !
                    ! DOC Fate:
                    !   - Labile fraction: Rapid bacterial uptake (hours-days)
                    !   - Semi-labile: Slower degradation (weeks-months)
                    !   - Refractory: Accumulates (years-millennia)
                    !
                    ! Ecological Role:
                    !   - Fuels microbial loop
                    !   - Carbon loss without grazing
                    !   - Nutrient regeneration (after remineralization)
                    !
                    ! Typical Rates:
                    !   - Percent of GPP: 5-30%
                    !   - Higher under nutrient stress
                    !   - Lower for healthy, growing cells
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton DOC excretion
                    vertdocexn(k) = vertdocexn(k) + ( &
                            +lossC * limitFacN * phyC & ! Stress-induced exudation
                            ) * recipbiostep

                    ! Diatom DOC excretion
                    vertdocexd(k) = vertdocexd(k) + ( &
                            +lossC_d * limitFacN_dia * DiaC &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore DOC excretion
                        vertdocexc(k) = vertdocexc(k) + ( &
                                +lossC_c * limitFacN_cocco * CoccoC &
                                ) * recipbiostep

                        ! Phaeocystis DOC excretion
                        vertdocexp(k) = vertdocexp(k) + ( &
                                +lossC_p * limitFacN_phaeo * PhaeoC &
                                ) * recipbiostep
                    end if

                    !===========================================================================
                    ! CALCIFICATION
                    !===========================================================================
                    ! Tracks calcium carbonate precipitation by phytoplankton
                    ! Key process for carbonate counter-pump
                    !
                    ! Variables:
                    !   vertcalcif(k)   : Calcification rate [mmolC m-3 day-1]
                    !   calcification   : CaCO3 precipitation flux [mmolC m-3 day-1]
                    !
                    ! Process:
                    !   Ca²⁺ + 2HCO3⁻ -> CaCO3 + CO2 + H2O
                    !   Consumes DIC (-1) and alkalinity (-2)
                    !   Releases CO2 (carbonate counter-pump)
                    !
                    ! Organisms:
                    !   - Coccolithophores (calcite plates)
                    !   - Foraminifera (calcite shells)
                    !   - Pteropods (aragonite shells - not in this model)
                    !
                    ! Environmental Controls:
                    !   - Temperature (warmer favors calcification)
                    !   - Carbonate saturation (Ω > 1 required)
                    !   - pH/CO2 (ocean acidification reduces rate)
                    !   - Nutrient availability (affects PIC:POC ratio)
                    !
                    ! Ecological Significance:
                    !   - Ballast for sinking (increases export)
                    !   - Protection from grazers
                    !   - CO2 source (opposite to photosynthesis)
                    !   - Sensitive to ocean acidification
                    !
                    ! Typical Rates:
                    !   - Surface (bloom): 0.1-2 mmolC m-3 day-1
                    !   - Background: 0.001-0.05 mmolC m-3 day-1
                    !---------------------------------------------------------------------------

                    vertcalcif(k) = vertcalcif(k) + ( &
                            +calcification & ! CaCO3 precipitation
                            ) * recipbiostep

                    !===========================================================================
                    ! PHYTOPLANKTON RESPIRATION
                    !===========================================================================
                    ! Tracks autotrophic respiration (maintenance costs)
                    ! Component of community respiration, complement to GPP
                    !
                    ! Variables:
                    !   vertrespn(k)    : Small phyto respiration [mmolC m-3 day-1]
                    !   PhyRespRate     : Respiration rate coefficient [day-1]
                    !   PhyC            : Small phytoplankton carbon [mmolC m-3]
                    !
                    ! Respiration Components:
                    !   - Basal metabolism: Maintenance of cellular machinery
                    !   - Biosynthesis: Costs of growth (protein synthesis, etc.)
                    !   - Active transport: Nutrient uptake against gradients
                    !   - Photorespiration: Oxygenase activity of Rubisco
                    !
                    ! Temperature Dependence:
                    !   - Q10 ≈ 2-3 (doubles per 10degC increase)
                    !   - Often higher than photosynthesis Q10
                    !   - Creates temperature-dependent growth efficiency
                    !
                    ! Relationship to NPP:
                    !   GPP = NPP + Respiration
                    !   Growth efficiency = NPP/GPP = 1 - (Resp/GPP)
                    !   Typical: 0.6-0.8 (60-80% efficiency)
                    !
                    ! Ecological Patterns:
                    !   - Increases with biomass (proportional)
                    !   - Higher in warm waters (temperature effect)
                    !   - Lower in nutrient-replete conditions (efficient growth)
                    !   - Higher during stress (maintenance costs)
                    !---------------------------------------------------------------------------

                    ! Small phytoplankton respiration
                    vertrespn(k) = vertrespn(k) + ( &
                            +PhyRespRate * PhyC & ! Maintenance respiration
                            ) * recipbiostep

                    ! Diatom respiration
                    vertrespd(k) = vertrespd(k) + ( &
                            +PhyRespRate_dia * DiaC &
                            ) * recipbiostep

                    if (enable_coccos) then
                        ! Coccolithophore respiration
                        vertrespc(k) = vertrespc(k) + ( &
                                +PhyRespRate_cocco * CoccoC &
                                ) * recipbiostep

                        ! Phaeocystis respiration
                        vertrespp(k) = vertrespp(k) + ( &
                                +PhyRespRate_phaeo * PhaeoC &
                                ) * recipbiostep
                    end if

                end if ! Diags

            end do ! Main vertikal loop ends

            !===============================================================================
            ! BENTHIC REMINERALIZATION
            !===============================================================================
            ! Simulates nutrient return from sediments to the bottom water layer through
            ! remineralization of organic matter and calcite dissolution.
            !
            ! This part calculates:
            !   1. Remineralization of benthic organic matter (N, C, Si)
            !   2. Dissolution of benthic calcite (CaCO3)
            !   3. Carbon isotope remineralization (13C, 14C - optional)
            !   4. Nutrient flux from sediments to bottom water
            !
            ! Key Features:
            !   - Single benthic layer (vertically integrated)
            !   - First-order decay kinetics for remineralization
            !   - Separate rates for N, C, Si, and calcite
            !   - Optional MEDUSA sediment flux forcing
            !   - Isotope fractionation during remineralization
            !   - Conservative tracer return to water column
            !
            ! Benthic Processes:
            !   - Aerobic respiration: Organic matter -> DIN + DIC
            !   - Silica dissolution: Biogenic SiO2 -> Si(OH)4
            !   - Calcite dissolution: CaCO3 -> Ca²⁺ + CO3²⁻ (affects DIC + Alk)
            !   - Denitrification: In anoxic sediments (not explicitly modeled)
            !
            ! Model Structure:
            !   - Benthos as single well-mixed compartment
            !   - Area-integrated concentrations [mmol m-2]
            !   - Decay rates calibrated to observations
            !   - Instant return to bottom water layer
            !
            ! Ecological/Biogeochemical Significance:
            !   - Closes nutrient cycles (prevents accumulation in benthos)
            !   - Critical for shallow shelf ecosystems
            !   - Supports benthic-pelagic coupling
            !   - Affects bottom water oxygen and pH
            !   - Important for coastal carbon budgets
            !
            ! Limitations:
            !   - No explicit burial (all material eventually remineralized)
            !   - No bioturbation or bioirrigation
            !   - No redox zonation in sediments
            !   - Simplified single-layer representation
            !
            ! References:
            !   - Schourup-Kristensen et al. (2013) - REcoM benthic module
            !   - MEDUSA forcing for data-constrained sediment fluxes
            !===============================================================================

            !===============================================================================
            ! REMINERALIZATION MODE SELECTION
            !===============================================================================
            ! Determines whether to use MEDUSA sediment flux forcing or internal
            ! benthic remineralization calculations.
            !
            ! Variables:
            !   use_MEDUSA      : Flag to enable MEDUSA sediment forcing [logical]
            !   sedflx_num      : Number of MEDUSA sediment flux entries [-]
            !   mype            : Processor ID for parallel output [integer]
            !
            ! MEDUSA Mode:
            !   - Uses observationally-constrained sediment fluxes
            !   - Overrides internal benthic calculations
            !   - Typically used for hindcast simulations
            !   - Requires external data file
            !
            ! Internal Mode:
            !   - Uses prognostic benthic compartment
            !   - First-order decay kinetics
            !   - Suitable for future projections
            !   - No external data required
            !-------------------------------------------------------------------------------

            if (use_MEDUSA .and. (sedflx_num /= 0)) then

                !===========================================================================
                ! MEDUSA SEDIMENT FLUX FORCING MODE
                !===========================================================================
                ! Use externally prescribed sediment nutrient fluxes
                ! This overrides all internal benthic calculations
                !---------------------------------------------------------------------------

                if (mype == 0) then
                    ! Print message only on master processor (parallel computing)
                    write(*, *) ' --> Sedimentary input of nutrients through MEDUSA'
                end if

                ! Note: MEDUSA fluxes are applied elsewhere in the code
                ! This section simply skips internal benthic calculations

            else

                !===========================================================================
                ! INTERNAL BENTHIC REMINERALIZATION (DEFAULT MODE)
                !===========================================================================
                ! Calculates nutrient return from benthic compartment using first-order
                ! decay kinetics for each element.
                !---------------------------------------------------------------------------

                !===========================================================================
                ! DISSOLVED INORGANIC NITROGEN (DIN) REMINERALIZATION
                !===========================================================================
                ! Aerobic remineralization of organic nitrogen in sediments releases
                ! bioavailable nitrogen (NH4+ and NO3-) back to the water column.
                !
                ! Variables:
                !   decayBenthos(1) : N remineralization flux [mmolN m-2 day-1]
                !   decayRateBenN   : Benthic N remineralization rate [day-1]
                !   LocBenthos(1)   : Benthic N pool (area-integrated) [mmolN m-2]
                !   dt_b            : REcoM time step [day]
                !
                ! Equation:
                !   Flux = decayRateBenN × LocBenthos(1)
                !   LocBenthos(1)_new = LocBenthos(1)_old - Flux × dt_b
                !
                ! Process Representation:
                !   Organic N (proteins, amino acids, etc.) -> NH4+ -> NO3-
                !   First-order decay: rate proportional to benthic N pool
                !
                ! Typical Values:
                !   decayRateBenN ~ 0.01-0.1 day-1 (10-100 day turnover)
                !   Higher rates in warm, oxygenated sediments
                !---------------------------------------------------------------------------

                ! Calculate N remineralization flux [mmolN m-2 day-1]
                decayBenthos(1) = decayRateBenN * LocBenthos(1)

                ! Update benthic N pool (remove remineralized N)
                LocBenthos(1) = LocBenthos(1) - decayBenthos(1) * dt_b

                !===========================================================================
                ! DISSOLVED INORGANIC CARBON (DIC) REMINERALIZATION
                !===========================================================================
                ! Aerobic respiration of organic carbon releases CO2 to bottom water.
                ! Affects carbonate system (pH, pCO2, saturation state).
                !
                ! Variables:
                !   decayBenthos(2) : C remineralization flux [mmolC m-2 day-1]
                !   decayRateBenC   : Benthic C remineralization rate [day-1]
                !   LocBenthos(2)   : Benthic C pool (area-integrated) [mmolC m-2]
                !
                ! Process Representation:
                !   Organic C (carbohydrates, lipids, etc.) + O2 -> CO2 + H2O
                !   Releases DIC and consumes O2 (not explicitly tracked)
                !
                ! Stoichiometry:
                !   C:N remineralization ratio typically near Redfield (6.6:1)
                !   Can vary with organic matter source and degradation state
                !
                ! Typical Values:
                !   decayRateBenC ~ 0.01-0.1 day-1
                !   May differ from N rate due to preferential degradation
                !---------------------------------------------------------------------------

                ! Calculate C remineralization flux [mmolC m-2 day-1]
                decayBenthos(2) = decayRateBenC * LocBenthos(2)

                ! Update benthic C pool
                LocBenthos(2) = LocBenthos(2) - decayBenthos(2) * dt_b

                !===========================================================================
                ! SILICATE (Si) DISSOLUTION
                !===========================================================================
                ! Dissolution of biogenic silica (diatom frustules) releases Si(OH)4.
                ! Temperature-dependent process (faster in warm water).
                !
                ! Variables:
                !   decayBenthos(3) : Si dissolution flux [mmolSi m-2 day-1]
                !   decayRateBenSi  : Benthic Si dissolution rate [day-1]
                !   LocBenthos(3)   : Benthic Si pool (area-integrated) [mmolSi m-2]
                !
                ! Process Representation:
                !   Biogenic SiO2 (opal) + H2O -> Si(OH)4 (silicic acid)
                !   Thermodynamically driven dissolution (opal undersaturated in seawater)
                !
                ! Typical Values:
                !   decayRateBenSi ~ 0.001-0.01 day-1 (100-1000 day turnover)
                !   Slower than organic matter (opal is recalcitrant)
                !   Temperature-dependent (see reminSiT calculation earlier)
                !
                ! Ecological Significance:
                !   - Returns Si to support diatom production
                !   - Critical in Si-limited regions
                !   - Permanent burial in deep sediments affects long-term Si cycle
                !---------------------------------------------------------------------------

                ! Calculate Si dissolution flux [mmolSi m-2 day-1]
                decayBenthos(3) = decayRateBenSi * LocBenthos(3)

                ! Update benthic Si pool
                LocBenthos(3) = LocBenthos(3) - decayBenthos(3) * dt_b

                !===========================================================================
                ! CALCITE (CaCO3) DISSOLUTION
                !===========================================================================
                ! Dissolution of calcium carbonate affects both DIC and alkalinity.
                ! Rate depends on saturation state (calculated in deepest water layer).
                !
                ! Variables:
                !   decayBenthos(4) : Calcite dissolution flux [mmolC m-2 day-1]
                !   calc_diss_ben   : Benthic calcite dissolution rate [day-1]
                !   LocBenthos(4)   : Benthic calcite pool [mmolC m-2]
                !
                ! Process Representation:
                !   CaCO3 + CO2 + H2O ⇌ Ca²⁺ + 2HCO3-
                !   Increases DIC (+1) and alkalinity (+2)
                !
                ! Rate Determination:
                !   calc_diss_ben taken from deepest water layer
                !   Uses either saturation-state or depth-dependent formulation
                !   (see calcite dissolution module earlier in code)
                !
                ! Carbonate Chemistry Effect:
                !   - Raises pH (consumes CO2)
                !   - Increases carbonate saturation
                !   - Buffers ocean acidification locally
                !
                ! Note: Changed from calc_diss to calc_diss_ben to allow independent
                !       control when using OmegaC_diss switch
                !---------------------------------------------------------------------------

                ! Calculate calcite dissolution flux [mmolC m-2 day-1]
                decayBenthos(4) = calc_diss_ben * LocBenthos(4)

                ! Update benthic calcite pool
                LocBenthos(4) = LocBenthos(4) - decayBenthos(4) * dt_b

                if (ciso) then

                    !=======================================================================
                    ! CARBON ISOTOPE REMINERALIZATION (OPTIONAL)
                    !=======================================================================
                    ! Tracks 13C and 14C through benthic remineralization processes.
                    ! Important for paleoclimate reconstructions and carbon cycle studies.
                    !-----------------------------------------------------------------------

                    !=======================================================================
                    ! Carbon-13 (13C) Remineralization
                    !=======================================================================
                    ! Organic matter remineralization with isotopic fractionation
                    !
                    ! Variables:
                    !   decayBenthos(5) : 13C remineralization flux [mmol13C m-2 day-1]
                    !   alpha_dcal_13   : 13C fractionation factor for remineralization [-]
                    !   LocBenthos(5)   : Benthic 13C pool [mmol13C m-2]
                    !
                    ! Note: Assumes same decay rate as total carbon (decayRateBenC)
                    !       Fractionation factor typically near 1.0 (minimal fractionation)
                    !-----------------------------------------------------------------------

                    ! Calculate 13C remineralization flux
                    ! Ignores isotopic fractionation during remineralization (alpha ≈ 1)
                    decayBenthos(5) = alpha_dcal_13 * decayRateBenC * LocBenthos(5)

                    ! Update benthic 13C pool
                    LocBenthos(5) = LocBenthos(5) - decayBenthos(5) * dt_b

                    !=======================================================================
                    ! Carbon-13 Calcite Dissolution
                    !=======================================================================
                    ! Dissolution of 13C-containing calcite
                    !
                    ! Variables:
                    !   decayBenthos(6) : 13C calcite dissolution flux [mmol13C m-2 day-1]
                    !   calc_diss_13    : 13C calcite dissolution rate [day-1]
                    !   LocBenthos(6)   : Benthic 13C-calcite pool [mmol13C m-2]
                    !-----------------------------------------------------------------------

                    ! Calculate 13C calcite dissolution flux
                    decayBenthos(6) = calc_diss_13 * LocBenthos(6)

                    ! Update benthic 13C-calcite pool
                    LocBenthos(6) = LocBenthos(6) - decayBenthos(6) * dt_b

                    if (ciso_14) then

                        if (ciso_organic_14) then

                            !===============================================================
                            ! Carbon-14 (14C) Remineralization
                            !===============================================================
                            ! Radiocarbon remineralization for age dating and residence times
                            !
                            ! Variables:
                            !   decayBenthos(7) : 14C remineralization flux [mmol14C m-2 day-1]
                            !   alpha_dcal_14   : 14C fractionation factor [-]
                            !   LocBenthos(7)   : Benthic 14C pool [mmol14C m-2]
                            !
                            ! Note: 14C also subject to radioactive decay (half-life ~5730 years)
                            !       This is typically handled separately in full carbon cycle
                            ! models
                            !---------------------------------------------------------------

                            ! Calculate 14C remineralization flux
                            ! Ignores isotopic fractionation during remineralization
                            decayBenthos(7) = alpha_dcal_14 * decayRateBenC * LocBenthos(7)

                            ! Update benthic 14C pool
                            LocBenthos(7) = LocBenthos(7) - decayBenthos(7) * dt_b

                            !===============================================================
                            ! 3.4 Carbon-14 Calcite Dissolution
                            !===============================================================
                            ! Dissolution of 14C-containing calcite
                            !
                            ! Variables:
                            !   decayBenthos(8) : 14C calcite dissolution flux [mmol14C m-2 day-1]
                            !   calc_diss_14    : 14C calcite dissolution rate [day-1]
                            !   LocBenthos(8)   : Benthic 14C-calcite pool [mmol14C m-2]
                            !---------------------------------------------------------------

                            ! Calculate 14C calcite dissolution flux
                            decayBenthos(8) = calc_diss_14 * LocBenthos(8)

                            ! Update benthic 14C-calcite pool
                            LocBenthos(8) = LocBenthos(8) - decayBenthos(8) * dt_b

                        else
                            !---------------------------------------------------------------
                            ! Alternative 14C handling
                            !---------------------------------------------------------------
                            ! When ciso_organic_14 = FALSE:
                            ! 14C defined as proportional to total DIC elsewhere in code
                            ! No separate remineralization calculation needed
                            ! (sms(idic_14) = sms(idic) × isotope_ratio)
                            !---------------------------------------------------------------

                        end if ! ciso_organic_14

                    end if ! ciso_14

                end if ! ciso

            end if ! use_MEDUSA

        end do ! Main time loop ends

    end subroutine REcoM_sms

    subroutine calculate_chlorophyll_degradation(maximum_photosynthesis_rate, &
            layer_available_radiation, &
            plastidic_chl_carbon_quota, &
            base_degradation_rate, &
            slope, &
            chlorophyll_loss_rate)

        use recom_declarations, only: wp
        use recom_config, only: tiny, one
        use ieee_arithmetic, only: ieee_is_nan

        implicit none

        real(kind=wp), intent(in) :: maximum_photosynthesis_rate
        real(kind=wp), intent(in) :: layer_available_radiation
        real(kind=wp), intent(in) :: plastidic_chl_carbon_quota
        real(kind=wp), intent(in) :: base_degradation_rate
        real(kind=wp), intent(in) :: slope
        real(kind=wp), intent(out) :: chlorophyll_loss_rate

        ! Keep a baseline degradation in darkness/invalid conditions.
        if (maximum_photosynthesis_rate < tiny .or. &
                ieee_is_nan(layer_available_radiation) .or. &
                ieee_is_nan(plastidic_chl_carbon_quota)) then
            chlorophyll_loss_rate = base_degradation_rate * 0.1d0
        else
            chlorophyll_loss_rate = base_degradation_rate * &
                    (real(one) - exp(-slope * plastidic_chl_carbon_quota * &
                    layer_available_radiation / &
                    maximum_photosynthesis_rate))

            chlorophyll_loss_rate = max((base_degradation_rate * 0.1d0), chlorophyll_loss_rate)
        end if
    end subroutine calculate_chlorophyll_degradation

    subroutine calculate_photosynthesis_rate(maximum_photosyntesis_rate, &
            layer_available_radiation, &
            chlorophyll_carbon_quota, &
            co2_limitation_factor, &
            slope, &
            photosynthesis_rate, &
            light_limitation_factor)

        use recom_declarations, only: wp
        use recom_config, only: tiny, zero, co2lim
        use ieee_arithmetic, only: ieee_is_nan

        implicit none

        real(kind=wp), intent(in) :: maximum_photosyntesis_rate, layer_available_radiation
        real(kind=wp), intent(in) :: chlorophyll_carbon_quota, co2_limitation_factor, slope
        real(kind=wp), intent(out) :: photosynthesis_rate, light_limitation_factor

        ! No photosynthesis in darkness or under invalid conditions
        if (maximum_photosyntesis_rate < tiny .or. &
                ieee_is_nan(layer_available_radiation) .or. &
                ieee_is_nan(chlorophyll_carbon_quota)) then
            photosynthesis_rate = zero
        else
            ! Calculate photosynthesis using exponential P-I curve
            ! Model saturates at high light (no photoinhibition)
            photosynthesis_rate = maximum_photosyntesis_rate * &
                    (1.0d0 - exp(-slope * chlorophyll_carbon_quota * &
                    layer_available_radiation / &
                    maximum_photosyntesis_rate))

            ! Store light limitation factor for diagnostics
            ! Ratio of actual to maximum rate indicates light limitation severity
            light_limitation_factor = photosynthesis_rate / maximum_photosyntesis_rate

            ! Apply CO2 limitation if ocean acidification sensitivity is enabled
            if (CO2lim) photosynthesis_rate = photosynthesis_rate * co2_limitation_factor
        end if

        ! Ensure non-negative values (numerical safety)
        if (photosynthesis_rate < tiny) photosynthesis_rate = zero
    end subroutine calculate_photosynthesis_rate

    !-------------------------------------------------------------------------------
    ! Function for calculating limiter
    !-------------------------------------------------------------------------------

    function recom_limiter(slope, qa, qb)
        use recom_declarations, only: wp
        use recom_config, only: REcoM_Geider_limiter

        implicit none

        real(kind=wp), intent(in) :: slope, qa, qb

        real(kind=wp) :: recom_limiter
        real(kind=wp) :: dq

        dq = qa - qb
        if (REcoM_Geider_limiter) then
            recom_limiter = max(min(-slope * dq, 1.d0), 0.d0)
        else
            recom_limiter = 1.d0 - exp(-slope * (abs(dq) - dq) ** 2)
        end if
        return
    end function recom_limiter

    !-------------------------------------------------------------------------------
    ! Function for iron chemistry
    !-------------------------------------------------------------------------------
    function iron_chemistry_2ligands(fet, l1t, l2t, k1, k2)
        use recom_declarations, only: wp

        implicit none

        real(kind=wp), intent(in) :: l1t, l2t, fet, k1, k2

        real(kind=wp) :: iron_chemistry_2ligands
        real(kind=wp) :: a3, a2, a1, a0, a, b, c, p, q, discr, rho, phi, amp, pi
        real(kind=wp) :: one3rd, one27th
        real(kind=wp) :: fe1, fe2, fe3

        ! coefficients of the 4th-order polynomial
        a3 = k1 * k2
        a2 = (k1 * k2 * (l1t + l2t - fet) + k1 + k2)
        a1 = (1 - (k1 + k2) * fet + k1 * l1t + k2 * l2t)
        a0 = -fet

        ! coefficients of the normalized polynomial
        a = a2 / a3
        b = a1 / a3
        c = a0 / a3

        ! some numbers that are used several times
        one3rd = 1.0 / 3.0
        one27th = 1.0 / 27.0

        ! now solve the polynomial stepwise
        p = b - a * a * one3rd
        q = c - a * b * one3rd + 2.0 * a * a * a * one27th
        discr = q * q / 4.0 + p * p * p * one27th

        rho = sqrt(-(p * p * p * one27th))
        phi = acos(-q / (2.0 * rho))
        amp = 2.0 * rho ** one3rd
        pi = 3.1415926535897931

        ! the equation has three real roots
        fe1 = amp * cos(phi * one3rd) - a * one3rd
        fe2 = amp * cos((phi + 2.0 * pi) * one3rd) - a * one3rd
        fe3 = amp * cos((phi + 4.0 * pi) * one3rd) - a * one3rd

        iron_chemistry_2ligands = max(fe1, fe2, fe3)

    end function iron_chemistry_2ligands

    !-------------------------------------------------------------------------------
    function iron_chemistry(Fe, totalLigand, ligandStabConst)
        use recom_declarations, only: wp
        use recom_config, only: tiny

        implicit none

        real(kind=wp), intent(in) :: Fe, totalLigand, ligandStabConst ! Input

        real(kind=wp) :: iron_chemistry
        real(kind=wp) :: FreeFe ! Output
        real(kind=wp) :: ligand, FeL, a, b, c, discrim

        ! Abbrevations
        a = ligandstabConst
        b = ligandstabConst * (Fe - totalLigand) + 1.d0
        c = -totalLigand
        discrim = b * b - 4.d0 * a * c

        if (abs(a) > tiny .and. discrim >= tiny) then
            ligand = (-b + sqrt(discrim)) / (2.d0 * a)
            FeL = totalLigand - ligand
            freeFe = Fe - FeL
        else ! No free iron
            freeFe = 0.d0
        end if

        iron_chemistry = freeFe

        return
    end function iron_chemistry

end module recom_sms_module
