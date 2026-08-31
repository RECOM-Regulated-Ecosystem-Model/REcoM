module recom_sms_update
    implicit none
    private

    public :: sms_update_tracer_scalars
    public :: sms_update_state

contains

    subroutine sms_update_tracer_scalars(n, k, state, sms, thick, Temp, Sali_depth, SurfSR, PAR, &
            kappa, DIN, DIC, Alk, PhyN, PhyC, PhyChl, DetN, DetC, HetN, HetC, DON, EOC, DiaN, &
            DiaC, DiaChl, DiaSi, DetSi, Si, Fe, PhyCalc, DetCalc, FreeFe, O2, CoccoN, CoccoC, &
            CoccoChl, PhaeoN, PhaeoC, PhaeoChl, Zoo2N, Zoo2C, DetZ2N, DetZ2C, DetZ2Si, DetZ2Calc, &
            MicZooN, MicZooC, recip_hetN_plus, REcoM_T_depth, REcoM_S_depth, REcoM_DIC_depth, &
            REcoM_Alk_depth, REcoM_Si_depth, REcoM_Phos_depth)

        use recom_declarations, only: wp, arrfunc, arrfunczoo2, calc_diss, chl2c, chl2c_cocco, &
                chl2c_dia, chl2c_phaeo, chl2c_plast, chl2c_plast_cocco, chl2c_plast_dia, &
                chl2c_plast_phaeo, chl2n, chl2n_cocco, chl2n_dia, chl2n_phaeo, chlave, kappastar, &
                kdzlower, kdzupper, lowerlight, o2func, parave, q10_mes, q10_mes_res, q10_mic, &
                q10_mic_res, qsic, qsin, quota, quota_cocco, quota_dia, quota_phaeo, recipdet, &
                recipdet2, recipquota, recipquota_cocco, recipquota_dia, recipquota_phaeo, &
                recipqzoo, recipqzoo2, recipqzoo3, reminsit, rtloc, rtref, temp_cocco, &
                temp_diatoms, temp_phaeo, temp_phyto, tiny_c, tiny_c_c, tiny_c_d, tiny_c_p, &
                tiny_n, tiny_n_c, tiny_n_d, tiny_n_p, tiny_si, vttemp_phyto, vttemp_diatoms, &
                vttemp_cocco, vttemp_phaeo, chl_lower, chl_upper

        use recom_config, only: a_chl, ae, beta_phaeo, c2k, chl2n_max, chl2n_max_c, chl2n_max_d, &
                chl2n_max_p, enable_3zoo2det, enable_coccos, expon_cocco, expon_d, expon_phy, &
                grazing_detritus, ialk, icchl, icocc, icocn, idchl, idetc, idetcal, idetn, &
                idetz2c, idetz2calc, idetz2n, idetz2si, idiac, idian, idiasi, idic, idin, idoc, &
                idon, ife, ihetc, ihetn, imiczooc, imiczoon, ioxy, ipchl, iphac, iphachl, iphan, &
                iphyc, iphyn, isi, izoo2c, izoo2n, k_o2_remin, k_w, ncmax, ncmax_c, ncmax_d, &
                ncmax_p, ncmin, ncmin_c, ncmin_d, ncmin_p, o2dep_remin, one, ord_cocco, ord_d, &
                ord_phy, recom_tref, res_het, sicmax, t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2, tiny, &
                tiny_chl, tiny_het, tmax_phaeo, topt_phaeo, uopt_phaeo, ciso, idetsi, reminsi, &
                zero, iphycal

        use recoM_ciso, only: alpha_dcal_13, alpha_dcal_14, calc_diss_13, calc_diss_14, ciso_14, &
                ciso_organic_14, detc_13, detc_14, detcalc_13, detcalc_14, diac_13, diac_14, &
                dic_13, dic_14, eoc_13, eoc_14, hetc_13, hetc_14, idetc_13, idetc_14, idetcal_13, &
                idetcal_14, idiac_13, idiac_14, idic_13, idic_14, idoc_13, idoc_14, ihetc_13, &
                ihetc_14, iphyc_13, iphyc_14, iphycal_13, iphycal_14, phyc_13, phyc_14, &
                phycalc_13, phycalc_14, quota_13, quota_14, quota_dia_13, quota_dia_14, &
                recipquota_13, recipquota_14, recipquota_dia_13, recipquota_dia_14, recipqzoo_13, &
                recipqzoo_14

        use recom_glovar, only: cosAI

        implicit none

        integer, intent(in) :: n, k

        real(kind=wp), intent(in) :: SurfSR
        real(kind=wp), intent(in), dimension(:) :: thick, Temp, Sali_depth
        real(kind=wp), intent(in), dimension(:, :) :: state, sms

        real(kind=wp), intent(inout) :: kappa
        real(kind=wp), intent(inout), dimension(:) :: PAR

        real(kind=wp), intent(inout) :: DIN, DIC, Alk, PhyN, PhyC, PhyChl, DetN, DetC, HetN, HetC, &
                DON, EOC, DiaN, DiaC, DiaChl, DiaSi, DetSi, Si, Fe, PhyCalc, DetCalc, FreeFe, O2, &
                CoccoN, CoccoC, CoccoChl, PhaeoN, PhaeoC, PhaeoChl, Zoo2N, Zoo2C, DetZ2N, DetZ2C, &
                DetZ2Si, DetZ2Calc, MicZooN, MicZooC, recip_hetN_plus

        real(kind=wp), intent(inout), dimension(1) :: REcoM_T_depth, REcoM_S_depth, REcoM_DIC_depth&
                , &
                REcoM_Alk_depth, REcoM_Si_depth, REcoM_Phos_depth

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
    end subroutine sms_update_tracer_scalars

    subroutine sms_update_state(k, dt_b, dt, MPI_COMM_FESOM, sms, state, &
            PhyN, PhyC, PhyChl, PhyCalc, &
            DiaN, DiaC, DiaChl, DiaSi, &
            CoccoN, CoccoC, CoccoChl, &
            PhaeoN, PhaeoC, PhaeoChl, &
            HetN, HetC, &
            Zoo2N, Zoo2C, &
            MicZooN, MicZooC, &
            DetN, DetC, DetSi, DetCalc, &
            DetZ2N, DetZ2C, DetZ2Si, DetZ2Calc, &
            DON, EOC, FreeFe)

        use recom_declarations, only: wp, arrfunc, o2func, &
                n_assim, n_assim_dia, n_assim_cocco, n_assim_phaeo, &
                cphot, cphot_dia, cphot_cocco, cphot_phaeo, &
                phyresprate, phyresprate_dia, phyresprate_cocco, phyresprate_phaeo, &
                limitfacn, limitfacn_dia, limitfacn_cocco, limitfacn_phaeo, &
                aggregationrate, grazeff, &
                grazingflux, grazingflux2, grazingflux3, &
                grazingflux_phy, grazingflux_phy2, grazingflux_phy3, &
                grazingflux_dia, grazingflux_dia2, grazingflux_dia3, &
                grazingflux_cocco, grazingflux_cocco2, grazingflux_cocco3, &
                grazingflux_phaeo, grazingflux_phaeo2, grazingflux_phaeo3, &
                grazingflux_het2, grazingflux_miczoo, grazingflux_miczoo2, &
                grazingflux_det, grazingflux_det2, grazingflux_detz2, grazingflux_detz22, &
                hetlossflux, zoo2lossflux, miczoolossflux, &
                hetrespflux, zoo2respflux, miczoorespflux, &
                mesfecalloss_c, mesfecalloss_n, zoo2fecalloss_c, zoo2fecalloss_n, &
                calc_diss, calc_diss2, calc_loss_agg, calc_loss_gra, &
                calc_loss_gra2, calc_loss_gra3, calcification, &
                chl2n, chl2n_dia, chl2n_cocco, chl2n_phaeo, &
                chlsynth, chlsynth_dia, chlsynth_cocco, chlsynth_phaeo, &
                kochl, kochl_dia, kochl_cocco, kochl_phaeo, &
                recipquota, recipquota_dia, recipquota_cocco, recipquota_phaeo, &
                recipdet, recipdet2, recipqzoo, recipqzoo2, recipqzoo3, &
                is_coccos, is_3zoo2det, si_assim, reminsit, qsin

        use recom_config, only: &
                idin, idic, ialk, iphyn, iphyc, ipchl, iphycal, &
                idian, idiac, idchl, idiasi, &
                icocn, icocc, icchl, &
                iphan, iphac, iphachl, &
                ihetn, ihetc, izoo2n, izoo2c, imiczoon, imiczooc, &
                idetn, idetc, idetsi, idetcal, &
                idetz2n, idetz2c, idetz2si, idetz2calc, &
                idon, idoc, isi, ife, ioxy, &
                grazing_detritus, enable_3zoo2det, enable_coccos, &
                lossn, lossn_d, lossn_c, lossn_p, &
                lossn_z, lossn_z2, lossn_z3, &
                lossc, lossc_d, lossc_c, lossc_p, &
                lossc_z, lossc_z2, lossc_z3, &
                reminn, reminc, rho_n, rho_c1, &
                fe2n, kscavfe, redo2c, calc_diss_guts, &
                grazeff2, grazeff3, ciso, tiny

        use recom_locvar, only: locriverdoc

        use recoM_ciso, only: ciso_14, ciso_organic_14, &
                calc_diss_13, calc_diss_14, &
                calc_loss_agg_13, calc_loss_agg_14, &
                calc_loss_gra_13, calc_loss_gra_14, &
                calcification_13, calcification_14, &
                idetc_13, idetc_14, idetcal_13, idetcal_14, &
                idiac_13, idiac_14, idic_13, idic_14, &
                idoc_13, idoc_14, ihetc_13, ihetc_14, &
                iphyc_13, iphyc_14, iphycal_13, iphycal_14, &
                phyc_13, phyc_14, diac_13, diac_14, &
                detc_13, detc_14, detcalc_13, detcalc_14, &
                hetc_13, hetc_14, eoc_13, eoc_14, &
                phycalc_13, phycalc_14, dic_13, dic_14, &
                recipquota_13, recipquota_14, &
                recipquota_dia_13, recipquota_dia_14, &
                recipqzoo_13, recipqzoo_14, &
                hetrespflux_13, hetrespflux_14, &
                r_iorg_13, r_iorg_14, &
                r_phyc_13, r_phyc_14, r_diac_13, r_diac_14

        implicit none

        ! ── Loop control & time
        ! ────────────────────────────────────────────
        integer, intent(in) :: k
        integer, intent(in) :: MPI_COMM_FESOM
        real(kind=wp), intent(in) :: dt_b
        real(kind=wp), intent(in) :: dt
        real(kind=wp), intent(inout) :: sms(:, :)
        real(kind=wp), intent(in) :: state(:, :)

        ! ── Small phytoplankton
        ! ────────────────────────────────────────────
        real(kind=wp), intent(in) :: PhyN, PhyC, PhyChl, PhyCalc
        ! ── Diatoms
        ! ────────────────────────────────────────────────────────
        real(kind=wp), intent(in) :: DiaN, DiaC, DiaChl, DiaSi
        ! ── Coccolithophores
        ! ───────────────────────────────────────────────
        real(kind=wp), intent(in) :: CoccoN, CoccoC, CoccoChl
        ! ── Phaeocystis
        ! ────────────────────────────────────────────────────
        real(kind=wp), intent(in) :: PhaeoN, PhaeoC, PhaeoChl
        ! ── Mesozooplankton
        ! ────────────────────────────────────────────────
        real(kind=wp), intent(in) :: HetN, HetC
        ! ── Macrozooplankton
        ! ──────────────────────────────────────────────
        real(kind=wp), intent(in) :: Zoo2N, Zoo2C
        ! ── Microzooplankton
        ! ──────────────────────────────────────────────
        real(kind=wp), intent(in) :: MicZooN, MicZooC
        ! ── Slow-sinking detritus
        ! ──────────────────────────────────────────
        real(kind=wp), intent(in) :: DetN, DetC, DetSi, DetCalc
        ! ── Fast-sinking detritus
        ! ─────────────────────────────────────────
        real(kind=wp), intent(in) :: DetZ2N, DetZ2C, DetZ2Si, DetZ2Calc
        ! ── Dissolved organics and free iron
        ! ──────────────────────────────
        real(kind=wp), intent(in) :: DON, EOC, FreeFe

        !===============================================================================
        ! 1. DISSOLVED INORGANIC NITROGEN (DIN)
        !===============================================================================
        ! Represents the pool of bioavailable nitrogen (nitrate + ammonium)
        !
        ! Variables:
        !   N_assim         : N assimilation rate for small phytoplankton [mmolN mmolC-1
        ! day-1]
        !   N_assim_Dia     : N assimilation rate for diatoms [mmolN mmolC-1 day-1]
        !   N_assim_Cocco   : N assimilation rate for coccolithophore [mmolN mmolC-1 day-1]
        !   N_assim_Phaeo   : N assimilation rate for Phaeocystis [mmolN mmolC-1 day-1]
        !   PhyC, DiaC      : Intracellular carbon concentration [mmolC m-3]
        !   CoccoC, PhaeoC  : Intracellular carbon concentration [mmolC m-3]
        !   rho_N           : Remineralization rate constant [day-1]
        !   arrFunc         : Arrhenius temperature dependency function [-]
        !   O2Func          : O2 dependency of organic matter remineralization [-]
        !   DON             : Dissolved organic nitrogen [mmolN m-3]
        !   dt_b            : REcoM time step [day]
        !
        ! Equation Reference: Schourup-Kristensen 2013, Eq. A2

        sms(k, idin) = ( &
        !---------------------------------------------------------------------------
        ! SINKS: Nitrogen Uptake (decreases DIN)
        !---------------------------------------------------------------------------
        ! Phytoplankton assimilation of NO3- and NH4+
                -N_assim * PhyC & ! Small phytoplankton
                - N_assim_Dia * DiaC & ! Diatoms
                - N_assim_Cocco * CoccoC * is_coccos & ! Coccolithophores
                - N_assim_Phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SOURCES: Remineralization (increases DIN)
        !---------------------------------------------------------------------------
        ! DON remineralization releases bioavailable nitrogen
                + rho_N * arrFunc * O2Func * DON & ! Temperature and O2 dependent
                ) * dt_b + sms(k, idin)

        !===============================================================================
        ! 2. DISSOLVED INORGANIC CARBON (DIC)
        !===============================================================================
        ! Represents the pool of inorganic carbon (CO2 + HCO3- + CO3-2)
        !
        ! Variables:
        !   Cphot           : Small phytoplankton photosynthesis rate [day-1]
        !   Cphot_Dia       : Diatom photosynthesis rate [day-1]
        !   Cphot_Cocco     : Coccolithophore photosynthesis rate [day-1]
        !   Cphot_Phaeo     : Phaeocystis photosynthesis rate [day-1]
        !   phyRespRate     : Small phytoplankton respiration rate [day-1]
        !   phyRespRate_Dia : Diatom respiration rate [day-1]
        !   phyRespRate_Cocco : Coccolithophore respiration rate [day-1]
        !   phyRespRate_Phaeo : Phaeocystis respiration rate [day-1]
        !   rho_C1          : Temperature-dependent DOC degradation rate [day-1]
        !   EOC             : Extracellular organic carbon [mmolC m-3]
        !   HetRespFlux     : Mesozooplankton respiration flux [mmolC m-3 day-1]
        !   Zoo2RespFlux    : Macrozooplankton respiration flux [mmolC m-3 day-1]
        !   MicZooRespFlux  : Microzooplankton respiration flux [mmolC m-3 day-1]
        !   calc_diss       : Slow-sinking calcite dissolution rate [day-1]
        !   calc_diss2      : Fast-sinking calcite dissolution rate [day-1]
        !   DetCalc         : Slow-sinking calcite detritus pool [mmolC m-3]
        !   DetZ2Calc       : Fast-sinking calcite detritus pool [mmolC m-3]
        !   calc_loss_gra   : Calcite loss via mesozooplankton grazing [mmolC m-3 day-1]
        !   calc_loss_gra2  : Calcite loss via macrozooplankton grazing [mmolC m-3 day-1]
        !   calc_loss_gra3  : Calcite loss via microzooplankton grazing [mmolC m-3 day-1]
        !   calc_diss_guts  : Calcite dissolution rate in zooplankton guts [-]
        !   calcification   : Rate of CaCO3 formation [mmolC m-3 day-1]

        sms(k, idic) = ( &
        !---------------------------------------------------------------------------
        ! SINKS: Carbon Fixation (decreases DIC)
        !---------------------------------------------------------------------------
        ! Photosynthetic uptake of CO2 by phytoplankton
                -Cphot * PhyC & ! Small phytoplankton
                - Cphot_Dia * DiaC & ! Diatoms
                - Cphot_Cocco * CoccoC * is_coccos & ! Coccolithophores
                - Cphot_Phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SINKS: Calcification (decreases DIC)
        !---------------------------------------------------------------------------
        ! CaCO3 formation: Ca2+ + CO3-2 -> CaCO3
                - calcification &
        !---------------------------------------------------------------------------
        ! SOURCES: Phytoplankton Respiration (increases DIC)
        !---------------------------------------------------------------------------
        ! Release of CO2 through autotrophic respiration
                + phyRespRate * PhyC & ! Small phytoplankton
                + phyRespRate_Dia * DiaC & ! Diatoms
                + phyRespRate_Cocco * CoccoC * is_coccos & ! Coccolithophores
                + phyRespRate_Phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SOURCES: DOC Remineralization (increases DIC)
        !---------------------------------------------------------------------------
        ! Microbial degradation of dissolved organic carbon
                + rho_C1 * arrFunc * O2Func * EOC & ! Temperature and O2 dependent
        !---------------------------------------------------------------------------
        ! SOURCES: Zooplankton Respiration (increases DIC)
        !---------------------------------------------------------------------------
        ! Release of CO2 through heterotrophic respiration
                + HetRespFlux & ! Mesozooplankton
                + Zoo2RespFlux * is_3zoo2det & ! Macrozooplankton
                + MicZooRespFlux * is_3zoo2det & ! Microzooplankton
        !---------------------------------------------------------------------------
        ! SOURCES: Calcite Dissolution (increases DIC)
        !---------------------------------------------------------------------------
        ! Reaction: CaCO3 + CO2 + H2O -> Ca2+ + 2HCO3-
                + calc_diss * DetCalc & ! Slow-sinking calcite
                + calc_loss_gra * calc_diss_guts & ! Mesozooplankton gut
                + calc_loss_gra2 * calc_diss_guts * is_3zoo2det & ! Macrozooplankton gut
                + calc_loss_gra3 * calc_diss_guts * is_3zoo2det & ! Microzooplankton gut
                + calc_diss2 * DetZ2Calc * is_3zoo2det & ! Fast-sinking calcite
                ) * dt_b + sms(k, idic)

        !  if((Latd(1)<-45.0) .and. ((state(k,idic)+sms(k,idic))>2500)) then
        !     !co2flux(1)=0.0
        !      print*,'ERROR: strange dic !'
        !      print*,'state(k,idic): ', state(k,idic)
        !      print*,'sms Cphot: ', -Cphot*PhyC
        !      print*,'sms resp: ', phyRespRate*PhyC
        !      print*,'sms Cphot dia: ', -Cphot_Dia*DiaC
        !      print*,'sms resp dia: ', phyRespRate_Dia * DiaC
        !      print*,'sms eoc: ', rho_C1* arrFunc *EOC
        !      print*,'sms het resp: ', HetRespFlux
        !      print*, 'sms co2: ',  dflux(1) * recipdzF(k) * max( 2-k, 0 )
        !      print*, 'sms calcdiss: ', calc_diss * DetCalc
        !      print*, 'sms calc_loss: ', calc_loss_gra * calc_diss_guts
        !      print*, 'sms calcification: ', -calcification
        !      stop
        !    endif

        !===============================================================================
        ! 3. ALKALINITY (Alk)
        !===============================================================================
        ! Total alkalinity affects ocean pH and CO2 uptake capacity
        ! Assumes constant Redfield N:P ratio
        !
        ! Key coefficient: 1.0625 = (1/16) + 1
        !   - Represents the change in alkalinity per mole of nitrogen
        !   - Includes both nitrate reduction (+ charge) and phosphate uptake

        sms(k, ialk) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Nutrient Uptake (increases alkalinity)
        !---------------------------------------------------------------------------
        ! Phytoplankton uptake of NO3- increases alkalinity
                +1.0625 * N_assim * PhyC & ! Small phytoplankton
                + 1.0625 * N_assim_Dia * DiaC & ! Diatoms
                + 1.0625 * N_assim_Cocco * CoccoC * is_coccos & ! Coccolithophores
                + 1.0625 * N_assim_Phaeo * PhaeoC * is_coccos & ! Phaeocystis

        !---------------------------------------------------------------------------
        ! SINKS: Remineralization (decreases alkalinity)
        !---------------------------------------------------------------------------
        ! DON remineralization releases H+ and decreases alkalinity
                - 1.0625 * rho_N * arrFunc * O2Func * DON &
        !---------------------------------------------------------------------------
        ! SOURCES: Calcite Dissolution (increases alkalinity)
        !---------------------------------------------------------------------------
        ! Reaction: CaCO3 + CO2 + H2O -> Ca2+ + 2HCO3-
        ! Increases alkalinity by 2 equivalents per mole CaCO3
                + 2.d0 * calc_diss * DetCalc & ! Slow-sinking calcite
                + 2.d0 * calc_loss_gra * calc_diss_guts & ! Mesozooplankton gut
        ! Macrozooplankton gut
                + 2.d0 * calc_loss_gra2 * calc_diss_guts * is_3zoo2det &
        ! Microzooplankton gut
                + 2.d0 * calc_loss_gra3 * calc_diss_guts * is_3zoo2det &
                + 2.d0 * calc_diss2 * DetZ2Calc * is_3zoo2det & ! Fast-sinking calcite

        !---------------------------------------------------------------------------
        ! SINKS: Calcification (decreases alkalinity)
        !---------------------------------------------------------------------------
        ! CaCO3 formation removes 2 equivalents of alkalinity
                - 2.d0 * calcification &
                ) * dt_b + sms(k, ialk)

        !===============================================================================
        ! SMALL PHYTOPLANKTON NITROGEN (PhyN)
        !===============================================================================
        ! Tracks the nitrogen content of small phytoplankton
        !
        ! Variables:
        !   N_assim            : N assimilation rate [day-1]
        !   lossN              : N loss rate [day-1]
        !   limitFacN          : Limiter function for N:C ratio regulation [-]
        !   aggregationRate    : Aggregation to detritus [day-1]
        !   grazingFlux_phy    : Mesozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_phy2   : Macrozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_phy3   : Microzooplankton grazing [mmolN m-3 day-1]

        sms(k, iphyn) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Nitrogen Assimilation
        !---------------------------------------------------------------------------
                +N_assim * PhyC &
        !---------------------------------------------------------------------------
        ! SINKS: Losses
        !---------------------------------------------------------------------------
                - lossN * limitFacN * PhyN & ! DON excretion (N:C regulated)
                - aggregationRate * PhyN & ! Aggregation to detritus
                - grazingFlux_phy & ! Mesozooplankton
                - grazingFlux_phy2 * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_phy3 * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, iphyn)

        !===============================================================================
        ! 5. SMALL PHYTOPLANKTON CARBON (PhyC)
        !===============================================================================
        ! Tracks the carbon content of small phytoplankton.
        !
        ! Variables:
        !   Cphot              : Gross photosynthesis rate [day-1]
        !   phyRespRate        : Autotrophic respiration rate [day-1]
        !   lossC              : C loss rate [day-1]
        !   recipQuota         : Reciprocal of N:C quota (for N->C conversion) [-]
        !
        ! Note: DOC excretion is downregulated by limitFacN when N:C ratio is too high

        sms(k, iphyc) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Net Photosynthesis
        !---------------------------------------------------------------------------
                +Cphot * PhyC & ! Gross photosynthesis
                - phyRespRate * PhyC & ! Autotrophic respiration
        !---------------------------------------------------------------------------
        ! SINKS: Losses
        !---------------------------------------------------------------------------
                - lossC * limitFacN * PhyC & ! DOC excretion (regulated)
                - aggregationRate * PhyC & ! Aggregation to detritus
                - grazingFlux_phy * recipQuota & ! Mesozooplankton (N->C)
                - grazingFlux_phy2 * recipQuota * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_phy3 * recipQuota * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, iphyc)

        !===============================================================================
        ! 6. PHYTOPLANKTON CHLOROPHYLL-A (PhyChl)
        !===============================================================================
        ! Tracks chlorophyll-a content for light harvesting and photoacclimation
        !
        ! Variables:
        !   chlSynth           : Chlorophyll synthesis rate [mgChl mmolC-1 day-1]
        !   KOchl              : Chlorophyll degradation rate constant [day-1]
        !   Chl2N              : Chl:N ratio = PhyChl/PhyN [mgChl mmolN-1]

        sms(k, ipchl) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Chlorophyll Synthesis
        !---------------------------------------------------------------------------
                +chlSynth * PhyC & ! Photoacclimation
        !---------------------------------------------------------------------------
        ! SINKS: Degradation and Losses
        !---------------------------------------------------------------------------

                - KOchl * PhyChl & ! Natural degradation
                - aggregationRate * PhyChl & ! Aggregation to detritus
                - grazingFlux_phy * Chl2N & ! Mesozooplankton
                - grazingFlux_phy2 * Chl2N * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_phy3 * Chl2N * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, ipchl)

        !===============================================================================
        ! 7. DETRITUS NITROGEN (DetN)
        !===============================================================================
        ! Tracks nitrogen content in slow-sinking organic particles.
        !
        ! Key Concepts:
        !   - Sloppy Feeding: Not all grazed material is assimilated
        !     Net detritus = Total grazing × (1 - grazing efficiency)
        !   - Four Configurations: Based on Grazing_detritus and enable_3zoo2det flags
        !
        ! Variables:
        !   grazEff, grazEff2, grazEff3         : Grazing efficiency (assimilation) [-]
        !   aggregationRate                     : Aggregation rate [day-1]
        !   hetLossFlux, miczooLossFlux         : Zooplankton mortality [mmolN m-3 day-1]
        !   reminN                              : Remineralization rate [day-1]
        !   arrFunc                             : Arrhenius temperature function [-]
        !   O2Func                              : Oxygen limitation function [-]
        !   grazingFlux_phy3, grazingFlux_dia3  : Grazing by microzooplankton [mmolN m-3
        ! day-1]
        !   grazingFlux_phy, grazingFlux_dia    : Grazing by mesozooplankton [mmolN m-3
        ! day-1]
        !   DetN                                : Detrital nitrogen concentration [mmolN
        ! m-3]
        !   dt_b                                : Time step [day]
        !-------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------
        ! Configuration 1: WITH Detritus Grazing + 3 Zooplankton Types
        !-------------------------------------------------------------------------------

        if (Grazing_detritus) then
            if (enable_3zoo2det) then
                sms(k, idetn) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Microzooplankton
                !-----------------------------------------------------------------------
                ! Net flux = Total grazing - Assimilated portion
                ! Small phytoplankton
                        +grazingFlux_phy3 - grazingFlux_phy3 * grazEff3 &
                        + grazingFlux_dia3 - grazingFlux_dia3 * grazEff3 & ! Diatoms

                ! Coccolithophores
                        + (grazingFlux_Cocco3 - grazingFlux_Cocco3 * grazEff3) * is_coccos &

                ! Phaeocystis
                        + (grazingFlux_Phaeo3 - grazingFlux_Phaeo3 * grazEff3) * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyN &
                        + aggregationRate * DiaN &
                        + aggregationRate * CoccoN * is_coccos &
                        + aggregationRate * PhaeoN * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + miczooLossFlux &
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption
                !-----------------------------------------------------------------------
                        - grazingFlux_Det * grazEff & ! Mesozooplankton
                        - grazingFlux_Det2 * grazEff2 & ! Macrozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetN & ! Bacterial decomposition
                        ) * dt_b + sms(k, idetn)
                !-------------------------------------------------------------------------------
                ! Configuration 2: WITH Detritus Grazing + 2 Zooplankton Types (Standard)
                !-------------------------------------------------------------------------------
            else
                sms(k, idetn) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton
                !-----------------------------------------------------------------------
                        +grazingFlux_phy - grazingFlux_phy * grazEff & ! Small phytoplankton
                        + grazingFlux_dia - grazingFlux_dia * grazEff & ! Diatoms
                ! Coccolithophores
                        + (grazingFlux_Cocco - grazingFlux_Cocco * grazEff) * is_coccos &
                ! Phaeocystis
                        + (grazingFlux_Phaeo - grazingFlux_Phaeo * grazEff) * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyN &
                        + aggregationRate * DiaN &
                        + aggregationRate * CoccoN * is_coccos &
                        + aggregationRate * PhaeoN * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + hetLossFlux &
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption
                !-----------------------------------------------------------------------
                        - grazingFlux_Det * grazEff & ! Mesozooplankton
                        - grazingFlux_Det2 * grazEff2 & ! Macrozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetN &
                        ) * dt_b + sms(k, idetn)
            end if
            !-------------------------------------------------------------------------------
            ! Configuration 3: WITHOUT Detritus Grazing + 3 Zooplankton Types
            !-------------------------------------------------------------------------------
        else
            if (enable_3zoo2det) then
                sms(k, idetn) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Direct Transfer from Grazing
                !-----------------------------------------------------------------------
                ! All grazed material enters detritus (no detritus grazing)
                        +grazingFlux_phy3 & ! Microzooplankton->small phyto
                        + grazingFlux_dia3 & ! Microzooplankton->diatoms
                        + grazingFlux_Cocco3 * is_coccos & ! Microzooplankton->coccoliths
                        + grazingFlux_Phaeo3 * is_coccos & ! Microzooplankton->Phaeocystis
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyN &
                        + aggregationRate * DiaN &
                        + aggregationRate * CoccoN * is_coccos &
                        + aggregationRate * PhaeoN * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + miczooLossFlux &
                !-----------------------------------------------------------------------
                ! SINKS: Generic Zooplankton Consumption
                !-----------------------------------------------------------------------
                        - grazingFlux * grazEff3 &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetN &
                        ) * dt_b + sms(k, idetn)

                !-------------------------------------------------------------------------------
                ! Configuration 4: WITHOUT Detritus Grazing + 2 Zooplankton Types
                !-------------------------------------------------------------------------------
            else
                sms(k, idetn) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Direct Transfer from Grazing
                !-----------------------------------------------------------------------
                        +grazingFlux_phy & ! Mesozooplankton->small phyto
                        + grazingFlux_dia & ! Mesozooplankton->diatoms
                        + grazingFlux_Cocco * is_coccos & ! Mesozooplankton->coccoliths
                        + grazingFlux_Phaeo * is_coccos & ! Mesozooplankton->Phaeocystis
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyN &
                        + aggregationRate * DiaN &
                        + aggregationRate * CoccoN * is_coccos &
                        + aggregationRate * PhaeoN * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + hetLossFlux &
                !-----------------------------------------------------------------------
                ! SINKS: Generic Zooplankton Consumption
                !-----------------------------------------------------------------------
                        - grazingFlux * grazEff &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetN &
                        ) * dt_b + sms(k, idetn)
            end if
        end if

        !===============================================================================
        ! 8. DETRITUS CARBON (DetC)
        !===============================================================================
        ! Tracks carbon content in slow-sinking organic particles
        !
        ! Key Concepts:
        !   - Stoichiometric Conversion: N-based fluxes -> C-based fluxes
        !     Uses reciprocal quotas (recipQuota = C:N ratio) for conversion
        !   - Sloppy Feeding: Net detritus = Total grazing × (1 - efficiency)
        !
        ! Variables:
        !   recipQuota, recipQuota_Dia, etc. : C:N quotas for phytoplankton [-]
        !   recipDet, recipDet2              : C:N ratios in detritus [-]
        !   recipQZoo, recipQZoo2, recipQZoo3: C:N ratios in zooplankton [-]
        !   reminC                              : C remineralization rate [day-1]
        !
        !===============================================================================
        ! KEY CONCEPTS:
        !===============================================================================
        ! 1. SLOPPY FEEDING: Not all grazed material is assimilated
        !    - Net detritus production = Total grazing × (1 - grazing efficiency)
        !    - Represents fecal pellets and inefficient consumption
        !
        ! 2. STOICHIOMETRIC CONVERSION: N-based fluxes -> C-based fluxes
        !    - recipQuota = C:N ratio of phytoplankton
        !    - recipDet = C:N ratio of detritus
        !    - recipQZoo = C:N ratio of zooplankton
        !
        ! 3. FOOD WEB CONFIGURATIONS:
        !    - Grazing_detritus ON: Zooplankton can feed on detritus (coprophagy)
        !    - Grazing_detritus OFF: Detritus only forms from grazing/mortality
        !    - enable_3zoo2det: Adds microzooplankton + fast-sinking detritus
        !
        ! 4. REMINERALIZATION: Temperature and oxygen dependent
        !    - arrFunc: Increases with temperature (Arrhenius kinetics)
        !    - O2Func: Decreases under low oxygen (anaerobic conditions)
        !===============================================================================
        !-------------------------------------------------------------------------------
        ! Configuration 1: WITH Detritus Grazing + 3 Zooplankton Types
        !-------------------------------------------------------------------------------
        if (Grazing_detritus) then
            if (enable_3zoo2det) then
                sms(k, idetc) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Microzooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy3 * recipQuota * (1.d0 - grazEff3) & ! Small phyto
                        + grazingFlux_Dia3 * recipQuota_Dia * (1.d0 - grazEff3) & ! Diatoms
                        + grazingFlux_Cocco3 * recipQuota_Cocco * (1.d0 - grazEff3) * &
                        is_coccos & ! Coccoliths
                        + grazingFlux_Phaeo3 * recipQuota_Phaeo * (1.d0 - grazEff3) * &
                        is_coccos & ! Phaeocystis
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation (C-basis)
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyC &
                        + aggregationRate * DiaC &
                        + aggregationRate * CoccoC * is_coccos &
                        + aggregationRate * PhaeoC * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + miczooLossFlux * recipQZoo3 & ! N->C conversion
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption (C-basis)
                !-----------------------------------------------------------------------
                        - grazingFlux_Det * recipDet * grazEff & ! Mesozooplankton
                        - grazingFlux_Det2 * recipDet * grazEff2 & ! Macrozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetC & ! Bacterial respiration
                        ) * dt_b + sms(k, idetc)

                !-------------------------------------------------------------------------------
                ! Configuration 2: WITH Detritus Grazing + 2 Zooplankton Types (Standard)
                !-------------------------------------------------------------------------------
            else
                sms(k, idetc) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy * recipQuota * (1.d0 - grazEff) &
                        + grazingFlux_Dia * recipQuota_Dia * (1.d0 - grazEff) &
                        + grazingFlux_Cocco * recipQuota_Cocco * (1.d0 - grazEff) * &
                        is_coccos &
                        + grazingFlux_Phaeo * recipQuota_Phaeo * (1.d0 - grazEff) * &
                        is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation (C-basis)
                !-----------------------------------------------------------------------
                        + aggregationRate * phyC &
                        + aggregationRate * DiaC &
                        + aggregationRate * CoccoC * is_coccos &
                        + aggregationRate * PhaeoC * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + hetLossFlux * recipQZoo &
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption (C-basis)
                !-----------------------------------------------------------------------
                        - grazingFlux_Det * recipDet * grazEff &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetC &
                        ) * dt_b + sms(k, idetc)

            end if

            !-------------------------------------------------------------------------------
            ! Configuration 3: WITHOUT Detritus Grazing + 3 Zooplankton Types
            !-------------------------------------------------------------------------------
        else
            if (enable_3zoo2det) then
                sms(k, idetc) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Microzooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy3 * recipQuota * (1.d0 - grazEff3) &
                        + grazingFlux_Dia3 * recipQuota_Dia * (1.d0 - grazEff3) &
                        + grazingFlux_Cocco3 * recipQuota_Cocco * (1.d0 - grazEff3) * &
                        is_coccos &
                        + grazingFlux_Phaeo3 * recipQuota_Phaeo * (1.d0 - grazEff3) * &
                        is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation (C-basis)
                !-----------------------------------------------------------------------
                        + aggregationRate * PhyC &
                        + aggregationRate * DiaC &
                        + aggregationRate * CoccoC * is_coccos &
                        + aggregationRate * PhaeoC * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + miczooLossFlux * recipQZoo3 &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetC &
                        ) * dt_b + sms(k, idetc)

                !-------------------------------------------------------------------------------
                ! Configuration 4: WITHOUT Detritus Grazing + 2 Zooplankton Types
                !-------------------------------------------------------------------------------
            else
                sms(k, idetc) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy * recipQuota * (1.d0 - grazEff) &
                        + grazingFlux_Dia * recipQuota_Dia * (1.d0 - grazEff) &
                        + grazingFlux_Cocco * recipQuota_Cocco * (1.d0 - grazEff) * &
                        is_coccos &
                        + grazingFlux_Phaeo * recipQuota_Phaeo * (1.d0 - grazEff) * &
                        is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Phytoplankton Aggregation (C-basis)
                !-----------------------------------------------------------------------
                        + aggregationRate * phyC &
                        + aggregationRate * DiaC &
                        + aggregationRate * CoccoC * is_coccos &
                        + aggregationRate * PhaeoC * is_coccos &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + hetLossFlux * recipQZoo &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetC &
                        ) * dt_b + sms(k, idetc)

            end if
        end if

        !===============================================================================
        ! 9. MESOZOOPLANKTON NITROGEN (HetN)
        !===============================================================================
        ! Primary herbivorous/omnivorous grazers that feed on phytoplankton and
        ! smaller prey items.
        !
        ! Variables:
        !   grazingFlux        : Total N grazing rate [mmolN m-3 day-1]
        !   grazEff            : Grazing/assimilation efficiency [-]
        !   grazingFlux_het2   : Predation by macrozooplankton [mmolN m-3 day-1]
        !   Mesfecalloss_n     : Fecal pellet production [mmolN m-3 day-1]
        !   hetLossFlux        : Mortality flux [mmolN m-3 day-1]
        !   lossN_z            : DON excretion rate [day-1]
        !-------------------------------------------------------------------------------

        sms(k, ihetn) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Grazing
        !---------------------------------------------------------------------------
                +grazingFlux * grazEff & ! Assimilated N
        !---------------------------------------------------------------------------
        ! SINKS: Predation, Mortality, Excretion
        !---------------------------------------------------------------------------
                - grazingFlux_het2 * is_3zoo2det & ! Predation by macrozooplankton
                - Mesfecalloss_n * is_3zoo2det & ! Fecal pellets
                - hetLossFlux & ! Mortality
                - lossN_z * HetN & ! DON excretion
                ) * dt_b + sms(k, ihetn)

        !===============================================================================
        ! 10. MESOZOOPLANKTON CARBON (HetC)
        !===============================================================================
        ! Carbon budget uses reciprocal quotas (C:N ratios) to convert N-based
        ! grazing rates to carbon equivalents.
        !
        ! Variables:
        !   recipQuota, recipQuota_Dia, etc. : C:N ratios of prey items [-]
        !   recipDet, recipDet2              : C:N ratios of detritus [-]
        !   recipQZoo, recipQZoo3            : C:N ratios of zooplankton [-]
        !   hetRespFlux                      : Respiration to CO2 [mmolC m-3 day-1]
        !   lossC_z                          : DOC excretion rate [day-1]
        !-------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------
        ! Configuration: Mesozooplankton CAN Graze on Detritus
        !-------------------------------------------------------------------------------
        if (Grazing_detritus) then
            sms(k, ihetc) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing (C-basis)
            !-----------------------------------------------------------------------
                    +grazingFlux_phy * recipQuota * grazEff & ! Small phytoplankton
                    + grazingFlux_Dia * recipQuota_Dia * grazEff & ! Diatoms
            ! Coccolithophores
                    + grazingFlux_Cocco * recipQuota_Cocco * grazEff * is_coccos &
            ! Phaeocystis
                    + grazingFlux_Phaeo * recipQuota_Phaeo * grazEff * is_coccos &
            ! Microzooplankton
                    + grazingFlux_miczoo * recipQZoo3 * grazEff * is_3zoo2det &
            ! Fast-sinking detritus
                    + grazingFlux_DetZ2 * recipDet2 * grazEff * is_3zoo2det &
                    + grazingFlux_Det * recipDet * grazEff & ! Slow-sinking detritus
            !-----------------------------------------------------------------------
            ! SINKS: Predation, Mortality, Respiration, Excretion
            !-----------------------------------------------------------------------
                    - grazingFlux_het2 * recipQZoo * is_3zoo2det & ! Predation
                    - Mesfecalloss_c * is_3zoo2det & ! Fecal pellets
                    - hetLossFlux * recipQZoo & ! Mortality
                    - lossC_z * HetC & ! DOC excretion
                    - hetRespFlux & ! Respiration to CO2
                    ) * dt_b + sms(k, ihetc)

            !-------------------------------------------------------------------------------
            ! Configuration: Mesozooplankton CANNOT Graze on Detritus
            !-------------------------------------------------------------------------------
        else
            sms(k, ihetc) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing (C-basis, herbivorous diet only)
            !-----------------------------------------------------------------------
                    +grazingFlux_phy * recipQuota * grazEff &
                    + grazingFlux_Dia * recipQuota_Dia * grazEff &
                    + grazingFlux_Cocco * recipQuota_Cocco * grazEff * is_coccos &
                    + grazingFlux_Phaeo * recipQuota_Phaeo * grazEff * is_coccos &
                    + grazingFlux_miczoo * recipQZoo3 * grazEff * is_3zoo2det &
            !-----------------------------------------------------------------------
            ! SINKS: Predation, Mortality, Respiration, Excretion
            !-----------------------------------------------------------------------
                    - grazingFlux_het2 * recipQZoo * is_3zoo2det &
                    - Mesfecalloss_c * is_3zoo2det &
                    - hetLossFlux * recipQZoo &
                    - lossC_z * HetC &
                    - hetRespFlux &
                    ) * dt_b + sms(k, ihetc)
        end if

        !===============================================================================
        ! 11. MACROZOOPLANKTON NITROGEN (Zoo2N)
        !===============================================================================
        ! Larger predatory zooplankton that feed on mesozooplankton, microzooplankton,
        ! and phytoplankton. Only active when enable_3zoo2det = .true.
        !
        ! Variables:
        !   grazingFlux2       : Total N grazing rate [mmolN m-3 day-1]
        !   grazEff2           : Grazing/assimilation efficiency [-]
        !   Zoo2LossFlux       : Mortality flux [mmolN m-3 day-1]
        !   lossN_z2           : DON excretion rate [day-1]
        !   Zoo2fecalloss_n    : Fecal pellet production [mmolN m-3 day-1]
        !-------------------------------------------------------------------------------

        if (enable_3zoo2det) then
            sms(k, izoo2n) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing
            !-----------------------------------------------------------------------
                    +grazingFlux2 * grazEff2 & ! Assimilated N
            !-----------------------------------------------------------------------
            ! SINKS: Mortality, Excretion, Fecal Pellets
            !-----------------------------------------------------------------------
                    - Zoo2LossFlux & ! Mortality
                    - lossN_z2 * Zoo2N & ! DON excretion
                    - Zoo2fecalloss_n & ! Fecal pellets
                    ) * dt_b + sms(k, izoo2n)

            !===============================================================================
            ! 12. MACROZOOPLANKTON CARBON (Zoo2C)
            !===============================================================================
            ! Carbon budget for macrozooplankton with stoichiometric conversions.
            !
            ! Variables:
            !   recipQuota, recipQuota_Dia, etc. : C:N ratios of prey [-]
            !   recipDet, recipDet2              : C:N ratios of detritus [-]
            !   recipQZoo, recipQZoo2, recipQZoo3: C:N ratios of zooplankton [-]
            !   Zoo2RespFlux                     : Respiration to CO2 [mmolC m-3 day-1]
            !   lossC_z2                         : DOC excretion rate [day-1]
            !-------------------------------------------------------------------------------

            !---------------------------------------------------------------------------
            ! Configuration: Macrozooplankton CAN Graze on Detritus
            !---------------------------------------------------------------------------

            if (Grazing_detritus) then
                sms(k, izoo2c) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Grazing (C-basis)
                !-------------------------------------------------------------------
                        +grazingFlux_phy2 * recipQuota * grazEff2 & ! Small phytoplankton
                        + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 & ! Diatoms
                ! Coccolithophores
                        + grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 * is_coccos &
                ! Phaeocystis
                        + grazingFlux_Phaeo2 * recipQuota_Phaeo * grazEff2 * is_coccos &
                ! Mesozooplankton (predation)
                        + grazingFlux_het2 * recipQZoo * grazEff2 &
                        + grazingFlux_miczoo2 * recipQZoo3 * grazEff2 & ! Microzooplankton
                        + grazingFlux_Det2 * recipDet * grazEff2 & ! Slow-sinking detritus
                ! Fast-sinking detritus
                        + grazingFlux_DetZ22 * recipDet2 * grazEff2 &
                !-------------------------------------------------------------------
                ! SINKS: Mortality, Respiration, Excretion, Fecal Pellets
                !-------------------------------------------------------------------
                        - zoo2LossFlux * recipQZoo2 & ! Mortality
                        - lossC_z2 * Zoo2C & ! DOC excretion
                        - Zoo2RespFlux & ! Respiration to CO2
                        - Zoo2fecalloss_c & ! Fecal pellets
                        ) * dt_b + sms(k, izoo2c)

                !---------------------------------------------------------------------------
                ! Configuration: Macrozooplankton CANNOT Graze on Detritus
                !---------------------------------------------------------------------------
            else
                sms(k, izoo2c) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Grazing (C-basis, no detritus feeding)
                !-------------------------------------------------------------------
                        +grazingFlux_phy2 * recipQuota * grazEff2 &
                        + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
                        + grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 * is_coccos &
                        + grazingFlux_Phaeo2 * recipQuota_Phaeo * grazEff2 * is_coccos &
                        + grazingFlux_het2 * recipQZoo * grazEff2 &
                        + grazingFlux_miczoo2 * recipQZoo3 * grazEff2 &
                !-------------------------------------------------------------------
                ! SINKS: Mortality, Respiration, Excretion, Fecal Pellets
                !-------------------------------------------------------------------
                        - zoo2LossFlux * recipQZoo2 &
                        - lossC_z2 * Zoo2C &
                        - Zoo2RespFlux &
                        - Zoo2fecalloss_c &
                        ) * dt_b + sms(k, izoo2c)

            end if

            !===============================================================================
            ! 13. MICROZOOPLANKTON NITROGEN (MicZooN)
            !===============================================================================
            ! Small heterotrophic protists that graze on phytoplankton and are prey for
            ! meso- and macrozooplankton.
            !
            ! Variables:
            !   grazingFlux3       : Total N grazing rate [mmolN m-3 day-1]
            !   grazEff3           : Grazing/assimilation efficiency [-]
            !   grazingFlux_miczoo : Predation by mesozooplankton [mmolN m-3 day-1]
            !   grazingFlux_miczoo2: Predation by macrozooplankton [mmolN m-3 day-1]
            !   MicZooLossFlux     : Mortality flux [mmolN m-3 day-1]
            !   lossN_z3           : DON excretion rate [day-1]
            !-------------------------------------------------------------------------------

            sms(k, imiczoon) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing
            !-----------------------------------------------------------------------
                    +grazingFlux3 * grazEff3 & ! Assimilated N
            !-----------------------------------------------------------------------
            ! SINKS: Predation, Mortality, Excretion
            !-----------------------------------------------------------------------
                    - grazingFlux_miczoo & ! Predation by mesozooplankton
                    - grazingFlux_miczoo2 & ! Predation by macrozooplankton
                    - MicZooLossFlux & ! Mortality
                    - lossN_z3 * MicZooN & ! DON excretion
                    ) * dt_b + sms(k, imiczoon)

            !===============================================================================
            ! 14. MICROZOOPLANKTON CARBON (MicZooC)
            !===============================================================================
            ! Carbon budget for microzooplankton with stoichiometric conversions.
            !
            ! Variables:
            !   recipQuota, recipQuota_Dia, etc. : C:N ratios of prey [-]
            !   recipQZoo3                       : C:N ratio of microzooplankton [-]
            !   MicZooRespFlux                   : Respiration to CO2 [mmolC m-3 day-1]
            !   lossC_z3                         : DOC excretion rate [day-1]
            !-------------------------------------------------------------------------------

            sms(k, imiczooc) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing (C-basis)
            !-----------------------------------------------------------------------
                    +grazingFlux_phy3 * recipQuota * grazEff3 & ! Small phytoplankton
                    + grazingFlux_Dia3 * recipQuota_Dia * grazEff3 & ! Diatoms
            ! Coccolithophores
                    + grazingFlux_Cocco3 * recipQuota_Cocco * grazEff3 * is_coccos &
            ! Phaeocystis
                    + grazingFlux_Phaeo3 * recipQuota_Phaeo * grazEff3 * is_coccos &
            !-----------------------------------------------------------------------
            ! SINKS: Predation, Mortality, Respiration, Excretion
            !-----------------------------------------------------------------------
                    - MicZooLossFlux * recipQZoo3 & ! Mortality
                    - grazingFlux_miczoo * recipQZoo3 & ! Predation by mesozooplankton
                    - grazingFlux_miczoo2 * recipQZoo3 & ! Predation by macrozooplankton
                    - lossC_z3 * MicZooC & ! DOC excretion
                    - MicZooRespFlux & ! Respiration to CO2
                    ) * dt_b + sms(k, imiczooc)

        end if

        !===============================================================================
        ! 15. FAST-SINKING DETRITUS NITROGEN (DetZ2N)
        !===============================================================================
        ! Particulate organic matter produced from zooplankton mortality, fecal pellets,
        ! and unassimilated grazing. Sinks faster than regular detritus.
        ! Only active when enable_3zoo2det = .true.
        !
        ! Variables:
        !   grazingFlux_phy, grazingFlux_phy2   : Grazing on small phyto [mmolN m-3 day-1]
        !   grazingFlux_dia, grazingFlux_dia2   : Grazing on diatoms [mmolN m-3 day-1]
        !   grazingFlux_het2                    : Predation on mesozooplankton [mmolN m-3
        ! day-1]
        !   grazingFlux_miczoo, grazingFlux_miczoo2 : Grazing on microzooplankton [mmolN
        ! m-3
        ! day-1]
        !   grazingFlux_DetZ2, grazingFlux_DetZ22   : Grazing on fast detritus [mmolN m-3
        ! day-1]
        !   Zoo2LossFlux, hetLossFlux           : Zooplankton mortality [mmolN m-3 day-1]
        !   Zoo2fecalloss_n, Mesfecalloss_n     : Fecal pellet production [mmolN m-3 day-1]
        !   reminN                              : Remineralization rate [day-1]
        !-------------------------------------------------------------------------------

        if (enable_3zoo2det) then

            !---------------------------------------------------------------------------
            ! Configuration: Zooplankton CAN Graze on Fast-Sinking Detritus
            !---------------------------------------------------------------------------
            if (Grazing_detritus) then
                sms(k, idetz2n) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Macrozooplankton
                !-----------------------------------------------------------------------
                        +grazingFlux_phy2 * (1.d0 - grazEff2) & ! Small phytoplankton
                        + grazingFlux_dia2 * (1.d0 - grazEff2) & ! Diatoms
                ! Coccoliths (meso)
                        + grazingFlux_Cocco * (1.d0 - grazEff) * is_coccos &
                ! Coccoliths (macro)
                        + grazingFlux_Cocco2 * (1.d0 - grazEff2) * is_coccos &
                ! Phaeocystis (meso)
                        + grazingFlux_Phaeo * (1.d0 - grazEff) * is_coccos &
                ! Phaeocystis (macro)
                        + grazingFlux_Phaeo2 * (1.d0 - grazEff2) * is_coccos &
                ! Mesozooplankton (predation)
                        + grazingFlux_het2 * (1.d0 - grazEff2) &
                        + grazingFlux_miczoo2 * (1.d0 - grazEff2) & ! Microzooplankton
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton
                !-----------------------------------------------------------------------
                        + grazingFlux_phy * (1.d0 - grazEff) & ! Small phytoplankton
                        + grazingFlux_dia * (1.d0 - grazEff) & ! Diatoms
                        + grazingFlux_miczoo * (1.d0 - grazEff) & ! Microzooplankton
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + Zoo2LossFlux & ! Macrozooplankton
                        + hetLossFlux & ! Mesozooplankton
                !-----------------------------------------------------------------------
                ! SOURCES: Fecal Pellet Production
                !-----------------------------------------------------------------------
                        + Zoo2fecalloss_n & ! Macrozooplankton
                        + Mesfecalloss_n & ! Mesozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption (Coprophagy)
                !-----------------------------------------------------------------------
                        - grazingFlux_DetZ2 * grazEff & ! Mesozooplankton
                        - grazingFlux_DetZ22 * grazEff2 & ! Macrozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetZ2N & ! Bacterial decomposition
                        ) * dt_b + sms(k, idetz2n)

                !---------------------------------------------------------------------------
                ! Configuration: Zooplankton CANNOT Graze on Fast-Sinking Detritus
                !---------------------------------------------------------------------------
            else
                sms(k, idetz2n) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Macrozooplankton
                !-----------------------------------------------------------------------
                        +grazingFlux_phy2 & ! All grazing enters detritus
                        + grazingFlux_dia2 &
                        + grazingFlux_Cocco * is_coccos &
                        + grazingFlux_Cocco2 * is_coccos &
                        + grazingFlux_Phaeo * is_coccos &
                        + grazingFlux_Phaeo2 * is_coccos &
                        + grazingFlux_het2 &
                        + grazingFlux_miczoo2 &
                        - grazingFlux2 * grazEff2 & ! Minus assimilated portion
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton
                !-----------------------------------------------------------------------
                        + grazingFlux_phy &
                        + grazingFlux_dia &
                        + grazingFlux_miczoo &
                        - grazingFlux * grazEff & ! Minus assimilated portion
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality
                !-----------------------------------------------------------------------
                        + Zoo2LossFlux &
                        + hetLossFlux &
                !-----------------------------------------------------------------------
                ! SOURCES: Fecal Pellet Production
                !-----------------------------------------------------------------------
                        + Zoo2fecalloss_n &
                        + Mesfecalloss_n &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminN * arrFunc * O2Func * DetZ2N &
                        ) * dt_b + sms(k, idetz2n)
            end if
            !===============================================================================
            ! 16. FAST-SINKING DETRITUS CARBON (DetZ2C)
            !===============================================================================
            ! Carbon budget for fast-sinking detritus with stoichiometric conversions.
            !
            ! Variables:
            !   recipQuota, recipQuota_Dia, etc. : C:N ratios of phytoplankton [-]
            !   recipQZoo, recipQZoo2, recipQZoo3: C:N ratios of zooplankton [-]
            !   recipDet2                        : C:N ratio of fast-sinking detritus [-]
            !   reminC                           : C remineralization rate [day-1]
            !-------------------------------------------------------------------------------

            !---------------------------------------------------------------------------
            ! Configuration: Zooplankton CAN Graze on Fast-Sinking Detritus
            !---------------------------------------------------------------------------
            if (Grazing_detritus) then
                sms(k, idetz2c) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Macrozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy2 * recipQuota * (1.d0 - grazEff2) & ! Small phyto
                        + grazingFlux_Dia2 * recipQuota_Dia * (1.d0 - grazEff2) & ! Diatoms
                        + grazingFlux_Cocco * recipQuota_Cocco * (1.d0 - grazEff) * &
                        is_coccos & ! Coccoliths (meso)
                        + grazingFlux_Cocco2 * recipQuota_Cocco * (1.d0 - grazEff2) * &
                        is_coccos & ! Coccoliths (macro)
                        + grazingFlux_Phaeo * recipQuota_Phaeo * (1.d0 - grazEff) * &
                        is_coccos & ! Phaeocystis (meso)
                        + grazingFlux_Phaeo2 * recipQuota_Phaeo * (1.d0 - grazEff2) * &
                        is_coccos & ! Phaeocystis (macro)
                ! Mesozooplankton
                        + grazingFlux_het2 * recipQZoo * (1.d0 - grazEff2) &
                ! Microzooplankton
                        + grazingFlux_miczoo2 * recipQZoo3 * (1.d0 - grazEff2) &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        + grazingFlux_phy * recipQuota * (1.d0 - grazEff) & ! Small phyto
                        + grazingFlux_Dia * recipQuota_Dia * (1.d0 - grazEff) & ! Diatoms
                ! Microzooplankton
                        + grazingFlux_miczoo * recipQZoo3 * (1.d0 - grazEff) &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + Zoo2LossFlux * recipQZoo2 & ! Macrozooplankton (N->C)
                        + hetLossFlux * recipQZoo & ! Mesozooplankton (N->C)
                !-----------------------------------------------------------------------
                ! SOURCES: Fecal Pellet Production (C-basis)
                !-----------------------------------------------------------------------
                        + Zoo2fecalloss_c & ! Macrozooplankton
                        + Mesfecalloss_c & ! Mesozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Detritus Consumption (C-basis)
                !-----------------------------------------------------------------------
                        - grazingFlux_DetZ2 * recipDet2 * grazEff & ! Mesozooplankton
                        - grazingFlux_DetZ22 * recipDet2 * grazEff2 & ! Macrozooplankton
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetZ2C & ! Bacterial respiration
                        ) * dt_b + sms(k, idetz2c)

                !---------------------------------------------------------------------------
                ! Configuration: Zooplankton CANNOT Graze on Fast-Sinking Detritus
                !---------------------------------------------------------------------------
            else
                sms(k, idetz2c) = ( &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Macrozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        +grazingFlux_phy2 * recipQuota * (1.d0 - grazEff2) &
                        + grazingFlux_Dia2 * recipQuota_Dia * (1.d0 - grazEff2) &
                        + grazingFlux_Cocco * recipQuota_Cocco * (1.d0 - grazEff) * &
                        is_coccos &
                        + grazingFlux_Cocco2 * recipQuota_Cocco * (1.d0 - grazEff2) * &
                        is_coccos &
                        + grazingFlux_Phaeo * recipQuota_Phaeo * (1.d0 - grazEff) * &
                        is_coccos &
                        + grazingFlux_Phaeo2 * recipQuota_Phaeo * (1.d0 - grazEff2) * &
                        is_coccos &
                        + grazingFlux_het2 * recipQZoo * (1.d0 - grazEff2) &
                        + grazingFlux_miczoo2 * recipQZoo3 * (1.d0 - grazEff2) &
                !-----------------------------------------------------------------------
                ! SOURCES: Sloppy Feeding by Mesozooplankton (C-basis)
                !-----------------------------------------------------------------------
                        + grazingFlux_phy * recipQuota * (1.d0 - grazEff) &
                        + grazingFlux_Dia * recipQuota_Dia * (1.d0 - grazEff) &
                        + grazingFlux_miczoo * recipQZoo3 * (1.d0 - grazEff) &
                !-----------------------------------------------------------------------
                ! SOURCES: Zooplankton Mortality (C-basis)
                !-----------------------------------------------------------------------
                        + Zoo2LossFlux * recipQZoo2 &
                        + hetLossFlux * recipQZoo &
                !-----------------------------------------------------------------------
                ! SOURCES: Fecal Pellet Production (C-basis)
                !-----------------------------------------------------------------------
                        + Zoo2fecalloss_c &
                        + Mesfecalloss_c &
                !-----------------------------------------------------------------------
                ! SINKS: Remineralization
                !-----------------------------------------------------------------------
                        - reminC * arrFunc * O2Func * DetZ2C &
                        ) * dt_b + sms(k, idetz2c)
            end if

            !===============================================================================
            ! 17. FAST-SINKING DETRITUS SILICA (DetZ2Si)
            !===============================================================================
            ! Biogenic silica from diatom frustules in fast-sinking detritus.
            !
            ! Variables:
            !   grazingFlux_dia, grazingFlux_dia2 : Grazing on diatoms [mmolN m-3 day-1]
            !   qSiN                              : Si:N ratio in diatoms [mmolSi mmolN-1]
            !   reminSiT                          : Temperature-dependent dissolution
            ! [day-1]
            !-------------------------------------------------------------------------------

            sms(k, idetz2si) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing on Diatoms
            !-----------------------------------------------------------------------
                    +grazingFlux_dia2 * qSiN & ! Macrozooplankton grazing
                    + grazingFlux_dia * qSiN & ! Mesozooplankton grazing
            !-----------------------------------------------------------------------
            ! SINKS: Dissolution
            !-----------------------------------------------------------------------
                    - reminSiT * DetZ2Si & ! Temperature-dependent
                    ) * dt_b + sms(k, idetz2si)

            !===============================================================================
            ! 18. FAST-SINKING DETRITUS CALCITE (DetZ2Calc)
            !===============================================================================
            ! Calcite particles from coccolithophore shells in fast-sinking detritus.
            !
            ! Variables:
            !   calc_loss_gra, calc_loss_gra2 : Calcite from grazing [mmolCaCO3 m-3 day-1]
            !   calc_diss_guts                : Gut dissolution fraction [-]
            !   calc_diss2                    : Water column dissolution rate [day-1]
            !-------------------------------------------------------------------------------

            sms(k, idetz2calc) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Grazing on Calcifying Phytoplankton
            !-----------------------------------------------------------------------
                    +calc_loss_gra2 * (1.d0 - calc_diss_guts) & ! Macrozooplankton (net)
                    + calc_loss_gra * (1.d0 - calc_diss_guts) & ! Mesozooplankton (net)
            !-----------------------------------------------------------------------
            ! SINKS: Dissolution in Water Column
            !-----------------------------------------------------------------------
                    - calc_diss2 * DetZ2Calc & ! CaCO3 dissolution
                    ) * dt_b + sms(k, idetz2calc)

        end if ! enable_3zoo2det

        !===============================================================================
        ! 19. DISSOLVED ORGANIC NITROGEN (DON)
        !===============================================================================
        ! Dissolved organic nitrogen pool from phytoplankton excretion, zooplankton
        ! metabolism, and detrital remineralization.
        !
        ! Variables:
        !   lossN, lossN_d, lossN_c, lossN_p : Phytoplankton DON excretion rates [day-1]
        !   limitFacN, limitFacN_Dia, etc.   : N:C ratio limiters (regulate excretion) [-]
        !   reminN                           : Detrital N remineralization rate [day-1]
        !   rho_N                            : DON remineralization rate [day-1]
        !   lossN_z, lossN_z2, lossN_z3      : Zooplankton DON excretion rates [day-1]
        !   arrFunc                          : Arrhenius temperature function [-]
        !   O2Func                           : Oxygen limitation function [-]
        !-------------------------------------------------------------------------------

        sms(k, idon) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Phytoplankton Excretion
        !---------------------------------------------------------------------------
                +lossN * limitFacN * phyN & ! Small phytoplankton
                + lossN_d * limitFacN_Dia * DiaN & ! Diatoms
                + lossN_c * limitFacN_Cocco * CoccoN * is_coccos & ! Coccolithophores
                + lossN_p * limitFacN_Phaeo * PhaeoN * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SOURCES: Detrital Remineralization
        !---------------------------------------------------------------------------
                + reminN * arrFunc * O2Func * DetN & ! Slow-sinking detritus
                + reminN * arrFunc * O2Func * DetZ2N * is_3zoo2det & ! Fast-sinking detritus
        !---------------------------------------------------------------------------
        ! SOURCES: Zooplankton Excretion
        !---------------------------------------------------------------------------
                + lossN_z * HetN & ! Mesozooplankton
                + lossN_z2 * Zoo2N * is_3zoo2det & ! Macrozooplankton
                + lossN_z3 * MicZooN * is_3zoo2det & ! Microzooplankton
        !---------------------------------------------------------------------------
        ! SINKS: Remineralization to NH4
        !---------------------------------------------------------------------------
                - rho_N * arrFunc * O2Func * DON & ! Bacterial remineralization
                ) * dt_b + sms(k, idon)

        !===============================================================================
        ! 20. EXTRACELLULAR ORGANIC CARBON (EOC / DOC)
        !===============================================================================
        ! Dissolved organic carbon pool from phytoplankton excretion, zooplankton
        ! metabolism, and detrital remineralization.
        !
        ! Variables:
        !   lossC, lossC_d, lossC_c, lossC_p : Phytoplankton DOC excretion rates [day-1]
        !   limitFacN, limitFacN_dia, etc.   : N:C ratio limiters (regulate excretion) [-]
        !   reminC                           : Detrital C remineralization rate [day-1]
        !   rho_c1                           : DOC remineralization rate [day-1]
        !   lossC_z, lossC_z2, lossC_z3      : Zooplankton DOC excretion rates [day-1]
        !-------------------------------------------------------------------------------

        sms(k, idoc) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Phytoplankton Excretion
        !---------------------------------------------------------------------------
                +lossC * limitFacN * phyC & ! Small phytoplankton
                + lossC_d * limitFacN_dia * DiaC & ! Diatoms
                + lossC_c * limitFacN_cocco * CoccoC * is_coccos & ! Coccolithophores
                + lossC_p * limitFacN_Phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SOURCES: Detrital Remineralization
        !---------------------------------------------------------------------------
                + reminC * arrFunc * O2Func * DetC & ! Slow-sinking detritus
                + reminC * arrFunc * O2Func * DetZ2C * is_3zoo2det & ! Fast-sinking detritus
        !---------------------------------------------------------------------------
        ! SOURCES: Zooplankton Excretion
        !---------------------------------------------------------------------------
                + lossC_z * HetC & ! Mesozooplankton
                + lossC_z2 * Zoo2C * is_3zoo2det & ! Macrozooplankton
                + lossC_z3 * MicZooC * is_3zoo2det & ! Microzooplankton
        !---------------------------------------------------------------------------
        ! SINKS: Remineralization to CO2
        !---------------------------------------------------------------------------
                - rho_c1 * arrFunc * O2Func * EOC & ! Bacterial respiration
                ) * dt_b + sms(k, idoc)

        !===============================================================================
        ! 21. DIATOM NITROGEN (DiaN)
        !===============================================================================
        ! Tracks nitrogen content of diatoms (large phytoplankton with silica frustules).
        !
        ! Variables:
        !   N_assim_dia                 : N assimilation rate [day-1]
        !   lossN_d                     : N loss rate [day-1]
        !   limitFacN_dia               : Limiter function for N:C ratio regulation [-]
        !   aggregationRate             : Aggregation to detritus [day-1]
        !   grazingFlux_Dia             : Mesozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_Dia2            : Macrozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_Dia3            : Microzooplankton grazing [mmolN m-3 day-1]
        !-------------------------------------------------------------------------------

        sms(k, idian) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Nitrogen Assimilation
        !---------------------------------------------------------------------------
                +N_assim_dia * DiaC &
        !---------------------------------------------------------------------------
        ! SINKS: DON Excretion
        !---------------------------------------------------------------------------
                - lossN_d * limitFacN_dia * DiaN &
        !---------------------------------------------------------------------------
        ! SINKS: Aggregation
        !---------------------------------------------------------------------------
                - aggregationRate * DiaN &
        !---------------------------------------------------------------------------
        ! SINKS: Grazing
        !---------------------------------------------------------------------------
                - grazingFlux_Dia & ! Mesozooplankton
                - grazingFlux_Dia2 * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_Dia3 * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, idian)

        !===============================================================================
        ! 22. DIATOM CARBON (DiaC)
        !===============================================================================
        ! Tracks carbon content of diatoms.
        !
        ! Variables:
        !   Cphot_dia                   : Gross photosynthesis rate [day-1]
        !   phyRespRate_dia             : Autotrophic respiration rate [day-1]
        !   lossC_d                     : C loss rate [day-1]
        !   recipQuota_dia              : Reciprocal of N:C quota (for N->C conversion) [-]
        !-------------------------------------------------------------------------------

        sms(k, idiac) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Net Photosynthesis
        !---------------------------------------------------------------------------
                +Cphot_dia * DiaC & ! Gross photosynthesis
        !---------------------------------------------------------------------------
        ! SINKS: DOC Excretion
        !---------------------------------------------------------------------------
                - lossC_d * limitFacN_dia * DiaC &
        !---------------------------------------------------------------------------
        ! SINKS: Respiration
        !---------------------------------------------------------------------------
                - phyRespRate_dia * DiaC &
        !---------------------------------------------------------------------------
        ! SINKS: Aggregation
        !---------------------------------------------------------------------------
                - aggregationRate * DiaC &
        !---------------------------------------------------------------------------
        ! SINKS: Grazing (C-basis)
        !---------------------------------------------------------------------------
                - grazingFlux_dia * recipQuota_dia & ! Mesozooplankton (N->C)
                - grazingFlux_dia2 * recipQuota_dia * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_dia3 * recipQuota_dia * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, idiac)

        !===============================================================================
        ! 23. DIATOM CHLOROPHYLL-A (DiaChl)
        !===============================================================================
        ! Tracks chlorophyll-a content for light harvesting and photoacclimation.
        !
        ! Variables:
        !   chlSynth_dia                : Chlorophyll synthesis rate [mgChl mmolC-1 day-1]
        !   KOchl_dia                   : Chlorophyll degradation rate [day-1]
        !   Chl2N_dia                   : Chl:N ratio = DiaChl/DiaN [mgChl mmolN-1]
        !-------------------------------------------------------------------------------

        sms(k, idchl) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Chlorophyll Synthesis
        !---------------------------------------------------------------------------
                +chlSynth_dia * DiaC & ! Photoacclimation
        !---------------------------------------------------------------------------
        ! SINKS: Photo-oxidation
        !---------------------------------------------------------------------------
                - KOchl_dia * DiaChl &
        !---------------------------------------------------------------------------
        ! SINKS: Aggregation
        !---------------------------------------------------------------------------
                - aggregationRate * DiaChl &
        !---------------------------------------------------------------------------
        ! SINKS: Grazing (Chl-basis)
        !---------------------------------------------------------------------------
                - grazingFlux_dia * Chl2N_dia & ! Mesozooplankton (N->Chl)
                - grazingFlux_dia2 * Chl2N_dia * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_dia3 * Chl2N_dia * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, idchl)

        !===============================================================================
        ! 24. DIATOM SILICA (DiaSi)
        !===============================================================================
        ! Tracks biogenic silica content in diatom frustules.
        !
        ! Variables:
        !   Si_assim                    : Silicic acid assimilation rate [day-1]
        !   qSiN                        : Si:N ratio in diatoms [mmolSi mmolN-1]
        !-------------------------------------------------------------------------------

        sms(k, idiasi) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Silicon Assimilation
        !---------------------------------------------------------------------------
                +Si_assim * DiaC &
        !---------------------------------------------------------------------------
        ! SINKS: Silicon Excretion
        !---------------------------------------------------------------------------
                - lossN_d * limitFacN_dia * DiaSi &
        !---------------------------------------------------------------------------
        ! SINKS: Aggregation
        !---------------------------------------------------------------------------
                - aggregationRate * DiaSi &
        !---------------------------------------------------------------------------
        ! SINKS: Grazing (Si-basis)
        !---------------------------------------------------------------------------
                - grazingFlux_dia * qSiN & ! Mesozooplankton (N->Si)
                - grazingFlux_dia2 * qSiN * is_3zoo2det & ! Macrozooplankton
                - grazingFlux_dia3 * qSiN * is_3zoo2det & ! Microzooplankton
                ) * dt_b + sms(k, idiasi)

        !===============================================================================
        ! 25. COCCOLITHOPHORE NITROGEN (CoccoN)
        !===============================================================================
        ! Tracks nitrogen content of coccolithophores (calcifying small phytoplankton).
        ! Only active when enable_coccos = .true.
        !
        ! Variables:
        !   N_assim_cocco               : N assimilation rate [day-1]
        !   lossN_c                     : N loss rate [day-1]
        !   limitFacN_cocco             : Limiter function for N:C ratio regulation [-]
        !   aggregationRate             : Aggregation to detritus [day-1]
        !   grazingFlux_Cocco           : Mesozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_Cocco2          : Macrozooplankton grazing [mmolN m-3 day-1]
        !   grazingFlux_Cocco3          : Microzooplankton grazing [mmolN m-3 day-1]
        !-------------------------------------------------------------------------------

        if (enable_coccos) then
            sms(k, icocn) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Nitrogen Assimilation
            !-----------------------------------------------------------------------
                    +N_assim_cocco * CoccoC &
            !-----------------------------------------------------------------------
            ! SINKS: DON Excretion
            !-----------------------------------------------------------------------
                    - lossN_c * limitFacN_cocco * CoccoN &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * CoccoN &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing
            !-----------------------------------------------------------------------
                    - grazingFlux_Cocco & ! Mesozooplankton
                    - grazingFlux_Cocco2 * is_3zoo2det & ! Macrozooplankton
                    - grazingFlux_Cocco3 * is_3zoo2det & ! Microzooplankton
                    ) * dt_b + sms(k, icocn)

            !===============================================================================
            ! 26. COCCOLITHOPHORE CARBON (CoccoC)
            !===============================================================================
            ! Tracks carbon content of coccolithophores.
            !
            ! Variables:
            !   Cphot_cocco                 : Gross photosynthesis rate [day-1]
            !   phyRespRate_cocco           : Autotrophic respiration rate [day-1]
            !   lossC_c                     : C loss rate [day-1]
            !   recipQuota_cocco            : Reciprocal of N:C quota (for N->C conversion)
            ! [-]
            !-------------------------------------------------------------------------------

            sms(k, icocc) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Net Photosynthesis
            !-----------------------------------------------------------------------
                    +Cphot_cocco * CoccoC & ! Gross photosynthesis
            !-----------------------------------------------------------------------
            ! SINKS: DOC Excretion
            !-----------------------------------------------------------------------
                    - lossC_c * limitFacN_cocco * CoccoC &
            !-----------------------------------------------------------------------
            ! SINKS: Respiration
            !-----------------------------------------------------------------------
                    - phyRespRate_cocco * CoccoC &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * CoccoC &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing (C-basis)
            !-----------------------------------------------------------------------
                    - grazingFlux_cocco * recipQuota_cocco & ! Mesozooplankton (N->C)
            ! Macrozooplankton
                    - grazingFlux_Cocco2 * recipQuota_cocco * is_3zoo2det &
            ! Microzooplankton
                    - grazingFlux_Cocco3 * recipQuota_cocco * is_3zoo2det &
                    ) * dt_b + sms(k, icocc)

            !---------------------------------------------------------------------------
            ! Error Check: Unrealistic CoccoC Growth
            !---------------------------------------------------------------------------
            if (sms(k, icocc) > 100) then
                print*, 'ERROR: Unrealistic CoccoC growth detected!'
                print*, 'k= ', k
                print*, 'dt= ', dt
                print*, 'dt_b= ', dt_b
                print*, 'state(k,icocc): ', state(k, icocc)
                print*, 'CoccoC: ', CoccoC
                print*, 'CoccoN: ', CoccoN
                print*, 'Cphot_cocco: ', Cphot_cocco * CoccoC
                print*, 'lossC_c: ', lossC_c
                print*, 'limitFacN_cocco: ', limitFacN_cocco
                print*, 'phyRespRate_cocco: ', phyRespRate_cocco
                print*, 'grazingFlux_cocco: ', grazingFlux_cocco
                print*, 'grazingFlux_Cocco2: ', grazingFlux_Cocco2
                print*, 'grazingFlux_Cocco3: ', grazingFlux_Cocco3
                print*, 'recipQuota_cocco: ', recipQuota_cocco
                call MPI_ABORT(MPI_COMM_FESOM, 1)
                stop
            end if

            !===============================================================================
            ! 27. COCCOLITHOPHORE CHLOROPHYLL-A (CoccoChl)
            !===============================================================================
            ! Tracks chlorophyll-a content for light harvesting and photoacclimation.
            !
            ! Variables:
            !   ChlSynth_cocco              : Chlorophyll synthesis rate [mgChl mmolC-1
            ! day-1]
            !   KOchl_cocco                 : Chlorophyll degradation rate [day-1]
            !   Chl2N_cocco                 : Chl:N ratio = CoccoChl/CoccoN [mgChl mmolN-1]
            !-------------------------------------------------------------------------------

            sms(k, icchl) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Chlorophyll Synthesis
            !-----------------------------------------------------------------------
                    +ChlSynth_cocco * CoccoC & ! Photoacclimation
            !-----------------------------------------------------------------------
            ! SINKS: Photo-oxidation
            !-----------------------------------------------------------------------
                    - KOchl_cocco * CoccoChl &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * CoccoChl &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing (Chl-basis)
            !-----------------------------------------------------------------------
                    - grazingFlux_cocco * Chl2N_cocco & ! Mesozooplankton (N->Chl)
                    - grazingFlux_Cocco2 * Chl2N_cocco * is_3zoo2det & ! Macrozooplankton
                    - grazingFlux_Cocco3 * Chl2N_cocco * is_3zoo2det & ! Microzooplankton
                    ) * dt_b + sms(k, icchl)

            !===============================================================================
            ! 28. PHAEOCYSTIS NITROGEN (PhaeoN)
            !===============================================================================
            ! Tracks nitrogen content of Phaeocystis (colony-forming phytoplankton).
            ! Only active when enable_coccos = .true.
            !
            ! Variables:
            !   N_assim_phaeo               : N assimilation rate [day-1]
            !   lossN_p                     : N loss rate [day-1]
            !   limitFacN_phaeo             : Limiter function for N:C ratio regulation [-]
            !   aggregationRate             : Aggregation to detritus [day-1]
            !   grazingFlux_phaeo           : Mesozooplankton grazing [mmolN m-3 day-1]
            !   grazingFlux_phaeo2          : Macrozooplankton grazing [mmolN m-3 day-1]
            !   grazingFlux_phaeo3          : Microzooplankton grazing [mmolN m-3 day-1]
            !-------------------------------------------------------------------------------

            sms(k, iphan) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Nitrogen Assimilation
            !-----------------------------------------------------------------------
                    +N_assim_phaeo * PhaeoC &
            !-----------------------------------------------------------------------
            ! SINKS: DON Excretion
            !-----------------------------------------------------------------------
                    - lossN_p * limitFacN_phaeo * PhaeoN &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * PhaeoN &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing
            !-----------------------------------------------------------------------
                    - grazingFlux_phaeo & ! Mesozooplankton
                    - grazingFlux_phaeo2 * is_3zoo2det & ! Macrozooplankton
                    - grazingFlux_phaeo3 * is_3zoo2det & ! Microzooplankton
                    ) * dt_b + sms(k, iphan)

            !===============================================================================
            ! 29. PHAEOCYSTIS CARBON (PhaeoC)
            !===============================================================================
            ! Tracks carbon content of Phaeocystis.
            !
            ! Variables:
            !   Cphot_phaeo                 : Gross photosynthesis rate [day-1]
            !   phyRespRate_phaeo           : Autotrophic respiration rate [day-1]
            !   lossC_p                     : C loss rate [day-1]
            !   recipQuota_phaeo            : Reciprocal of N:C quota (for N->C conversion)
            ! [-]
            !-------------------------------------------------------------------------------

            sms(k, iphac) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Net Photosynthesis
            !-----------------------------------------------------------------------
                    +Cphot_phaeo * PhaeoC & ! Gross photosynthesis
            !-----------------------------------------------------------------------
            ! SINKS: DOC Excretion
            !-----------------------------------------------------------------------
                    - lossC_p * limitFacN_phaeo * PhaeoC &
            !-----------------------------------------------------------------------
            ! SINKS: Respiration
            !-----------------------------------------------------------------------
                    - phyRespRate_phaeo * PhaeoC &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * PhaeoC &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing (C-basis)
            !-----------------------------------------------------------------------
                    - grazingFlux_phaeo * recipQuota_phaeo & ! Mesozooplankton (N->C)
            ! Macrozooplankton
                    - grazingFlux_phaeo2 * recipQuota_phaeo * is_3zoo2det &
            ! Microzooplankton
                    - grazingFlux_phaeo3 * recipQuota_phaeo * is_3zoo2det &
                    ) * dt_b + sms(k, iphac)

            !===============================================================================
            ! 30. PHAEOCYSTIS CHLOROPHYLL-A (PhaeoChl)
            !===============================================================================
            ! Tracks chlorophyll-a content for light harvesting and photoacclimation.
            !
            ! Variables:
            !   chlSynth_phaeo              : Chlorophyll synthesis rate [mgChl mmolC-1
            ! day-1]
            !   KOchl_phaeo                 : Chlorophyll degradation rate [day-1]
            !   Chl2N_phaeo                 : Chl:N ratio = PhaeoChl/PhaeoN [mgChl mmolN-1]
            !-------------------------------------------------------------------------------

            sms(k, iphachl) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Chlorophyll Synthesis
            !-----------------------------------------------------------------------
                    +chlSynth_phaeo * PhaeoC & ! Photoacclimation
            !-----------------------------------------------------------------------
            ! SINKS: Photo-oxidation
            !-----------------------------------------------------------------------
                    - KOchl_phaeo * PhaeoChl &
            !-----------------------------------------------------------------------
            ! SINKS: Aggregation
            !-----------------------------------------------------------------------
                    - aggregationRate * PhaeoChl &
            !-----------------------------------------------------------------------
            ! SINKS: Grazing (Chl-basis)
            !-----------------------------------------------------------------------
                    - grazingFlux_phaeo * Chl2N_phaeo & ! Mesozooplankton (N->Chl)
                    - grazingFlux_phaeo2 * Chl2N_phaeo * is_3zoo2det & ! Macrozooplankton
                    - grazingFlux_phaeo3 * Chl2N_phaeo * is_3zoo2det & ! Microzooplankton
                    ) * dt_b + sms(k, iphachl)

        end if ! enable_coccos

        !===============================================================================
        ! 31. DETRITAL SILICA (DetSi)
        !===============================================================================
        ! Biogenic silica from diatom frustules in slow-sinking detritus.
        !
        ! Variables:
        !   aggregationRate             : Diatom aggregation rate [day-1]
        !   lossN_d                     : Diatom mortality/excretion rate [day-1]
        !   grazingFlux_dia, grazingFlux_dia3 : Grazing on diatoms [mmolN m-3 day-1]
        !   qSiN                        : Si:N ratio in diatoms [mmolSi mmolN-1]
        !   reminSiT                    : Temperature-dependent dissolution [day-1]
        !-------------------------------------------------------------------------------

        sms(k, idetsi) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Diatom Aggregation
        !---------------------------------------------------------------------------
                +aggregationRate * DiaSi &
        !---------------------------------------------------------------------------
        ! SOURCES: Diatom Excretion/Mortality
        !---------------------------------------------------------------------------
                + lossN_d * limitFacN_dia * DiaSi &
        !---------------------------------------------------------------------------
        ! SOURCES: Grazing on Diatoms (Si-basis)
        !---------------------------------------------------------------------------
                + grazingFlux_dia3 * qSiN * is_3zoo2det & ! Microzooplankton
        ! Mesozooplankton (when 3zoo disabled)
                + grazingFlux_dia * qSiN * (1.0 - is_3zoo2det) &
        !---------------------------------------------------------------------------
        ! SINKS: Dissolution
        !---------------------------------------------------------------------------
                - reminSiT * DetSi & ! Temperature-dependent
                ) * dt_b + sms(k, idetsi)

        !===============================================================================
        ! 32. DISSOLVED SILICATE (DSi)
        !===============================================================================
        ! Dissolved silicate available for diatom uptake.
        ! Based on Schourup 2013 Eq. A3
        !
        ! Variables:
        !   Si_assim                    : Silicic acid assimilation by diatoms [day-1]
        !   reminSiT                    : Temperature-dependent dissolution [day-1]
        !   DetSi, DetZ2Si              : Detrital silica pools [mmolSi m-3]
        !-------------------------------------------------------------------------------

        sms(k, isi) = ( &
        !---------------------------------------------------------------------------
        ! SINKS: Biological Uptake
        !---------------------------------------------------------------------------
                -Si_assim * DiaC &
        !---------------------------------------------------------------------------
        ! SOURCES: Remineralization
        !---------------------------------------------------------------------------
                + reminSiT * DetSi & ! Slow-sinking detritus
                + reminSiT * DetZ2Si * is_3zoo2det & ! Fast-sinking detritus
                ) * dt_b + sms(k, isi)

        !===============================================================================
        ! 33. DISSOLVED IRON (Fe)
        !===============================================================================
        ! Tracks dissolved iron, a limiting micronutrient for phytoplankton growth.
        !
        ! Key Concept: Iron cycling is coupled to nitrogen via the Fe:N ratio (Fe2N)
        !              All N-based fluxes are converted to Fe equivalents
        !
        ! Variables:
        !   Fe2N                        : Intracellular Fe:N ratio [μmol Fe mmol N-1]
        !                                 Note: Fe2N = Fe2C × 6.625 (Redfield conversion)
        !   N_assim, N_assim_dia, etc.  : N assimilation rates [mmol N mmol C-1 day-1]
        !   lossN, lossN_d, etc.        : N excretion rates [day-1]
        !   limitFacN, etc.             : Nutrient limitation factors [-]
        !   reminN                      : Temperature-dependent remineralization [day-1]
        !   kScavFe                     : Iron scavenging rate [m3 mmol C-1 day-1]
        !   FreeFe                      : Free dissolved iron [μmol Fe m-3]
        !-------------------------------------------------------------------------------

        sms(k, ife) = ( &
        !---------------------------------------------------------------------------
        ! Iron Uptake/Release Coupled to Nitrogen Cycling (via Fe:N ratio)
        !---------------------------------------------------------------------------
                Fe2N * ( &
        !-----------------------------------------------------------------------
        ! SINKS: Phytoplankton Assimilation
        !-----------------------------------------------------------------------
                -N_assim * PhyC & ! Small phytoplankton
                - N_assim_dia * DiaC & ! Diatoms
                - N_assim_cocco * CoccoC * is_coccos & ! Coccolithophores
                - N_assim_phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !-----------------------------------------------------------------------
        ! SOURCES: Phytoplankton Excretion
        !-----------------------------------------------------------------------
                + lossN * limitFacN * PhyN & ! Small phytoplankton
                + lossN_d * limitFacN_dia * DiaN & ! Diatoms
                + lossN_c * limitFacN_cocco * CoccoN * is_coccos & ! Coccolithophores
                + lossN_p * limitFacN_phaeo * PhaeoN * is_coccos & ! Phaeocystis
        !-----------------------------------------------------------------------
        ! SOURCES: Detrital Remineralization
        !-----------------------------------------------------------------------
                + reminN * arrFunc * O2Func * DetN & ! Slow-sinking detritus
                + reminN * arrFunc * O2Func * DetZ2N * is_3zoo2det & ! Fast-sinking detritus
        !-----------------------------------------------------------------------
        ! SOURCES: Zooplankton Excretion
        !-----------------------------------------------------------------------
                + lossN_z * HetN & ! Mesozooplankton
                + lossN_z2 * Zoo2N * is_3zoo2det & ! Macrozooplankton
                + lossN_z3 * MicZooN * is_3zoo2det & ! Microzooplankton
                ) &
        !---------------------------------------------------------------------------
        ! SINKS: Abiotic Iron Scavenging onto Particles
        !---------------------------------------------------------------------------
                - kScavFe * DetC * FreeFe & ! Slow-sinking detritus
                - kScavFe * DetZ2C * FreeFe * is_3zoo2det & ! Fast-sinking detritus
                ) * dt_b + sms(k, ife)

        !===============================================================================
        ! 34. PHYTOPLANKTON CALCITE (PhyCalc)
        !===============================================================================
        ! Tracks calcium carbonate in living phytoplankton (coccoliths).
        !
        ! Variables:
        !   calcification               : CaCO3 production rate [mmolCaCO3 m-3 day-1]
        !   calc_loss_agg               : Calcite loss to aggregation [mmolCaCO3 m-3 day-1]
        !   calc_loss_gra, calc_loss_gra2, calc_loss_gra3 : Calcite loss to grazing
        ! [mmolCaCO3 m-3 day-1]
        !   lossC, lossC_c              : C excretion rates [day-1]
        !   phyRespRate, phyRespRate_cocco : Respiration rates [day-1]
        !   limitFacN, limitFacN_cocco  : N:C ratio limiters [-]
        !-------------------------------------------------------------------------------

        if (enable_coccos) then
            !---------------------------------------------------------------------------
            ! Configuration: Coccolithophore-Specific Calcite Dynamics
            !---------------------------------------------------------------------------
            sms(k, iphycal) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Calcification
            !-----------------------------------------------------------------------
                    +calcification & ! New CaCO3 production
            !-----------------------------------------------------------------------
            ! SINKS: Losses from Living Cells
            !-----------------------------------------------------------------------
                    - lossC_c * limitFacN_cocco * PhyCalc & ! Excretion/exudation
                    - phyRespRate_cocco * PhyCalc & ! Respiration-associated loss
                    - calc_loss_agg & ! Aggregation/sinking
                    - calc_loss_gra & ! Mesozooplankton grazing
                    - calc_loss_gra2 * is_3zoo2det & ! Macrozooplankton grazing
                    - calc_loss_gra3 * is_3zoo2det & ! Microzooplankton grazing
                    ) * dt_b + sms(k, iphycal)
        else
            !---------------------------------------------------------------------------
            ! Configuration: Generic Phytoplankton Calcite (Small Calcifiers)
            !---------------------------------------------------------------------------
            sms(k, iphycal) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Calcification
            !-----------------------------------------------------------------------
                    +calcification &
            !-----------------------------------------------------------------------
            ! SINKS: Losses from Living Cells
            !-----------------------------------------------------------------------
                    - lossC * limitFacN * PhyCalc &
                    - phyRespRate * PhyCalc &
                    - calc_loss_agg &
                    - calc_loss_gra &
                    - calc_loss_gra2 * is_3zoo2det &
                    - calc_loss_gra3 * is_3zoo2det &
                    ) * dt_b + sms(k, iphycal)
        end if

        !===============================================================================
        ! 35. DETRITAL CALCITE (DetCalc)
        !===============================================================================
        ! Tracks calcium carbonate in slow-sinking organic particles.
        !
        ! Variables:
        !   calc_loss_agg               : Calcite from aggregation [mmolCaCO3 m-3 day-1]
        !   calc_loss_gra, calc_loss_gra3 : Calcite from grazing [mmolCaCO3 m-3 day-1]
        !   calc_diss_guts              : Gut dissolution fraction [-]
        !   calc_diss                   : Water column dissolution rate [day-1]
        !-------------------------------------------------------------------------------

        if (enable_coccos) then
            if (enable_3zoo2det) then
                !-----------------------------------------------------------------------
                ! Configuration: Coccolithophore Calcite with 3-Zooplankton Model
                !-----------------------------------------------------------------------
                sms(k, idetcal) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Transfer from Living Cells
                !-------------------------------------------------------------------
                        +lossC_c * limitFacN_cocco * PhyCalc & ! Excretion
                        + phyRespRate_cocco * PhyCalc & ! Respiration products
                        + calc_loss_agg & ! Aggregation products
                        + calc_loss_gra3 & ! Microzooplankton grazing
                !-------------------------------------------------------------------
                ! SINKS: Dissolution
                !-------------------------------------------------------------------
                        - calc_loss_gra3 * calc_diss_guts & ! Gut dissolution
                        - calc_diss * DetCalc & ! Water column dissolution
                        ) * dt_b + sms(k, idetcal)
            else
                !-----------------------------------------------------------------------
                ! Configuration: Coccolithophore Calcite with Standard Zooplankton
                !-----------------------------------------------------------------------
                sms(k, idetcal) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Transfer from Living Cells
                !-------------------------------------------------------------------
                        +lossC_c * limitFacN_cocco * PhyCalc &
                        + phyRespRate_cocco * PhyCalc &
                        + calc_loss_agg &
                        + calc_loss_gra &
                !-------------------------------------------------------------------
                ! SINKS: Dissolution
                !-------------------------------------------------------------------
                        - calc_loss_gra * calc_diss_guts &
                        - calc_diss * DetCalc &
                        ) * dt_b + sms(k, idetcal)
            end if
        else
            if (enable_3zoo2det) then
                !-----------------------------------------------------------------------
                ! Configuration: Generic Phytoplankton Calcite with 3-Zooplankton
                !-----------------------------------------------------------------------
                sms(k, idetcal) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Transfer from Living Cells
                !-------------------------------------------------------------------
                        +lossC * limitFacN * PhyCalc &
                        + phyRespRate * PhyCalc &
                        + calc_loss_agg &
                        + calc_loss_gra3 &
                !-------------------------------------------------------------------
                ! SINKS: Dissolution
                !-------------------------------------------------------------------
                        - calc_loss_gra3 * calc_diss_guts &
                        - calc_diss * DetCalc &
                        ) * dt_b + sms(k, idetcal)
            else
                !-----------------------------------------------------------------------
                ! Configuration: Generic Phytoplankton Calcite with Standard Zooplankton
                !-----------------------------------------------------------------------
                sms(k, idetcal) = ( &
                !-------------------------------------------------------------------
                ! SOURCES: Transfer from Living Cells
                !-------------------------------------------------------------------
                        +lossC * limitFacN * PhyCalc &
                        + phyRespRate * PhyCalc &
                        + calc_loss_agg &
                        + calc_loss_gra &
                !-------------------------------------------------------------------
                ! SINKS: Dissolution
                !-------------------------------------------------------------------
                        - calc_loss_gra * calc_diss_guts &
                        - calc_diss * DetCalc &
                        ) * dt_b + sms(k, idetcal)
            end if
        end if

        !===============================================================================
        ! 36. DISSOLVED OXYGEN (O2)
        !===============================================================================
        ! Tracks oxygen production (photosynthesis) and consumption (respiration,
        ! remineralization).
        !
        ! Variables:
        !   Cphot, Cphot_dia, etc.      : Gross photosynthesis rates [day-1]
        !   phyRespRate, phyRespRate_dia, etc. : Autotrophic respiration rates [day-1]
        !   rho_C1                      : DOC remineralization rate [day-1]
        !   hetRespFlux, Zoo2RespFlux, MicZooRespFlux : Zooplankton respiration [mmolC m-3
        ! day-1]
        !   redO2C                      : O2:C stoichiometric ratio (Redfield) [-]
        !                                 Typically ~170/122 = 1.39 mol O2/mol C
        !-------------------------------------------------------------------------------

        sms(k, ioxy) = ( &
        !---------------------------------------------------------------------------
        ! SOURCES: Photosynthetic Oxygen Production
        !---------------------------------------------------------------------------
                +Cphot * phyC & ! Small phytoplankton
                + Cphot_dia * diaC & ! Diatoms
                + Cphot_cocco * CoccoC * is_coccos & ! Coccolithophores
                + Cphot_phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SINKS: Autotrophic Respiration
        !---------------------------------------------------------------------------
                - phyRespRate * phyC & ! Small phytoplankton
                - phyRespRate_dia * diaC & ! Diatoms
                - phyRespRate_cocco * CoccoC * is_coccos & ! Coccolithophores
                - phyRespRate_phaeo * PhaeoC * is_coccos & ! Phaeocystis
        !---------------------------------------------------------------------------
        ! SINKS: Heterotrophic Respiration and Remineralization
        !---------------------------------------------------------------------------
                - rho_C1 * arrFunc * O2Func * EOC & ! DOC remineralization
                - hetRespFlux & ! Mesozooplankton
                - Zoo2RespFlux * is_3zoo2det & ! Macrozooplankton
                - MicZooRespFlux * is_3zoo2det & ! Microzooplankton
                ) * redO2C * dt_b + sms(k, ioxy)
        ! Note: redO2C converts carbon-based rates to oxygen equivalents
        !       using the Redfield ratio (typically ~170/122 = 1.39 mol O2/mol C)

        if (ciso) then

            !===========================================================================
            ! 1. CARBON-13 (13C) BUDGETS
            !===========================================================================
            ! Calculates 13C budgets for all carbon pools with isotope fractionation.
            ! Parallel structure to total carbon budgets.
            !---------------------------------------------------------------------------

            !===========================================================================
            ! DISSOLVED INORGANIC CARBON (DIC_13)
            !===========================================================================
            ! Source-minus-sink budget for 13C in dissolved inorganic carbon pool.
            !
            ! SOURCES (+):
            !   - Phytoplankton respiration (returns 13C to DIC)
            !   - Diatom respiration
            !   - DOC remineralization (aerobic respiration)
            !   - Heterotroph respiration
            !   - Calcite dissolution (releases 13C from CaCO3)
            !   - Calcite dissolution in grazer guts
            !
            ! SINKS (-):
            !   - Phytoplankton photosynthesis (fixes 13C into organic matter)
            !   - Diatom photosynthesis
            !   - Calcification (removes 13C for CaCO3 formation)
            !
            ! Variables:
            !   Cphot, Cphot_Dia       : Photosynthesis rates [day-1]
            !   PhyC_13, DiaC_13       : Phytoplankton 13C pools [mmol13C m-3]
            !   phyRespRate            : Phytoplankton respiration rates [day-1]
            !   rho_C1                 : DOC remineralization rate [day-1]
            !   arrFunc                : Temperature function [-]
            !   EOC_13                 : Dissolved organic 13C [mmol13C m-3]
            !   HetRespFlux_13         : Heterotroph respiration flux [mmol13C m-3 day-1]
            !   calc_diss_13           : Calcite dissolution rate [day-1]
            !   DetCalc_13             : Detrital calcite 13C [mmol13C m-3]
            !   calc_loss_gra_13       : Calcite grazing flux [mmol13C m-3 day-1]
            !   calc_diss_guts         : Gut dissolution fraction [-]
            !   calcification_13       : Calcification flux [mmol13C m-3 day-1]
            !   dt_b                   : Biogeochemistry time step [day]
            !
            ! Note: Photosynthesis preferentially takes up 12C, leaving DIC enriched in 13C
            !       Respiration returns carbon with original isotopic composition
            !---------------------------------------------------------------------------

            sms(k, idic_13) = ( &
            !-----------------------------------------------------------------------
            ! SINKS: Carbon fixation (removes 13C-DIC)
            !-----------------------------------------------------------------------
            ! Fixation flux uses the isotope ratio of the source DIC pool
            ! (r_phyc_13/r_diac_13) applied to the bulk photosynthesis flux, rather
            ! than the phyto 13C pool itself -- corr fix eb707e35 "changes wrt. C
            ! isotopes"
                    -Cphot * r_phyc_13 * PhyC & ! Small phyto photosynthesis
                    - Cphot_Dia * r_diac_13 * DiaC & ! Diatom photosynthesis
            !
            !-----------------------------------------------------------------------
            ! SOURCES: Respiration and remineralization (returns 13C to DIC)
            !-----------------------------------------------------------------------
                    + phyRespRate * PhyC_13 & ! Small phyto respiration
                    + phyRespRate_Dia * DiaC_13 & ! Diatom respiration
                    + rho_C1 * arrFunc * EOC_13 & ! DOC remineralization
                    + HetRespFlux_13 & ! Heterotroph respiration
            !
            !-----------------------------------------------------------------------
            ! SOURCES: Calcite dissolution (releases 13C from CaCO3)
            !-----------------------------------------------------------------------
                    + calc_diss_13 * DetCalc_13 & ! Water column dissolution
                    + calc_loss_gra_13 * calc_diss_guts & ! Gut dissolution
            !
            !-----------------------------------------------------------------------
            ! SINKS: Calcification (removes 13C for CaCO3 formation)
            !-----------------------------------------------------------------------
                    - calcification_13 & ! CaCO3 precipitation
            !
                    ) * dt_b + sms(k, idic_13)

            !===========================================================================
            ! SMALL PHYTOPLANKTON ORGANIC CARBON (PhyC_13)
            !===========================================================================
            ! 13C budget for small phytoplankton biomass.
            !
            ! SOURCES (+):
            !   - Photosynthesis (fixes 13C-DIC into biomass)
            !
            ! SINKS (-):
            !   - Nutrient-stress mortality (lysis)
            !   - Respiration (maintenance costs)
            !   - Aggregation (particle formation)
            !   - Grazing by zooplankton
            !
            ! Variables:
            !   lossC              : Mortality rate constant [day-1]
            !   limitFacN          : Nutrient limitation factor [0-1]
            !   aggregationRate    : Aggregation rate [day-1]
            !   grazingFlux_phy    : Grazing flux on small phyto [mmolN m-3 day-1]
            !   recipQuota_13      : 13C:N ratio [mmol13C mmolN-1]
            !
            ! Note: Grazing uses recipQuota_13 to convert N-based flux to 13C flux
            !---------------------------------------------------------------------------

            sms(k, iphyc_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Photosynthetic production
            !-----------------------------------------------------------------------
            ! corr fix eb707e35: fixation uses r_phyc_13 * PhyC (source-pool ratio),
            ! not Cphot * PhyC_13
                    +Cphot * r_phyc_13 * PhyC & ! 13C fixation
            !
            !-----------------------------------------------------------------------
            ! SINKS: Mortality, respiration, and losses
            !-----------------------------------------------------------------------
                    - lossC * limitFacN * PhyC_13 & ! Nutrient-stress mortality
                    - phyRespRate * PhyC_13 & ! Respiration
                    - aggregationRate * PhyC_13 & ! Aggregation loss
                    - grazingFlux_phy * recipQuota_13 & ! Grazing loss (N->C conversion)
            !
                    ) * dt_b + sms(k, iphyc_13)

            !===========================================================================
            ! DETRITAL ORGANIC CARBON (DetC_13)
            !===========================================================================
            ! 13C budget for dead organic matter (detritus pool).
            !
            ! SOURCES (+):
            !   - Unassimilated grazing (sloppy feeding + egestion)
            !   - Phytoplankton aggregation
            !   - Diatom aggregation
            !   - Heterotroph mortality
            !
            ! SINKS (-):
            !   - Remineralization (aerobic respiration)
            !   - Assimilated grazing (efficient consumption)
            !
            ! Variables:
            !   grazEff       : Grazing efficiency (fraction assimilated) [-]
            !   hetLossFlux   : Heterotroph mortality flux [mmolN m-3 day-1]
            !   recipQZoo_13  : Heterotroph 13C:N ratio [mmol13C mmolN-1]
            !   reminC        : Detritus remineralization rate [day-1]
            !
            ! Grazing Partitioning:
            !   Total ingestion = Assimilated + Unassimilated
            !   Assimilated: Goes to heterotroph biomass (grazEff × flux)
            !   Unassimilated: Goes to detritus ((1-grazEff) × flux)
            !---------------------------------------------------------------------------

            sms(k, idetc_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Unassimilated grazing (sloppy feeding)
            !-----------------------------------------------------------------------
                    +grazingFlux_phy * recipQuota_13 & ! Total small phyto grazing
            ! Minus assimilated portion
                    - grazingFlux_phy * recipQuota_13 * grazEff &
                    + grazingFlux_Dia * recipQuota_dia_13 & ! Total diatom grazing
            ! Minus assimilated portion
                    - grazingFlux_Dia * recipQuota_dia_13 * grazEff &
            !
            !-----------------------------------------------------------------------
            ! SOURCES: Aggregation and mortality
            !-----------------------------------------------------------------------
                    + aggregationRate * phyC_13 & ! Small phyto aggregation
                    + aggregationRate * DiaC_13 & ! Diatom aggregation
                    + hetLossFlux * recipQZoo_13 & ! Heterotroph mortality
            !
            !-----------------------------------------------------------------------
            ! SINKS: Remineralization
            !-----------------------------------------------------------------------
                    - reminC * arrFunc * DetC_13 & ! Aerobic respiration
            !
                    ) * dt_b + sms(k, idetc_13)

            !===========================================================================
            ! HETEROTROPH ORGANIC CARBON (HetC_13)
            !===========================================================================
            ! 13C budget for zooplankton biomass.
            !
            ! SOURCES (+):
            !   - Assimilated grazing on phytoplankton
            !   - Assimilated grazing on diatoms
            !
            ! SINKS (-):
            !   - Mortality (density-dependent)
            !   - Non-predatory losses (diseases, senescence)
            !   - Respiration (metabolic costs)
            !
            ! Variables:
            !   lossC_z       : Non-predatory loss rate [day-1]
            !   hetRespFlux_13: Heterotroph respiration flux [mmol13C m-3 day-1]
            !
            ! Note: Grazing efficiency (grazEff) determines assimilation fraction
            !       Typical values: 0.6-0.8 (60-80% assimilated)
            !---------------------------------------------------------------------------

            sms(k, ihetc_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Assimilated food
            !-----------------------------------------------------------------------
                    +grazingFlux_phy * recipQuota_13 * grazEff & ! Small phyto consumption
                    + grazingFlux_Dia * recipQuota_dia_13 * grazEff & ! Diatom consumption
            !
            !-----------------------------------------------------------------------
            ! SINKS: Mortality and respiration
            !-----------------------------------------------------------------------
                    - hetLossFlux * recipQZoo_13 & ! Mortality flux
                    - lossC_z * HetC_13 & ! Non-predatory losses
                    - hetRespFlux_13 & ! Respiration
            !
                    ) * dt_b + sms(k, ihetc_13)

            !===========================================================================
            ! DISSOLVED ORGANIC CARBON (EOC_13)
            !===========================================================================
            ! 13C budget for dissolved organic carbon pool.
            !
            ! SOURCES (+):
            !   - Phytoplankton exudation (nutrient-stress losses)
            !   - Diatom exudation
            !   - Detritus remineralization (solubilization)
            !   - Heterotroph exudation (sloppy feeding, excretion)
            !   - River input (terrestrial DOC)
            !
            ! SINKS (-):
            !   - Remineralization (microbial respiration)
            !
            ! Variables:
            !   lossC, lossC_d    : Exudation rate constants [day-1]
            !   limitFacN         : Nutrient limitation factors [0-1]
            !   LocRiverDOC       : River DOC input flux [mmolC m-3 day-1]
            !   r_iorg_13         : River 13C:12C ratio (isotopic signature) [-]
            !
            ! DOC Pool Characteristics:
            !   - Labile fraction: Days to weeks turnover
            !   - Semi-labile fraction: Months to years turnover
            !   - Model uses bulk DOC with single remineralization rate
            !
            ! River Isotope Signature:
            !   - Terrestrial organic matter typically depleted in 13C
            !   - δ13C ≈ -27‰ for C3 plants, -13‰ for C4 plants
            !   - Marine phytoplankton: δ13C ≈ -20 to -22‰
            !---------------------------------------------------------------------------

            sms(k, idoc_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Exudation and solubilization
            !-----------------------------------------------------------------------
                    +lossC * limitFacN * phyC_13 & ! Small phyto exudation
                    + lossC_d * limitFacN_dia * DiaC_13 & ! Diatom exudation
                    + reminC * arrFunc * DetC_13 & ! Detritus solubilization
                    + lossC_z * HetC_13 & ! Heterotroph exudation
            !
            !-----------------------------------------------------------------------
            ! SOURCES: River input (terrestrial DOC)
            !-----------------------------------------------------------------------
                    + LocRiverDOC * r_iorg_13 & ! River 13C input
            !
            !-----------------------------------------------------------------------
            ! SINKS: Remineralization
            !-----------------------------------------------------------------------
                    - rho_c1 * arrFunc * EOC_13 & ! Microbial respiration
            !
                    ) * dt_b + sms(k, idoc_13)

            !===========================================================================
            ! DIATOM ORGANIC CARBON (DiaC_13)
            !===========================================================================
            ! 13C budget for diatom biomass (large phytoplankton).
            !
            ! Structure identical to small phytoplankton (section 1.2)
            ! with diatom-specific parameters and fluxes.
            !---------------------------------------------------------------------------

            sms(k, idiac_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Photosynthetic production
            !-----------------------------------------------------------------------
            ! corr fix eb707e35: fixation uses r_diac_13 * DiaC (source-pool ratio),
            ! not Cphot_dia * DiaC_13
                    +Cphot_dia * r_diac_13 * DiaC & ! 13C fixation
            !
            !-----------------------------------------------------------------------
            ! SINKS: Mortality, respiration, and losses
            !-----------------------------------------------------------------------
                    - lossC_d * limitFacN_dia * DiaC_13 & ! Nutrient-stress mortality
                    - phyRespRate_dia * DiaC_13 & ! Respiration
                    - aggregationRate * DiaC_13 & ! Aggregation loss
                    - grazingFlux_dia * recipQuota_dia_13 & ! Grazing loss
            !
                    ) * dt_b + sms(k, idiac_13)

            !===========================================================================
            ! PHYTOPLANKTON CALCITE (PhyCalc_13)
            !===========================================================================
            ! 13C budget for calcium carbonate associated with living phytoplankton.
            !
            ! SOURCES (+):
            !   - Calcification (CaCO3 precipitation on cells)
            !
            ! SINKS (-):
            !   - Nutrient-stress mortality (CaCO3 to detritus)
            !   - Cell death respiration (CaCO3 to detritus)
            !   - Aggregation (CaCO3 incorporated in aggregates)
            !   - Grazing (CaCO3 consumed with cells)
            !
            ! Variables:
            !   calcification_13  : 13C calcification flux [mmol13C m-3 day-1]
            !   phyCalc_13        : Phytoplankton calcite 13C [mmol13C m-3]
            !   calc_loss_agg_13  : Aggregation loss flux [mmol13C m-3 day-1]
            !   calc_loss_gra_13  : Grazing loss flux [mmol13C m-3 day-1]
            !
            ! Calcite Isotope Fractionation:
            !   - Small enrichment in 13C relative to DIC (~+1‰)
            !   - Temperature-dependent fractionation
            !   - Important for paleoclimate proxies (foraminifera, coccoliths)
            !---------------------------------------------------------------------------

            sms(k, iphycal_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Calcification
            !-----------------------------------------------------------------------
                    +calcification_13 & ! CaCO3 precipitation
            !
            !-----------------------------------------------------------------------
            ! SINKS: Mortality and losses
            !-----------------------------------------------------------------------
                    - lossC * limitFacN * phyCalc_13 & ! Mortality to detritus
                    - phyRespRate * phyCalc_13 & ! Death to detritus
                    - calc_loss_agg_13 & ! Aggregation
                    - calc_loss_gra_13 & ! Grazing
            !
                    ) * dt_b + sms(k, iphycal_13)

            !===========================================================================
            ! DETRITAL CALCITE (DetCalc_13)
            !===========================================================================
            ! 13C budget for calcium carbonate in detritus/particles.
            !
            ! SOURCES (+):
            !   - Phytoplankton mortality (CaCO3 from dead cells)
            !   - Cell death (respiratory loss to detritus)
            !   - Aggregation (CaCO3 in aggregates)
            !   - Grazing (CaCO3 in fecal pellets)
            !
            ! SINKS (-):
            !   - Dissolution in water column
            !   - Dissolution in grazer guts
            !
            ! Variables:
            !   calc_diss_guts : Fraction dissolved in guts (typically 0.1-0.5) [-]
            !
            ! Calcite Dissolution:
            !   - Thermodynamically driven (saturation state dependent)
            !   - Faster in undersaturated waters (deep ocean, high CO2)
            !   - Gut dissolution: Acidic environment accelerates dissolution
            !   - Returns 13C to DIC pool (source for DIC_13 budget)
            !---------------------------------------------------------------------------

            sms(k, idetcal_13) = ( &
            !-----------------------------------------------------------------------
            ! SOURCES: Mortality and transfer from living cells
            !-----------------------------------------------------------------------
                    +lossC * limitFacN * phyCalc_13 & ! Mortality
                    + phyRespRate * phyCalc_13 & ! Death
                    + calc_loss_agg_13 & ! Aggregation
                    + calc_loss_gra_13 & ! Grazing (to fecal pellets)
            !
            !-----------------------------------------------------------------------
            ! SINKS: Dissolution
            !-----------------------------------------------------------------------
                    - calc_loss_gra_13 * calc_diss_guts & ! Gut dissolution
                    - calc_diss_13 * DetCalc_13 & ! Water column dissolution
            !
                    ) * dt_b + sms(k, idetcal_13)

            if (ciso_14) then

                if (ciso_organic_14) then

                    !===================================================================
                    ! CARBON-14 (14C) BUDGETS
                    !===================================================================
                    ! Calculates 14C budgets for radiocarbon applications.
                    ! Structure identical to 13C budgets (sections 1.1-1.8).
                    !
                    ! Additional Consideration:
                    !   - Radioactive decay (half-life 5,730 years)
                    !   - Decay term handled separately in forcing module
                    !   - Bomb radiocarbon (anthropogenic 14C pulse)
                    !
                    ! Applications:
                    !   - Ocean ventilation age (Δ14C)
                    !   - Carbon residence time
                    !   - Mixing timescales
                    !   - Model validation (WOCE/CLIVAR 14C data)
                    !
                    ! Notation:
                    !   All variables end in _14 (e.g., PhyC_14, DiaC_14)
                    !   Structure parallels 13C budgets exactly
                    !-------------------------------------------------------------------

                    !===================================================================
                    ! DIC_14
                    !===================================================================
                    sms(k, idic_14) = ( &
                    ! corr fix eb707e35: r_phyc_14/r_diac_14 * bulk C, not
                    ! Cphot * PhyC_14
                            -Cphot * r_phyc_14 * PhyC &
                            + phyRespRate * PhyC_14 &
                            - Cphot_Dia * r_diac_14 * DiaC &
                            + phyRespRate_Dia * DiaC_14 &
                            + rho_C1 * arrFunc * EOC_14 &
                            + HetRespFlux_14 &
                            + calc_diss_14 * DetCalc_14 &
                            + calc_loss_gra_14 * calc_diss_guts &
                            - calcification_14 &
                            ) * dt_b + sms(k, idic_14)

                    !===================================================================
                    ! PhyC_14
                    !===================================================================
                    sms(k, iphyc_14) = ( &
                    ! corr fix eb707e35: r_phyc_14 * PhyC (source-pool ratio)
                            +Cphot * r_phyc_14 * PhyC &
                            - lossC * limitFacN * PhyC_14 &
                            - phyRespRate * PhyC_14 &
                            - aggregationRate * PhyC_14 &
                            - grazingFlux_phy * recipQuota_14 &
                            ) * dt_b + sms(k, iphyc_14)

                    !===================================================================
                    ! DetC_14
                    !===================================================================
                    sms(k, idetc_14) = ( &
                            +grazingFlux_phy * recipQuota_14 &
                            - grazingFlux_phy * recipQuota_14 * grazEff &
                            + grazingFlux_Dia * recipQuota_dia_14 &
                            - grazingFlux_Dia * recipQuota_dia_14 * grazEff &
                            + aggregationRate * phyC_14 &
                            + aggregationRate * DiaC_14 &
                            + hetLossFlux * recipQZoo_14 &
                            - reminC * arrFunc * DetC_14 &
                            ) * dt_b + sms(k, idetc_14)

                    !===================================================================
                    ! HetC_14
                    !===================================================================
                    sms(k, ihetc_14) = ( &
                            +grazingFlux_phy * recipQuota_14 * grazEff &
                            + grazingFlux_Dia * recipQuota_dia_14 * grazEff &
                            - hetLossFlux * recipQZoo_14 &
                            - lossC_z * HetC_14 &
                            - hetRespFlux_14 &
                            ) * dt_b + sms(k, ihetc_14)

                    !===================================================================
                    ! EOC_14
                    !===================================================================
                    sms(k, idoc_14) = ( &
                            +lossC * limitFacN * phyC_14 &
                            + lossC_d * limitFacN_dia * DiaC_14 &
                            + reminC * arrFunc * DetC_14 &
                            + lossC_z * HetC_14 &
                            - rho_c1 * arrFunc * EOC_14 &
                            + LocRiverDOC * r_iorg_14 &
                            ) * dt_b + sms(k, idoc_14)

                    !===================================================================
                    ! DiaC_14
                    !===================================================================
                    sms(k, idiac_14) = ( &
                    ! corr fix eb707e35: r_diac_14 * DiaC (source-pool ratio)
                            +Cphot_dia * r_diac_14 * DiaC &
                            - lossC_d * limitFacN_dia * DiaC_14 &
                            - phyRespRate_dia * DiaC_14 &
                            - aggregationRate * DiaC_14 &
                            - grazingFlux_dia * recipQuota_dia_14 &
                            ) * dt_b + sms(k, idiac_14)

                    !===================================================================
                    ! PhyCalc_14
                    !===================================================================
                    sms(k, iphycal_14) = ( &
                            +calcification_14 &
                            - lossC * limitFacN * phyCalc_14 &
                            - phyRespRate * phyCalc_14 &
                            - calc_loss_agg_14 &
                            - calc_loss_gra_14 &
                            ) * dt_b + sms(k, iphycal_14)

                    !===================================================================
                    ! DetCalc_14
                    !===================================================================
                    sms(k, idetcal_14) = ( &
                            +lossC * limitFacN * phyCalc_14 &
                            + phyRespRate * phyCalc_14 &
                            + calc_loss_agg_14 &
                            + calc_loss_gra_14 &
                            - calc_loss_gra_14 * calc_diss_guts &
                            - calc_diss_14 * DetCalc_14 &
                            ) * dt_b + sms(k, idetcal_14)

                else

                    !===================================================================
                    ! ABIOTIC DIC_14 (SIMPLIFIED MODE)
                    !===================================================================
                    ! "Abiotic" 14C tracking without explicit organic pools
                    ! DIC_14 tracks total carbon with radioactive decay only
                    !
                    ! Use Case:
                    !   - Simplified radiocarbon tracer
                    !   - Tracks ventilation/mixing without biology
                    !   - Computationally efficient
                    !   - Decay handled in forcing module (recom_forcing)
                    !
                    ! Equation:
                    !   DIC_14 changes identically to DIC (conservative tracer)
                    !   Plus: Radioactive decay (handled separately)
                    !
                    ! Limitation:
                    !   - No isotope fractionation during biological processes
                    !   - Cannot capture biological isotope signals
                    !   - Suitable only for physical circulation studies
                    !-------------------------------------------------------------------

                    sms(k, idic_14) = sms(k, idic)

                end if ! ciso_organic_14

            end if ! ciso_14

        end if ! ciso

    end subroutine sms_update_state
end module recom_sms_update
