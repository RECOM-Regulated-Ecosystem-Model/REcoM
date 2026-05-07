module recom_sms_init
    implicit none
    private
    public :: sms_initialize_variables

contains

    subroutine sms_initialize_variables(n, k, state, sms, thick, Temp, Sali_depth, SurfSR, PAR, &
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
                temp_phaeo, temp_phyto, tiny_c, tiny_c_c, tiny_c_d, tiny_c_p, tiny_n, tiny_n_c, &
                tiny_n_d, tiny_n_p, tiny_si, vttemp_phyto, vttemp_diatoms, vttemp_cocco, &
                vttemp_phaeo, temp_diatoms, chl_lower, chl_upper

        use recom_config, only: a_chl, ae, beta_phaeo, c2k, chl2n_max, chl2n_max_c, chl2n_max_d, &
                chl2n_max_p, enable_3zoo2det, enable_coccos, expon_cocco, expon_d, expon_phy, &
                grazing_detritus, ialk, icchl, icocc, icocn, idchl, idetc, idetcal, idetn, &
                idetz2c, idetz2calc, idetz2n, idetz2si, idiac, idian, idiasi, idic, idin, &
                idoc, idon, ife, ihetc, ihetn, imiczooc, imiczoon, ioxy, ipchl, iphac, iphachl, &
                iphan, iphyc, iphyn, isi, izoo2c, izoo2n, k_o2_remin, k_w, ncmax, ncmax_c, &
                ncmax_d, ncmax_p, ncmin, ncmin_c, ncmin_d, ncmin_p, o2dep_remin, one, ord_cocco, &
                recom_tref, res_het, sicmax, t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2, tiny, tiny_chl, &
                tiny_het, tmax_phaeo, topt_phaeo, uopt_phaeo, ciso, idetsi, reminsi, zero, &
                iphycal, ord_d, ord_phy

        use recoM_ciso, only: alpha_dcal_13, alpha_dcal_14, calc_diss_13, calc_diss_14, ciso_14, &
                ciso_organic_14, detc_13, detc_14, detcalc_13, detcalc_14, diac_13, diac_14, &
                dic_13, dic_14, eoc_13, eoc_14, hetc_13, hetc_14, idetc_13, idetc_14, idetcal_13, &
                idetcal_14, idiac_13, idiac_14, idic_13, idic_14, idoc_13, idoc_14, ihetc_13, &
                ihetc_14, iphyc_13, iphyc_14, iphycal_13, iphycal_14, phyc_13, phyc_14, &
                phycalc_13, phycalc_14, quota_13, quota_14, quota_dia_13, quota_dia_14, &
                recipquota_13, recipquota_14, recipquota_dia_13, recipquota_dia_14, &
                recipqzoo_13, recipqzoo_14

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
    end subroutine sms_initialize_variables
end module recom_sms_init
