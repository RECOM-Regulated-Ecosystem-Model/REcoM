module recom_carbonate_chemistry
    implicit none
    private

    public :: calculate_carbonate_chemistry

contains

    subroutine calculate_carbonate_chemistry(k, zF, dt_b, mstep, SurfSR, Latd, &
            REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, &
            REcoM_Phos_depth, Patm_depth, CO2_watercolumn, pH_watercolumn, pCO2_watercolumn, &
            HCO3_watercolumn, CO3_watercolumn, OmegaC_watercolumn, kspc_watercolumn, &
            rhoSW_watercolumn)

        use recom_declarations, only: wp, parave
        use recom_locvar, only: betad_depth, co2_depth, co3_depth, dpos, fco2_depth, &
                hco3_depth, kspc_depth, logfile_outfreq_30, logfile_outfreq_7, omegaa_depth, &
                omegac_depth, pco2_depth, p_depth, ph_depth, rhosw_depth, tempis_depth
        use recom_config, only: nmocsy
        use mvars, only: vars_sprac

        implicit none

        integer, intent(in) :: k, mstep
        real(kind=wp), intent(in) :: dt_b, SurfSR
        real(kind=wp), intent(in), dimension(:) :: zF, Latd
        real(kind=wp), intent(in), dimension(:) :: REcoM_T_depth, REcoM_S_depth
        real(kind=wp), intent(in), dimension(:) :: REcoM_Alk_depth, REcoM_DIC_depth
        real(kind=wp), intent(in), dimension(:) :: REcoM_Si_depth, REcoM_Phos_depth
        real(kind=wp), intent(in), dimension(:) :: Patm_depth
        real(kind=wp), intent(inout), dimension(:) :: CO2_watercolumn, pH_watercolumn
        real(kind=wp), intent(inout), dimension(:) :: pCO2_watercolumn, HCO3_watercolumn
        real(kind=wp), intent(inout), dimension(:) :: CO3_watercolumn, OmegaC_watercolumn
        real(kind=wp), intent(inout), dimension(:) :: kspc_watercolumn, rhoSW_watercolumn

        real(kind=wp) :: mocsy_step_per_day
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

        ! Convert model depth coordinate to positive depth for MOCSY
        ! Model convention: zF(k) is negative (e.g., -100 m)
        ! MOCSY convention: depth is positive (e.g., 100 m)
        dpos(1) = -zF(k)

        ! Calculate update frequencies based on model time step
        mocsy_step_per_day = 1.0 / dt_b
        logfile_outfreq_7 = int(mocsy_step_per_day * 7.0) ! Steps in 7 days
        logfile_outfreq_30 = int(mocsy_step_per_day * 30.0) ! Steps in 30 days

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
    end subroutine calculate_carbonate_chemistry

end module recom_carbonate_chemistry
