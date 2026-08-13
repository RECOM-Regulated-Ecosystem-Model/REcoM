module recom_forcing_module
    implicit none
    private

    public :: recom_forcing

contains

    !===============================================================================
    ! REcoM_Forcing
    !
    !   Cavity convention (nzmin > 1):
    !     CO2 and O2 surface fluxes are set to zero; all other BGC progresses
    !     normally through the water column below the ice shelf.
    !===============================================================================
    subroutine REcoM_Forcing(n, nzmin, Nn, state, SurfSW, Loc_slp, Temp, Sali, Sali_depth, &
            CO2_watercolumn, pH_watercolumn, pCO2_watercolumn, &
            HCO3_watercolumn, CO3_watercolumn, OmegaC_watercolumn, &
            kspc_watercolumn, rhoSW_watercolumn, PAR, MPI_COMM_FESOM, &
            mype, myDim_nod2D, eDim_nod2D, nl, hnode, zbar_3d_n, &
            geo_coord_nod2D, daynew, ndpyr, dt, kappa, mstep, rad)

        use recom_declarations, only: wp, tiny_si, tiny_n, tiny_c, tiny_n_d, tiny_n_c, tiny_n_p, &
                tiny_c_p, tiny_c_d, tiny_c_c, vertNPPn, vertGPPn, locGPPn, locNPPn, vertNNAn, &
                vertChldegn, vertNPPd, locNPPd, locchldegn, vertGPPd, locGPPd, vertNNAd, locNNAd, &
                vertChldegd, locChldegd, vertNPPc, vertGPPc, locNPPc, locGPPc, vertNNAc, &
                locNNAc, locChldegc, vertNPPp, locNPPp, vertGPPp, vertNNAp, vertchldegp, locNNAp, &
                locGPPp, locChldegp, vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d, &
                vertgrazmeso_p, vertgrazmeso_det, vertgrazmeso_mic, vertgrazmeso_det2, &
                vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d, vertgrazmacro_c, &
                vertgrazmacro_mes, vertgrazmacro_det, vertgrazmacro_mic, vertgrazmacro_det2, &
                vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d, vertgrazmicro_c, &
                locgrazmacro_c, locgrazmacro_d, locgrazmacro_det, &
                locgrazmacro_mes, locgrazmacro_mic, locgrazmacro_n, locgrazmacro_p, &
                locgrazmeso_c, locgrazmeso_d, locgrazmeso_det, locgrazmeso_det2, locgrazmeso_mic, &
                locgrazmeso_n, locgrazmeso_p, locgrazmeso_tot, locgrazmicro_c, locgrazmicro_d, &
                locgrazmicro_n, locgrazmicro_p, locgrazmicro_tot, vertgrazmicro_p, &
                locNNAn, vertChldegc, vertgrazmeso_c, locgrazmacro_tot, locgrazmacro_det2, &
                vertgrazmacro_p

        use recom_config, only: bgc_num, chl2n_max, chl2n_max_c, chl2n_max_d, chl2n_max_p, ciso, &
                diags, enable_3zoo2det, enable_coccos, grazing_detritus, ialk, icchl, icocc, &
                idchl, idiac, idian, idiasi, idic, idin, imiczooc, imiczoon, ioxy, idicremin, ipchl, &
                iphac, iphachl, iphan, iphyc, iphyn, isi, ncmax, ncmax_c, ncmax_d, ncmax_p, nmocsy, &
                one, pa2atm, recom_debug, secondsperday, sicmax, tiny, tiny_chl, icocn

        use recom_ciso, only: alpha_aq_13, alpha_aq_14, alpha_dic_13, alpha_dic_14, alpha_k_13, &
                alpha_k_14, alpha_p_13, alpha_p_14, alpha_p_dia_13, alpha_p_dia_14, ciso_14, &
                ciso_organic_14, co2flux_13, co2flux_14, co2flux_seaicemask_13, &
                co2sat, idiac_13, idiac_14, idic_13, idic_14, iphyc_13, iphyc_14, kwco2, r_atm_13, &
                r_atm_14, r_co2s_13, r_co2s_14, r_diac_13, r_diac_14, r_dic_13, r_dic_14, &
                r_phyc_13, r_phyc_14, co2flux_seaicemask_14, recom_ciso_airsea, recom_ciso_photo

        use recom_locvar, only: betad, co2, co2ex, dpco2surf, fco2, hco3, k0, kw660, loc_ice_conc, &
                locatmco2, o2ex, o2flux_seaicemask, oflux, omegaa, omegac, p, pco2surf, ph, rhosw, &
                tempis, uloc, co2flux, co3, co2flux_seaicemask, dflux

        use recom_extra, only: cobeta, depth_calculations
        use recom_sms_module, only: recom_sms
        use gasx, only: scco2, flxco2, pistonvel, o2flux

        implicit none

        integer, intent(in) :: daynew, ndpyr, mype, myDim_nod2D, eDim_nod2D, nl, mstep

        ! Node index, top active level, bottom active level
        integer, intent(in) :: n, nzmin, Nn
        integer, intent(in) :: MPI_COMM_FESOM

        real(kind=wp), intent(in) :: rad
        real(kind=wp), intent(in) :: Sali ! Salinity of current surface layer
        real(kind=wp), intent(in) :: SurfSW ! Shortwave at surface [W/m2]
        real(kind=wp), intent(in) :: Loc_slp ! Sea-level pressure [Pa]
        real(kind=wp), intent(inout) :: kappa, dt

        real(kind=wp), intent(in), dimension(nl - 1) :: Temp ! Potential temperature [deg C]
        real(kind=wp), intent(in), dimension(nl - 1) :: Sali_depth ! Salinity profile

        ! Carbonate chemistry profiles (mocsy + DISS), updated in place
        real(kind=wp), intent(inout), dimension(nl - 1) :: CO2_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: pH_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: pCO2_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: HCO3_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: CO3_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: OmegaC_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: kspc_watercolumn
        real(kind=wp), intent(inout), dimension(nl - 1) :: rhoSW_watercolumn

        ! PAR profile [W/m2], computed here and returned to caller
        real(kind=wp), intent(inout), dimension(nl - 1) :: PAR

        real(kind=WP), intent(in), dimension(:, :) :: hnode, zbar_3d_n
        real(kind=WP), intent(in), dimension(:, :) :: geo_coord_nod2D
        real(kind=wp), intent(inout), dimension(nl - 1, bgc_num) :: state

        !---------------------------------------------------------------------------
        ! Mocsy / CO2-flux inputs (size-1 arrays - mocsy uses array interfaces)
        !   Note: Sali shadows the dummy argument Sali_depth for the surface
        !   value only, avoiding a name collision with the profile array.
        !---------------------------------------------------------------------------
        real(kind=wp) :: REcoM_DIC(1)  ! DIC at surface [mol/m3]
        real(kind=wp) :: REcoM_Alk(1)  ! Alkalinity at surface [mol/m3]
        real(kind=wp) :: REcoM_Si(1)   ! Silicate at surface [mol/m3]
        real(kind=wp) :: REcoM_Phos(1) ! Phosphate at surface [mol/m3] (from DIN/Redfield)
        real(kind=wp) :: Latd(1)       ! latitude [degrees]
        real(kind=wp) :: Lond(1)       ! longitude [degrees]
        real(kind=wp) :: REcoM_T(1)    ! Surface temperature, clamped for Lueker K1/K2 [deg C]
        real(kind=wp) :: REcoM_S(1)    ! clamped for Lueker K1/K2
        real(kind=wp) :: Patm(1)       ! Atmospheric pressure [atm]
        real(kind=wp) :: Latr          ! Latitude [radians], scalar scratch

        !---------------------------------------------------------------------------
        ! O2-flux inputs
        !---------------------------------------------------------------------------
        real(kind=wp) :: ppo(1)        ! Sea-level pressure normalised to 1 atm [-]
        real(kind=wp) :: REcoM_O2(1)   ! Surface O2 [mol/m3]

        !---------------------------------------------------------------------------
        ! SMS tendency [mmol/m3/timestep]
        !---------------------------------------------------------------------------
        real(kind=wp), dimension(nl - 1, bgc_num) :: sms

        !---------------------------------------------------------------------------
        ! Local geometry
        !---------------------------------------------------------------------------
        real(kind=wp), dimension(nl)     :: zF         ! Depth at flux points [m]
        real(kind=wp), dimension(nl, 6)  :: SinkVel    ! Sinking velocities [m/d]
        real(kind=wp), dimension(nl - 1) :: thick      ! Layer thickness [m]
        real(kind=wp), dimension(nl - 1) :: recipthick ! 1/thick [1/m]

        !---------------------------------------------------------------------------
        ! Threshold concentrations derived from config parameters
        !   tiny_X = minimum meaningful concentration for species X, used as lower
        !   bound after the SMS update and as a near-zero guard in ratio calculations.
        !---------------------------------------------------------------------------
        tiny_N   = tiny_chl / chl2N_max    ! [mmol N/m3]  small phytoplankton N
        tiny_N_d = tiny_chl / chl2N_max_d  ! [mmol N/m3]  diatom N
        tiny_C   = tiny_N   / NCmax        ! [mmol C/m3]  small phytoplankton C
        tiny_C_d = tiny_N_d / NCmax_d      ! [mmol C/m3]  diatom C
        tiny_Si  = tiny_C_d / SiCmax       ! [mmol Si/m3] diatom Si

        if (enable_coccos) then
            tiny_N_c = tiny_chl / chl2N_max_c  ! [mmol N/m3]  coccolithophore N
            tiny_C_c = tiny_N_c / NCmax_c      ! [mmol C/m3]  coccolithophore C
            tiny_N_p = tiny_chl / chl2N_max_p  ! [mmol N/m3]  phaeocystis N
            tiny_C_p = tiny_N_p / NCmax_p      ! [mmol C/m3]  phaeocystis C
        end if

        !---------------------------------------------------------------------------
        ! Grid geometry for this column
        !---------------------------------------------------------------------------
        call Cobeta(daynew, ndpyr, myDim_nod2D, geo_coord_nod2D)
        call Depth_calculations(n, nzmin, Nn, SinkVel, zF, thick, recipthick, nl, hnode, zbar_3d_n)

        !===========================================================================
        ! Surface gas exchange
        !   Cavity nodes (nzmin > 1): fluxes are forced to zero - there is no
        !   atmosphere - ocean interface below an ice shelf.
        !===========================================================================

        !---------------------------------------------------------------------------
        ! Prepare surface-layer scalars for mocsy
        !   Lueker (2000) K1/K2 valid range: T in [2, 35] deg C, S in [19, 43].
        !   We clamp to [2, 40] deg C and [21, 43] to avoid extrapolation errors,
        !   including in near-freezing, low-salinity, high-ice-cover cells.
        !---------------------------------------------------------------------------
        REcoM_DIC  = max(tiny * 1e-3, state(one, idic) * 1e-3) ! mmol/m3 -> mol/m3
        REcoM_Alk  = max(tiny * 1e-3, state(one, ialk) * 1e-3)
        REcoM_Si   = max(tiny * 1e-3, state(one, isi)  * 1e-3)
        REcoM_Phos = max(tiny * 1e-3, state(one, idin) * 1e-3) / 16.d0  ! N->P Redfield
        REcoM_O2   = max(tiny * 1e-3, state(one, ioxy) * 1e-3)

        !!---- minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
        REcoM_T = max(2.d0, Temp(nzmin))
        !!---- maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
        REcoM_T = min(40.d0, REcoM_T(1))

        !!---- minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes
        !!trouble in regions with S between 19 and 21 and ice conc above 97%
        REcoM_S = max(21.d0, Sali)
        !!---- maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes
        !!trouble   REcoM_S    = min(REcoM_S, 43.d0)  !!!!!!!!

        Patm = Loc_slp / Pa2atm ! Pa -> atm

        !!---- lon
        Lond = geo_coord_nod2D(1, n) / rad ! rad -> degrees
        !!---- lat
        Latr = geo_coord_nod2D(2, n)       ! radians (scratch)
        Latd = geo_coord_nod2D(2, n) / rad ! rad -> degrees

        !!---- calculate piston velocity kw660, which is an input to the flxco2 calculation
        !!---- pistonvel already scaled for ice-free area
        !!---- compute piston velolicty kw660 (at 25 C) from wind speed
        !!---- BUT without Schmidt number temperature correction (Sc differs each gas)
        !! ULoc: wind speed at 10-m height
        !! Loc_ice_conc: modeled sea-ice cover: fraction of grid cell, varying between 0.0 (no ice)
        !! and 1.0 (full cover)
        !! kw660: piston velocity at 25°C [m/s], uncorrected by the Schmidt number for different
        !! temperatures

        !---------------------------------------------------------------------------
        ! Piston velocity kw660 at 25 deg C [m/s]
        !   Computed unconditionally; set to zero inside REcoM_sms for cavities.
        !   Ice-cover fraction is already folded in by pistonvel - do not reapply.
        !---------------------------------------------------------------------------

        call pistonvel(ULoc, Loc_ice_conc, Nmocsy, kw660)

        ! Guard against unphysical DIC (indicates upstream tracer corruption)
        if (REcoM_DIC(1) > 10.d0) then  ! > 10 mol/m3 = 10000 mmol/m3
            print *, 'FATAL: DIC out of range at n=', n
            print *, '  DIC   [mol/m3] =', REcoM_DIC
            print *, '  Alk   [mol/m3] =', REcoM_Alk
            print *, '  T     [deg C]  =', REcoM_T
            print *, '  S              =', REcoM_S
            print *, '  Si    [mol/m3] =', REcoM_Si
            print *, '  Phos  [mol/m3] =', REcoM_Phos
            print *, '  pCO2           =', pco2surf
            print *, '  rhoSW          =', rhoSW
            print *, '  tempis         =', tempis
            print *, '  kw660          =', kw660
            print *, '  AtmCO2         =', LocAtmCO2
            print *, '  Patm   [atm]   =', Patm
            print *, '  ULoc   [m/s]   =', ULoc
            print *, '  ice_conc       =', Loc_ice_conc
            print *, '  thick(nzmin)   =', thick(nzmin)
            print *, '  Nmocsy         =', Nmocsy
            print *, '  Lond, Latd     =', Lond, Latd
            stop
        end if

        if (nzmin == 1) then
            call flxco2(co2flux, co2ex, dpco2surf, ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, &
                    BetaD, rhoSW, p, tempis, K0, REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, &
                    REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Nmocsy, Lond, Latd, &
                    optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', &
                    optGAS='Pinsitu', optS='Sprc')

            ! Sanity check on computed flux magnitude

            if (abs(co2flux(1)) > 1.e10) then
                print *, 'FATAL: CO2 flux out of range at n=', n
                print *, '  co2flux        =', co2flux
                print *, '  pco2surf       =', pco2surf
                print *, '  co2            =', co2
                print *, '  rhoSW          =', rhoSW
                print *, '  T     [deg C]  =', REcoM_T
                print *, '  tempis         =', tempis
                print *, '  S              =', REcoM_S
                print *, '  Alk   [mol/m3] =', REcoM_Alk
                print *, '  DIC   [mol/m3] =', REcoM_DIC
                print *, '  Si    [mol/m3] =', REcoM_Si
                print *, '  Phos  [mol/m3] =', REcoM_Phos
                print *, '  kw660          =', kw660
                print *, '  AtmCO2         =', LocAtmCO2
                print *, '  Patm   [atm]   =', Patm
                print *, '  thick(nzmin)   =', thick(nzmin)
                print *, '  Nmocsy         =', Nmocsy
                print *, '  Lond, Latd     =', Lond, Latd
                print *, '  ULoc   [m/s]   =', ULoc
                print *, '  ice_conc       =', Loc_ice_conc
                stop
            end if

            ! Ice fraction already in piston velocity - do not re-apply here
            dflux              = co2flux * 1.e3 * SecondsPerDay  ! mol/m2/s -> mmol/m2/d
            co2flux_seaicemask = co2flux * 1.e3                  ! mol/m2/s -> mmol/m2/s
        else  ! cavity - no atmosphere above
            dflux              = 0.0_WP
            co2flux_seaicemask = 0.0_WP
            pco2surf = 0.0_WP
            dpco2surf = 0.0_WP
            ph = 0.0_WP
            kw660 = 0.0_WP
            K0 = 0.0_WP
            omegaC = 0.0_WP
            omegaA = 0.0_WP
        end if
 
        !---------------------------------------------------------------------------
        ! O2 air-sea flux (open ocean only)
        !---------------------------------------------------------------------------
        ppo = Loc_slp / Pa2atm  ! sea-level pressure normalised to 1 atm

        if (nzmin == 1) then
            call o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, Nmocsy, o2ex)
            oflux              = o2ex * 1.e3 * SecondsPerDay  ! mol/m2/s -> mmol/m2/d
            o2flux_seaicemask  = o2ex * 1.e3                  ! mol/m2/s -> mmol/m2/s
        else
            oflux              = 0.0_WP
            o2flux_seaicemask  = 0.0_WP
        end if

        !===========================================================================
        ! Water-column BGC <E2><80><94> sources minus sinks
        !===========================================================================
        if (recom_debug .and. mype == 0) print *, achar(27) // '[36m' // '     --> REcoM_sms' // &
                achar(27) // '[0m'

        call REcoM_sms(n, nzmin, Nn, state, thick, SurfSW, sms, Temp, Sali_depth, &
                CO2_watercolumn, & ! MOCSY [mol/m3]
                pH_watercolumn, & ! MOCSY on total scale
                pCO2_watercolumn, & ! MOCSY [uatm]
                HCO3_watercolumn, & ! MOCSY [mol/m3]
                CO3_watercolumn, & ! DISS [mol/m3]
                OmegaC_watercolumn, & ! DISS calcite saturation state
                kspc_watercolumn, & ! DISS stoichiometric solubility product [mol^2/kg^2]
                rhoSW_watercolumn, & ! DISS in-situ density of seawater [kg/m3]
                Loc_slp, zF, PAR, Latd, daynew, dt, kappa, mstep, MPI_COMM_FESOM, mype, &
                myDim_nod2D, eDim_nod2D, nl, geo_coord_nod2D)

        !---------------------------------------------------------------------------
        ! Apply SMS tendency and enforce lower concentration bounds
        !---------------------------------------------------------------------------
        state(nzmin:nn, :) = max(tiny, state(nzmin:nn, :) + sms(nzmin:nn, :))

        ! Per-species lower bounds (tighter than the generic 'tiny')
        state(nzmin:nn, ipchl)  = max(tiny_chl, state(nzmin:nn, ipchl))
        state(nzmin:nn, iphyn)  = max(tiny_N,   state(nzmin:nn, iphyn))
        state(nzmin:nn, iphyc)  = max(tiny_C,   state(nzmin:nn, iphyc))
        state(nzmin:nn, idchl)  = max(tiny_chl, state(nzmin:nn, idchl))
        state(nzmin:nn, idian)  = max(tiny_N_d, state(nzmin:nn, idian))
        state(nzmin:nn, idiac)  = max(tiny_C_d, state(nzmin:nn, idiac))
        state(nzmin:nn, idiasi) = max(tiny_Si,  state(nzmin:nn, idiasi))

        if (enable_coccos) then
            state(nzmin:nn, icchl) = max(tiny_chl, state(nzmin:nn, icchl))
            state(nzmin:nn, icocn) = max(tiny_N_c, state(nzmin:nn, icocn))
            state(nzmin:nn, icocc) = max(tiny_C_c, state(nzmin:nn, icocc))

            state(nzmin:nn, iphachl) = max(tiny_chl, state(nzmin:nn, iphachl))
            state(nzmin:nn, iphan)   = max(tiny_N_p, state(nzmin:nn, iphan))
            state(nzmin:nn, iphac)   = max(tiny_C_p, state(nzmin:nn, iphac))
        end if

        if (enable_3zoo2det) then
            state(nzmin:nn, imiczoon) = max(tiny, state(nzmin:nn, imiczoon))
            state(nzmin:nn, imiczooc) = max(tiny, state(nzmin:nn, imiczooc))
        end if

        state(nzmin:nn, idicremin) = max(tiny, state(nzmin:nn, idicremin))

        if (recom_debug .and. mype == 0) print *, achar(27) // '[36m' // '     --> ciso after' // &
                ' REcoM_Forcing' // achar(27) // '[0m'

        if (ciso) then
            if (nzmin == 1) then
            ! Calculate carbon-isotopic fractionation, radioactive decay is calculated in
            ! oce_ale_tracer.F90

            ! Fractionation due to air-sea exchange and chemical speciation of CO2
            ! -> alpha_aq, alpha_dic. CO3 is taken from mocsy
            call recom_ciso_airsea(recom_t(1), co3(1), recom_dic(1))

            ! Isotopic ratios of dissolved CO2, also needed to calculate biogenic fractionation
            r_dic_13 = max(tiny * 1e-3, state(1, idic_13) * 1e-3) / recom_dic(1)
            r_co2s_13 = alpha_aq_13 / alpha_dic_13 * r_dic_13
            ! Calculate air-sea fluxes of 13|14CO2 in mmol / m**2 / s
            kwco2 = kw660(1) * (660 / scco2(REcoM_T(1))) ** 0.5 ! Piston velocity (via mocsy)

            ! Saturation concentration of CO2 (via mocsy)
            co2sat = co2flux(1) / (kwco2 + tiny) + co2(1)
            ! co2flux_13   = kwco2 * alpha_k_13 * (alpha_aq_13 * r_atm_13 * co2sat - r_co2s_13 *
            ! co2(1))
            ! co2flux_13   = alpha_k_13 * alpha_aq_13 * kwco2 * (r_atm_13 * co2sat - r_dic_13 *
            ! co2(1) / alpha_dic_13)
            ! Fractionation factors were determined for freshwater, include a correction for
            ! enhanced fractionation in seawater
            co2flux_13 = (alpha_k_13 * alpha_aq_13 - 0.0002) * kwco2 &
                    * (r_atm_13 * co2sat - r_dic_13 * co2(1) / alpha_dic_13)
            co2flux_seaicemask_13 = co2flux_13 * 1.e3

            ! Biogenic fractionation due to photosynthesis of plankton
            ! phyc_13|14 and diac_13|14 are only used in REcoM_sms to calculate DIC_13|14,
            ! DOC_13|14 and DetC_13|14

            call recom_ciso_photo(co2(1)) ! -> alpha_p
            r_phyc_13 = r_co2s_13 / alpha_p_13
            r_diac_13 = r_co2s_13 / alpha_p_dia_13
            state(nzmin:nn, iphyc_13) = max((tiny_C * r_phyc_13), (state(nzmin:nn, iphyc) * r_phyc_13))
            state(nzmin:nn, idiac_13) = max((tiny_C_d * r_diac_13), (state(nzmin:nn, idiac) * r_diac_13))

            ! The same for radiocarbon, fractionation factors have been already derived above
            if (ciso_14) then
                ! Air-sea exchange
                r_dic_14 = max(tiny * 1e-3, state(1, idic_14) * 1e-3) / recom_dic(1)
                r_co2s_14 = alpha_aq_14 / alpha_dic_14 * r_dic_14
                ! co2flux_14 = kwco2 * alpha_k_14 * (alpha_aq_14 * r_atm_14 * co2sat - r_co2s_14
                ! * co2(1))
                ! Fractionation factors were determined for freshwater, include a correction for
                ! enhanced fractionation seawater
                co2flux_14 = (alpha_k_14 * alpha_aq_14 - 0.0004) * kwco2 &
                        * (r_atm_14 * co2sat - r_dic_14 * co2(1) / alpha_dic_14)
                co2flux_seaicemask_14 = co2flux_14 * 1.e3
                ! Biogenic fractionation
                if (ciso_organic_14) then
                    r_phyc_14 = r_co2s_14 / alpha_p_14
                    r_diac_14 = r_co2s_14 / alpha_p_dia_14
                    state(nzmin:nn, iphyc_14) = max((tiny_C * r_phyc_14), (state(nzmin:nn, iphyc) &
                            * r_phyc_14))
                    state(nzmin:nn, idiac_14) = max((tiny_C_d * r_diac_14), (state(nzmin:nn, idiac) &
                            * r_diac_14))
                end if
            end if
            else
                co2flux_13 = 0.0_WP
                co2flux_seaicemask_13 = 0.0_WP
                if (ciso_14) then
                    co2flux_14 = 0.0_WP
                    co2flux_seaicemask_14 = 0.0_WP
                endif
            endif 
            ! Radiocarbon
        end if
        ! ciso

        !-------------------------------------------------------------------------------
        ! Diagnostics
        if (Diags) then

            ! logical, optional                 :: lNPPn

            ! if (present(lNPPn))then
            !     locNPPn = sum(diags3Dloc(nzmin:nn,idiags) * thick(nzmin:nn))
            ! endif
            locNPPn = sum(vertNPPn(nzmin:nn) * thick(nzmin:nn))
            locGPPn = sum(vertGPPn(nzmin:nn) * thick(nzmin:nn))
            locNNAn = sum(vertNNAn(nzmin:nn) * thick(nzmin:nn))
            locChldegn = sum(vertChldegn(nzmin:nn) * thick(nzmin:nn))

            locNPPd = sum(vertNPPd(nzmin:nn) * thick(nzmin:nn))
            locGPPd = sum(vertGPPd(nzmin:nn) * thick(nzmin:nn))
            locNNAd = sum(vertNNAd(nzmin:nn) * thick(nzmin:nn))
            locChldegd = sum(vertChldegd(nzmin:nn) * thick(nzmin:nn))

            if (enable_coccos) then
                locNPPc = sum(vertNPPc(nzmin:nn) * thick(nzmin:nn))
                locGPPc = sum(vertGPPc(nzmin:nn) * thick(nzmin:nn))
                locNNAc = sum(vertNNAc(nzmin:nn) * thick(nzmin:nn))
                locChldegc = sum(vertChldegc(nzmin:nn) * thick(nzmin:nn))

                locNPPp = sum(vertNPPp(nzmin:nn) * thick(nzmin:nn))
                locGPPp = sum(vertGPPp(nzmin:nn) * thick(nzmin:nn))
                locNNAp = sum(vertNNAp(nzmin:nn) * thick(nzmin:nn))
                locChldegp = sum(vertChldegp(nzmin:nn) * thick(nzmin:nn))
            end if

            ! only for the case if grazing detritus is used, as probably only
            ! needed for tuning which uses detritus grazing
            if (Grazing_detritus) then
                ! Mesozooplankton
                locgrazmeso_tot = sum(vertgrazmeso_tot(nzmin:nn) * thick(nzmin:nn))
                locgrazmeso_n = sum(vertgrazmeso_n(nzmin:nn) * thick(nzmin:nn))
                locgrazmeso_d = sum(vertgrazmeso_d(nzmin:nn) * thick(nzmin:nn))
                if (enable_coccos) then
                    locgrazmeso_c = sum(vertgrazmeso_c(nzmin:nn) * thick(nzmin:nn))
                    locgrazmeso_p = sum(vertgrazmeso_p(nzmin:nn) * thick(nzmin:nn))
                end if
                locgrazmeso_det = sum(vertgrazmeso_det(nzmin:nn) * thick(nzmin:nn))
                if (enable_3zoo2det) then
                    locgrazmeso_mic = sum(vertgrazmeso_mic(nzmin:nn) * thick(nzmin:nn))
                    locgrazmeso_det2 = sum(vertgrazmeso_det2(nzmin:nn) * thick(nzmin:nn))
                end if

                if (enable_3zoo2det) then
                    ! Macrozooplankton
                    locgrazmacro_tot = sum(vertgrazmacro_tot(nzmin:nn) * thick(nzmin:nn))
                    locgrazmacro_n = sum(vertgrazmacro_n(nzmin:nn) * thick(nzmin:nn))
                    locgrazmacro_d = sum(vertgrazmacro_d(nzmin:nn) * thick(nzmin:nn))
                    if (enable_coccos) then
                        locgrazmacro_c = sum(vertgrazmacro_c(nzmin:nn) * thick(nzmin:nn))
                        locgrazmacro_p = sum(vertgrazmacro_p(nzmin:nn) * thick(nzmin:nn))
                    end if
                    locgrazmacro_mes = sum(vertgrazmacro_mes(nzmin:nn) * thick(nzmin:nn))
                    locgrazmacro_det = sum(vertgrazmacro_det(nzmin:nn) * thick(nzmin:nn))
                    locgrazmacro_mic = sum(vertgrazmacro_mic(nzmin:nn) * thick(nzmin:nn))
                    locgrazmacro_det2 = sum(vertgrazmacro_det2(nzmin:nn) * thick(nzmin:nn))

                    ! Microzooplankton
                    locgrazmicro_tot = sum(vertgrazmicro_tot(nzmin:nn) * thick(nzmin:nn))
                    locgrazmicro_n = sum(vertgrazmicro_n(nzmin:nn) * thick(nzmin:nn))
                    locgrazmicro_d = sum(vertgrazmicro_d(nzmin:nn) * thick(nzmin:nn))
                    if (enable_coccos) then
                        locgrazmicro_c = sum(vertgrazmicro_c(nzmin:nn) * thick(nzmin:nn))
                        locgrazmicro_p = sum(vertgrazmicro_p(nzmin:nn) * thick(nzmin:nn))
                    end if

                end if
            end if

        end if
    end subroutine REcoM_Forcing

end module recom_forcing_module
