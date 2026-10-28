# REcoM units reference

This file lists the physical units of the variables, parameters and fields used in REcoM.
Units were taken from the declarations and comments in the source and, where a comment was
missing or contradicted the code, worked out from how the variable is computed. Those cases
are listed in [Known inconsistencies in source comments](#known-inconsistencies-in-source-comments).

Notation: `mmol m-3` = millimole per cubic metre, `d` = day, `-` = dimensionless,
`mgChl` = milligram chlorophyll a. `ppm` = mole fraction x 10^6.

## Contents

1. [Conventions](#conventions)
2. [Prognostic BGC tracers](#prognostic-bgc-tracers)
3. [Carbon isotope tracers (ciso)](#carbon-isotope-tracers-ciso)
4. [Physical inputs and forcing](#physical-inputs-and-forcing)
5. [Carbonate chemistry (mocsy)](#carbonate-chemistry-mocsy)
6. [Air-sea gas exchange](#air-sea-gas-exchange)
7. [Atmospheric CO2 box model](#atmospheric-co2-box-model)
8. [Surface boundary inputs](#surface-boundary-inputs)
9. [Benthos and sediment](#benthos-and-sediment)
10. [Sinking and ballasting](#sinking-and-ballasting)
11. [Internal SMS variables](#internal-sms-variables)
12. [Namelist parameters](#namelist-parameters)
13. [Diagnostics](#diagnostics)
14. [Unit conversions in the code](#unit-conversions-in-the-code)
15. [Known inconsistencies in source comments](#known-inconsistencies-in-source-comments)

---

## Conventions

- **Concentrations** of all BGC tracers are in `mmol m-3`, except chlorophyll (`mgChl m-3`)
  and dissolved iron (`umol m-3`).
- **Biological rates** inside `REcoM_sms` are in `d-1`, because `dt_b` is in days.
- **The SMS array** (`sms`) holds the concentration change per REcoM sub-step, in
  `mmol m-3` per `dt_b`. Every update has the form `sms(k, i) = (...) * dt_b + sms(k, i)`.
- **Time steps.** The ocean time step `dt` is in `s`. Then `dt_d = dt / SecondsPerDay` [d] and
  `dt_b = dt_d / biostep` [d].
- **Fluxes exchanged with FESOM** are in `mmol m-2 s-1`. Fluxes written for output are mostly in
  `mmol m-2 d-1`.
- **Vertical coordinate.** `zF` and `zbar_3d_n` are in `m` and negative downward. mocsy receives
  positive depth `dpos = -zF(k)` [m].
- **Sign of air-sea fluxes.** Positive means into the ocean (downward).
- **mocsy** works in `mol m-3`, because it is called with `optCON='mol/m3'`. REcoM converts with
  a factor of `1e-3` before calling it.

## Prognostic BGC tracers

Tracer layout: `T, S | BGC (bgc_num) | [age] | [transit]`. BGC tracers start at slot 3.

| Index | Tracer ID | Variable | Description | Unit |
|---|---|---|---|---|
| `idin` | 1001 | `DIN` | Dissolved inorganic nitrogen | mmolN m-3 |
| `idic` | 1002 | `DIC` | Dissolved inorganic carbon | mmolC m-3 |
| `ialk` | 1003 | `Alk` | Total alkalinity | mmol m-3 (meq m-3) |
| `iphyn` | 1004 | `PhyN` | Small phytoplankton N | mmolN m-3 |
| `iphyc` | 1005 | `PhyC` | Small phytoplankton C | mmolC m-3 |
| `ipchl` | 1006 | `PhyChl` | Small phytoplankton chlorophyll | mgChl m-3 |
| `idetn` | 1007 | `DetN` | Detritus N (slow sinking) | mmolN m-3 |
| `idetc` | 1008 | `DetC` | Detritus C (slow sinking) | mmolC m-3 |
| `ihetn` | 1009 | `HetN` | Mesozooplankton (heterotrophs) N | mmolN m-3 |
| `ihetc` | 1010 | `HetC` | Mesozooplankton (heterotrophs) C | mmolC m-3 |
| `idon` | 1011 | `DON` | Dissolved organic N | mmolN m-3 |
| `idoc` | 1012 | `EOC` | Dissolved (extracellular) organic C | mmolC m-3 |
| `idian` | 1013 | `DiaN` | Diatom N | mmolN m-3 |
| `idiac` | 1014 | `DiaC` | Diatom C | mmolC m-3 |
| `idchl` | 1015 | `DiaChl` | Diatom chlorophyll | mgChl m-3 |
| `idiasi` | 1016 | `DiaSi` | Diatom silica (frustule) | mmolSi m-3 |
| `idetsi` | 1017 | `DetSi` | Detrital silica | mmolSi m-3 |
| `isi` | 1018 | `Si` | Dissolved silicate | mmolSi m-3 |
| `ife` | 1019 | `Fe` | Dissolved iron | **umolFe m-3** (= nmol L-1) |
| `iphycal` | 1020 | `PhyCalc` | Phytoplankton calcite | mmolC m-3 (as CaCO3) |
| `idetcal` | 1021 | `DetCalc` | Detrital calcite | mmolC m-3 (as CaCO3) |
| `ioxy` | 1022 | `O2` | Dissolved oxygen | mmolO2 m-3 |
| `idicremin` | 1037 | `DICremin` | DIC from remineralisation (diagnostic tracer) | mmolC m-3 |

Only present with `enable_3zoo2det = .true.`:

| Index | Variable | Description | Unit |
|---|---|---|---|
| `izoo2n` | `Zoo2N` | Macrozooplankton (krill) N | mmolN m-3 |
| `izoo2c` | `Zoo2C` | Macrozooplankton (krill) C | mmolC m-3 |
| `idetz2n` | `DetZ2N` | Fast-sinking detritus N | mmolN m-3 |
| `idetz2c` | `DetZ2C` | Fast-sinking detritus C | mmolC m-3 |
| `idetz2si` | `DetZ2Si` | Fast-sinking detritus Si | mmolSi m-3 |
| `idetz2calc` | `DetZ2Calc` | Fast-sinking detritus calcite | mmolC m-3 |
| `imiczoon` | `MicZooN` | Microzooplankton N | mmolN m-3 |
| `imiczooc` | `MicZooC` | Microzooplankton C | mmolC m-3 |

Only present with `enable_coccos = .true.`:

| Index | Variable | Description | Unit |
|---|---|---|---|
| `icocn` | `CoccoN` | Coccolithophore N | mmolN m-3 |
| `icocc` | `CoccoC` | Coccolithophore C | mmolC m-3 |
| `icchl` | `CoccoChl` | Coccolithophore chlorophyll | mgChl m-3 |
| `iphan` | `PhaeoN` | Phaeocystis N | mmolN m-3 |
| `iphac` | `PhaeoC` | Phaeocystis C | mmolC m-3 |
| `iphachl` | `PhaeoChl` | Phaeocystis chlorophyll | mgChl m-3 |

Lower bounds applied after each update:

| Variable | Definition | Unit |
|---|---|---|
| `tiny` | 2.23e-16 | same as the tracer |
| `tiny_chl` | 0.00001 | mgChl m-3 |
| `tiny_N`, `tiny_N_d`, `tiny_N_c`, `tiny_N_p` | `tiny_chl / Chl2N_max*` | mmolN m-3 |
| `tiny_C`, `tiny_C_d`, `tiny_C_c`, `tiny_C_p` | `tiny_N* / NCmax*` | mmolC m-3 |
| `tiny_Si` | `tiny_C_d * SiCmax` | mmolSi m-3 |
| `tiny_het` | 1e-5, added to `HetN` | mmolN m-3 |

## Carbon isotope tracers (ciso)

These are active with `ciso = .true.`. 14C needs `ciso_14`, and organic 14C needs
`ciso_organic_14`. Isotope tracers are scaled to 12C, so they carry the same unit as their
parent pool.

| Variable | Description | Unit |
|---|---|---|
| `DIC_13`, `DIC_14` | Dissolved inorganic 13C / 14C | mmolC m-3 |
| `PhyC_13|14`, `DiaC_13|14` | Phytoplankton / diatom 13C / 14C | mmolC m-3 |
| `DetC_13|14`, `HetC_13|14`, `EOC_13|14` | Detritus / heterotroph / DOC 13C / 14C | mmolC m-3 |
| `PhyCalc_13|14`, `DetCalc_13|14` | Calcite 13C / 14C | mmolC m-3 |
| `r_dic_13|14`, `r_co2s_13|14`, `r_phyc_*`, `r_diac_*`, `r_atm_13|14` | Isotope ratios | - |
| `alpha_*_13|14` | Fractionation factors | - |
| `delta_co2_13`, `delta_co2_14`, `big_delta_co2_14` | Atmospheric d13C, d14C, D14C | permil |
| `delta_dic_13_init`, `delta_dic_14_init`, `big_delta_dic_14_init` | Initial DIC isotope profiles | permil |
| `lambda_14` | 14C decay constant (3.8561e-12) | s-1 |
| `co2sat` | Dissolved CO2 at saturation | mol m-3 |
| `kwco2` | Schmidt-corrected piston velocity for CO2 | m s-1 |
| `co2flux_13|14` | Air-sea 13CO2 / 14CO2 flux | mol m-2 s-1 |
| `co2flux_seaicemask_13|14`, `GloCO2flux_seaicemask_13|14` | Same, converted | mmol m-2 s-1 |
| `AtmCO2_13`, `AtmCO2_14`, `LocAtmCO2_13|14` | Prescribed atmospheric 13CO2 / 14CO2 | ppm (see note on `AtmCO2`) |
| `GloPCO2surf_13|14` | Surface ocean p13CO2 / p14CO2 | uatm |

## Physical inputs and forcing

| Variable | Where | Description | Unit |
|---|---|---|---|
| `dt` | `recom` | Ocean time step | s |
| `dt_d`, `dt_b` | `REcoM_sms` | Time step, REcoM sub-step | d |
| `Temp` | `recom` / `REcoM_Forcing` | Potential temperature profile | degC |
| `Sali`, `Sali_depth` | `recom` / `REcoM_Forcing` | Practical salinity (surface, profile) | psu (PSS-78) |
| `REcoM_T`, `REcoM_T_depth` | forcing / sms_init | Temperature clamped to 2-40 for mocsy | degC |
| `REcoM_S`, `REcoM_S_depth` | forcing / sms_init | Salinity clamped to 21-43 for mocsy | psu |
| `hnode`, `thick` | mesh / `Depth_calculations` | Layer thickness | m |
| `recipthick` | `Depth_calculations` | 1 / thick | m-1 |
| `zF`, `zbar_3d_n` | `Depth_calculations` | Depth of layer interfaces (negative) | m |
| `Z_3d_n`, `zr` | `recom` | Depth of layer mid-points (negative) | m |
| `areasvol` | mesh | Area per level and node | m2 |
| `ocean_area` | mesh | Total ocean surface area (excludes cavities) | m2 |
| `geo_coord_nod2D` | mesh | Longitude, latitude | rad |
| `Latd`, `Lond` | `REcoM_Forcing` | Latitude, longitude | deg |
| `rad` | FESOM | pi / 180 | rad deg-1 |
| `shortwave` | FESOM | Shortwave radiation at the surface | W m-2 |
| `SW`, `SurfSW`, `SurfSR` | `recom` / sms | PAR at the surface: `parFrac * shortwave * (1 - ice)` | W m-2 |
| `PAR`, `PARave`, `PAR3D` | sms / GloVar | PAR in the water column | W m-2 |
| `ice_data_values`, `Loc_ice_conc` | FESOM | Sea-ice concentration (0-1) | - |
| `u_wind`, `v_wind`, `ULoc` | FESOM | 10 m wind components, wind speed | m s-1 |
| `press_air`, `Loc_slp` | FESOM | Sea-level pressure | **Pa** |
| `Patm`, `ppo` | `REcoM_Forcing` | Sea-level pressure, `Loc_slp / Pa2atm` | atm |
| `cosAI` | `Cobeta` | Cosine of the refracted angle of incidence | - |
| `daynew`, `ndpyr` | FESOM | Day of year, days per year | d |
| `mstep` | FESOM | Time-step counter | - |

## Carbonate chemistry (mocsy)

mocsy is called with `optCON='mol/m3'`, `optT='Tpot'`, `optP='m'`, `optS='Sprc'`.

**Inputs**

| Variable | Description | Unit |
|---|---|---|
| `REcoM_DIC(_depth)` | DIC | mol m-3 |
| `REcoM_Alk(_depth)` | Alkalinity | mol m-3 |
| `REcoM_Si(_depth)` | Silicate | mol m-3 |
| `REcoM_Phos(_depth)` | Phosphate, from `DIN / 16` (Redfield N:P) | mol m-3 |
| `REcoM_O2` | Surface oxygen | mol m-3 |
| `dpos` | Positive depth | m |

**Outputs**

| Variable | 3D storage | Description | Unit |
|---|---|---|---|
| `ph`, `ph_depth`, `pH_watercolumn` | `pH3D` | pH, total scale | - |
| `pco2`, `pco2_depth`, `pCO2_watercolumn` | `pCO23D` | Seawater pCO2 | uatm |
| `fco2`, `fco2_depth` | – | Seawater fCO2 | uatm |
| `pCO2surf` | `GloPCO2surf` | Surface pCO2 | uatm |
| `dpCO2surf`, `dpco2` | `GlodPCO2surf` | Ocean minus atmosphere pCO2 | uatm |
| `co2`, `co2_depth`, `CO2_watercolumn` | `CO23D` | Aqueous CO2 | mol m-3 |
| `hco3`, `hco3_depth`, `HCO3_watercolumn` | `HCO33D` | Bicarbonate | mol m-3 |
| `co3`, `co3_depth`, `CO3_watercolumn` | `CO33D` | Carbonate ion | mol m-3 |
| `OmegaC`, `OmegaA` (and `_depth`) | `OmegaC3D` | Calcite / aragonite saturation state | - |
| `kspc_depth`, `kspc_watercolumn` | `kspc3D` | Stoichiometric solubility product of calcite | mol2 kg-2 |
| `BetaD`, `BetaD_depth` | – | Revelle factor | - |
| `rhoSW`, `rhoSW_depth`, `rhoSW_watercolumn` | `rhoSW3D` | In-situ seawater density | kg m-3 |
| `p`, `p_depth` | – | Pressure | dbar |
| `tempis`, `tempis_depth` | – | In-situ temperature | degC |
| `K0` | `alphaCO2` | CO2 solubility | mol L-1 atm-1 (per source comment; see inconsistency 20) |
| `Hplus` | – | H+ concentration | mol kg-1 |
| `h_depth` | – | H+ concentration from pH | mol kg-1 |

**Calcite dissolution (DISS)**

| Variable | Description | Unit |
|---|---|---|
| `Ca` | Calcium ion, `(0.02128 / 40.078) * S / 1.80655` | mol kg-1 |
| `CO3_sat` | CO3 at calcite saturation, `kspc / Ca * rhoSW` | mol m-3 |

## Air-sea gas exchange

| Variable | Description | Unit |
|---|---|---|
| `kw660`, `PistonVelocity` | Piston velocity at Sc = 660, ice fraction included | m s-1 |
| `co2flux` | Air-to-sea CO2 flux (mocsy output) | mol m-2 s-1 |
| `co2ex` | Rate of change of surface CO2 from gas exchange | mol m-3 s-1 |
| `co2flux_seaicemask`, `GloCO2flux_seaicemask` | Air-to-sea CO2 flux, `co2flux * 1e3` | mmol m-2 s-1 |
| `dflux`, `GloCO2flux` | Air-to-sea CO2 flux, `co2flux * 1e3 * 86400` | mmol m-2 d-1 |
| `o2ex` | Air-to-sea O2 flux (mocsy output) | mol m-2 s-1 |
| `o2flux_seaicemask`, `GloO2flux_seaicemask` | Air-to-sea O2 flux, `o2ex * 1e3` | mmol m-2 s-1 |
| `oflux`, `GloO2flux` | Air-to-sea O2 flux, `o2ex * 1e3 * 86400` | mmol m-2 d-1 |

At cavity nodes (`nzmin > 1`) the fluxes that enter the boundary condition are set to 0. The
output fields are filled with NaN.

## Atmospheric CO2 box model

These variables are used with `use_atbox = .true.`.

| Variable | Description | Unit |
|---|---|---|
| `x_co2atm`, `x_co2atm_13`, `x_co2atm_14` | Atmospheric CO2 mixing ratio | ppm |
| `LocAtmCO2` | Atmospheric CO2 at the node | ppm |
| `AtmCO2(12)` | Prescribed monthly atmospheric CO2 | ppm (see note on `AtmCO2`) |
| `CO2_for_spinup` | Initial atmospheric CO2 (278) | ppm |
| `mol_allatm` | Molar inventory of the dry atmosphere (1.7726e20) | mol |
| `total_co2flux`, `total_co2flux_13|14` | Globally integrated air-sea flux | mol s-1 |
| `cosmic_14` | Cosmogenic 14C production | mol s-1 |
| `cosmic_14_init` | Initial cosmogenic 14C production flux | atoms cm-2 s-1 |
| `production_rate_to_flux_14` | Conversion from mol s-1 to atoms cm-2 s-1 | atoms cm-2 mol-1 |

## Surface boundary inputs

| Variable | Description | Unit |
|---|---|---|
| `GloFeDust` | Iron dissolved from dust, monthly | umol m-2 s-1 |
| `FeDust`, `AtmFeInput` | Iron dust input including ice and `dust_sol` | umol m-2 s-1 |
| `GloNDust` | Nitrogen from dust, 10-year mean | mmol m-2 s-1 |
| `NDust`, `AtmNInput` | Nitrogen dust input including ice | mmol m-2 s-1 |
| `relax_alk` | Alkalinity restoring flux, global mean removed | mmol m-2 s-1 |
| `Alk_surf` | Surface alkalinity climatology used for restoring | mmol m-3 |
| `virtual_din`, `virtual_dic`, `virtual_alk`, `virtual_dsi`, `virtual_dfe`, `virtual_oxy` | Virtual freshwater tracer fluxes (linfs) | mmol m-2 s-1 |
| `RiverDIN2D`, `RiverDON2D`, `RiverDOC2D`, `RiverDSi2D`, `RiverDIC2D`, `RiverAlk2D`, `RiverFe` | River input | not stated in these sources |
| `ErosionTON2D`, `ErosionTOC2D`, `ErosionTSi2D` | Erosion input | not stated in these sources |
| `RiverFeConc` | Mean DFe concentration in rivers (100) | not stated in source |

## Benthos and sediment

| Variable | Description | Unit |
|---|---|---|
| `Benthos(:, 1:benthos_num)`, `LocBenthos` | Benthic N, C, Si, calcite (plus 13C / 14C with ciso) | mmol m-2 |
| `Benthos_tr` | Per-tracer buffer of `Benthos` | mmol m-2 |
| `decayRateBenN`, `decayRateBenC`, `decayRateBenSi` | Benthic decay rate constants | d-1 |
| `decayBenthos` | Benthic release, `decayRateBen* * LocBenthos` | mmol m-2 d-1 |
| `GlodecayBenthos` | Benthic release passed to FESOM | mmol m-2 s-1 |
| `SinkFlx`, `SinkFlx_tr` | Sinking flux into the bottom layer (MEDUSA) | mmol m-2 s-1 |
| `wFluxDet`, `wFluxPhy`, `wFluxDia`, `wFluxCocco`, `wFluxPhaeo` | Sinking fluxes into the benthos | mmol m-2 d-1 |
| `Vben_*` | Sinking speed into the benthos, converted to m s-1 in `recom_sinking` | m d-1 |
| `Ironflux` | Iron flux from sediment | umolFe m-2 d-1 |
| `Fe2N_benthos` | Fe:N ratio of benthic release (0.15) | umolFe mmolN-1 |
| `q_NC_Denit` | N:C ratio of denitrification (0.86) | mmolN mmolC-1 |
| `GloSed` | Yearly input to bottom water from MEDUSA | not stated in source |
| `lb_tscale` | Time scale to balance burial loss | not stated in source |

## Sinking and ballasting

| Variable | Description | Unit |
|---|---|---|
| `VDet`, `VDet_zoo2`, `VPhy`, `VDia`, `VCocco`, `VPhaeo` | Base sinking speeds (20, 200, 0, 0, 0, 0) | m d-1 |
| `wF`, `SinkVel` | Sinking speed at layer interfaces | m d-1 |
| `Wvel_flux` | Sinking speed used in the advection scheme (negative down) | m s-1 |
| `Vdet_a` | Increase of detritus sinking speed with depth (0.0288) | d-1 (m d-1 per m) |
| `Vcalc` | Depth dependence of `calc_diss` (0.0216) | d-1 |
| `Sink_Vel` | Sinking speed used for `calc_diss` | m d-1 |
| `w_ref1`, `w_ref2` | Reference sinking speeds of small and large detritus (10, 200) | **m d-1** |
| `depth_scaling1`, `depth_scaling2` | Increase of sinking speed with depth (0.015, 0) | **d-1** (m d-1 per m) |
| `max_sinking_velocity` | Upper limit on ballasted sinking speed (250) | m d-1 |
| `Sinkingvel1`, `Sinkingvel2`, `Sinkvel1_tr`, `Sinkvel2_tr` | Sinking speeds of particle classes 1 and 2 | m s-1 (after conversion) |
| `rho_POC`, `rho_PON`, `rho_CaCO3`, `rho_opal` | Particle component densities | kg m-3 |
| `rho_ref_part`, `rho_ref_water` | Reference particle and seawater densities | kg m-3 |
| `rho_particle1`, `rho_particle2` | Particle densities | kg m-3 |
| `rho_seawater` | In-situ seawater density | kg m-3 |
| `sa` | Absolute salinity | g kg-1 |
| `visc_ref_water`, `seawater_visc_3D` | Dynamic viscosity of seawater | kg m-1 s-1 |
| `scaling_density*_3D`, `scaling_visc_3D` | Ballasting scaling factors | - |

## Internal SMS variables

**Quotas**

| Variable | Description | Unit |
|---|---|---|
| `quota`, `quota_dia`, `quota_cocco`, `quota_phaeo` | N:C | mmolN mmolC-1 |
| `recipQuota*` | C:N | mmolC mmolN-1 |
| `Chl2C*`, `CHL2C_plast*` | Chl:C (total, plastidic) | mgChl mmolC-1 |
| `Chl2N*` | Chl:N | mgChl mmolN-1 |
| `qSiC` | Si:C of diatoms | mmolSi mmolC-1 |
| `qSiN` | Si:N of diatoms | mmolSi mmolN-1 |
| `recipQZoo`, `recipQZoo2`, `recipQZoo3` | Zooplankton C:N | mmolC mmolN-1 |
| `recipDet`, `recipDet2` | Detritus C:N | mmolC mmolN-1 |
| `recip_hetN_plus` | 1 / (HetN + tiny_het) | m3 mmolN-1 |

**Temperature and oxygen functions**

| Variable | Description | Unit |
|---|---|---|
| `rTloc`, `rTref` | 1 / T and 1 / Tref | K-1 |
| `arrFunc`, `arrFuncZoo2` | Arrhenius temperature factors | - |
| `Temp_phyto`, `Temp_diatoms`, `Temp_cocco` | Exponential temperature factors | - |
| `Temp_phaeo` | Blanchard temperature function (scaled by `uopt_phaeo`) | d-1 |
| `q10_mes`, `q10_mic`, `q10_mes_res`, `q10_mic_res` | Q10 factors | - |
| `O2Func` | Oxygen limitation of remineralisation | - |
| `reminSiT` | Temperature-dependent opal dissolution rate | d-1 |

**Light**

| Variable | Description | Unit |
|---|---|---|
| `kappa`, `kappar` | Light attenuation coefficient | m-1 |
| `kappastar` | Attenuation corrected for angle of incidence | m-1 |
| `kdzUpper`, `kdzLower` | Optical depth at layer borders | - |
| `chl_upper`, `chl_lower`, `Chlave` | Chlorophyll at layer borders, layer mean | mgChl m-3 |
| `Upperlight`, `Lowerlight` | PAR at layer borders | W m-2 |

**Phytoplankton physiology**

| Variable | Description | Unit |
|---|---|---|
| `pMax*` | Maximum C-specific photosynthesis rate | d-1 |
| `Cphot*` | C-specific photosynthesis rate | d-1 |
| `phyRespRate*` | Phytoplankton respiration rate | d-1 |
| `N_assim*` | C-specific N assimilation rate | mmolN mmolC-1 d-1 |
| `Si_assim` | C-specific Si assimilation rate | mmolSi mmolC-1 d-1 |
| `ChlSynth*` | Chlorophyll synthesis term | mgChl mmolN-1 |
| `KOchl*` | Chlorophyll degradation rate | d-1 |
| `qlimitFac`, `limitFacN*`, `limitFacSi`, `FeLimitFac` | Limitation factors | - |
| `PhyCO2`, `DiaCO2`, `CoccoCO2`, `PhaeoCO2` | CO2 dependence factors | - |
| `FreeFe` | Free (unbound) iron | umolFe m-3 |

**Grazing, losses and calcite**

| Variable | Description | Unit |
|---|---|---|
| `varpz*` | Variable grazing preferences | - |
| `f*N`, `food`, `food2`, `food3` | Available food | mmolN m-3 |
| `foodsq*` | Squared food | mmolN2 m-6 |
| `grazingFlux*` | N grazing fluxes | mmolN m-3 d-1 |
| `grazingFluxcarbon_mes`, `grazingFluxcarbonzoo2` | C grazing fluxes | mmolC m-3 d-1 |
| `grazEff` | Grazing efficiency | - |
| `HetRespFlux`, `Zoo2RespFlux`, `MicZooRespFlux` | Zooplankton respiration | mmolC m-3 d-1 |
| `HetLossFlux`, `Zoo2LossFlux`, `MicZooLossFlux` | Zooplankton mortality | mmolN m-3 d-1 |
| `Zoo2fecalloss_n`, `Mesfecalloss_n` | Fecal pellet production, N | mmolN m-3 d-1 |
| `Zoo2fecalloss_c`, `Mesfecalloss_c` | Fecal pellet production, C | mmolC m-3 d-1 |
| `res_zoo2_a`, `res_zoo2_f` | Krill respiration modifiers | - |
| `AggregationRate` | Aggregation rate | d-1 |
| `calcification` | Calcite production | mmolC m-3 d-1 |
| `calc_loss_agg`, `calc_loss_gra*` | Calcite loss by aggregation and grazing | mmolC m-3 d-1 |
| `calc_diss`, `calc_diss2`, `calc_diss_ben` | Calcite dissolution rates | d-1 |
| `PICPOCtemp`, `PICPOCCO2`, `PICPOCN`, `calc_prod_final` | PIC:POC modifiers | - |

## Namelist parameters

**`pavariables` (selected)**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `VDet`, `VDet_zoo2`, `VPhy`, `VDia`, `VCocco`, `VPhaeo` | 20, 200, 0, 0, 0, 0 | Sinking speeds | m d-1 |
| `biostep` | 1 | REcoM sub-steps per ocean step | - |
| `lb_tscale` | 1 | Loopback time scale | not stated |

**`paArrhenius`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `recom_Tref` | 288.15 | Reference temperature | K |
| `C2K` | 273.15 | degC to K offset | K |
| `Ae` | 4500 | Arrhenius slope (activation energy / R) | K |
| `Tmax_phaeo`, `Topt_phaeo` | 16, 7.5272 | Blanchard maximum and optimum temperature | degC |
| `uopt_phaeo` | 0.7328 | Growth at optimum temperature | d-1 |
| `beta_phaeo` | 0.7829 | Blanchard shape | - |
| `ord_*`, `expon_*` | – | Exponential temperature function, `exp(ord + expon * T)` | - and degC-1 |
| `reminSi` | 0.02 | Minimum opal dissolution rate | d-1 |
| `k_o2_remin` | 15 | Half saturation of O2-limited remineralisation | mmolO2 m-3 |

**`palimiter_function`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `NMinSlope`, `SiMinSlope` | 50, 1000 | Limiter slopes | - |
| `NCmin*` | 0.04 | Minimum N:C | mmolN mmolC-1 |
| `SiCmin` | 0.04 | Minimum Si:C | mmolSi mmolC-1 |
| `k_Fe`, `k_Fe_d`, `k_Fe_c`, `k_Fe_p` | 0.04, 0.12, 0.04, 0.09 | Fe half saturation | umolFe m-3 |
| `k_si` | 4 | Si half saturation | mmolSi m-3 |
| `P_cm`, `P_cm_d` | 3.0, 3.5 | Maximum C-specific photosynthesis rate | d-1 |

**`palight_calculations` and `paphotosynthesis`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `k_w` | 0.04 | Attenuation by water | m-1 |
| `a_chl` | 0.03 | Chlorophyll-specific attenuation | m2 mgChl-1 |
| `alfa`, `alfa_d`, `alfa_c`, `alfa_p` | 0.14, 0.19, 0.10, 0.10 | Initial slope of the P-I curve | mmolC m2 mgChl-1 W-1 d-1 |
| `parFrac` | 0.43 | PAR fraction of shortwave | - |

**`paassimilation`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `V_cm_fact*` | 0.7 | Scaling of maximum N uptake | - |
| `NMaxSlope`, `SiMaxSlope` | 1000 | Limiter slopes | - |
| `NCmax`, `NCmax_d`, `NCmax_c`, `NCmax_p` | 0.2, 0.2, 0.15, 0.1 | Maximum N:C | mmolN mmolC-1 |
| `SiCmax` | 0.8 | Maximum Si:C | mmolSi mmolC-1 |
| `NCuptakeRatio*` | 0.2 | Maximum N:C uptake ratio | mmolN mmolC-1 |
| `SiCUptakeRatio` | 0.2 | Maximum Si:C uptake ratio | mmolSi mmolC-1 |
| `k_din`, `k_din_d`, `k_din_c`, `k_din_p` | 0.55, 1.0, 0.55, 0.55 | DIN half saturation | mmolN m-3 |
| `Chl2N_max`, `Chl2N_max_d`, `Chl2N_max_c`, `Chl2N_max_p` | 3.15, 4.2, 3.5, 3.5 | Maximum Chl:N | mgChl mmolN-1 |
| `res_phy*` | 0.01, 0.01, 0.0075, 0.008 | Maintenance respiration | d-1 |
| `biosynth` | 2.33 | Cost of biosynthesis | mmolC mmolN-1 |
| `biosynthSi` | 0 | Cost of Si assimilation | mmolC mmolSi-1 |

**`pairon_chem` and `pairon`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `totalligand` | 1 | Total ligand concentration | umol m-3 |
| `ligandStabConst` | 100 | Ligand stability constant | m3 umol-1 |
| `Fe2N` | 0.033 | Intracellular Fe:N | umolFe mmolN-1 |
| `Fe2N_benthos` | 0.15 | Fe:N of benthic release | umolFe mmolN-1 |
| `kScavFe` | 0.07 | Scavenging onto detritus | m3 mmolC-1 d-1 |
| `kScavFe2` | 0 | Background scavenging | d-1 |
| `dust_sol` | 0.02 | Dust iron solubility | - |

**`pazooplankton`, `pasecondzooplankton`, `pathirdzooplankton`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `graz_max`, `graz_max2`, `graz_max3` | 2.4, 0.1, 0.46 | Maximum specific grazing rate | d-1 |
| `epsilonr`, `epsilon2`, `epsilon3` | 0.35, 0.0144, 0.64 | Half saturation of squared food | mmolN2 m-6 |
| `res_het`, `res_zoo2`, `res_miczoo` | 0.01, 0.0107, 0.01 | Respiration rate | d-1 |
| `loss_het`, `loss_zoo2`, `loss_miczoo` | 0.05, 0.003, 0.01 | Quadratic mortality | m3 mmolN-1 d-1 |
| `Redfield` | 6.625 | C:N Redfield ratio | mmolC mmolN-1 |
| `pz*` | – | Maximum grazing preferences | - |
| `s*Nsq*` | 0 | Half saturation of squared prey ('sigmoid' form only) | mmolN2 m-6 |
| `fecal_rate_n*`, `fecal_rate_c*` | – | Fecal pellet fractions | - |
| `t1_zoo2`, `t2_zoo2`, `t3_zoo2`, `t4_zoo2` | – | Krill temperature function constants | K |
| `tiny_het` | 1e-5 | Regularisation of `HetN` | mmolN m-3 |

**Remaining process parameters**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `agg_PD`, `agg_PP` | 0.165, 0.015 | Aggregation parameters | m3 mmolN-1 d-1 |
| `rho_N` | 0.11 | DON remineralisation rate | d-1 |
| `rho_C1` | 0.1 | DOC remineralisation rate | d-1 |
| `rho_c2` | 0.1 | TEP-C degradation rate | d-1 |
| `lossN*`, `lossC*` | 0.05, 0.10 | Phytoplankton exudation | d-1 |
| `lossN_z*`, `lossC_z*` | – | Zooplankton excretion | d-1 |
| `deg_Chl*` | 0.25, 0.25, 0.20, 0.25 | Chlorophyll degradation | d-1 |
| `gfin`, `grazEff2`, `grazEff3` | 0.3, 0.8, 0.8 | Grazing efficiency | - |
| `reminN`, `reminC` | 0.165, 0.15 | Detritus remineralisation | d-1 |
| `calc_prod_ratio` | 0.02 | Calcite production ratio | - |
| `calc_diss_guts` | 0 | Gut dissolution fraction | - |
| `calc_diss_rate`, `calc_diss_rate2` | 0.005714 | Calcite dissolution rate | d-1 |
| `calc_diss_omegac` | 0.197 | Omega-dependent dissolution rate | d-1 |
| `calc_diss_exp` | 1 | Dissolution order | - |
| `surf_relax_Alk` | 3.2e-7 | Alkalinity restoring velocity (10 m per year) | m s-1 |

**`paco2lim` and `paco2_flux_param`**

| Parameter | Default | Description | Unit |
|---|---|---|---|
| `Cunits` | 976.5625 | Conversion from mol m-3 to umol kg-1 | (umol kg-1) / (mol m-3) |
| `a_co2_*` | – | CO2 limitation coefficient | - |
| `b_co2_*` | – | CO2 limitation coefficient | mol kg-1 |
| `c_co2_*`, `d_co2_*` | – | CO2 limitation coefficients | kg mol-1 |
| `permil` | 9.76e-7 | Conversion from mmol m-3 to mol kg-1 (1e-3 / 1024.5) | (mol kg-1) / (mmol m-3) |
| `permeg` | 1e-6 | Conversion from uatm to atm | atm uatm-1 |
| `Xacc` | 1e-12 | pH iteration accuracy | - |
| `CO2_for_spinup` | 278 | Spin-up atmospheric CO2 | ppm |

**Constants**

| Constant | Value | Unit |
|---|---|---|
| `SecondsPerDay` | 86400 | s d-1 |
| `Pa2atm` | 101325 | Pa atm-1 |
| `redO2C` | 1.453 | molO2 molC-1 |
| `nWater` | 1.33 | - |

## Diagnostics

The 3D diagnostics are per-layer rates. The 2D diagnostics are the same rates integrated over
the water column (`sum(vert* * thick)`).

| Variable | Description | Unit |
|---|---|---|
| `vertNPP*`, `NPPn3D`, `NPPd3D`, `NPPc3D`, `NPPp3D` | Net primary production | mmolC m-3 d-1 |
| `NPPn`, `NPPd`, `NPPc`, `NPPp` | Vertically integrated NPP | mmolC m-2 d-1 |
| `vertGPP*` / `GPPn`, `GPPd`, `GPPc`, `GPPp` | Gross primary production | mmolC m-3 d-1 / mmolC m-2 d-1 |
| `vertNNA*` / `NNAn`, `NNAd`, `NNAc`, `NNAp` | Net N assimilation | mmolN m-3 d-1 / mmolN m-2 d-1 |
| `vertChldeg*` / `Chldegn`, `Chldegd`, `Chldegc`, `Chldegp` | Chlorophyll degradation rate | d-1 / m d-1 (see the notes at the end) |
| `vertgrazmeso_*`, `vertgrazmacro_*`, `vertgrazmicro_*` / `grazmeso_*`, `grazmacro_*`, `grazmicro_*` | Grazing (C units, efficiency applied) | mmolC m-3 d-1 / mmolC m-2 d-1 |
| `respmeso`, `respmacro`, `respmicro` | Zooplankton respiration | mmolC m-3 d-1 |
| `respn`, `respd`, `respc`, `respp` | Phytoplankton respiration | mmolC m-3 d-1 |
| `calcif`, `calcdiss` | Calcification and calcite dissolution | mmolC m-3 d-1 |
| `aggn`, `aggd`, `aggc`, `aggp` | Aggregation | mmolC m-3 d-1 |
| `docexn`, `docexd`, `docexc`, `docexp` | DOC exudation | mmolC m-3 d-1 |
| `TTemp_*`, `TPhyCO2`, `TDiaCO2`, `TCoccoCO2`, `TPhaeoCO2`, `TqlimitFac_*`, `TCphotLigLim_*` | Growth limitation factors | - |
| `TCphot_*` | Photosynthesis rate | d-1 |
| `TSi_assimDia` | Diatom Si assimilation | mmolSi mmolC-1 d-1 |
| `tra_recom_sms` | Layer-integrated SMS tendency for `ltra_diag` tracers | mmol m-2 s-1 |
| `PAR3D` | PAR in the water column | W m-2 |
| `DenitBen` | Benthic denitrification | not stated in source |

## Unit conversions in the code

| Location | Expression | Conversion |
|---|---|---|
| `REcoM_Forcing`, `sms_initialize_variables` | `state * 1e-3` | mmol m-3 to mol m-3 (mocsy inputs) |
| `REcoM_Forcing`, `sms_initialize_variables` | `DIN * 1e-3 / 16` | DIN to phosphate (Redfield N:P) |
| `REcoM_Forcing` | `co2flux * 1e3 * SecondsPerDay` | mol m-2 s-1 to mmol m-2 d-1 |
| `REcoM_Forcing` | `co2flux * 1e3` | mol m-2 s-1 to mmol m-2 s-1 |
| `REcoM_Forcing` | `Loc_slp / Pa2atm` | Pa to atm |
| `REcoM_Forcing` | `geo_coord_nod2D / rad` | rad to deg |
| `recom_init` | iron initial field `* 1e9` | mol L-1 to umol m-3 |
| `REcoM_sms` | `dt / SecondsPerDay` | s to d |
| `recom` | `decayBenthos / SecondsPerDay` | mmol m-2 d-1 to mmol m-2 s-1 |
| `recom_sinking` | `-Vsink / SecondsPerDay` | m d-1 to m s-1 |
| `recom_atbox` | `0.001 * GloCO2flux_seaicemask` | mmol m-2 s-1 to mol m-2 s-1 |
| `recom_atbox` | `flux / mol_allatm * dt * 1e6` | mol to ppm |
| `REcoM_sms_init` | `Temp + C2K` | degC to K |

---

## Known inconsistencies in source comments

The units above follow the code. In these places the comment in the source says something
else, so the comment should be corrected.

| # | File | Variable | Comment says | Code implies |
|---|---|---|---|---|
| 1 | `recom_modules.F90` (GloVar) | `GloCO2flux_seaicemask`, `GloO2flux_seaicemask` | mmol/m2/day | mmol m-2 s-1 (`co2flux * 1e3`; `recom_main.F90` also says /s) |
| 2 | `recom_modules.F90` (GloVar) | `GlodPCO2surf` | mmol/m2/day | uatm (it stores `dpco2surf`) |
| 3 | `recom_modules.F90` (GloVar) | `GloHplus` | mol/kg H+ | holds **pH** (`GloHplus(n) = ph(1)`), although it is initialised with 1e-8 (an H+ value) |
| 4 | `recom_modules.F90` (GloVar) | `AtmNInput` | umol/m2/s | mmol m-2 s-1, like its source `GloNDust` |
| 5 | `recom_modules.F90` (GloVar), `recom_main.F90` | `rhoSW3D`, `rhoSW_watercolumn` | mol/m3 | kg m-3 (mocsy density) |
| 6 | `recom_main.F90` | `Loc_slp` | atm | Pa (defaults to `Pa2atm` = 101325, then `press_air`) |
| 7 | `recom_sms_init.F90`, `recom_sms_update.F90`, `recom_sms.F90` | `Fe`, `FreeFe`, `k_Fe` | mmolFe m-3 | umolFe m-3 (initial field converted to umol m-3, `Fe2N` in umol mmolN-1, init log prints umol m-3; the iron chemistry block in `recom_sms.F90` already says umol m-3) |
| 8 | `recom_modules.F90` (config) | `w_ref1`, `w_ref2` | m s-1 | m d-1 (divided by `SecondsPerDay` before use) |
| 9 | `recom_modules.F90` (config) | `depth_scaling1`, `depth_scaling2` | s-1 | d-1 (added to a speed in m d-1 before conversion) |
| 10 | `recom_modules.F90` (declarations, GloVar) | `decayBenthos`, `GlodecayBenthos` | 1/day | `decayBenthos` = mmol m-2 d-1 (`decayRateBen* * LocBenthos`); `GlodecayBenthos` = mmol m-2 s-1 |
| 11 | `recom_modules.F90` (config) | `graz_max`, `graz_max2`, `graz_max3` | mmol N/(m3 day) | d-1 (multiplied by `HetN` in the grazing flux; `recom_sms.F90` also says d-1) |
| 12 | `recom_modules.F90` (declarations) | `HetLossFlux`, `Zoo2LossFlux`, `MicZooLossFlux` | (mmol N)2/(m6 day) | mmolN m-3 d-1 |
| 13 | `recom_modules.F90` (config) | `loss_het`, `loss_zoo2`, `loss_miczoo` | "1/day, temperature dependent N degradation of EON" | m3 mmolN-1 d-1, quadratic mortality (the description is also wrong) |
| 14 | `recom_modules.F90` (declarations) | `PARave`, `Upperlight`, `Lowerlight` | [?] | W m-2 |
| 15 | `recom_modules.F90` (config) | `AtmCO2` | uatm | used as `LocAtmCO2` and as `x_co2atm` (ppm mixing ratio); the mocsy flux routine converts to pCO2 using `Patm` |
| 16 | `recom_sms.F90` | `alfa` | mmolC (mgChl)-1 m2 umol-1 ... | mmolC m2 mgChl-1 W-1 d-1 (`PARave` is in W m-2) |
| 17 | `recom_forcing.F90` | `sms` | mmol/m3/timestep | consistent with the code; `recom_sms.F90` says mmol m-3 day-1, which is wrong because each term is multiplied by `dt_b` |
| 18 | `recom_modules.F90` (declarations) | `Zoo2fecalloss_c` | (mmol N)/(m3 day) | mmolC m-3 d-1 |
| 19 | `recom_modules.F90` (declarations) | `recipDet`, `recipDet2` | "Quota between second zoo C and N" | detritus C:N (unit is correct, description is copied) |
| 20 | `recom_modules.F90` (GloVar) | `alphaCO2` (`K0`) | mol L-1 atm-1 | mocsy documents K0 in mol kg-1 atm-1; check against the mocsy version in use |

Two further notes:

- `Chldegn`, `Chldegd`, `Chldegc`, `Chldegp` integrate a rate in d-1 over layer thickness, which
  gives m d-1. That is probably not a meaningful 2D diagnostic.
- `FeLimitFac` is commented as "Mumol/m3, half sat constant for iron" but is used as a
  dimensionless limitation factor.
