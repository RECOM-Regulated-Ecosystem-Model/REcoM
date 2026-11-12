# Changelog

All notable changes to REcoM are listed here. Versions follow
[semantic versioning](https://semver.org) within the v3 series:

- **minor** (v3.x.0): new options, changes in results, and changes that need
  an update outside REcoM (marked **breaking** below);
- **patch** (v3.x.y): bit-identical with the previous release.

FESOM pins REcoM by release tag.

## [v3.2.0] - #56

Changes since `r3.1.3` (2026-07-16).

### Breaking (FESOM or namelists must change)

- `recom()` takes `stress_atmoce_x`, `stress_atmoce_y` after `v_wind`, used to
  derive the 10 m wind from the wind stress in ECHAM-coupled runs (#33).
  Needs the matching FESOM call.
- The DIC remineralization tracer (ID 1037, #19) is optional with
  `use_DICremin`, default `.false.` (#33). Setups that carry tracer 1037 in
  `namelist.tra` and count it in `bgc_num` must set `use_DICremin = .true.`.
- Unused configuration flags and PFT-specific parameters removed (#38); a
  `namelist.recom` that still sets them fails to read.
- `__oasis` is defined for coupled builds and unsupported couplings are
  rejected at build time (#35).
- `fe_2ligands = .true.` without `fe_compl_nica = .true.` stops at start-up
  (#51); before, iron scavenging was silently off in that combination.

### Changed results

- Deposition into the benthos uses the same sinking speed as the water
  column (ballasted detritus 1 and 2); phytoplankton that do not sink are no
  longer deposited; deposition is capped at 0.9 of the layer content per
  step (#53). Removes REcoM's carbon source from re-flooring deposited
  phytoplankton (~8.7e11 mmol C/day in the 2p3z2d test).
- Sign of the limited face value in the DST3 sinking flux corrected (#55):
  sharper sinking profiles, export and remineralization depths change.
- The sinking increment is divided by `hnode` instead of `hnode_new` (#54);
  the sinking operator now conserves mass to rounding under zstar.
- O2 air-sea flux uses the surface temperature limited to -2..40 deg C
  instead of the 2..40 deg C clamp of the carbonate constants (#44).
- Tendencies of all bio sub-steps are kept; before, `biostep > 1` kept only
  the last sub-step (#50). No change with `biostep = 1`.
- Carbonate chemistry scheduled by the step of the year, not of the run, so
  restarts reproduce continuous runs (#43).
- Fe scavenging onto calcite and opal (`kScavFe_calc`, `kScavFe_si`, default
  0) and background scavenging loss `kScavFe2` (#33).
- Depth dependence of calcite dissolution (#33).
- Variable grazing efficiency and two-ligand Fe chemistry (#33, #36).
- Cavity-safe sinking, surface state, benthos deposition and release (#25,
  #27); nzmin propagated through the BGC forcing chain (#20).
- Virtual BGC tracer fluxes for linear free surface (#29).
- Coupled CO2 with OpenIFS (#22) and sea-level pressure for gas exchange in
  the OpenIFS-coupled build (#34).
- Carbon isotopes: atmospheric box model, isotope tracer initialisation and
  C-isotope corrections ported from int_recom (#33).
- BGC-only loops exclude transient tracers (#18).

### Fixed without change in results

- One `recom_abort()` for all error exits: `MPI_Abort` on `MPI_COMM_WORLD`
  with `ierror`, then `error stop` (#45).
- Two O(N^2) loops in the node loop removed (`Cobeta`, benthos reset) (#47).
- GNU builds use `-fdefault-real-8 -fdefault-double-8 -finit-local-zero`
  everywhere, as Intel uses `-r8 -init=zero` (#52).
- NNA diagnostics use the group loss rates `lossN_d/_c/_p` (#49).
- Per-column diagnostic allocation matches its use (#37).
- Duplicated benthos accumulation that broke the build removed (#24).
- Unused `sms_update_tracer_scalars` (#48), redundant `DIC_PI` (#32) and
  unused benthos-flux declarations (#46) removed; comments corrected (#46).

### Added

- `recom_log` module for tagged, coloured terminal output (#40).
- Tracer budget diagnostic of `recom_sms` as a rate in mmol m-2 s-1 (#39).
- CMIP7 output (#41).
- `UNITS.md` with the units of REcoM variables (#42).

### Known limitations

- The CO2 air-sea flux (solubility, Schmidt number, carbonate chemistry) is
  evaluated at the surface temperature clamped to 2 deg C in colder water.
- The ballasting density scaling has no lower limit: organic-rich detritus
  can get negative (upward) speeds at depth, where in-situ seawater is denser
  than `rho_POC`; the speed cap turns large negative speeds into +max.
- mocsy is evaluated at the depth of the layer top, not its mid-depth.
- Carbon isotopes (`ciso`) are not mass-conserving outside the base
  configuration.

### Not itemised

PRs #23 (Yings branch), #28 (Rdev) and #31 (Fixes) bundle several changes;
their content is included in the sections above where known.
