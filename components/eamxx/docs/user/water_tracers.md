# Water tracers

Water tracers follow a labeled portion of the model's water through the
atmosphere. They are built from two processes, `water_tracers` and
`water_tracer_transfers`. Only vapor transport (HOMME advection and SHOC
diffusion) and the bulk water phase changes of SHOC and P3 are covered so far.
Generic tracers move in proportion to the donor phase's tracer-to-bulk ratio
in every transfer (no fractionation). Isotope tracers also fractionate, with
equilibrium fractionation factors, in the phase changes through the vapor (see
"Isotope fractionation" below).

## Tracers file

The tracers are defined in a YAML file with a top-level `tracers` map, set with
`water_tracers::water_tracers_file`:

``` {.yaml .copy}
tracers:
  wtrc_t:    {type: generic, init: ratio, r_init: 1}
  wtrc_half: {type: generic, init: ratio, r_init: 0.5}
```

Each entry has a `type` (`generic` or `isotope`), an `init` (`ratio`: start from
`r_init` times the bulk water, or `zero`), and `r_init`. Generic entries also
accept `ocean_flux` and `r_ocn`; isotope entries need a `species` (and a build
with `EAMXX_ENABLE_WATER_ISOTOPES=ON`). The key is the tracer's base name.

Both kinds accept `r_land` (default: `r_ocn`; for isotopes, the species' ocean
ratio), the tracer-to-bulk ratio of the evaporation from land and sea ice. See
"Surface evaporation" below.

## Fields

Each tracer has one field per bulk water phase, named `<base>_{qv,qc,qr,qi}`
(kg tracer per kg moist air, like the bulk fields), for example `wtrc_t_qv`.
The tracer vapor is advected by HOMME and diffused by SHOC. The other phases
only change through `water_tracer_transfers`.

## Surface evaporation

For each tracer, `water_tracers` computes the field `<base>_qv_surf_flux`
(kg/m2/s, positive upward; for example `wtrc_t_qv_surf_flux`) every step, in the
group `tracer_surface_fluxes`. `shoc` applies these fluxes in the same implicit
solve as the surface evaporation of bulk water (`surf_evap`, the merged ocean,
land and sea-ice flux), so a tracer with R = 1 stays equal to bulk water. SHOC
does not know the tracers: it applies each flux to the tracer whose field is the
flux's name without `_surf_flux`. The driver turns this on (SHOC's parameter
`apply_tracer_surface_fluxes`) when both `water_tracers` and `shoc` are in the
tree, and requires `water_tracers` to come before `shoc`, since the flux uses
the `surf_evap` of the start of the step.

With S = `surf_evap`, E the ocean-only bulk evaporation (computed from `ocnfrac`,
`ocn_ssq`, `ocn_re`, `ocn_ustar` and the bottom level) and L = S - E the rest:

* the ocean part follows `ocean_flux`/`r_ocn` as before (for `bulk_proportional`,
  `r_ocn*E`);
* the land and sea-ice part evaporates at `r_land` (`r_land*L`, for L > 0);
* dew (E < 0 or L < 0) removes the tracer at its bottom-level ratio to `qv`.

For `bulk_proportional` tracers, when both parts evaporate, the flux is formed as
`r_land*S + (r_ocn - r_land)*E`, so with `r_ocn = r_land = 1` it is exactly `S`.
`bulk_formula` tracers and isotopes compute the ocean part from the tracer's own
bottom-level concentration, so they do not close exactly against `surf_evap`
(they do, to roundoff, for R = 1); their land part is the same.

For isotopes, the ocean part also carries two fractionation factors: the
equilibrium factor alpha_l at the SST, and the kinetic factor alpkn of
evaporation into the boundary layer (a port of iCAM's `wiso_kmol`: Brutsaert
(1975) resistances with the Merlivat & Jouzel (1979) kinetic modifier, a
smooth or rough regime from the surface Reynolds number). Two options of
`water_tracers` control alpkn:

* `kinetic_fractionation_formulation`: `brutsaert_merlivat_jouzel_1979`
  (default) or `none` (alpkn = 1);
* `diffusivity_formulation`: the molecular diffusivity ratios it uses,
  `merlivat_1978` (default) or `cappa_2003`.

alpkn lies in (0, 1] for ustar times the bottom-level height above about
1e-5 m2/s; there is no low-wind guard, as in the Fortran (debug builds check
the range). The land and sea-ice part never fractionates.

## The two processes

* `water_tracers` registers the fields, sets their initial conditions and
  computes the surface evaporation flux of each tracer (next section). It never
  changes the bulk water or a tracer after initialization.
* `water_tracer_transfers` replays what SHOC (condensation/evaporation) and P3
  (all process rates, sedimentation, clipping) did to the bulk water on the
  tracers, in the order the bulk model applied it. It also computes
  `<base>_precip_{liq,ice}_surf_mass` (kg/m2, accumulated like the bulk
  precipitation) and the diagnostics `<base>_residual_{shoc,p3}`
  (bulk water minus replayed tracer water, per phase), and once for all
  tracers `water_transfer_stage_T` (see "Stage temperatures") and
  `water_transfer_alpha_extrapolations` (see "Isotope fractionation").

Parameters of `water_tracer_transfers`:

* `replay_shoc`, `replay_p3` (default true): which sources to replay.
* `water_tracers_file`: copied from `water_tracers` when unset.
* `liquid_vapor_formulation`, `ice_vapor_formulation`,
  `standard_ratio_formulation`, `ocean_enrichment_formulation`,
  `kinetic_fractionation_formulation`, `diffusivity_formulation`: the isotope
  options, copied from `water_tracers` when unset (and an error if set to a
  different value), so both processes use the same fractionation factors.
* `alpha_override` (test only, not in the namelist defaults; default 0, off):
  if > 0, every isotope uses this value as both alpha_l and alpha_i in the
  replay (so Bergeron's factor is 1). With 1, an isotope follows the bulk
  water like a generic tracer, bit for bit.

## Isotope fractionation

Only tracers of `type: isotope` fractionate; generic tracers never do. For an
isotope, each replayed transfer through the vapor removes tracer from its donor
phase with a fractionation factor beta = R_sink/R_donor, from the equilibrium
factors alpha_l = R_liquid/R_vapor and alpha_i = R_ice/R_vapor of the species
(the formulations set by `liquid_vapor_formulation` and
`ice_vapor_formulation`):

| Transfer | beta | Temperature |
|---|---|---|
| SHOC net condensation (vapor to cloud liquid) | alpha_l | `ShocAfterTransport` |
| SHOC net evaporation (cloud liquid to vapor) | 1/alpha_l | `ShocAfterTransport` |
| P3 vapor deposition and ice nucleation | alpha_i | `P3Part2` |
| P3 rain evaporation | 1/alpha_l | `P3Part2` |
| P3 Bergeron process (cloud liquid to ice) | alpha_i/alpha_l | `P3Part2` |
| P3 sublimation, and every transfer between condensates (autoconversion, accretion, riming, collection, freezing, melting, shedding) | 1 | |
| P3 clips and homogeneous freezing, sedimentation, the checkpoint vapor repairs, surface fluxes | never fractionate | |

The Bergeron factor treats the transfer as going through the vapor (the liquid
evaporates and the vapor deposits on the ice). This is a deliberate departure
from the EAM Fortran, which has no Bergeron branch and so uses 1. Sublimation
uses 1, as in the Fortran: ice does not exchange isotopes with the vapor while
it sublimates.

Each alpha is evaluated once per level and tracer, at the temperature at the
start of its stage (`water_transfer_stage_T`, below), on the levels where the
transfers that use it moved water.

**Rayleigh step.** All the sinks of one donor phase act together, as a
Rayleigh distillation of the donor's water: if the sinks remove A of the
donor's bulk water m (for P3, its reservoir plus its inflows in the same
call), the fraction of the donor left is f = (m - A)/m, and the donor's tracer
t loses t*(1 - f^beta_eff), with beta_eff the sinks' amount-weighted mean
beta. That tracer is split among the sinks in proportion to amount*beta. This
cannot move more tracer than the donor has, or leave tracer in a donor whose
bulk water is all gone, and it does not depend on the order of the sinks or on
splitting the step. Where every beta of a donor is 1 (always for generic
tracers) the transfer is the linear one, amount*t/m, exactly as before. After
P3's fixed-point passes over the donor ratios, a guard scales down the
outflows of any fractionating donor that would exceed its tracer plus its
final inflows (this conserves the tracer).

**Extrapolated fits.** The fits of alpha are used as they are outside their
fitted temperature range. Liquid alpha counts as extrapolated below 253.15 K,
and supercooled liquid in P3 can be as cold as about 233 K; cirrus ice is
colder than the 233-240 K lower limit of the ice fits. The diagnostic
`water_transfer_alpha_extrapolations` (per column, accumulated over the atm
step) counts the alpha evaluations at such temperatures, and their total over
the run is logged once, at the end of the run. A temperature below 50 K, or
NaN, stops the run.

**Not represented.** Kinetic effects (of evaporation into undersaturated air
and of deposition in supersaturation) are not included. Neither is any
isotopic exchange at zero net mass change: existing cloud liquid and rain do
not re-equilibrate with the vapor (only the water SHOC or P3 move between
phases carries a fractionation). SHOC exports its net condensation per level,
so condensation and evaporation within a level in the same call are not
seen separately.

## Stage temperatures

`water_transfer_stage_T` (K; layout column, `water_transfer_stage`, level) is
the temperature at the start of each replayed stage, at which the isotopes'
fractionation factors are evaluated. Its three components are, in order:

1. `ShocAfterTransport`: SHOC after its turbulent transport, before its
   condensation;
2. `P3Entry`: P3's entry state;
3. `P3Part2`: P3 after its first clipping of small masses, the state its
   process rates are computed from.

SHOC and P3 export no temperature. Both heat in proportion to the water they
move, so `water_tracer_transfers` rebuilds the stage temperatures from `T_mid`
and the exports, undoing P3's latent heating (`p3::water_transfer_phases` in
`p3_water_transfers.hpp` gives the phases, and so the latent heat, of each
transfer) and then SHOC's condensation heating. This needs `T_mid` to be
changed only by SHOC and P3 between SHOC and `water_tracer_transfers`; the
driver checks it (see below). A component whose stage is not replayed
(`replay_shoc` or `replay_p3` false) holds EAMxx's fill value (about 3.4e33,
marked as `_FillValue` in output files).

`ShocAfterTransport` includes SHOC's energy fixer, which corrects SHOC's final
temperature and is not exported separately.

## Enabling

Neither process is in a default `atm_procs_list`. Put `water_tracers` before
`homme` (and so before `shoc`) and `water_tracer_transfers` directly after `p3`, inside the
(sub-cycled) group that holds `shoc` and `p3`:

``` {.shell .copy}
./atmchange atm_procs_list=sc_import,water_tracers,homme,physics,sc_export
./atmchange physics::atm_procs_list=mac_aero_mic,rrtmgp
./atmchange mac_aero_mic::atm_procs_list=shoc,cld_fraction,spa,p3,water_tracer_transfers
./atmchange water_tracers::water_tracers_file=/path/to/tracers.yaml
```

Builds with `EAMXX_ENABLE_WATER_TRACERS=ON` make SHOC and P3 export the bulk
water transfers that `water_tracer_transfers` replays (in every run of such a
build, whether or not the tracer processes are used). At startup the
atmosphere driver copies the tracers file to `water_tracer_transfers` and checks
that:

* `water_tracers` exists, and `shoc` and `p3` (when replayed) are in the same
  group as `water_tracer_transfers`, in the order `shoc`, `p3`,
  `water_tracer_transfers`;
* no process between them has its own `number_of_subcycles` other than 1
  (sub-cycle the enclosing group instead; the exports only keep the last call);
* `shoc::check_flux_state_consistency` is false when SHOC is replayed;
* `water_tracers` comes before `shoc`, if both are in the tree, and
  `shoc::check_flux_state_consistency` is then false (the fixer zeroes
  `surf_evap` during the step, so bulk water and the tracers would get different
  surface fluxes);
* no process between SHOC/P3 and the replay updates `qv`, `qc`, `qr`, `qi`, the
  tracer groups or `T_mid` (the stage temperatures are rebuilt from `T_mid`).

A warning is logged if `zm`, `gw`, `nudging` or `iop_forcing` is in the tree.

## Limits

* Only SHOC and P3 are replayed. Other bulk water changes (ZM, GW, nudging,
  IOP forcing, the SHOC flux-consistency fixer) are not followed.
* Only equilibrium fractionation, in the transfers listed above (see "Not
  represented"). The surface flux of land and sea ice uses `r_land` (default
  `r_ocn`); the coupler's land isotope fluxes are not imported.
* A generic tracer with a uniform ratio R keeps R through every transfer; the
  tests check R = 1 against the bulk fields to roundoff.

## Tests

``` {.shell}
ctest -R "water_tracer|wtt|shoc_water|p3_water" --output-on-failure
```

The multi-process tests are in
`tests/multi-process/physics_only/shoc_water_tracer_transfers` and
`shoc_cld_p3_water_tracer_transfers`. The latter includes a link test
(`shoc_cld_p3_wtt_link`, nothing configured by hand) and negative tests
(`wtt_link_neg_*`) for mis-ordered or mis-grouped processes, a missing
`water_tracers`, sub-cycled P3, and a process in between that updates bulk water
or `T_mid` (a test-only dummy process, since no linked physics does). The P3-only
run checks the rebuilt `P3Entry` temperature against the previous step's `T_mid`
(nothing else changes it there); the SHOC stage is only checked for consistent
bookkeeping, since no output holds SHOC's temperature after transport.

In builds with `EAMXX_ENABLE_WATER_ISOTOPES=ON`, the `*_water_iso*` tests run
HDO and H218O next to a generic R = 1 tracer (`check_water_isotope_transfers.py`):
column budgets (tracer + precipitation - surface flux) to roundoff, no negative
tracer beyond roundoff (also with `alpha_override` 1.5 and 0.7), new SHOC cloud
with R_qc/R_qv equal to the Rayleigh result for alpha_l, condensing cloud
enriched relative to the vapor, and ice deposited from vapor near alpha_i.
With `alpha_override: 1.0` the isotopes must equal the generic tracer bit for
bit. The unit test `water_tracer_transfers_functions` covers the Rayleigh step
itself.
