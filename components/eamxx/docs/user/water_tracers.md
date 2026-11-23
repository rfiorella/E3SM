# Water tracers

Water tracers follow a labeled portion of the model's water through the
atmosphere. Generic tracers move in proportion to the donor phase's
tracer-to-bulk ratio in every transfer (no fractionation). Isotope tracers also
fractionate, with equilibrium fractionation factors (and, optionally, kinetic
ones), in the phase changes through the vapor and in ocean evaporation. Only
vapor transport (HOMME advection and SHOC diffusion), surface evaporation and
the bulk water phase changes of SHOC and P3 are covered so far.

## Overview

Two processes do the work, with SHOC and P3 in between:

1. `water_tracers` registers the tracer fields, sets their initial conditions
   and computes, every step, the surface flux of each tracer. It never changes
   the bulk water or a tracer after initialization.
2. `shoc` applies the tracer surface fluxes in the same implicit solve as the
   bulk surface evaporation, and diffuses the tracers like bulk water (HOMME
   advects them).
3. `shoc` and `p3` export the bulk water transfers they perform (the amount
   moved between phases by each process, per level), in builds with
   `EAMXX_ENABLE_WATER_TRACERS=ON`. The exports are read-only: the bulk model
   does not know the tracers.
4. `water_tracer_transfers` replays those transfers on the tracers, in the order
   the bulk model applied them. Isotopes are fractionated in the replay: the
   tracer moved by a transfer is the donor's tracer-to-bulk ratio times the
   amount of bulk water moved, times a fractionation factor.

The atmosphere driver links the processes at startup and checks the process
order (see "Enabling").

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

`water_tracers` also computes `<base>_qv_surf_flux` (see "Surface evaporation").
`water_tracer_transfers` computes `<base>_precip_{liq,ice}_surf_mass` (kg/m2,
accumulated like the bulk precipitation) and the diagnostics
`<base>_residual_{shoc,p3}` (bulk water minus replayed tracer water, per phase),
and once for all tracers `water_transfer_stage_T` (see "Stage temperatures") and
`water_transfer_alpha_extrapolations` (see "Isotope fractionation").

## Configuration

Parameters of `water_tracers`:

* `water_tracers_file`: the tracers file above.
* The isotope options (only used if the tracers file defines isotope entries):
    * `liquid_vapor_formulation`: `horita_wesolowski_1994` (default) or
      `majoube_1971`;
    * `ice_vapor_formulation`: `merlivat_nief_1967` (default) or `isocam3`;
    * `standard_ratio_formulation`: `normalized` (default; every standard ratio
      is 1) or `natural_abundance`;
    * `ocean_enrichment_formulation`: `modern` (default) or `LGM`;
    * `kinetic_fractionation_formulation`: kinetic factor of the ocean
      evaporation, `brutsaert_merlivat_jouzel_1979` (default) or `none`;
    * `diffusivity_formulation`: the molecular diffusivity ratios of the kinetic
      factors, `merlivat_1978` (default) or `cappa_2003`.

Parameters of `water_tracer_transfers`:

* `replay_shoc`, `replay_p3` (default true): which sources to replay.
* `water_tracers_file`: copied from `water_tracers` when unset.
* `liquid_vapor_formulation`, `ice_vapor_formulation`,
  `standard_ratio_formulation`, `ocean_enrichment_formulation`,
  `kinetic_fractionation_formulation`, `diffusivity_formulation`: the isotope
  options, copied from `water_tracers` when unset (and an error if set to a
  different value), so both processes use the same fractionation factors.
* `kinetic_fractionation` (default false): kinetic isotope fractionation in
  the replay (see "Kinetic fractionation" below). It affects the replay only;
  the kinetic factor of ocean evaporation is
  `water_tracers::kinetic_fractionation_formulation` (on by default).
* `kinetic_ice_supersaturation` (`empirical`, default, or `model`): the
  supersaturation over ice of the kinetic factor of deposition, with
  `kinetic_fractionation`.
* `alpha_override` (test only, not in the namelist defaults; default 0, off):
  if > 0, every isotope uses this value as both alpha_l and alpha_i in the
  replay (so Bergeron's factor is 1), and no kinetic factor. With 1, an
  isotope follows the bulk water like a generic tracer, bit for bit.

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

* the ocean part follows `ocean_flux`/`r_ocn` (for `bulk_proportional`,
  `r_ocn*E`);
* the land and sea-ice part evaporates at `r_land` (`r_land*L`, for L > 0);
* dew (E < 0 or L < 0) removes the tracer at its bottom-level ratio to `qv`.

For `bulk_proportional` tracers, when both parts evaporate, the flux is formed as
`r_land*S + (r_ocn - r_land)*E`, so with `r_ocn = r_land = 1` it is exactly `S`.
`bulk_formula` tracers and isotopes compute the ocean part from the tracer's own
bottom-level concentration, so they do not close exactly against `surf_evap`
(they do, to roundoff, for R = 1); their land part is the same.

For isotopes, the ocean part also carries two fractionation factors, described
under "Isotope fractionation": the equilibrium factor alpha_l at the SST and the
kinetic factor alpkn. The land and sea-ice part never fractionates.

## The replay

`water_tracer_transfers` replays what SHOC (condensation/evaporation) and P3
(all process rates, sedimentation, clipping) did to the bulk water on the
tracers, in the order the bulk model applied it. Each transfer moves tracer out
of its donor phase in proportion to the donor's tracer-to-bulk ratio (for
isotopes, with the fractionation factors below). It must be placed directly
after `p3` in the group that holds `shoc` and `p3`, because the exports only
hold the last call of each.

## Isotope fractionation

Only tracers of `type: isotope` fractionate; generic tracers never do.

### Equilibrium fractionation

For an isotope, each replayed transfer through the vapor removes tracer from its
donor phase with a fractionation factor beta = R_sink/R_donor, from the
equilibrium factors alpha_l = R_liquid/R_vapor and alpha_i = R_ice/R_vapor of
the species (the formulations set by `liquid_vapor_formulation` and
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
start of its stage (`water_transfer_stage_T`, see "Stage temperatures"), on the
levels where the transfers that use it moved water.

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

### Kinetic fractionation

With `water_tracer_transfers::kinetic_fractionation: true` (off by default),
evaporation of liquid into subsaturated air and deposition on ice from
supersaturated air use effective fractionation factors, ports of iCAM's
`wiso_akel` and `wiso_akci` (`share/util/water_isotopes.F90`), in the same
R_condensed/R_vapor convention, with D_i/D the species' molecular diffusivity
ratio (`diffusivity_formulation`):

* liquid evaporating at relative humidity RH over liquid (Stewart 1975):
  akel = alpha_l*h/(alpha_l*(D/D_i)^0.58*(h - 1) + 1), with the effective
  humidity h = min(1, 0.25*min(1, RH) + 0.75). akel >= alpha_l, equal at
  RH >= 1, and grows as RH falls: the evaporate is more depleted than at
  equilibrium.
* vapor depositing on ice at supersaturation ratio S over ice (Jouzel &
  Merlivat 1984): akci = alpha_i*S/(alpha_i*(D/D_i)*(S - 1) + 1), with S
  clipped to [1, 2], **only below 253.15 K**: a hard step, as in the Fortran
  (akci = alpha_i at and above 253.15 K). akci <= alpha_i and falls as S
  grows: cold ice takes in less of the heavy isotopes than at equilibrium.
  S is set by `kinetic_ice_supersaturation`: `empirical` (the Fortran default,
  `wiso_ssatf`: S = min(max(1 - 0.002*(T - 273.16), 1), 2), 1.04 at 253.15 K,
  1.15 at 200 K) or `model` (P3's RH over ice, see below).

Every species uses its own diffusivity ratio, as in the Fortran (`difrm` of
`water_isotopes.F90`), here and in the ocean's kinetic factor:

| `diffusivity_formulation` | HDO | H218O | H217O | HTO |
|---|---|---|---|---|
| `merlivat_1978` (default; the Fortran's active row, Merlivat 1978 + Schoenemann et al. 2014) | 0.9757 | 0.9727 | 0.9858 | 0.9679 |
| `cappa_2003` (the Fortran's Cappa et al. 2003 + Merlivat 1978 row) | 0.9839 | 0.9691 | 0.9856 | 0.9679 |

(HDO's is per hydrogen atom, D/H.) Only the equilibrium factors of H217O and
HTO derive from those of H218O and HDO (exponents 0.529 and 2); their kinetic
part uses their own ratios. H217O's 0.9858 is 0.9727^0.518, the diffusion
exponent of Barkan & Luz (2007), so kinetic fractionation produces 17O-excess.

Where they apply (isotopes only; everything else in the table above is
unchanged):

| Transfer | beta with kinetic fractionation | Humidity |
|---|---|---|
| SHOC net evaporation (cloud liquid to vapor) | 1/akel | RH over liquid of SHOC's vapor after transport, at `ShocAfterTransport` |
| SHOC net condensation | alpha_l (equilibrium) | |
| P3 rain evaporation | 1/akel | RH over liquid of P3's vapor after its part1 clips, at `P3Part2` |
| P3 vapor deposition and ice nucleation | akci | S over ice: `empirical` (of `P3Part2`), or `model` (P3's vapor after its part1 clips) |
| P3 Bergeron process | akci/akel | both of the above (the liquid evaporates at RH < 1, the vapor deposits at S > 1) |

The relative humidities are those of EAMxx's RH diagnostic (Murphy & Koop
2005 saturation vapor pressure at `p_dry_mid`, on the wet basis of the vapor),
which is also P3's own (P3 forms it on its dry basis). The process then also
requires `p_dry_mid` and `pseudo_density_dry`. The humidity is that of the
state before the transfer (P3's state before its process rates, SHOC's after
its transport), not after it: P3's limiters bring the vapor close to saturation
by the end of its call.

The kinetic factor of cloud evaporation rests on two proxies. SHOC exports only
its net condensation per level, so the sign of the net transfer stands in for
"evaporation that no condensation pairs with": net condensation stays at
equilibrium and net evaporation is taken as kinetic. And SHOC gives no humidity
for the sub-grid part of the level where the cloud evaporates, so the
grid-mean RH stands in for it. Treat the cloud-evaporation kinetics as the
first suspect if results look wrong.

### Ocean evaporation

The ocean part of an isotope's surface flux carries the equilibrium factor
alpha_l at the SST and the kinetic factor alpkn of evaporation into the boundary
layer (a port of iCAM's `wiso_kmol`: Brutsaert (1975) resistances with the
Merlivat & Jouzel (1979) kinetic modifier, a smooth or rough regime from the
surface Reynolds number). Two options of `water_tracers` control alpkn:

* `kinetic_fractionation_formulation`: `brutsaert_merlivat_jouzel_1979`
  (default) or `none` (alpkn = 1);
* `diffusivity_formulation`: the molecular diffusivity ratios it uses, one
  per species, `merlivat_1978` (default) or `cappa_2003` (values in the table
  above).

alpkn lies in (0, 1] for ustar times the bottom-level height above about
1e-5 m2/s; there is no low-wind guard, as in the Fortran (debug builds check
the range). The land and sea-ice part never fractionates.

### Not represented

Isotopic exchange at zero net mass change: existing cloud liquid and rain do
not re-equilibrate with the vapor (only the water SHOC or P3 move between
phases carries a fractionation), with or without kinetic fractionation.
Rain-vapor and cloud-vapor exchange are not implemented yet. SHOC exports its
net condensation per level, so condensation and evaporation within a level in
the same call are not seen separately.

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
atmosphere driver copies the tracers file and the six isotope options to
`water_tracer_transfers` (an error if the latter sets one to a different value)
and checks that:

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
* Equilibrium fractionation in the transfers listed above, and optionally
  kinetic fractionation of net evaporation and cold deposition (see "Not
  represented" for what is missing). The surface flux of land and sea ice uses
  `r_land` (default `r_ocn`); the coupler's land isotope fluxes are not
  imported.
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
bit. The `*_iso_kinetic*` runs turn `kinetic_fractionation` on: the same
budget and positivity checks; in one SHOC call from R = 1, cloud that SHOC's
net transfer evaporates must keep the Rayleigh result for 1/akel (computed in
the check script); ice deposited from vapor must be near akci; and against the
run with kinetic fractionation off, rain evaporating below cloud must be more
enriched with a lower deuterium excess, and cold deposited ice less enriched.
The unit test `water_tracer_transfers_functions` covers the Rayleigh step
itself, and `water_isotopes_kinetic` compares akel and akci with reference
values of the Fortran functions.

The driver's copy of the isotope options is tested by `shoc_cld_p3_wtt_link_iso`
(the options set only on `water_tracers`: `liquid_vapor_formulation:
majoube_1971`, `ice_vapor_formulation: isocam3`, `diffusivity_formulation:
cappa_2003`, `kinetic_fractionation_formulation: none`, with
`kinetic_fractionation: true` on the transfers) and
`shoc_cld_p3_wtt_link_iso_manual` (the same, set by hand on both processes). The
two outputs must be bit for bit equal. A missing copy of the first three options
would change the replay's results; the other three options
(`standard_ratio_formulation`, `ocean_enrichment_formulation`,
`kinetic_fractionation_formulation`) do not enter the replay's results, so a
missing copy of those shows in no output. `wtt_link_neg_iso_*` check that an
option set to different values on the two processes stops the run.
