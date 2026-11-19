# Water tracers

Water tracers follow a labeled portion of the model's water through the
atmosphere. They are built from two processes, `water_tracers` and
`water_tracer_transfers`. Only vapor transport (HOMME advection and SHOC
diffusion) and the bulk water phase changes of SHOC and P3 are covered so far.
There is no fractionation yet: every transfer moves tracer mass in proportion
to the donor phase's tracer-to-bulk ratio (alpha = 1).

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
  tracers `water_transfer_stage_T` (next section).

Parameters of `water_tracer_transfers`:

* `replay_shoc`, `replay_p3` (default true): which sources to replay.
* `water_tracers_file`: copied from `water_tracers` when unset.

## Stage temperatures

`water_transfer_stage_T` (K; layout column, `water_transfer_stage`, level) is
the temperature at the start of each replayed stage, for the fractionation
factors to come (nothing uses it yet). Its three components are, in order:

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
* No fractionation yet. The surface flux of land and sea ice uses `r_land`
  (default `r_ocn`); the coupler's land isotope fluxes are not imported.
* A tracer with a uniform ratio R keeps R through every transfer; the tests
  check R = 1 against the bulk fields to roundoff.

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
