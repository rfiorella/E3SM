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

## Fields

Each tracer has one field per bulk water phase, named `<base>_{qv,qc,qr,qi}`
(kg tracer per kg moist air, like the bulk fields), for example `wtrc_t_qv`.
The tracer vapor is advected by HOMME and diffused by SHOC. The other phases
only change through `water_tracer_transfers`.

## The two processes

* `water_tracers` registers the fields, sets their initial conditions and
  computes the ocean evaporation flux. That flux is diagnostic only: SHOC does
  not apply a surface flux to the tracers, so runs with surface evaporation
  do not conserve the tracer mass yet.
* `water_tracer_transfers` replays what SHOC (condensation/evaporation) and P3
  (all process rates, sedimentation, clipping) did to the bulk water on the
  tracers, in the order the bulk model applied it. It also computes
  `<base>_precip_{liq,ice}_surf_mass` (kg/m2, accumulated like the bulk
  precipitation) and the diagnostics `<base>_residual_{shoc,p3}`
  (bulk water minus replayed tracer water, per phase).

Parameters of `water_tracer_transfers`:

* `replay_shoc`, `replay_p3` (default true): which sources to replay.
* `water_tracers_file`: copied from `water_tracers` when unset.

## Enabling

Neither process is in a default `atm_procs_list`. Put `water_tracers` before
`homme` and `water_tracer_transfers` directly after `p3`, inside the
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
* no process between SHOC/P3 and the replay updates `qv`, `qc`, `qr`, `qi` or the
  tracer groups.

A warning is logged if `zm`, `gw`, `nudging` or `iop_forcing` is in the tree.

## Limits

* Only SHOC and P3 are replayed. Other bulk water changes (ZM, GW, nudging,
  IOP forcing, the SHOC flux-consistency fixer) are not followed.
* No fractionation yet, and no surface evaporation applied to the tracers.
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
(a test-only dummy process, since no linked physics does).
