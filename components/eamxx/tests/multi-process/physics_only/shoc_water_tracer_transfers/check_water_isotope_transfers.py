#!/usr/bin/env python3

"""
Check the water_tracer_transfers replay with fractionating isotope tracers.

Isotopes do not keep a uniform ratio to bulk water, so the checks of
check_water_tracer_transfers.py do not apply. Instead, in every output slice
after the first (one atm step each):

  - Conservation. For bulk water, for the reference tracer (--reference, a
    generic R = 1 tracer) and for each isotope (--tracer), per column,
      M(t) - M(t-1) + P(t) - F(t)*dt
    with M = sum over levels and phases of q*dp/g [kg/m2], P the surface
    precipitation accumulated over the step (with --p3-header) and F the
    surface flux (surf_evap for bulk water, <name>_qv_surf_flux for tracers),
    constant over the step. The replay moves tracer only between phases, so
    this is roundoff, apart from the vapor repairs of the checkpoints (which
    move tracer with the bulk water's repairs, so they would show in bulk
    water's budget too). Relative to M(t-1) it must be below --rtol-budget plus
    bulk water's.
  - Positivity: every phase of every isotope is >= -(--neg-rtol)*qtot, qtot
    the bulk total water at the point. Donors that fractionate cannot go below
    0 (the Rayleigh step and its guard), but the transfers that move tracer at
    the donor's ratio (e.g. ice's sinks) and the sedimentation sweep leave
    roundoff of either sign where the bulk water goes to 0, as for generic
    tracers.

With --alpha-check (the run uses the equilibrium alphas, not alpha_override),
the physics of the fractionation is checked too, against alpha_l and alpha_i
computed here from their published fits, independently of the code (the
default formulations: Horita & Wesolowski 1994 for liquid; Merlivat & Nief
1967 for HDO and Majoube 1971 for H218O for ice):
  - SHOC (--stages-header, and SHOC is the last process to change the
    tracers, i.e. no --p3-header): where SHOC's last call made new cloud from
    vapor (bulk cloud liquid after transport below 1e-8 of the net
    condensation), R_qc/R_qv must be the result of one Rayleigh step for
    alpha_l at the stage temperature ShocAfterTransport, to --alpha-rtol
    (to first order in the condensed fraction that is alpha_l). Over the run
    there must be at least --min-new-cloud such points.
  - SHOC: where it condensed (net condensation > 0, qc > 1e-8), the cloud
    must be enriched relative to the vapor, R_qc > R_qv, at more than 95% of
    the points.
  - P3 (--p3-header): where P3's last call deposited (and nucleated) more ice
    from vapor than half of the ice it started with plus the ice it gained
    otherwise, the normalized enrichment e = (R_qi/R_qv - 1)/(alpha_i - 1) at
    P3's part2 temperature must be near 1 (equilibrium deposition): median
    within 0.1 of 1, every point in [0.5, 1.5], at least 10 points over the
    run. Not exact, since the ice also holds older ice and ice sedimented
    from above (deposited at other temperatures). Without fractionation e
    would be 0, with an inverted factor about -1.

--tracer takes NAME=SPECIES (HDO or H218O).
"""

import argparse, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4] / "scripts"))
from utils import expect, ensure_netcdf4  # pylint: disable=wrong-import-position

ensure_netcdf4()

from netCDF4 import Dataset  # pylint: disable=wrong-import-position
import numpy as np            # pylint: disable=wrong-import-position

PHASES = ["qv", "qc", "qr", "qi"]
GRAVIT = 9.80616  # m/s2, physics::Constants::gravit

###############################################################################
def alpha_liquid(species, T):
###############################################################################
    # R_liquid/R_vapor, Horita & Wesolowski (1994), their fits of 10^3 ln(alpha)
    # for D/H and 18O/16O (T in K)
    if species == "HDO":
        ln = 1158.8e-9*T**3 - 1620.1e-6*T**2 + 794.84e-3*T - 161.04 + 2.9992e9/T**3
    else:
        ln = -7.685 + 6.7123e3/T - 1.6664e6/T**2 + 0.35041e9/T**3
    return np.exp(ln * 1e-3)

###############################################################################
def alpha_ice(species, T):
###############################################################################
    # R_ice/R_vapor: Merlivat & Nief (1967) for D/H, Majoube (1971) for 18O/16O
    if species == "HDO":
        ln = 16289e3/T**2 - 94.5
    else:
        ln = 11.839e3/T - 28.224
    return np.exp(ln * 1e-3)

###############################################################################
def read_enum(header, enum_re):
###############################################################################
    text = pathlib.Path(header).read_text()
    m = re.search(enum_re + r"\s*\{(.*?)\}", text, re.S)
    expect(m is not None, f"Could not find '{enum_re}' in {header}")
    body = re.sub(r"//[^\n]*", "", m.group(1))
    names = [n.strip() for n in body.split(",") if n.strip()]
    expect(names[-1] == "Count", f"Enum '{enum_re}' in {header} must end with Count")
    return names[:-1]

###############################################################################
def components(ds, name, t, names):
###############################################################################
    var = ds[name]
    expect(len(var.dimensions) == 4 and var.dimensions[0] == "time" and var.dimensions[-1] == "lev",
           f"Unexpected dims for {name}: {var.dimensions}")
    expect(var.shape[2] == len(names), f"{name} has {var.shape[2]} components, header has {len(names)}")
    data = np.asarray(var[t])
    return {n: data[:, i, :] for i, n in enumerate(names)}

###############################################################################
def column_mass(ds, t, prefix):
###############################################################################
    # sum over levels and phases of q*dp/g [kg/m2], per column
    dp = np.asarray(ds["pseudo_density"][t])
    q = sum(np.asarray(ds[f"{prefix}{p}"][t]) for p in PHASES)
    return np.sum(q * dp, axis=1) / GRAVIT

###############################################################################
def budget(ds, t, prefix, flux_name, with_precip, dt):
###############################################################################
    m0, m1 = column_mass(ds, t-1, prefix), column_mass(ds, t, prefix)
    res = m1 - m0 - np.asarray(ds[flux_name][t]) * dt
    if with_precip:
        res += np.asarray(ds[f"{prefix}precip_liq_surf_mass"][t]) + np.asarray(ds[f"{prefix}precip_ice_surf_mass"][t])
    return res, m0

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True, help="output file of the test")
    parser.add_argument("--tracer", nargs="+", required=True,
                        help="isotope tracers as NAME=SPECIES (HDO or H218O)")
    parser.add_argument("--reference", required=True, help="a generic tracer with R = 1")
    parser.add_argument("--stages-header", help="path to shoc_water_transfers.hpp (SHOC is replayed)")
    parser.add_argument("--p3-header", help="path to p3_water_transfers.hpp (P3 is replayed)")
    parser.add_argument("--stage-T-header", required=True,
                        help="path to eamxx_water_tracer_transfers_functions.hpp (TransferStage)")
    parser.add_argument("--dt", type=float, required=True, help="atm time step [s]")
    parser.add_argument("--dt-p3", type=float, help="dt of P3's last call of a step [s] (with --p3-header)")
    parser.add_argument("--rtol-budget", type=float, default=1e-12,
                        help="tolerance of the column budgets, relative to the column's water mass")
    parser.add_argument("--neg-rtol", type=float, default=1e-12,
                        help="tolerance of the positivity check, relative to the point's total water")
    parser.add_argument("--alpha-check", action="store_true",
                        help="the run uses the equilibrium alphas: check the fractionation's physics")
    parser.add_argument("--alpha-rtol", type=float, default=1e-10,
                        help="tolerance of R_qc/R_qv against the Rayleigh result for new SHOC cloud")
    parser.add_argument("--min-new-cloud", type=int, default=1,
                        help="new SHOC cloud points needed over the run (with --alpha-check, SHOC only)")
    args = parser.parse_args()
    expect(args.stages_header or args.p3_header, "Need --stages-header and/or --p3-header")
    expect(not (args.alpha_check and args.p3_header) or args.dt_p3, "The P3 alpha check needs --dt-p3")

    tracers = {}
    for entry in args.tracer:
        name, species = entry.split("=")
        expect(species in ["HDO", "H218O"], f"Unsupported species {species}")
        tracers[name] = species

    stage_T_names = read_enum(args.stage_T_header, r"struct\s+TransferStage\s*\{\s*enum\s*:\s*int")
    if args.stages_header:
        stages = read_enum(args.stages_header, r"struct\s+WaterStage\s*\{\s*enum\s*:\s*int")
    if args.p3_header:
        transfers = read_enum(args.p3_header, r"struct\s+WaterTransfer\s*\{\s*enum\s*:\s*int")
        entry_names = read_enum(args.p3_header, r"struct\s+WaterEntry\s*\{\s*enum\s*:\s*int")

    ds = Dataset(args.file, "r")
    nt = ds.dimensions["time"].size
    expect(nt >= 2, "Need the t0 slice and at least one step")

    success = True
    # Over the run, per tracer: new SHOC cloud points, and e of deposited ice
    new_cloud = {name: 0 for name in tracers}
    e_ice = {name: [] for name in tracers}
    for t in range(1, nt):
        print(f" t={t}")
        q = {p: np.asarray(ds[p][t]) for p in PHASES}
        qtot = sum(q.values())

        # ---- Conservation ----
        with_precip = args.p3_header is not None
        res_b, m_b = budget(ds, t, "", "surf_evap", with_precip, args.dt)
        rel_b = np.abs(res_b) / m_b
        print(f"  bulk water: max |budget residual|/M = {np.max(rel_b):.3e}")
        for name in [args.reference] + list(tracers):
            res, m = budget(ds, t, f"{name}_", f"{name}_qv_surf_flux", with_precip, args.dt)
            rel = np.abs(res) / m
            nbad = int(np.sum(rel > args.rtol_budget + rel_b))
            print(f"  {name}: max |budget residual|/M = {np.max(rel):.3e}, columns over tol = {nbad}")
            success &= nbad == 0

        # ---- Positivity (the reference tracer for comparison) ----
        for name in list(tracers) + [args.reference]:
            for p in PHASES:
                x = np.asarray(ds[f"{name}_{p}"][t])
                nneg = int(np.sum(x < 0))
                nbad = int(np.sum(x < -args.neg_rtol * qtot))
                if nneg > 0 or name in tracers:
                    print(f"  {name}_{p}: min = {np.min(x):.3e}, points < 0 = {nneg}, "
                          f"min/qtot = {np.min(x / qtot):.3e}, points below tol = {nbad}")
                if name in tracers:
                    success &= nbad == 0

        if "water_transfer_alpha_extrapolations" in ds.variables:
            n = np.asarray(ds["water_transfer_alpha_extrapolations"][t])
            print(f"  alpha evaluations extrapolating a fit over the step: {int(np.sum(n))} "
                  f"(in {int(np.sum(n > 0))} columns)")

        if not args.alpha_check:
            continue
        T_stage = components(ds, "water_transfer_stage_T", t, stage_T_names)

        # ---- SHOC: condensation ----
        if args.stages_header:
            st = components(ds, "shoc_water_stages", t, stages)
            qv_tr, qc_tr = st["VaporAfterTransport"], st["LiquidAfterTransport"]
            cond = st["LiquidAfterCondensation"] - qc_tr
            new = (cond > 0) & (qc_tr <= 1e-8 * cond)
            cloudy = (cond > 0) & (q["qc"] > 1e-8)
            T = T_stage["ShocAfterTransport"]
            for name, sp in tracers.items():
                R_v = np.asarray(ds[f"{name}_qv"][t]) / q["qv"]
                R_c = np.asarray(ds[f"{name}_qc"][t]) / np.where(q["qc"] > 0, q["qc"], 1)
                if args.p3_header is None and np.any(new):
                    # One Rayleigh step removing x = cond/qv of the vapor
                    # removes F = 1 - (1-x)^alpha of its tracer, so the cloud's
                    # ratio over the vapor's is F*(1-x)/(x*(1-F)) (log1p/expm1:
                    # x is as small as 1e-12)
                    alpha = alpha_liquid(sp, T[new])
                    x = cond[new] / qv_tr[new]
                    F = -np.expm1(alpha * np.log1p(-x))
                    expected = F * (1 - x) / (x * (1 - F))
                    ratio = R_c[new] / R_v[new]
                    err = np.abs(ratio / expected - 1)
                    nbad = int(np.sum(err > args.alpha_rtol))
                    print(f"  {name} new SHOC cloud ({int(np.sum(new))} points): "
                          f"max |(R_qc/R_qv)/Rayleigh - 1| = {np.max(err):.3e}, points over tol = {nbad}; "
                          f"max |(R_qc/R_qv)/alpha_l - 1| = {np.max(np.abs(ratio/alpha - 1)):.3e}, "
                          f"alpha_l - 1 in [{np.min(alpha)-1:.3e}, {np.max(alpha)-1:.3e}], "
                          f"x in [{np.min(x):.1e}, {np.max(x):.1e}]")
                    success &= nbad == 0
                    new_cloud[name] += int(np.sum(new))
                frac = float(np.mean(R_c[cloudy] > R_v[cloudy])) if np.any(cloudy) else 0.0
                print(f"  {name} condensing cloudy points ({int(np.sum(cloudy))}): "
                      f"fraction with R_qc > R_qv = {frac:.4f}")
                success &= frac > 0.95

        # ---- P3: ice deposited from vapor ----
        if args.p3_header:
            rate = components(ds, "p3_water_transfers", t, transfers)
            entry = components(ds, "p3_water_entry", t, entry_names)
            dep = (rate["Part2VapdepV2I"] + rate["Part2NucleatV2I"]) * args.dt_p3
            other = sum(rate[n] for n in ["Part2BergC2I", "Part2CollectC2I", "Part2FreezeC2I",
                                          "Part2CollectR2I", "Part2FreezeR2I",
                                          "HomFreezeC2I", "HomFreezeR2I"]) * args.dt_p3
            dom = (dep > 0.5 * (entry["Ice"] + other)) & (q["qi"] > 1e-10)
            T = T_stage["P3Part2"][dom]
            for name, sp in tracers.items():
                R_i = np.asarray(ds[f"{name}_qi"][t])[dom] / q["qi"][dom]
                R_v = np.asarray(ds[f"{name}_qv"][t])[dom] / q["qv"][dom]
                e_ice[name] += list((R_i / R_v - 1) / (alpha_ice(sp, T) - 1))

    if args.alpha_check:
        print(" over the run")
        for name in tracers:
            if args.stages_header and args.p3_header is None:
                print(f"  {name}: new SHOC cloud points = {new_cloud[name]} (need >= {args.min_new_cloud})")
                success &= new_cloud[name] >= args.min_new_cloud
            if args.p3_header:
                e = np.array(e_ice[name])
                ok = len(e) >= 10 and abs(np.median(e) - 1) < 0.1 and np.all((e > 0.5) & (e < 1.5))
                print(f"  {name}: deposition-dominated ice points = {len(e)}, "
                      f"e = (R_qi/R_qv - 1)/(alpha_i - 1): "
                      + (f"min {np.min(e):.3f}, median {np.median(e):.3f}, max {np.max(e):.3f}" if len(e) else "none")
                      + ("" if ok else "  <-- needs >= 10 points, |median - 1| < 0.1, all in [0.5, 1.5]"))
                success &= ok

    print("Water isotope transfers check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
