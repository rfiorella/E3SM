#!/usr/bin/env python3

"""
Check the SHOC and P3 replays of the water_tracer_transfers process.

Every tracer given as NAME=R starts with ratio R to bulk water in every phase
(init: ratio). SHOC diffuses tracers and bulk water in the same solve, and the
replay moves the tracer at its donor phase's ratio, so a uniform ratio is
preserved. In every output slice, each phase of each tracer must be R times
its bulk phase:

  NAME_qv = R*qv, NAME_qc = R*qc   (to roundoff, relative to the total water)
  NAME_qr = R*qr, NAME_qi = R*qi   (to roundoff; with --exact-qr-qi, exactly
                                    for R = 1: SHOC alone only transports
                                    them, in the same solve on the same data)

Per replay (--shoc and/or --p3-header), after the first step the
NAME_residual_<proc> field (bulk water minus tracer after that replay, per
phase) must equal the bulk-minus-tracer difference in the output, to roundoff
(SHOC: only if SHOC is the last process; P3: before P3's final vapor
correction, so only the condensate phases). For R = 1 the residual is the
closure error, so it must itself be roundoff.

With --p3-header, the tracer's surface precipitation must also be R times
P3's, NAME_precip_{liq,ice}_surf_mass = R*precip_{liq,ice}_surf_mass, to
roundoff relative to the lowest level's total water per area (see
check_precip; needs pseudo_density in the output).

The checks prove nothing if nothing moved water, so every slice after the
first must show, for SHOC, condensation and evaporation (shoc_water_stages)
and, with SHOC alone, turbulent transport (a change of qv+qc, qr and qi since
the previous slice; with no surface fluxes only transport changes them); for
P3, transfers in at least --min-p3-transfers components of
p3_water_transfers, rain and ice sedimentation, and surface precipitation.

Enum orders are read from shoc_water_transfers.hpp (WaterStage),
p3_water_transfers.hpp (WaterTransfer, SedFlux) and
eamxx_water_tracers_registry.hpp (WaterPhase), so this script cannot drift
from the code.
"""

import argparse, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4] / "scripts"))
from utils import expect, ensure_netcdf4  # pylint: disable=wrong-import-position

ensure_netcdf4()

from netCDF4 import Dataset  # pylint: disable=wrong-import-position
import numpy as np            # pylint: disable=wrong-import-position

BULK = {"Vapor": "qv", "CloudLiquid": "qc", "Rain": "qr", "Ice": "qi"}
EXACT = ["qr", "qi"]
PRECIP = ["precip_liq_surf_mass", "precip_ice_surf_mass"]
GRAVIT = 9.80616  # m/s2, physics::Constants::gravit

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
    # (time, ncol, <component>, lev)
    expect(len(var.dimensions) == 4 and var.dimensions[0] == "time" and var.dimensions[-1] == "lev",
           f"Unexpected dims for {name}: {var.dimensions}")
    expect(var.shape[2] == len(names), f"{name} has {var.shape[2]} components, header has {len(names)}")
    data = np.asarray(var[t])
    return {n: data[:, i, :] for i, n in enumerate(names)}

###############################################################################
def report(label, a, b, scale, rtol, exact):
###############################################################################
    diff = np.abs(a - b)
    rel = diff / (scale + 1e-300)
    if exact:
        nbad = int(np.sum(diff != 0))
        print(f"    {label}: max |diff| = {np.max(diff):.3e} kg/kg, points differing = {nbad}")
    else:
        nbad = int(np.sum(rel > rtol))
        print(f"    {label}: max |diff| = {np.max(diff):.3e} kg/kg, "
              f"max |diff|/qtot = {np.max(rel):.3e}, points over tol = {nbad}")
    return nbad == 0

###############################################################################
def check_precip(ds, t, name, ratio, qtot, rtol):
###############################################################################
    # The tracer's surface precipitation leaves the lowest level at that
    # level's tracer ratio, so its error is the ratio's error times the flux.
    # Where little rain or ice is left of much larger transfers, the ratio is
    # only accurate to roundoff relative to the level's total water, like the
    # phases above. Hence the scale: the lowest level's total water per area,
    # qtot*dp/g [kg/m2] (plus the precipitation itself). The error relative
    # to the precipitation alone is printed too.
    dp = np.asarray(ds["pseudo_density"][t])
    level_mass = qtot[:, -1] * dp[:, -1] / GRAVIT
    ok = True
    for p in PRECIP:
        bulk = ratio * np.asarray(ds[p][t])
        diff = np.abs(np.asarray(ds[f"{name}_{p}"][t]) - bulk)
        rel = diff / (np.abs(bulk) + ratio*level_mass + 1e-300)
        rel_p = np.max(diff[bulk != 0] / np.abs(bulk[bulk != 0])) if np.any(bulk != 0) else 0
        nbad = int(np.sum(rel > rtol))
        print(f"    {name}_{p} vs {ratio}*{p}: max |diff| = {np.max(diff):.3e} kg/m2, "
              f"max |diff|/(|precip| + qtot*dp/g) = {np.max(rel):.3e}, "
              f"max |diff|/|precip| = {rel_p:.3e}, points over tol = {nbad}")
        ok &= nbad == 0
    return ok

###############################################################################
def check_nonvacuous_shoc(ds, t, stage, transport):
###############################################################################
    # Condensation and evaporation in SHOC's last call of the step
    cond = stage["LiquidAfterCondensation"] - stage["LiquidAfterTransport"]
    ncond, nevap = int(np.sum(cond > 0)), int(np.sum(cond < 0))
    print(f"    condensation points = {ncond} (max {np.max(cond):.3e} kg/kg), "
          f"evaporation points = {nevap} (min {np.min(cond):.3e} kg/kg)")
    ok = ncond > 0 and nevap > 0

    # Transport over the step: with no surface fluxes, SHOC changes qv+qc,
    # qr and qi only by transport
    if transport:
        q = {p: np.asarray(ds[p][t]) for p in BULK.values()}
        q0 = {p: np.asarray(ds[p][t-1]) for p in BULK.values()}
        dqt = np.max(np.abs((q["qv"] + q["qc"]) - (q0["qv"] + q0["qc"])))
        dqr = np.max(np.abs(q["qr"] - q0["qr"]))
        dqi = np.max(np.abs(q["qi"] - q0["qi"]))
        print(f"    max |transport| over the step of qv+qc = {dqt:.3e}, qr = {dqr:.3e}, qi = {dqi:.3e} kg/kg")
        ok &= dqt > 0 and dqr > 0 and dqi > 0
    if not ok:
        print("    ERROR: this test needs SHOC condensation, evaporation and transport")
    return ok

###############################################################################
def check_nonvacuous_p3(ds, t, transfers, sed, min_transfers):
###############################################################################
    # P3's last call of the step: transfers, sedimentation, precipitation
    rate = components(ds, "p3_water_transfers", t, transfers)
    active = [n for n in transfers if np.any(rate[n] != 0)]
    print(f"    nonzero P3 transfers: {len(active)}/{len(transfers)}: {', '.join(active)}")
    ok = len(active) >= min_transfers
    flux = components(ds, "p3_sed_flux", t, sed)
    for s in ["Rain", "Ice"]:
        print(f"    max {s} sedimentation flux = {np.max(np.abs(flux[s])):.3e} kg/m2/s")
        ok &= bool(np.any(flux[s] != 0))
    for p in PRECIP:
        print(f"    max {p} = {np.max(np.asarray(ds[p][t])):.3e} kg/m2")
        ok &= bool(np.any(np.asarray(ds[p][t]) > 0))
    if not ok:
        print(f"    ERROR: this test needs >= {min_transfers} P3 transfer components, "
              "rain and ice sedimentation, and surface precipitation")
    return ok

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True, help="output file of the test")
    parser.add_argument("--phases-header", required=True, help="path to eamxx_water_tracers_registry.hpp")
    parser.add_argument("--stages-header",
                        help="path to shoc_water_transfers.hpp: check the SHOC replay")
    parser.add_argument("--p3-header",
                        help="path to p3_water_transfers.hpp: check the P3 replay")
    parser.add_argument("--tracer", nargs="+", required=True,
                        help="water tracers as NAME=R, R being the tracer's uniform ratio to bulk water")
    parser.add_argument("--rtol", type=float, default=1e-12,
                        help="tolerance of the roundoff checks, relative to the total water at a point")
    parser.add_argument("--exact-qr-qi", action="store_true",
                        help="require R = 1 tracers' qr and qi to match bulk bit for bit (SHOC alone)")
    parser.add_argument("--min-p3-transfers", type=int, default=1,
                        help="P3 transfer components that must be nonzero in every step")
    args = parser.parse_args()
    expect(args.stages_header or args.p3_header, "Need --stages-header and/or --p3-header")

    phases = read_enum(args.phases_header, r"enum\s+class\s+WaterPhase")
    expect(sorted(phases) == sorted(BULK), f"Unexpected WaterPhase {phases}")
    if args.stages_header:
        stages = read_enum(args.stages_header, r"struct\s+WaterStage\s*\{\s*enum\s*:\s*int")
    if args.p3_header:
        transfers = read_enum(args.p3_header, r"struct\s+WaterTransfer\s*\{\s*enum\s*:\s*int")
        sed = read_enum(args.p3_header, r"struct\s+SedFlux\s*\{\s*enum\s*:\s*int")
        expect(sed == ["CloudLiquid", "Rain", "Ice"], f"Unexpected SedFlux {sed}")
    exact = EXACT if args.exact_qr_qi else []

    tracers = {}
    for entry in args.tracer:
        name, ratio = entry.split("=")
        tracers[name] = float(ratio)

    ds = Dataset(args.file, "r")
    nt = ds.dimensions["time"].size
    expect(nt >= 2, "Need the t0 slice and at least one step")

    success = True
    for t in range(nt):
        print(f" t={t}" + (" (initial state)" if t == 0 else ""))
        q = {p: np.asarray(ds[p][t]) for p in BULK.values()}
        qtot = q["qv"] + q["qc"] + q["qr"] + q["qi"]

        if t > 0 and args.stages_header:
            success &= check_nonvacuous_shoc(ds, t, components(ds, "shoc_water_stages", t, stages),
                                             transport=not args.p3_header)
        if t > 0 and args.p3_header:
            success &= check_nonvacuous_p3(ds, t, transfers, sed, args.min_p3_transfers)

        for name, ratio in tracers.items():
            print(f"  {name} (R = {ratio})")
            tr = {p: np.asarray(ds[f"{name}_{p}"][t]) for p in BULK.values()}
            for p in BULK.values():
                success &= report(f"{name}_{p} vs {ratio}*{p}", tr[p], ratio*q[p],
                                  ratio*qtot, args.rtol, p in exact and ratio == 1)
            if t == 0:
                continue
            if args.p3_header:
                success &= check_precip(ds, t, name, ratio, qtot, args.rtol)
            # The residual fields must be the bulk-minus-tracer difference
            # (and for R = 1, roundoff themselves). SHOC's only matches the
            # output if SHOC is the last process; P3's is from before P3's
            # final vapor correction, so its vapor is not in the output.
            residuals = []
            if args.stages_header and not args.p3_header:
                residuals.append(("residual_shoc", list(BULK.items())))
            if args.p3_header:
                residuals.append(("residual_p3", [(ph, p) for ph, p in BULK.items() if p != "qv"]))
            for field, checked in residuals:
                res = components(ds, f"{name}_{field}", t, phases)
                for phase, p in checked:
                    success &= report(f"{name}_{field}[{phase}] vs {p}-{name}_{p}", res[phase],
                                      q[p] - tr[p], qtot, args.rtol, False)
            # For R = 1 every residual is the closure error. P3's vapor
            # residual also holds any repair of P3's qv postcondition.
            if ratio == 1:
                names = (["residual_shoc"] if args.stages_header else []) + \
                        (["residual_p3"] if args.p3_header else [])
                for field in names:
                    res = components(ds, f"{name}_{field}", t, phases)
                    for phase, p in BULK.items():
                        success &= report(f"{name}_{field}[{phase}] vs 0", res[phase],
                                          np.zeros_like(qtot), qtot, args.rtol,
                                          p in exact and field == "residual_shoc")

    print("Water tracer transfers replay check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
