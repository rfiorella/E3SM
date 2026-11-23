#!/usr/bin/env python3

"""
Check the SHOC replay of the water_tracer_transfers process.

Every tracer given as NAME=R starts with ratio R to bulk water in every phase
(init: ratio). SHOC diffuses tracers and bulk water in the same solve, and the
replay moves the tracer at its donor phase's ratio, so a uniform ratio is
preserved. In every output slice, each phase of each tracer must be R times
its bulk phase:

  NAME_qv = R*qv, NAME_qc = R*qc   (to roundoff, relative to the total water)
  NAME_qr = R*qr, NAME_qi = R*qi   (SHOC only transports them: exactly for
                                    R = 1, same solve on the same data; to
                                    roundoff otherwise, since R*q can be
                                    subnormal, where scaling is not exact)

and, after the first step, the NAME_residual_shoc field (bulk water minus
tracer after the replay, per phase) must equal the bulk-minus-tracer
difference in the output, to roundoff. For R = 1 that residual is the closure
error, so it must itself be roundoff.

The checks prove nothing if SHOC moved nothing, so every slice after the
first must show condensation and evaporation (shoc_water_stages) and turbulent
transport (a change of qv+qc, qr and qi since the previous slice; with no
surface fluxes only transport changes them).

Enum orders are read from shoc_water_transfers.hpp (WaterStage) and
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
def check_nonvacuous(ds, t, stage):
###############################################################################
    # Condensation and evaporation in SHOC's last call of the step
    cond = stage["LiquidAfterCondensation"] - stage["LiquidAfterTransport"]
    ncond, nevap = int(np.sum(cond > 0)), int(np.sum(cond < 0))
    print(f"    condensation points = {ncond} (max {np.max(cond):.3e} kg/kg), "
          f"evaporation points = {nevap} (min {np.min(cond):.3e} kg/kg)")
    ok = ncond > 0 and nevap > 0

    # Transport over the step: with no surface fluxes, SHOC changes qv+qc,
    # qr and qi only by transport
    q = {p: np.asarray(ds[p][t]) for p in BULK.values()}
    q0 = {p: np.asarray(ds[p][t-1]) for p in BULK.values()}
    dqt = np.max(np.abs((q["qv"] + q["qc"]) - (q0["qv"] + q0["qc"])))
    dqr = np.max(np.abs(q["qr"] - q0["qr"]))
    dqi = np.max(np.abs(q["qi"] - q0["qi"]))
    print(f"    max |transport| over the step of qv+qc = {dqt:.3e}, qr = {dqr:.3e}, qi = {dqi:.3e} kg/kg")
    ok &= dqt > 0 and dqr > 0 and dqi > 0
    if not ok:
        print("    ERROR: this test needs condensation, evaporation and transport")
    return ok

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True, help="output file of the test")
    parser.add_argument("--stages-header", required=True, help="path to shoc_water_transfers.hpp")
    parser.add_argument("--phases-header", required=True, help="path to eamxx_water_tracers_registry.hpp")
    parser.add_argument("--tracer", nargs="+", required=True,
                        help="water tracers as NAME=R, R being the tracer's uniform ratio to bulk water")
    parser.add_argument("--rtol", type=float, default=1e-12,
                        help="tolerance of the roundoff checks, relative to the total water at a point")
    args = parser.parse_args()

    stages = read_enum(args.stages_header, r"struct\s+WaterStage\s*\{\s*enum\s*:\s*int")
    phases = read_enum(args.phases_header, r"enum\s+class\s+WaterPhase")
    expect(sorted(phases) == sorted(BULK), f"Unexpected WaterPhase {phases}")

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

        if t > 0:
            success &= check_nonvacuous(ds, t, components(ds, "shoc_water_stages", t, stages))

        for name, ratio in tracers.items():
            print(f"  {name} (R = {ratio})")
            tr = {p: np.asarray(ds[f"{name}_{p}"][t]) for p in BULK.values()}
            for p in BULK.values():
                success &= report(f"{name}_{p} vs {ratio}*{p}", tr[p], ratio*q[p],
                                  ratio*qtot, args.rtol, p in EXACT and ratio == 1)
            if t == 0:
                continue
            # The residual field must be the bulk-minus-tracer difference
            # (and for R = 1, roundoff itself)
            res = components(ds, f"{name}_residual_shoc", t, phases)
            for phase, p in BULK.items():
                success &= report(f"{name}_residual_shoc[{phase}] vs {p}-{name}_{p}", res[phase],
                                  q[p] - tr[p], qtot, args.rtol, False)
                if ratio == 1:
                    success &= report(f"{name}_residual_shoc[{phase}] vs 0", res[phase],
                                      np.zeros_like(qtot), qtot, args.rtol, p in EXACT)

    print("Water tracer transfers SHOC replay check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
