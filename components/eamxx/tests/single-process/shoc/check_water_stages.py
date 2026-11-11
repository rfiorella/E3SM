#!/usr/bin/env python3

"""
Check SHOC's export_water_transfers field, shoc_water_stages.

Bookkeeping (default): in every output slice (each follows a SHOC call), the
stages must reproduce the qv and qc SHOC returned:

  qc = LiquidAfterCondensation                                  (exactly)
  qv = VaporAfterTransport + LiquidAfterTransport - LiquidAfterCondensation
                                                                (to roundoff)

where qv is clipped to [0, 0.2] like SHOC's qv postcondition repair.

Transport (--tracer NAME): a water tracer with R = 1 in every phase is in the
same implicit solve as bulk water, so after the first SHOC call (slice 1;
slice 0 is the initial state) its liquid must equal LiquidAfterTransport and
its vapor VaporAfterTransport. Its rain and ice must equal bulk qr and qi.
Later slices are not checked: SHOC then resets bulk qc to its condensed
value, which no process applies to the tracer yet.

The stage components are read from the WaterStage enum in
shoc_water_transfers.hpp, so this script cannot drift from SHOC.
"""

import argparse, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3] / "scripts"))
from utils import expect, ensure_netcdf4  # pylint: disable=wrong-import-position

ensure_netcdf4()

from netCDF4 import Dataset  # pylint: disable=wrong-import-position
import numpy as np            # pylint: disable=wrong-import-position

QV_MIN, QV_MAX = 0.0, 0.2  # kg/kg, SHOC's repairable qv postcondition

###############################################################################
def read_enum(header, struct_name):
###############################################################################
    text = pathlib.Path(header).read_text()
    m = re.search(r"struct\s+" + struct_name + r"\s*\{\s*enum\s*:\s*int\s*\{(.*?)\}", text, re.S)
    expect(m is not None, f"Could not find struct {struct_name} in {header}")
    body = re.sub(r"//[^\n]*", "", m.group(1))
    names = [n.strip() for n in body.split(",") if n.strip()]
    expect(names[-1] == "Count", f"{struct_name} must end with Count")
    return names[:-1]

###############################################################################
def read_stages(ds, t, names):
###############################################################################
    var = ds["shoc_water_stages"]
    # (time, ncol, <component>, lev)
    expect(len(var.dimensions) == 4 and var.dimensions[0] == "time" and var.dimensions[-1] == "lev",
           f"Unexpected dims for shoc_water_stages: {var.dimensions}")
    expect(var.shape[2] == len(names),
           f"shoc_water_stages has {var.shape[2]} components, header has {len(names)}")
    data = var[t]
    return {n: np.asarray(data[:, i, :]) for i, n in enumerate(names)}

###############################################################################
def report_exact(label, a, b):
###############################################################################
    diff = np.abs(a - b)
    nbad = int(np.sum(diff != 0))
    print(f"  {label}: max |diff| = {np.max(diff):.3e} kg/kg, points differing = {nbad}")
    return nbad == 0

###############################################################################
def report_close(label, a, b, scale, rtol):
###############################################################################
    diff = np.abs(a - b)
    tot = scale + 1e-300
    nbad = int(np.sum(diff > rtol * tot))
    print(f"  {label}: max |diff|/scale = {np.max(diff / tot):.3e}, "
          f"max |diff| = {np.max(diff):.3e} kg/kg, points over tol = {nbad}")
    return nbad == 0

###############################################################################
def check_condensation(stage):
###############################################################################
    # A check that passes because SHOC condensed nothing proves nothing
    cond = stage["LiquidAfterCondensation"] - stage["LiquidAfterTransport"]
    ncond, nevap = int(np.sum(cond > 0)), int(np.sum(cond < 0))
    print(f"  condensation points = {ncond} (max {np.max(cond):.3e} kg/kg), "
          f"evaporation points = {nevap} (min {np.min(cond):.3e} kg/kg)")
    ok = ncond > 0 and nevap > 0
    if not ok:
        print("  ERROR: this test needs both condensation and evaporation somewhere")
    return ok

###############################################################################
def bookkeeping(ds, names, rtol):
###############################################################################
    success = True
    for t in range(ds.dimensions["time"].size):
        print(f" t={t}")
        stage = read_stages(ds, t, names)
        vat = stage["VaporAfterTransport"]
        lat = stage["LiquidAfterTransport"]
        lac = stage["LiquidAfterCondensation"]
        qv = np.asarray(ds["qv"][t])
        qc = np.asarray(ds["qc"][t])

        success &= report_exact("qc vs LiquidAfterCondensation", qc, lac)

        pred = vat + lat - lac
        clipped = (pred < QV_MIN) | (pred > QV_MAX)
        print(f"  points where the qv postcondition repair clips: {int(np.sum(clipped))}")
        scale = np.abs(vat) + np.abs(lat) + np.abs(lac)
        success &= report_close("qv vs VaporAfterTransport+LiquidAfterTransport-LiquidAfterCondensation",
                                qv, np.clip(pred, QV_MIN, QV_MAX), scale, rtol)
        success &= check_condensation(stage)
    return success

###############################################################################
def transport(ds, names, tracer, rtol):
###############################################################################
    expect(ds.dimensions["time"].size >= 2, "Need the t0 slice and one after the first SHOC call")
    tr = {p: f"{tracer}_{p}" for p in ["qv", "qc", "qr", "qi"]}

    success = True
    print(" t=0 (initial state, init: ratio with r_init = 1)")
    for p, name in tr.items():
        success &= report_exact(f"{name} vs {p}", np.asarray(ds[name][0]), np.asarray(ds[p][0]))

    print(" t=1 (after the first SHOC call)")
    stage = read_stages(ds, 1, names)
    vat = stage["VaporAfterTransport"]
    lat = stage["LiquidAfterTransport"]
    qv0, qc0 = np.asarray(ds["qv"][0]), np.asarray(ds["qc"][0])

    # Same solve, same input, same column of the right-hand side
    success &= report_exact(f"{tr['qc']} vs LiquidAfterTransport", np.asarray(ds[tr["qc"]][1]), lat)
    for p in ["qr", "qi"]:
        success &= report_exact(f"{tr[p]} vs {p}", np.asarray(ds[tr[p]][1]), np.asarray(ds[p][1]))
    # The bulk solve diffuses qw = qv + qc, the tracer qv alone: equal by linearity
    success &= report_close(f"{tr['qv']} vs VaporAfterTransport", np.asarray(ds[tr["qv"]][1]), vat,
                            np.abs(vat) + np.abs(qv0) + np.abs(qc0), rtol)

    # A check that passes because SHOC moved nothing proves nothing
    dqv, dqc = np.max(np.abs(vat - qv0)), np.max(np.abs(lat - qc0))
    print(f"  max |transport| of qv = {dqv:.3e} kg/kg, of qc = {dqc:.3e} kg/kg")
    if not (dqv > 0 and dqc > 0):
        print("  ERROR: this test needs nonzero turbulent transport of qv and qc")
        success = False
    success &= check_condensation(stage)

    # Not checked (see above); printed to show the size of what a tracer
    # process must still replay
    for t in range(2, ds.dimensions["time"].size):
        dq = np.max(np.abs(np.asarray(ds[tr["qc"]][t]) - np.asarray(ds["qc"][t])))
        print(f" t={t} (not checked) max |{tr['qc']} - qc| = {dq:.3e} kg/kg")
    return success

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True, help="output file with shoc_water_stages")
    parser.add_argument("--header", required=True, help="path to shoc_water_transfers.hpp")
    parser.add_argument("--tracer", default=None,
                        help="water tracer base name: run the transport check instead")
    parser.add_argument("--rtol", type=float, default=1e-12,
                        help="tolerance for the roundoff checks, relative to the magnitude of the terms")
    args = parser.parse_args()

    names = read_enum(args.header, "WaterStage")
    expect(names == ["VaporAfterTransport", "LiquidAfterTransport", "LiquidAfterCondensation"],
           f"Unexpected WaterStage {names}")

    ds = Dataset(args.file, "r")
    if args.tracer is None:
        success = bookkeeping(ds, names, args.rtol)
        what = "bookkeeping"
    else:
        success = transport(ds, names, args.tracer, args.rtol)
        what = "transport"

    print(f"SHOC water stages {what} check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
