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

With --stage-T-header, water_transfer_stage_T (the temperature at the start
of each replayed stage, rebuilt from T_mid and the exports; needs T_mid and,
with P3, pseudo_density_dry in the output) is checked, to --T-rtol relative
to T_mid. Stages that are not replayed must hold the fill value. The
heating of each p3_water_transfers component is taken from its name (C2V:
cloud liquid to vapor, so -LatVap/cp per kg/kg), not from the code's table.
  - P3, independent (--T-from-previous-slice; only valid if nothing but P3
    changes T_mid, as in a P3-only run without subcycling): P3Entry must be
    the previous slice's T_mid. This checks the rebuild against P3's own
    temperature, so a wrong heating table or a heating P3 applies but does
    not export (incl. the clips and homogeneous freezing, which no unit test
    reaches) shows here.
  - P3, bookkeeping: P3Entry = T_mid - sum_k L_k*rate_k*dt/cp, and P3Part2 -
    P3Entry = the part1 clips' heating times pseudo_density/pseudo_density_dry
    (P3's dry basis). Only checks that the code's table matches the names
    and that the stages are assembled as documented.
  - SHOC, bookkeeping only: ShocAfterTransport = T_after_SHOC -
    (LatVap/cp)*(LiquidAfterCondensation - LiquidAfterTransport), with
    T_after_SHOC = P3Entry with P3, else T_mid. No SHOC output holds the
    temperature after transport, so nothing independent is available; the
    physics rests on reading the code: shoc_main sets dse from
    thl/inv_exner + (LatVap/cp)*ql (shoc_update_host_dse_impl.hpp), thl is
    only read by the condensation (shoc_assumed_pdf, shoc_main_impl.hpp), and
    SHOCPostprocess sets T_mid from dse (eamxx_shoc_process_interface.hpp).

Enum orders are read from shoc_water_transfers.hpp (WaterStage),
p3_water_transfers.hpp (WaterTransfer, SedFlux),
eamxx_water_tracers_registry.hpp (WaterPhase) and
eamxx_water_tracer_transfers_functions.hpp (TransferStage), so this script
cannot drift from the code.
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
LATVAP = 2501000.0  # J/kg, physics::Constants::LatVap
LATICE = 333700.0   # J/kg, physics::Constants::LatIce
CP     = 1004.64    # J/kg/K, physics::Constants::CP
# Latent enthalpy relative to vapor [J/kg] of the phase each letter of a
# p3_water_transfers component name stands for (C2V: cloud liquid to vapor)
ENTHALPY = {"V": 0.0, "C": LATVAP, "R": LATVAP, "I": LATVAP + LATICE}
# constants::fill_value<Real>: float max / 1e5, computed in float
FILL = float(np.float32(np.finfo(np.float32).max) / np.float32(1e5))

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
def report_T(label, a, b, T, rtol):
###############################################################################
    diff = np.abs(a - b)
    rel = diff / T
    nbad = int(np.sum(rel > rtol))
    print(f"    {label}: max |diff| = {np.max(diff):.3e} K, max |diff|/T = {np.max(rel):.3e}, "
          f"points over tol = {nbad}")
    return nbad == 0

###############################################################################
def transfer_heating(names):
###############################################################################
    # Heating per unit mass moved by each p3_water_transfers component, from
    # its name: (H(to) - H(from))/cp [K per kg/kg]
    heating = {}
    for n in names:
        m = re.search(r"([VCRI])2([VCRI])$", n)
        expect(m is not None, f"Cannot tell the phases of p3_water_transfers component {n}")
        heating[n] = (ENTHALPY[m.group(2)] - ENTHALPY[m.group(1)]) / CP
    return heating

###############################################################################
def check_stage_temperatures(ds, t, names, transfers, stages, dt, from_previous, rtol):
###############################################################################
    # transfers/stages: the WaterTransfer/WaterStage names if P3/SHOC is
    # replayed, else None
    T_stage = components(ds, "water_transfer_stage_T", t, names)
    T = np.asarray(ds["T_mid"][t])
    ok = True

    # Stages that are not replayed hold the fill value
    filled = ([] if stages else ["ShocAfterTransport"]) + ([] if transfers else ["P3Entry", "P3Part2"])
    for s in filled:
        nbad = int(np.sum(T_stage[s] != FILL))
        print(f"    water_transfer_stage_T[{s}] (not replayed): points not the fill value = {nbad}")
        ok &= nbad == 0

    if transfers:
        rate = components(ds, "p3_water_transfers", t, transfers)
        heat = transfer_heating(transfers)
        heating = sum(heat[n] * rate[n] * dt for n in transfers)
        part1 = sum(heat[n] * rate[n] * dt for n in transfers if n.startswith("Part1"))
        wet_to_dry = np.asarray(ds["pseudo_density"][t]) / np.asarray(ds["pseudo_density_dry"][t])
        print(f"    P3 heating over its last call: max |T_mid - P3Entry| = "
              f"{np.max(np.abs(T - T_stage['P3Entry'])):.3e} K; part1 clips: "
              f"max |P3Part2 - P3Entry| = {np.max(np.abs(T_stage['P3Part2'] - T_stage['P3Entry'])):.3e} K "
              f"at {int(np.sum(part1 != 0))} points")
        if from_previous:
            ok &= report_T("P3Entry vs previous slice's T_mid (independent)", T_stage["P3Entry"],
                           np.asarray(ds["T_mid"][t-1]), T, rtol)
        ok &= report_T("P3Entry vs T_mid - P3 heating (bookkeeping)", T_stage["P3Entry"],
                       T - heating, T, rtol)
        ok &= report_T("P3Part2 - P3Entry vs part1 clip heating (bookkeeping)",
                       T_stage["P3Part2"] - T_stage["P3Entry"], part1 * wet_to_dry, T, rtol)

    if stages:
        stage = components(ds, "shoc_water_stages", t, stages)
        cond = stage["LiquidAfterCondensation"] - stage["LiquidAfterTransport"]
        T_after_shoc = T_stage["P3Entry"] if transfers else T
        print(f"    SHOC condensation heating: max |T_after_SHOC - ShocAfterTransport| = "
              f"{np.max(np.abs(T_after_shoc - T_stage['ShocAfterTransport'])):.3e} K")
        ok &= report_T("ShocAfterTransport vs T_after_SHOC - (LatVap/cp)*cond (bookkeeping)",
                       T_stage["ShocAfterTransport"], T_after_shoc - LATVAP / CP * cond, T, rtol)

    if not ok:
        print("    ERROR: water_transfer_stage_T is not the rebuilt stage temperature")
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
    parser.add_argument("--stage-T-header",
                        help="path to eamxx_water_tracer_transfers_functions.hpp: check water_transfer_stage_T")
    parser.add_argument("--dt", type=float,
                        help="dt of the last P3 call of a step [s] (needed by the stage-T checks with P3)")
    parser.add_argument("--T-from-previous-slice", action="store_true",
                        help="nothing but P3 changes T_mid: P3Entry must be the previous slice's T_mid")
    parser.add_argument("--T-rtol", type=float, default=1e-13,
                        help="tolerance of the stage-T checks, relative to T_mid")
    args = parser.parse_args()
    expect(args.stages_header or args.p3_header, "Need --stages-header and/or --p3-header")
    expect(not (args.stage_T_header and args.p3_header) or args.dt,
           "The stage-T checks with P3 need --dt")
    expect(not args.T_from_previous_slice or (args.stage_T_header and args.p3_header),
           "--T-from-previous-slice needs --stage-T-header and --p3-header")

    phases = read_enum(args.phases_header, r"enum\s+class\s+WaterPhase")
    expect(sorted(phases) == sorted(BULK), f"Unexpected WaterPhase {phases}")
    if args.stages_header:
        stages = read_enum(args.stages_header, r"struct\s+WaterStage\s*\{\s*enum\s*:\s*int")
    if args.p3_header:
        transfers = read_enum(args.p3_header, r"struct\s+WaterTransfer\s*\{\s*enum\s*:\s*int")
        sed = read_enum(args.p3_header, r"struct\s+SedFlux\s*\{\s*enum\s*:\s*int")
        expect(sed == ["CloudLiquid", "Rain", "Ice"], f"Unexpected SedFlux {sed}")
    if args.stage_T_header:
        stage_T_names = read_enum(args.stage_T_header, r"struct\s+TransferStage\s*\{\s*enum\s*:\s*int")
        expect(sorted(stage_T_names) == sorted(["ShocAfterTransport", "P3Entry", "P3Part2"]),
               f"Unexpected TransferStage {stage_T_names}")
    exact = EXACT if args.exact_qr_qi else []

    tracers = {}
    for entry in args.tracer:
        name, ratio = entry.split("=")
        tracers[name] = float(ratio)

    ds = Dataset(args.file, "r")
    # Read the fill value of stages that are not replayed as it is
    ds.set_auto_mask(False)
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
        if t > 0 and args.stage_T_header:
            print("  water_transfer_stage_T")
            success &= check_stage_temperatures(
                ds, t, stage_T_names,
                transfers if args.p3_header else None, stages if args.stages_header else None,
                args.dt, args.T_from_previous_slice, args.T_rtol)

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
