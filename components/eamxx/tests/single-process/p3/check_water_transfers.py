#!/usr/bin/env python3

"""
Check that P3's water transfer export fields close the bulk water budget.

For every P3 call (one per output time slice here, since P3 runs once per
step with no subcycling), the state P3 returns must equal the state at P3
entry plus every exported transfer:

  q_x(final) = p3_water_entry_x
             + dt * sum over transfers of (+rate into x, -rate out of x)
             + dt * g * (F_x(k-1) - F_x(k)) / pseudo_density   (sedimentation)

for x in qv, qc, qr, qi, with F = p3_sed_flux (downward flux out of the bottom
of level k; level 0 is the model top). Also checks that the flux out of the
lowest level, times dt, is the surface precipitation P3 accumulated.

The transfer components and their direction are read from the WaterTransfer
enum in p3_water_transfers.hpp, so this script cannot drift from P3.
"""

import argparse, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3] / "scripts"))
from utils import expect, ensure_netcdf4  # pylint: disable=wrong-import-position

ensure_netcdf4()

from netCDF4 import Dataset  # pylint: disable=wrong-import-position
import numpy as np            # pylint: disable=wrong-import-position

GRAVIT = 9.80616  # m/s2, physics::Constants::gravit
SPECIES = ["qv", "qc", "qr", "qi"]
LETTER  = {"V": "qv", "C": "qc", "R": "qr", "I": "qi"}

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
def species_axis(var, name):
###############################################################################
    # Variables are (time, ncol, <component>, lev); return the component axis
    dims = var.dimensions
    expect(len(dims) == 4 and dims[0] == "time" and dims[-1] == "lev",
           f"Unexpected dims for {name}: {dims}")
    return 2

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True, help="P3 standalone output file")
    parser.add_argument("--header", required=True, help="path to p3_water_transfers.hpp")
    parser.add_argument("--dt", type=float, required=True, help="P3 time step [s]")
    parser.add_argument("--rtol", type=float, default=1e-12,
                        help="tolerance, relative to the magnitude of the terms at each point")
    args = parser.parse_args()

    transfers = read_enum(args.header, "WaterTransfer")
    sed_species = read_enum(args.header, "SedFlux")
    entry_species = read_enum(args.header, "WaterEntry")
    expect(entry_species == ["Vapor", "CloudLiquid", "Rain", "Ice"], f"Unexpected WaterEntry {entry_species}")
    expect(sed_species == ["CloudLiquid", "Rain", "Ice"], f"Unexpected SedFlux {sed_species}")

    ds = Dataset(args.file, "r")
    nt = ds.dimensions["time"].size
    expect(nt >= 1, "Need at least one output slice")

    dt = args.dt
    rates = ds["p3_water_transfers"]
    species_axis(rates, "p3_water_transfers")
    expect(rates.shape[2] == len(transfers),
           f"p3_water_transfers has {rates.shape[2]} components, header has {len(transfers)}")

    success = True
    # The output skips t0, so every slice follows a P3 call
    for t in range(nt):
        entry = ds["p3_water_entry"][t]   # (ncol, 4, lev)
        sed   = ds["p3_sed_flux"][t]      # (ncol, 3, lev)
        rate  = rates[t]                  # (ncol, ntr, lev)
        dp    = ds["pseudo_density"][t]   # (ncol, lev)

        change = {s: np.zeros_like(dp) for s in SPECIES}
        scale  = {s: np.zeros_like(dp) for s in SPECIES}
        for itr, name in enumerate(transfers):
            m = re.search(r"([VCRI])2([VCRI])$", name)
            expect(m is not None, f"Cannot parse direction of transfer {name}")
            src, dst = LETTER[m.group(1)], LETTER[m.group(2)]
            amount = rate[:, itr, :] * dt
            change[src] -= amount
            change[dst] += amount
            scale[src] += np.abs(amount)
            scale[dst] += np.abs(amount)

        for ised, s in zip(range(3), ["qc", "qr", "qi"]):
            out_mass = sed[:, ised, :] * dt              # kg/m2 out of the bottom of k
            in_mass  = np.zeros_like(out_mass)
            in_mass[:, 1:] = out_mass[:, :-1]            # level 0 is the top
            div = (in_mass - out_mass) * GRAVIT / dp
            change[s] += div
            scale[s] += np.abs(in_mass * GRAVIT / dp) + np.abs(out_mass * GRAVIT / dp)

        for i, s in enumerate(SPECIES):
            q0 = entry[:, i, :]
            q1 = ds[s][t]
            resid = q1 - (q0 + change[s])
            tot = np.abs(q0) + np.abs(q1) + scale[s] + 1e-300
            worst = np.max(np.abs(resid) / tot)
            nbad = int(np.sum(np.abs(resid) > args.rtol * tot))
            print(f"  t={t} {s}: max |residual|/scale = {worst:.3e}, "
                  f"max |residual| = {np.max(np.abs(resid)):.3e} kg/kg, "
                  f"max |change| = {np.max(np.abs(change[s])):.3e} kg/kg, points over tol = {nbad}")
            if nbad > 0:
                success = False

        # Surface precipitation: flux out of the lowest level
        liq = (sed[:, 0, -1] + sed[:, 1, -1]) * dt
        ice = sed[:, 2, -1] * dt
        for name, exported in [("precip_liq_surf_mass", liq), ("precip_ice_surf_mass", ice)]:
            accum = ds[name][t]
            resid = np.abs(accum - exported)
            tol = args.rtol * (np.abs(accum) + np.abs(exported)) + 1e-300
            nbad = int(np.sum(resid > tol))
            print(f"  t={t} {name}: max |residual| = {np.max(resid):.3e} kg/m2, "
                  f"max value = {np.max(np.abs(accum)):.3e} kg/m2, points over tol = {nbad}")
            if nbad > 0:
                success = False

        # The test must exercise every stage that can move water here; a budget
        # that closes because nothing happened proves nothing.
        active = [name for itr, name in enumerate(transfers) if np.any(rate[:, itr, :] != 0)]
        print(f"  t={t} nonzero transfers: {len(active)}/{len(transfers)}: {', '.join(active)}")
        for s, ised in [("qr", 1), ("qi", 2)]:
            if not np.any(sed[:, ised, :] != 0):
                print(f"  t={t} ERROR: no {s} sedimentation in this test")
                success = False

    print("Water transfer budget check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
