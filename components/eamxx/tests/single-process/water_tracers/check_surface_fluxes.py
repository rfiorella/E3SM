#!/usr/bin/env python3

"""
Check the tracer surface fluxes of the water_tracers standalone test.

The test has constant initial conditions: a half-ocean column that evaporates
(E_ocn > 0) with a merged evaporation surf_evap = S > E_ocn, so the rest,
L = S - E_ocn, evaporates from land and sea ice. With wtrc_bp (R_ocn = R_land = 1,
bulk_proportional) the flux must be S bit for bit. The other tracers have
R_ocn = 0.5 and R_land = 2: with the bulk-proportional wtrc_rp the flux is
R_ocn*E_ocn + R_land*L = 2*S - 1.5*E_ocn, and wtrc_rl (bulk_formula, whose tracer
starts at R_ocn*qv, so its ocean part is R_ocn*E_ocn too) must agree with it to
roundoff. wtrc_bf (bulk_formula, R = 1) follows the bulk water to roundoff, but
not bit for bit: its ocean part comes from the tracer's own concentration.
E_ocn is recovered from wtrc_rp and checked against the bulk formula
f*rho*u*r*(ssq - qv) with rho = p/(Rd*Tv), to the accuracy of that rho.
"""

import argparse, pathlib, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3] / "scripts"))
from utils import expect, ensure_netcdf4  # pylint: disable=wrong-import-position

ensure_netcdf4()

from netCDF4 import Dataset  # pylint: disable=wrong-import-position
import numpy as np            # pylint: disable=wrong-import-position

RD = 287.04  # J/kg/K; only used for the approximate density below

###############################################################################
def main():
###############################################################################
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--file", required=True)
    parser.add_argument("--surf-evap", type=float, required=True)
    parser.add_argument("--ocnfrac", type=float, required=True)
    parser.add_argument("--ssq", type=float, required=True)
    parser.add_argument("--re", type=float, required=True)
    parser.add_argument("--ustar", type=float, required=True)
    args = parser.parse_args()

    ds = Dataset(args.file, "r")
    success = True
    for t in range(ds.dimensions["time"].size):
        get = lambda n: np.asarray(ds[n][t], dtype=np.float64)
        bp, bf, rp, rl = (get(f"wtrc_{n}_qv_surf_flux") for n in ("bp", "bf", "rp", "rl"))
        evap = get("surf_evap")
        s = evap  # the IC as stored in the output (to the output precision)
        tol = 1e-12 if ds["surf_evap"].dtype == np.float64 else 1e-5

        # E_ocn from the R = (0.5, 2) proportional tracer: rp = 2*S - 1.5*E
        e_ocn = (2*s - rp) / 1.5
        qv = get("qv")[:, -1]
        rho = get("p_mid")[:, -1] / (RD * get("T_mid")[:, -1] * (1 + 0.608*qv))
        e_ref = args.ocnfrac * rho * args.ustar * args.re * (args.ssq - qv)

        print(f" t={t}: surf_evap = {evap.min():.6e}, E_ocn = {e_ocn.min():.6e}..{e_ocn.max():.6e}, "
              f"E_ocn/E_ref - 1 = {np.max(np.abs(e_ocn/e_ref - 1)):.2e}")
        checks = {
            "surf_evap is the IC":                 bool(np.all(evap == args.surf_evap)) if tol < 1e-6 else bool(np.all(np.abs(evap/args.surf_evap - 1) < tol)),
            "0 < E_ocn < surf_evap":               bool(np.all((e_ocn > 0) & (e_ocn < s))),
            "E_ocn matches the bulk formula (1%)": bool(np.all(np.abs(e_ocn/e_ref - 1) < 1e-2)),
            "wtrc_bp flux == surf_evap (exactly)": bool(np.all(bp == s)),
            "wtrc_bf flux == surf_evap (roundoff)":  bool(np.all(np.abs(bf/s - 1) < tol)),
            "wtrc_rl flux == wtrc_rp flux (roundoff)": bool(np.all(np.abs(rl/rp - 1) < tol)),
            "wtrc_rp differs from surf_evap":      bool(np.all(np.abs(rp/s - 1) > 0.1)),
        }
        for label, ok in checks.items():
            print(f"    {'ok  ' if ok else 'FAIL'} {label}")
            success &= ok

    print("Water tracer surface flux check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
