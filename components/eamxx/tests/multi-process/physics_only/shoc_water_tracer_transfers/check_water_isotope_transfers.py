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

With --kinetic, the run has water_tracer_transfers::kinetic_fractionation on.
The checks above then use the effective factors where the replay does: the
P3 check normalizes by alpha_eff - 1, alpha_eff = akci below 253.15 K (with
--kinetic-ice-S, the supersaturation the run used: the empirical ssatf(T), or
P3's RH over ice of the bulk vapor after the part1 clips), and skips points
where |alpha_eff - 1| < 0.2*|alpha_i - 1|; with P3's RH over ice every point
must be in [0.5, 2] (S varies between the P3 calls whose deposits the ice
holds, and akci uses the last call's). SHOC's net condensation stays at
equilibrium, so its checks are unchanged. akel and akci are computed here
from the Fortran formulas (wiso_akel, wiso_akci) and the relative humidity of
EAMxx's RH diagnostic (Murphy & Koop 2005 at p_dry_mid, wet basis), which the
output must then hold (p_dry_mid, pseudo_density_dry). Further:
  - --shoc-evap-rayleigh (one SHOC call from a state with R = 1 in every
    phase, and SHOC the last process): where SHOC's net transfer evaporated
    cloud, leaving a fraction f in [0.01, 0.99] of the cloud after transport,
    the cloud left must have R_qc = f^(beta - 1), the Rayleigh result for
    beta = 1/akel (1/alpha_l without --kinetic) at ShocAfterTransport and the
    relative humidity over liquid of SHOC's vapor after transport, to
    --evap-rtol, at >= --min-evap-points points. With --kinetic, at least
    --min-evap-kinetic of them must also be off the equilibrium result by more
    than 100 times --evap-rtol (the kinetic factor is resolved; SHOC
    evaporates cloud close to saturation, so the effect is small).
  - --equilibrium-file F (with --p3-header): F is the same run with
    kinetic_fractionation off. The bulk water must be the same bit for bit (the
    tracers do not feed back), and the kinetic effects must have their sign:
     - rain evaporating below cloud into subsaturated air (P3's last call
       evaporated more than 10% of the rain available to it, no cloud liquid,
       relative humidity over liquid below 0.95), where that rain is mostly
       the level's own (the rain falling in from above over the call is less
       than 10 times the rain available: otherwise it carries the history of
       the levels above, e.g. melted snow, whose kinetic deposition has the
       other sign): the rain left must be more enriched than in F, in both
       isotopes, and have a lower deuterium excess (d = delta_D -
       8*delta_18O, the standard ratios being 1), at > 98% of at least 10
       points;
     - cold ice deposited from vapor (the P3 points above, below 253.15 K):
       less enriched relative to the vapor than in F, e_kinetic < e_F (e as
       above, with the equilibrium alpha_i), at > 90% of at least 5 points.
    An inverted kinetic factor fails these.

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
EP_2 = 18.016 / 28.966  # physics::Constants::ep_2 = MWH2O/MWdry
TKINI = 253.15  # K, below it deposition on ice is kinetic (wiso_akci)
# D_i/D of Merlivat (1978), the default diffusivity_formulation
DIFRM = {"HDO": 0.9757, "H218O": 0.9727}

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
def svp_murphy_koop(T, ice):
###############################################################################
    # Saturation vapor pressure [Pa] of Murphy & Koop (2005), eqs. (7) (ice,
    # below 273.15 K) and (10) (liquid), as physics::Functions::MurphyKoop_svp
    logt = np.log(T)
    liq = np.exp(54.842763 - 6763.22 / T - 4.210 * logt + 0.000367 * T
                 + np.tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * logt + 0.014025 * T))
    if not ice:
        return liq
    ice_svp = np.exp(9.550426 - 5723.265 / T + 3.53068 * logt - 0.00728332 * T)
    return np.where(T < 273.15, ice_svp, liq)

###############################################################################
def saturation_ratio(qv, T, p_dry, dp, dp_dry, ice):
###############################################################################
    # qv [kg/kg, wet] over qv_sat_wet, as share/diagnostics/relative_humidity.cpp
    qv_sat = EP_2 * svp_murphy_koop(T, ice) / np.maximum(p_dry, 1e-3) * dp_dry / dp
    return qv / qv_sat

###############################################################################
def akel(species, alpha, rh):
###############################################################################
    # wiso_akel (Stewart 1975): effective liquid/vapor factor of evaporation
    h = np.minimum(0.25 * np.clip(rh, 0, 1) + 1 - 0.25, 1)
    dondi = (1 / DIFRM[species]) ** 0.58
    return alpha * h / (alpha * dondi * (h - 1) + 1)

###############################################################################
def ssatf(T):
###############################################################################
    # wiso_ssatf: empirical supersaturation over ice
    return np.minimum(np.maximum(1 - 0.002 * (T - 273.16), 1), 2)

###############################################################################
def akci(species, alpha, T, S):
###############################################################################
    # wiso_akci (Jouzel & Merlivat 1984): effective ice/vapor factor of
    # deposition, below TKINI only
    sat = np.clip(S, 1, 2)
    dondi = 1 / DIFRM[species]
    return np.where(T < TKINI, alpha * sat / (alpha * dondi * (sat - 1) + 1), alpha)

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
    parser.add_argument("--kinetic", action="store_true",
                        help="the run has water_tracer_transfers::kinetic_fractionation on")
    parser.add_argument("--kinetic-ice-S", choices=["empirical", "model"], default="empirical",
                        help="water_tracer_transfers::kinetic_ice_supersaturation of the run (with --kinetic)")
    parser.add_argument("--shoc-evap-rayleigh", action="store_true",
                        help="one SHOC call from R = 1: check the Rayleigh result of SHOC's net evaporation")
    parser.add_argument("--evap-rtol", type=float, default=1e-12,
                        help="tolerance of R_qc against the Rayleigh result of SHOC's net evaporation")
    parser.add_argument("--min-evap-points", type=int, default=10,
                        help="SHOC net evaporation points needed (with --shoc-evap-rayleigh)")
    parser.add_argument("--min-evap-kinetic", type=int, default=1,
                        help="of them, points off the equilibrium result needed (with --kinetic)")
    parser.add_argument("--equilibrium-file",
                        help="the same run with kinetic_fractionation off: check the kinetic effects' sign")
    args = parser.parse_args()
    expect(args.stages_header or args.p3_header, "Need --stages-header and/or --p3-header")
    expect(not (args.alpha_check and args.p3_header) or args.dt_p3, "The P3 alpha check needs --dt-p3")
    expect(not args.shoc_evap_rayleigh or (args.stages_header and args.p3_header is None),
           "--shoc-evap-rayleigh needs --stages-header and no --p3-header")
    expect(not args.equilibrium_file or (args.kinetic and args.p3_header and args.dt_p3),
           "--equilibrium-file needs --kinetic, --p3-header and --dt-p3")

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
        sed_names = read_enum(args.p3_header, r"struct\s+SedFlux\s*\{\s*enum\s*:\s*int")

    ds = Dataset(args.file, "r")
    nt = ds.dimensions["time"].size
    expect(nt >= 2, "Need the t0 slice and at least one step")
    ds_eq = Dataset(args.equilibrium_file, "r") if args.equilibrium_file else None
    if ds_eq is not None:
        expect(ds_eq.dimensions["time"].size == nt, "The equilibrium file has a different number of slices")

    success = True
    # Over the run, per tracer: new SHOC cloud points, and e of deposited ice
    new_cloud = {name: 0 for name in tracers}
    e_ice = {name: [] for name in tracers}
    # --shoc-evap-rayleigh: evaporating points, and those off equilibrium
    evap_pts = {name: 0 for name in tracers}
    evap_kin = {name: 0 for name in tracers}
    # --equilibrium-file: rain evaporation points (R higher, d lower than in
    # the equilibrium run) and cold deposition points (e lower)
    rain_cmp = {"n": 0, "R_higher": {name: 0 for name in tracers}, "d_lower": 0}
    cold_e = {name: [] for name in tracers}
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

            # ---- SHOC: net evaporation of cloud, one call from R = 1 ----
            if args.shoc_evap_rayleigh:
                f_left = st["LiquidAfterCondensation"] / np.where(qc_tr > 0, qc_tr, 1)
                ev = (cond < 0) & (qc_tr > 0) & (f_left >= 0.01) & (f_left <= 0.99)
                if args.kinetic:
                    rh_l = saturation_ratio(qv_tr[ev], T[ev], np.asarray(ds["p_dry_mid"][t])[ev],
                                            np.asarray(ds["pseudo_density"][t])[ev],
                                            np.asarray(ds["pseudo_density_dry"][t])[ev], False)
                for name, sp in tracers.items():
                    if not np.any(ev):
                        continue
                    alpha = alpha_liquid(sp, T[ev])
                    beta_eq = 1 / alpha
                    beta = 1 / akel(sp, alpha, rh_l) if args.kinetic else beta_eq
                    log_f = np.log(f_left[ev])
                    R_c = np.asarray(ds[f"{name}_qc"][t])[ev] / q["qc"][ev]
                    err = np.abs(R_c / np.exp((beta - 1) * log_f) - 1)
                    err_eq = np.abs(R_c / np.exp((beta_eq - 1) * log_f) - 1)
                    nbad = int(np.sum(err > args.evap_rtol))
                    nkin = int(np.sum(err_eq > 100 * args.evap_rtol))
                    msg = (f"  {name} SHOC net evaporation ({int(np.sum(ev))} points): "
                           f"max |R_qc/Rayleigh - 1| = {np.max(err):.3e}, points over tol = {nbad}; "
                           f"|R_qc/Rayleigh(1/alpha_l) - 1| in [{np.min(err_eq):.3e}, {np.max(err_eq):.3e}]")
                    if args.kinetic:
                        msg += (f", RH_l in [{np.min(rh_l):.3f}, {np.max(rh_l):.3f}], "
                                f"beta/beta_eq - 1 in [{np.min(beta/beta_eq)-1:.3e}, {np.max(beta/beta_eq)-1:.3e}], "
                                f"points off equilibrium = {nkin}")
                    print(msg)
                    success &= nbad == 0
                    evap_pts[name] += int(np.sum(ev))
                    evap_kin[name] += nkin

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
            if args.kinetic:
                # P3's state before its part2 rates: the vapor after the part1 clips
                qv_p2 = entry["Vapor"] + (rate["Part1ClipC2V"] + rate["Part1ClipR2V"]
                                          + rate["Part1ClipI2V"]) * args.dt_p3
                sat_in = [np.asarray(ds[n][t]) for n in ["p_dry_mid", "pseudo_density", "pseudo_density_dry"]]
                S_i = (saturation_ratio(qv_p2[dom], T, *[x[dom] for x in sat_in], True)
                       if args.kinetic_ice_S == "model" else ssatf(T))
            for name, sp in tracers.items():
                R_i = np.asarray(ds[f"{name}_qi"][t])[dom] / q["qi"][dom]
                R_v = np.asarray(ds[f"{name}_qv"][t])[dom] / q["qv"][dom]
                a_i = alpha_ice(sp, T)
                a_eff = akci(sp, a_i, T, S_i) if args.kinetic else a_i
                # (where akci is close to 1 the normalization is ill-posed)
                ok = np.abs(a_eff - 1) > 0.2 * np.abs(a_i - 1)
                e_ice[name] += list(((R_i / R_v - 1) / (a_eff - 1))[ok])

            # ---- Against the equilibrium run: the kinetic effects' sign ----
            if ds_eq is not None:
                for p in PHASES:
                    expect(np.array_equal(q[p], np.asarray(ds_eq[p][t])),
                           f"Bulk {p} differs from the equilibrium run at t={t}: not the same run")
                # Rain evaporating below cloud into subsaturated air
                rain_in = entry["Rain"] + (rate["Part2AccretC2R"] + rate["Part2AutoconvC2R"] + rate["Part2ShedC2R"]
                                           + rate["Part2MeltI2R"] + rate["Part1ClipI2R"]) * args.dt_p3
                evap = rate["Part2EvapR2V"] * args.dt_p3
                rh_l = saturation_ratio(qv_p2, T_stage["P3Part2"], *sat_in, False)
                # rain falling in from above over the call [kg/kg]
                sed_rain = components(ds, "p3_sed_flux", t, sed_names)["Rain"]
                rain_from_above = np.zeros_like(sed_rain)
                rain_from_above[:, 1:] = sed_rain[:, :-1] * GRAVIT * args.dt_p3 / sat_in[1][:, 1:]
                rain = ((evap > 0.1 * rain_in) & (entry["CloudLiquid"] < 1e-8) & (q["qc"] < 1e-8)
                        & (q["qr"] > 1e-9) & (rh_l < 0.95) & (rain_from_above < 10 * rain_in))
                n = int(np.sum(rain))
                if n > 0:
                    delta, delta_eq = {}, {}
                    for name, sp in tracers.items():
                        R = np.asarray(ds[f"{name}_qr"][t])[rain] / q["qr"][rain]
                        R_eq = np.asarray(ds_eq[f"{name}_qr"][t])[rain] / q["qr"][rain]
                        rain_cmp["R_higher"][name] += int(np.sum(R > R_eq))
                        delta[sp], delta_eq[sp] = R - 1, R_eq - 1
                        print(f"  {name} rain evaporating below cloud ({n} points, RH_l in "
                              f"[{np.min(rh_l[rain]):.2f}, {np.max(rh_l[rain]):.2f}]): R_qr/R_qr(eq) - 1 in "
                              f"[{np.min(R/R_eq)-1:.3e}, {np.max(R/R_eq)-1:.3e}], median {np.median(R/R_eq)-1:.3e}")
                    d = delta["HDO"] - 8 * delta["H218O"]
                    d_eq = delta_eq["HDO"] - 8 * delta_eq["H218O"]
                    rain_cmp["d_lower"] += int(np.sum(d < d_eq))
                    rain_cmp["n"] += n
                    print(f"  rain deuterium excess - equilibrium run's: [{np.min(d-d_eq):.3e}, "
                          f"{np.max(d-d_eq):.3e}], median {np.median(d-d_eq):.3e}")
                    # (the vapor there, for information)
                    dv, dv_eq = {}, {}
                    for name, sp in tracers.items():
                        dv[sp] = np.asarray(ds[f"{name}_qv"][t])[rain] / q["qv"][rain] - 1
                        dv_eq[sp] = np.asarray(ds_eq[f"{name}_qv"][t])[rain] / q["qv"][rain] - 1
                    ddv = (dv["HDO"] - 8 * dv["H218O"]) - (dv_eq["HDO"] - 8 * dv_eq["H218O"])
                    print(f"  (vapor there: deuterium excess - equilibrium run's in [{np.min(ddv):.3e}, "
                          f"{np.max(ddv):.3e}], higher at {int(np.sum(ddv > 0))} of {n})")
                # Cold ice deposited from vapor
                cold = T < TKINI
                if np.any(cold):
                    for name, sp in tracers.items():
                        a_i = alpha_ice(sp, T[cold])
                        e_k = (np.asarray(ds[f"{name}_qi"][t])[dom][cold] / np.asarray(ds[f"{name}_qv"][t])[dom][cold]
                               * q["qv"][dom][cold] / q["qi"][dom][cold] - 1) / (a_i - 1)
                        e_e = (np.asarray(ds_eq[f"{name}_qi"][t])[dom][cold] / np.asarray(ds_eq[f"{name}_qv"][t])[dom][cold]
                               * q["qv"][dom][cold] / q["qi"][dom][cold] - 1) / (a_i - 1)
                        cold_e[name] += list(zip(e_k, e_e))

    if args.alpha_check:
        print(" over the run")
        for name in tracers:
            if args.stages_header and args.p3_header is None:
                print(f"  {name}: new SHOC cloud points = {new_cloud[name]} (need >= {args.min_new_cloud})")
                success &= new_cloud[name] >= args.min_new_cloud
            if args.shoc_evap_rayleigh:
                ok = (evap_pts[name] >= args.min_evap_points and
                      (not args.kinetic or evap_kin[name] >= args.min_evap_kinetic))
                print(f"  {name}: SHOC net evaporation points = {evap_pts[name]}"
                      + (f", off equilibrium = {evap_kin[name]}" if args.kinetic else "")
                      + ("" if ok else f"  <-- needs >= {args.min_evap_points} (and >= {args.min_evap_kinetic} "
                                       "off equilibrium with --kinetic)"))
                success &= ok
            if args.p3_header:
                e = np.array(e_ice[name])
                # (with P3's RH over ice, S varies between the P3 calls whose
                # deposits the ice holds, but akci uses the last call's: wider
                # bounds per point)
                e_max = 2.0 if args.kinetic and args.kinetic_ice_S == "model" else 1.5
                ok = len(e) >= 10 and abs(np.median(e) - 1) < 0.1 and np.all((e > 0.5) & (e < e_max))
                print(f"  {name}: deposition-dominated ice points = {len(e)}, "
                      f"e = (R_qi/R_qv - 1)/({'alpha_eff' if args.kinetic else 'alpha_i'} - 1): "
                      + (f"min {np.min(e):.3f}, median {np.median(e):.3f}, max {np.max(e):.3f}" if len(e) else "none")
                      + ("" if ok else f"  <-- needs >= 10 points, |median - 1| < 0.1, all in [0.5, {e_max}]"))
                success &= ok

    if ds_eq is not None:
        print(" against the equilibrium run, over the run")
        n = rain_cmp["n"]
        for name in tracers:
            frac = rain_cmp["R_higher"][name] / n if n else 0.0
            ok = n >= 10 and frac > 0.98
            print(f"  {name}: rain evaporating below cloud: R_qr above the equilibrium run's at {frac:.4f} "
                  f"of {n} points" + ("" if ok else "  <-- needs >= 10 points, > 0.98"))
            success &= ok
        frac = rain_cmp["d_lower"] / n if n else 0.0
        ok = n >= 10 and frac > 0.98
        print(f"  rain deuterium excess below the equilibrium run's at {frac:.4f} of {n} points"
              + ("" if ok else "  <-- needs >= 10 points, > 0.98"))
        success &= ok
        for name in tracers:
            pairs = np.array(cold_e[name]).reshape(-1, 2)
            m = len(pairs)
            frac = float(np.mean(pairs[:, 0] < pairs[:, 1])) if m else 0.0
            ok = m >= 5 and frac > 0.9
            print(f"  {name}: cold deposition-dominated ice ({m} points): e below the equilibrium run's at "
                  f"{frac:.4f}" + (f", median e {np.median(pairs[:, 0]):.3f} vs {np.median(pairs[:, 1]):.3f}" if m else "")
                  + ("" if ok else "  <-- needs >= 5 points, > 0.9"))
            success &= ok

    print("Water isotope transfers check: " + ("SUCCESS" if success else "FAIL"))
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
