#!/usr/bin/env python3
"""
Validate tracer ratio preservation in EAMxx water tracer implementation.

This script verifies that passive tracer transport is correct by checking
that tracer ratios are preserved to machine precision throughout the model
integration. A scaling factor (typically 0.5) is applied to tracer 1's
surface fluxes, and this script verifies that all water reservoir ratios
converge to the same scaling factor.

Usage:
    python3 validate_tracer_scaling.py <output_nc> [rtol] [atol]

Arguments:
    output_nc: Path to model output NetCDF file
    rtol: Relative tolerance (default: 1e-12 for double precision)
    atol: Absolute tolerance (default: 1e-15 for double precision)

Exit codes:
    0: All reservoirs passed validation
    1: One or more reservoirs failed validation
    2: File not found or other error
"""

import sys
import numpy as np

def validate_tracer_scaling(output_file, rtol=1e-12, atol=1e-15, expected_ratio=0.5):
    """
    Validate that tracer ratios match expected value to within tolerance.

    Args:
        output_file: Path to NetCDF output file
        rtol: Relative tolerance
        atol: Absolute tolerance
        expected_ratio: Expected tracer 1 / tracer 0 ratio

    Returns:
        True if all reservoirs pass, False otherwise
    """
    try:
        import netCDF4 as nc
    except ImportError:
        print("ERROR: netCDF4 module not available. Install with: pip install netCDF4")
        return False

    try:
        ds = nc.Dataset(output_file, 'r')
    except FileNotFoundError:
        print(f"ERROR: Output file not found: {output_file}")
        return False
    except Exception as e:
        print(f"ERROR: Failed to open NetCDF file: {e}")
        return False

    # Water reservoir variables to validate
    # 3D fields: (time, tracer, col, lev)
    reservoirs_3d = ['qv', 'qc', 'qi', 'qr']
    # 2D surface flux fields: (time, tracer, col)
    reservoirs_2d = ['precip_liq_surf_mass', 'precip_ice_surf_mass']

    all_passed = True
    results = []

    print("\n" + "="*80)
    print("EAMxx Water Tracer Ratio Validation")
    print("="*80)
    print(f"Output file: {output_file}")
    print(f"Expected ratio: {expected_ratio}")
    print(f"Tolerance: rtol={rtol}, atol={atol}")
    print("="*80)

    # Validate 3D fields
    for var_name in reservoirs_3d:
        if var_name not in ds.variables:
            print(f"WARNING: Variable {var_name} not found in output, skipping")
            continue

        var = ds.variables[var_name]

        # Get final timestep data
        bulk = var[-1, 0, ...]      # Final time, tracer 0 (bulk)
        tracer = var[-1, 1, ...]    # Final time, tracer 1

        # Compute ratios only where bulk water is significant
        mask = bulk > 1e-20
        n_valid = np.sum(mask)

        if n_valid == 0:
            print(f"{var_name:30s}: No valid points (all bulk < 1e-20), SKIP")
            continue

        ratio = tracer[mask] / bulk[mask]
        mean_ratio = np.mean(ratio)
        max_error = np.max(np.abs(ratio - expected_ratio))
        rel_error = max_error / expected_ratio

        # Check against tolerance
        passed = (max_error <= atol + rtol * expected_ratio)
        status = "PASS" if passed else "FAIL"

        print(f"{var_name:30s}: mean ratio = {mean_ratio:.15f}, "
              f"max rel error = {rel_error:.3e}, {status}")

        results.append({
            'variable': var_name,
            'mean_ratio': mean_ratio,
            'max_error': max_error,
            'rel_error': rel_error,
            'passed': passed,
            'n_valid': n_valid
        })

        if not passed:
            all_passed = False

    # Validate 2D fields
    for var_name in reservoirs_2d:
        if var_name not in ds.variables:
            print(f"WARNING: Variable {var_name} not found in output, skipping")
            continue

        var = ds.variables[var_name]

        # Get final timestep data
        bulk = var[-1, 0, ...]      # Final time, tracer 0
        tracer = var[-1, 1, ...]    # Final time, tracer 1

        mask = bulk > 1e-20
        n_valid = np.sum(mask)

        if n_valid == 0:
            print(f"{var_name:30s}: No valid points (all bulk < 1e-20), SKIP")
            continue

        ratio = tracer[mask] / bulk[mask]
        mean_ratio = np.mean(ratio)
        max_error = np.max(np.abs(ratio - expected_ratio))
        rel_error = max_error / expected_ratio

        passed = (max_error <= atol + rtol * expected_ratio)
        status = "PASS" if passed else "FAIL"

        print(f"{var_name:30s}: mean ratio = {mean_ratio:.15f}, "
              f"max rel error = {rel_error:.3e}, {status}")

        results.append({
            'variable': var_name,
            'mean_ratio': mean_ratio,
            'max_error': max_error,
            'rel_error': rel_error,
            'passed': passed,
            'n_valid': n_valid
        })

        if not passed:
            all_passed = False

    ds.close()

    print("="*80)
    if all_passed:
        print("RESULT: All reservoirs passed tracer ratio validation")
        print("="*80 + "\n")
        return True
    else:
        print("RESULT: One or more reservoirs FAILED validation")
        print("\nFailed variables:")
        for res in results:
            if not res['passed']:
                print(f"  - {res['variable']}: rel_error = {res['rel_error']:.3e}")
        print("\nDebugging suggestions:")
        print("  1. Check tracer initialization (tracer 1 should = 0.5 * tracer 0)")
        print("  2. Check surface flux scaling (tracer 1 flux should = 0.5 * tracer 0 flux)")
        print("  3. Look for processes that may not preserve tracer dimension")
        print("  4. Verify no hardcoded slot-0 access without tracer loops")
        print("  5. Check field copies preserve tracer dimension")
        print("="*80 + "\n")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    output_file = sys.argv[1]
    rtol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-12
    atol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-15

    success = validate_tracer_scaling(output_file, rtol, atol)
    sys.exit(0 if success else 1)
