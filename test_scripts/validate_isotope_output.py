#!/usr/bin/env python3
"""
EAMv2 Water Isotope Output Validation Script

This script validates that isotope tracer output from EAM simulations
contains the expected variables and that they have physically reasonable values.

Usage:
    python validate_isotope_output.py <case_directory> <tracer_package>

Author: Auto-generated test suite
Date: 2026-04-07
"""

import sys
import os
import glob
import argparse
from pathlib import Path

try:
    import netCDF4 as nc
    import numpy as np
    HAS_NETCDF = True
except ImportError:
    HAS_NETCDF = False
    print("WARNING: netCDF4 not available, using basic validation only")


# Define expected tracers for each package
EXPECTED_TRACERS = {
    'h2o_hdo': {
        'vapor': ['Q_H2O', 'Q_HDO'],
        'liquid': ['CLDLIQ_H2O', 'CLDLIQ_HDO'],
        'ice': ['CLDICE_H2O', 'CLDICE_HDO'],
        'description': 'H2O + HDO (14 tracers)'
    },
    'h2o_h218o': {
        'vapor': ['Q_H2O', 'Q_H218O'],
        'liquid': ['CLDLIQ_H2O', 'CLDLIQ_H218O'],
        'ice': ['CLDICE_H2O', 'CLDICE_H218O'],
        'description': 'H2O + H218O (14 tracers)'
    },
    'h2o_hdo_h218o': {
        'vapor': ['Q_H2O', 'Q_HDO', 'Q_H218O'],
        'liquid': ['CLDLIQ_H2O', 'CLDLIQ_HDO', 'CLDLIQ_H218O'],
        'ice': ['CLDICE_H2O', 'CLDICE_HDO', 'CLDICE_H218O'],
        'description': 'H2O + HDO + H218O (21 tracers)'
    },
    'all_stable_wiso': {
        'vapor': ['Q_H2O', 'Q_H216O', 'Q_HDO', 'Q_H218O', 'Q_H217O'],
        'liquid': ['CLDLIQ_H2O', 'CLDLIQ_H216O', 'CLDLIQ_HDO', 'CLDLIQ_H218O', 'CLDLIQ_H217O'],
        'ice': ['CLDICE_H2O', 'CLDICE_H216O', 'CLDICE_HDO', 'CLDICE_H218O', 'CLDICE_H217O'],
        'description': 'All stable isotopes (35 tracers)'
    }
}


def find_history_files(case_dir):
    """Find EAM history files in the case directory"""
    run_dir = os.path.join(case_dir, 'run')
    
    if not os.path.exists(run_dir):
        print(f"ERROR: Run directory not found: {run_dir}")
        return []
    
    hist_files = sorted(glob.glob(os.path.join(run_dir, '*.eam.h0.*.nc')))
    return hist_files


def check_run_completion(case_dir):
    """Check if the model run completed successfully"""
    run_dir = os.path.join(case_dir, 'run')
    cpl_logs = glob.glob(os.path.join(run_dir, 'cpl.log.*'))
    
    if not cpl_logs:
        print("ERROR: No coupler log files found")
        return False
    
    # Check the most recent log
    latest_log = sorted(cpl_logs)[-1]
    
    with open(latest_log, 'r') as f:
        content = f.read()
        if 'SUCCESSFUL TERMINATION' in content:
            print("✓ Model run completed successfully")
            return True
        else:
            print("✗ Model run did not complete successfully")
            return False


def validate_with_netcdf(hist_file, tracer_package):
    """Validate output using netCDF4 library (detailed validation)"""
    
    if not HAS_NETCDF:
        print("WARNING: Cannot perform detailed validation without netCDF4")
        return False
    
    print(f"\nValidating: {os.path.basename(hist_file)}")
    
    try:
        ds = nc.Dataset(hist_file, 'r')
    except Exception as e:
        print(f"ERROR: Cannot open history file: {e}")
        return False
    
    expected = EXPECTED_TRACERS.get(tracer_package)
    if not expected:
        print(f"ERROR: Unknown tracer package: {tracer_package}")
        ds.close()
        return False
    
    print(f"Package: {expected['description']}")
    
    all_passed = True
    
    # Check each category of tracers
    for category in ['vapor', 'liquid', 'ice']:
        if category not in expected:
            continue
            
        print(f"\n{category.upper()} tracers:")
        
        for tracer in expected[category]:
            if tracer not in ds.variables:
                print(f"  ✗ MISSING: {tracer}")
                all_passed = False
                continue
            
            # Get the variable
            var = ds.variables[tracer]
            data = var[:]
            
            # Check for all zeros
            if np.all(data == 0):
                print(f"  ✗ ALL ZEROS: {tracer}")
                all_passed = False
                continue
            
            # Check for NaN or Inf
            if np.any(np.isnan(data)):
                print(f"  ✗ CONTAINS NaN: {tracer}")
                all_passed = False
                continue
            
            if np.any(np.isinf(data)):
                print(f"  ✗ CONTAINS Inf: {tracer}")
                all_passed = False
                continue
            
            # Check for negative values (physical constraint)
            if np.any(data < 0):
                neg_count = np.sum(data < 0)
                neg_percent = 100.0 * neg_count / data.size
                print(f"  ⚠ NEGATIVE VALUES: {tracer} ({neg_percent:.2f}% of points)")
                # This is a warning, not a failure for this basic test
            
            # Get some statistics
            nonzero = data[data != 0]
            if nonzero.size > 0:
                min_val = np.min(nonzero)
                max_val = np.max(nonzero)
                mean_val = np.mean(nonzero)
                print(f"  ✓ ACTIVE: {tracer} (min={min_val:.2e}, max={max_val:.2e}, mean={mean_val:.2e})")
            else:
                print(f"  ✓ PRESENT: {tracer} (all zeros - may be OK for some fields)")
    
    ds.close()
    
    return all_passed


def validate_basic(hist_file, tracer_package):
    """Basic validation without netCDF4 (just check variable presence)"""
    
    import subprocess
    
    print(f"\nValidating: {os.path.basename(hist_file)}")
    
    expected = EXPECTED_TRACERS.get(tracer_package)
    if not expected:
        print(f"ERROR: Unknown tracer package: {tracer_package}")
        return False
    
    print(f"Package: {expected['description']}")
    
    try:
        # Use ncdump to get variable list
        result = subprocess.run(['ncdump', '-h', hist_file], 
                              capture_output=True, text=True, check=True)
        header = result.stdout
    except Exception as e:
        print(f"ERROR: Cannot run ncdump: {e}")
        return False
    
    all_passed = True
    
    # Check for presence of key tracers (vapor only for basic check)
    print("\nVAPOR tracers:")
    for tracer in expected.get('vapor', []):
        if f"float {tracer}(" in header or f"double {tracer}(" in header:
            print(f"  ✓ PRESENT: {tracer}")
        else:
            print(f"  ✗ MISSING: {tracer}")
            all_passed = False
    
    return all_passed


def main():
    parser = argparse.ArgumentParser(
        description='Validate EAMv2 water isotope output')
    parser.add_argument('case_dir', help='Path to case directory')
    parser.add_argument('tracer_package', 
                       choices=list(EXPECTED_TRACERS.keys()),
                       help='Water tracer package name')
    parser.add_argument('--basic', action='store_true',
                       help='Use basic validation (no netCDF4 required)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("EAMv2 Water Isotope Output Validation")
    print("="*80)
    
    # Check if case directory exists
    if not os.path.exists(args.case_dir):
        print(f"ERROR: Case directory not found: {args.case_dir}")
        return 1
    
    print(f"\nCase directory: {args.case_dir}")
    print(f"Tracer package: {args.tracer_package}")
    
    # Check run completion
    print("\n" + "-"*80)
    print("Checking run completion...")
    print("-"*80)
    
    if not check_run_completion(args.case_dir):
        print("\nERROR: Run did not complete successfully")
        return 1
    
    # Find history files
    print("\n" + "-"*80)
    print("Finding history files...")
    print("-"*80)
    
    hist_files = find_history_files(args.case_dir)
    
    if not hist_files:
        print("ERROR: No history files found")
        return 1
    
    print(f"Found {len(hist_files)} history file(s)")
    
    # Validate the most recent history file
    print("\n" + "-"*80)
    print("Validating isotope tracers...")
    print("-"*80)
    
    latest_file = hist_files[-1]
    
    if args.basic or not HAS_NETCDF:
        success = validate_basic(latest_file, args.tracer_package)
    else:
        success = validate_with_netcdf(latest_file, args.tracer_package)
    
    # Print summary
    print("\n" + "="*80)
    if success:
        print("VALIDATION PASSED")
        print("="*80)
        return 0
    else:
        print("VALIDATION FAILED")
        print("="*80)
        return 1


if __name__ == '__main__':
    sys.exit(main())
