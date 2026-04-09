# Test Script Fixes - 2026-04-07

## Overview

This document summarizes the fixes applied to resolve the "Case directory does not exist" error and improve the overall robustness of the EAMv2 water isotope test suite.

## Issues Fixed

### 1. **Primary Issue: Case Directory Path Errors**

**Problem:** The `create_test_case()` function was returning error messages mixed with stdout, causing `case_dir` variable to contain invalid paths like `--` or error text instead of actual directory paths.

**Root Cause:**
- Log messages were going to stdout instead of stderr
- `create_newcase` output was being captured in command substitution
- No validation of the returned path
- Error handling failed to create proper fallback directory

**Fix Applied:**
- Redirected all log messages to stderr using `>&2`
- Redirected `create_newcase` output to a dedicated log file
- Added comprehensive validation of returned case directory path
- Added checks for directory existence and required case-specific files (xmlquery, xmlchange are created in the case directory by create_newcase)
- Improved error messages with log file references

**Note:** xmlquery and xmlchange are case-specific scripts created by `create_newcase` in each case directory, not global CIME scripts.

### 2. **Removed Unsupported Flag**

**Problem:** E3SM's CIME no longer supports `--run-unsupported` flag.

**Fix Applied:**
- Removed `--run-unsupported` from `create_newcase` command

### 3. **Path Validation Issues**

**Problem:** Functions tried to `cd` to invalid directories without proper validation.

**Fix Applied:**
- Added directory existence checks in all functions before `cd`
- Added validation for required scripts (xmlquery, xmlchange, case.setup, case.build)
- Added better error messages indicating what's missing and where

### 4. **Improved Error Handling in `run_single_test()`**

**Problem:** Invalid case_dir values weren't properly detected and handled.

**Fix Applied:**
- Validate case_dir is not empty
- Check for error text in case_dir string
- Verify it's an absolute path
- Confirm directory exists after creation
- Create placeholder directory for reporting even on failures

## Enhancements Made

### A. **Configurable Paths in Docker Wrapper**

**Enhancement:** Made all paths configurable instead of hardcoded.

**Features Added:**
- Auto-detection of repository path from script location
- Command-line options: `--input-data`, `--output`, `--repo`, `--image`
- Environment variable support: `E3SM_INPUT_DATA`, `E3SM_OUTPUT`, `E3SM_REPO_PATH`, `E3SM_DOCKER_IMAGE`
- `--dry-run` option to validate configuration without running
- Smarter directory creation (try without sudo first, fall back to sudo if needed)

**Benefits:**
- Users don't need to edit the script for different paths
- Works on different systems without modification
- Easy to use with different data directories
- Better for CI/CD and automated testing

### B. **Enhanced Environment Validation**

**Enhancement:** Added comprehensive environment checks before running tests.

**Checks Added:**
- Validate repository directory exists
- Check CIME scripts directory exists
- Verify required CIME scripts are present (create_newcase, xmlquery, xmlchange)
- Confirm EAM component directory exists
- Test directory creation with proper error handling
- Better error messages indicating what's wrong and how to fix it

**Benefits:**
- Fail fast with clear error messages
- Prevent cryptic failures later in the process
- Help users diagnose configuration problems

### C. **Improved Logging and Debugging**

**Enhancement:** Better separation of output streams and log files.

**Improvements:**
- Create dedicated log files for each operation (e.g., `create_newcase_*.log`)
- Show last 10 lines of failed operation logs in error messages
- Better use of stderr vs stdout for proper command substitution
- More informative success/failure messages

**Benefits:**
- Easier to diagnose failures
- Logs are preserved for later analysis
- Error messages point directly to relevant log files

### D. **Documentation Updates**

**Enhancement:** Updated all documentation with new features and troubleshooting.

**Updates Made:**
- Added new command-line options to README
- Updated Quick Start with wrapper script examples
- Added troubleshooting section for common errors
- Updated QUICKSTART.md with configuration examples
- Documented environment variables

## Files Modified

1. **test_scripts/test_isotopes.sh**
   - Fixed `create_test_case()` function
   - Enhanced `check_environment()` function
   - Improved `run_single_test()` validation
   - Added better error handling throughout

2. **test_scripts/run_isotope_tests_docker.sh**
   - Made paths configurable
   - Added command-line argument parsing
   - Added `--dry-run` option
   - Improved directory creation logic
   - Auto-detect repository path
   - Added help text with all options

3. **test_scripts/README.md**
   - Added wrapper script usage section
   - Updated Quick Start with new methods
   - Added troubleshooting for path errors
   - Documented configuration options

4. **test_scripts/QUICKSTART.md**
   - Updated with new command-line options
   - Added environment variable examples
   - Updated directory setup section

5. **test_scripts/FIXES_APPLIED.md** (this file)
   - Documented all changes for future reference

## Testing Recommendations

### 1. Test with Default Paths
```bash
./test_scripts/run_isotope_tests_docker.sh --dry-run
```

### 2. Test with Custom Paths
```bash
./test_scripts/run_isotope_tests_docker.sh \
  --input-data ~/e3sm_data/input \
  --output ~/e3sm_data/output \
  --dry-run
```

### 3. Test Environment Variables
```bash
E3SM_INPUT_DATA=~/data/input E3SM_OUTPUT=~/data/output \
  ./test_scripts/run_isotope_tests_docker.sh --dry-run
```

### 4. Run a Single Test Interactively
```bash
./test_scripts/run_isotope_tests_docker.sh --interactive

# Then inside container:
cd /home/e3smuser/EAMv2-wiso
bash test_scripts/test_isotopes.sh
```

## What to Watch For

1. **First Run:** Create_newcase logs should show successful case creation
2. **Path Validation:** All paths should be absolute and exist
3. **Error Messages:** Should be clear and actionable, pointing to specific log files
4. **Directory Creation:** Should work without sudo when possible
5. **Case Creation:** Should complete successfully with proper validation

## Migration Notes

If you have existing scripts or workflows that use the old hardcoded paths:

1. **No changes needed** if you're using default paths (`/data/inputdata`, `/data/eamv2`)
2. **Update your scripts** if you were editing the wrapper script for custom paths:
   ```bash
   # Old way (editing script):
   # REPO_PATH="/custom/path"  # Edit the script
   # ./run_isotope_tests_docker.sh
   
   # New way (command line):
   ./run_isotope_tests_docker.sh --repo /custom/path
   
   # Or (environment variable):
   export E3SM_REPO_PATH=/custom/path
   ./run_isotope_tests_docker.sh
   ```

## Future Improvements

Potential enhancements for future versions:

1. Add `--test-subset` option to run specific tracer packages only
2. Add `--skip-build` option to reuse existing builds
3. Add `--continue-on-failure` option for debugging
4. Create a configuration file for persistent settings
5. Add progress indicators for long-running operations
6. Implement parallel test execution for multiple packages
7. Add email notifications on completion
8. Create HTML report generation

## Summary

These fixes address the immediate "Case directory does not exist" error while also making the test suite more robust, configurable, and user-friendly. The wrapper script now auto-detects paths and provides clear error messages, making it easier to use across different systems and configurations.
