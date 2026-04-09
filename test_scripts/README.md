# EAMv2 Water Isotope Test Suite

Automated testing framework for validating water isotope functionality in EAMv2.

## Overview

This test suite verifies that the water isotope/tracer implementation in EAMv2 is working correctly by:
1. Building the model with different isotope tracer packages
2. Running short (5-day) simulations with the Fi2000climo compset
3. Validating that isotope tracer fields are present and active in the output

## Test Configurations

The suite tests four different water isotope tracer packages:

| Package | Isotopes | Tracers | Description |
|---------|----------|---------|-------------|
| `h2o_hdo` | H₂O, HDO | 14 | Deuterium tracking |
| `h2o_h218o` | H₂O, H₂¹⁸O | 14 | Oxygen-18 tracking |
| `h2o_hdo_h218o` | H₂O, HDO, H₂¹⁸O | 21 | Combined deuterium and oxygen-18 |
| `all_stable_wiso` | All stable isotopes | 35 | H₂O, H₂¹⁶O, HDO, H₂¹⁸O, H₂¹⁷O |

Each configuration runs a 5-day simulation at ne30pg2 resolution using the Fi2000climo compset.

## Requirements

### Docker Container
- **Container:** `rfiorella/model-containers:e3sm-openmpi-dev-latest`
- Contains all necessary build tools and libraries

### Host System
- Docker installed and running
- Sufficient disk space:
  - Input data: ~100-500 GB (downloaded on first run, reused subsequently)
  - Output data: ~20-40 GB per test (4 tests total)
  - Build artifacts: ~10-20 GB per test

### Volume Mounts
Three directories must be mounted:
1. **EAMv2-wiso repository** → `/home/e3smuser/EAMv2-wiso`
2. **Input data** → `/home/e3smuser/inputdata` (mount `/data/inputdata`)
3. **Output directory** → `/home/e3smuser/output` (mount `/data/eamv2`)

## Quick Start

### Method 1: Using the Wrapper Script (Recommended)

The easiest way to run the tests is using the provided wrapper script, which handles Docker configuration automatically:

```bash
# Navigate to repository
cd /path/to/EAMv2-wiso

# Run tests with default settings
./test_scripts/run_isotope_tests_docker.sh

# Or with custom paths
./test_scripts/run_isotope_tests_docker.sh \
  --input-data ~/e3sm_data/input \
  --output ~/e3sm_data/output

# Check configuration without running
./test_scripts/run_isotope_tests_docker.sh --dry-run

# Start interactive shell for debugging
./test_scripts/run_isotope_tests_docker.sh --interactive
```

**Configuration Options:**
- `--input-data PATH` - Set input data directory (default: `/data/inputdata`)
- `--output PATH` - Set output directory (default: `/data/eamv2`)
- `--repo PATH` - Set repository path (default: auto-detected)
- `--image IMAGE` - Set Docker image (default: `rfiorella/model-containers:e3sm-openmpi-dev-latest`)
- `--dry-run` - Show configuration without running
- `--interactive` - Start interactive shell instead of running tests

**Environment Variables:**
You can also configure paths using environment variables:
```bash
export E3SM_INPUT_DATA=~/e3sm_data/input
export E3SM_OUTPUT=~/e3sm_data/output
./test_scripts/run_isotope_tests_docker.sh
```

### Method 2: Manual Docker Command

If you prefer manual control or need custom Docker options:

```bash
# Create directories (if not using defaults)
mkdir -p /data/inputdata /data/eamv2

# Run tests directly
docker run -it \
  -v $(pwd):/home/e3smuser/EAMv2-wiso \
  -v /data/inputdata:/home/e3smuser/inputdata \
  -v /data/eamv2:/home/e3smuser/output \
  --hostname docker \
  rfiorella/model-containers:e3sm-openmpi-dev-latest \
  bash /home/e3smuser/EAMv2-wiso/test_scripts/test_isotopes.sh
```

**Note:** The first run will download input data (~100-500 GB) which may take several hours. Subsequent runs will reuse the cached data.

## Alternative: Interactive Mode

For debugging or manual testing:

```bash
# Start container interactively
docker run -it \
  -v /Users/rfiorella/code/E3SM/EAMv2-wiso:/home/e3smuser/EAMv2-wiso \
  -v /data/inputdata:/home/e3smuser/inputdata \
  -v /data/eamv2:/home/e3smuser/output \
  --hostname docker \
  rfiorella/model-containers:e3sm-openmpi-dev-latest \
  bash

# Inside container, run the test suite
bash /home/e3smuser/EAMv2-wiso/test_scripts/test_isotopes.sh

# Or run individual components manually
cd /home/e3smuser/EAMv2-wiso/cime/scripts
./create_newcase --case /home/e3smuser/output/test_case --compset Fi2000climo --res ne30pg2_oECv3 --mach docker
# ... continue with manual configuration
```

## Test Workflow

Each test follows this sequence:

1. **Create Case**
   - Uses CIME `create_newcase` script
   - Compset: Fi2000climo (2000 climatology with water isotopes)
   - Resolution: ne30pg2_oECv3 (standard E3SM low-res grid)

2. **Configure Case**
   - Sets `CAM_CONFIG_OPTS` with appropriate `-water_tracer <package>`
   - Enables `FLDS_WISO=TRUE` for coupler isotope field passing
   - Sets 5-day run duration
   - Configures namelist options (trace_water, wisotope, etc.)

3. **Setup Case**
   - Runs `case.setup` to prepare build environment

4. **Check Input Data**
   - Downloads required input data if not already present
   - First run takes longest due to downloads

5. **Build Model**
   - Compiles EAM with isotope tracers
   - First build: ~45-60 minutes
   - Subsequent builds: ~8-15 minutes (with caching)

6. **Run Simulation**
   - Executes 5-day simulation
   - Runtime: ~10-20 minutes per test

7. **Validate Output**
   - Checks for successful completion
   - Verifies presence of isotope tracer variables
   - Confirms tracers have non-zero values

## Output Structure

Test results are saved to `/data/eamv2/isotope_tests_YYYYMMDD_HHMMSS/`:

```
isotope_tests_20260407_143000/
├── isotope_test_h2o_hdo/           # Test case for h2o_hdo
│   ├── run/                        # Model output
│   │   ├── *.eam.h0.*.nc          # History files with isotope data
│   │   ├── atm.log.*              # Atmosphere component log
│   │   └── cpl.log.*              # Coupler log
│   ├── bld/                        # Build artifacts
│   ├── build.log                   # Build log
│   ├── run.log                     # Run log
│   └── TEST_STATUS                 # PASS or FAIL
├── isotope_test_h2o_h218o/         # Test case for h2o_h218o
├── isotope_test_h2o_hdo_h218o/     # Test case for h2o_hdo_h218o
├── isotope_test_all_stable_wiso/   # Test case for all_stable_wiso
├── test_summary.txt                # Overall test summary
└── test_suite.log                  # Combined test suite log
```

## Validation Criteria

### Phase 1: Basic Smoke Tests (Current)

✓ **Build Success:** Model compiles with isotope tracers  
✓ **Run Completion:** Simulation completes without crashes  
✓ **Tracer Presence:** Expected isotope variables exist in output  
✓ **Tracer Activity:** Variables contain non-zero values

### Phase 2: Physical Validation (Future)

- δD and δ18O isotopic ratios are physically reasonable
- Meteoric water line relationship (δD ≈ 8×δ18O + 10)
- Spatial patterns (polar depletion, tropical enrichment)
- d-excess parameter validation

### Phase 3: Mass Conservation (Future)

- Total water = sum of isotope tracers
- No unphysical negative values
- Global water budget closure
- Conservation across component boundaries

### Phase 4: Performance (Future)

- Runtime comparison with/without isotopes
- Memory usage profiling
- Scalability testing

### Phase 5: Scientific Validation (Future)

- Comparison with GNIP precipitation observations
- Validation against published isotope studies
- Long-term climate simulations (1-5 years)

## Interpreting Results

### Successful Test

```
===============================================================================
Test 1/4: h2o_hdo (14 tracers)
===============================================================================

[INFO] Creating test case for h2o_hdo
[SUCCESS] Case created: /home/e3smuser/output/isotope_tests_*/isotope_test_h2o_hdo
[INFO] Configuring case for h2o_hdo
[SUCCESS] Case configured for h2o_hdo
[INFO] Building case for h2o_hdo (this may take 30-60 minutes)...
[SUCCESS] Build complete in 45 minutes
[INFO] Running 5-day simulation for h2o_hdo...
[SUCCESS] Run complete in 12 minutes
[INFO] Validating output for h2o_hdo
[SUCCESS] Model run completed successfully
[INFO] Found 6 history file(s)
[SUCCESS] Found tracer: Q_H2O
[SUCCESS] Found tracer: Q_HDO
[SUCCESS] All expected isotope tracers found
[SUCCESS] Test PASSED for h2o_hdo (57 minutes)
```

### Failed Test

Check the following if a test fails:

1. **Build Failure:**
   - Check `bld/atm.bldlog.*` for compilation errors
   - Verify water_tracer package is recognized
   - Look for missing dependencies

2. **Run Failure:**
   - Check `run/atm.log.*` for runtime errors
   - Review `run/cpl.log.*` for coupling issues
   - Verify input data was downloaded correctly

3. **Validation Failure:**
   - Check if model completed (look for "SUCCESSFUL TERMINATION")
   - Verify history files exist in `run/` directory
   - Use `validate_isotope_output.py` for detailed diagnostics

## Manual Validation

For detailed validation of a specific test case:

```bash
# Inside container or with ncdump available
cd /data/eamv2/isotope_tests_YYYYMMDD_HHMMSS/isotope_test_h2o_hdo

# Check run completion
grep "SUCCESSFUL TERMINATION" run/cpl.log.*

# List isotope variables
ncdump -h run/*.eam.h0.*.nc | grep -E "Q_H2O|Q_HDO|Q_H218O"

# Use validation script
python /home/e3smuser/EAMv2-wiso/test_scripts/validate_isotope_output.py \
    /data/eamv2/isotope_tests_YYYYMMDD_HHMMSS/isotope_test_h2o_hdo \
    h2o_hdo
```

## Troubleshooting

### Issue: "Case directory does not exist" error

**Solution:** This usually indicates a problem during case creation. Check:
1. The create_newcase log: `/path/to/output/isotope_tests_*/create_newcase_*.log`
2. Ensure CIME scripts directory exists: `ls -la /home/e3smuser/EAMv2-wiso/cime/scripts`
3. Verify repository is properly mounted in Docker
4. Run with `--dry-run` to check configuration

### Issue: "cd: --: invalid option" or path errors

**Solution:** This was fixed in the latest version. Make sure you're using the updated scripts:
```bash
cd /path/to/EAMv2-wiso
git pull  # Get latest fixes
./test_scripts/run_isotope_tests_docker.sh --dry-run  # Verify config
```

### Issue: Input data download fails or times out

**Solution:** Input data download can take several hours. If it fails:
```bash
# Download manually using CIME tools
cd /home/e3smuser/output/isotope_test_h2o_hdo
./check_input_data --download
```

### Issue: Permission denied when creating directories

**Solution:** The wrapper script will try to create directories automatically. If this fails:
```bash
# Use custom paths in your home directory
./test_scripts/run_isotope_tests_docker.sh \
  --input-data ~/e3sm_data/input \
  --output ~/e3sm_data/output

# Or create directories manually
mkdir -p /data/inputdata /data/eamv2
# May need sudo if creating in /data
```

### Issue: Build fails with "water_tracer not recognized"

**Solution:** Ensure CAM_CONFIG_OPTS includes the `-water_tracer` flag:
```bash
./xmlquery CAM_CONFIG_OPTS
# Should show: ... -water_tracer h2o_hdo
```

### Issue: Out of memory during build or run

**Solution:** Increase Docker memory limit:
```bash
# In Docker Desktop: Settings → Resources → Memory
# Or via command line with --memory flag
docker run --memory=16g ...
```

### Issue: NODENAME_REGEX doesn't match

**Solution:** Ensure `--hostname docker` is set in docker run command, matching config_machines.xml

### Issue: Permission denied errors in mounted volumes

**Solution:** Check directory permissions on host:
```bash
sudo chown -R $USER:$USER /data/inputdata /data/eamv2
```

## Configuration Files

Key configuration files in the repository:

- **`components/eam/bld/configure`** - Build-time water_tracer configuration
- **`components/eam/cime_config/config_compsets.xml`** - WISO compset definitions
- **`components/eam/bld/config_files/definition.xml`** - Valid water_tracer packages
- **`/code/model-containers/cime_files/config_machines.xml`** - Docker machine config

## Isotope Physics Code

Core implementation files:

- **`share/util/water_isotopes.F90`** - Fractionation calculations
- **`share/util/water_types.F90`** - Water phase type definitions
- **`components/eam/src/physics/cam/water_tracers.F90`** - Main tracer implementation
- **`components/eam/src/physics/cam/water_tracer_vars.F90`** - Configuration variables

## References

### Water Isotope Science
- Noone, D. (2003). *Water isotope tracing in climate models*
- Nusbaumer, J. (2011). *CAM5 isotope implementation*
- Jouzel, J. & Merlivat, L. (1984). *Deuterium and oxygen 18 in precipitation*
- Majoube, M. (1971). *Fractionation factors for O18 and D*

### E3SM Documentation
- E3SM User Guide: https://e3sm.org/
- CIME Documentation: https://esmci.github.io/cime/

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review logs in test case directories
3. Consult E3SM documentation and forums
4. Contact: rfiorella@lanl.gov

## License

This test suite is part of the EAMv2-wiso project and follows the same license as E3SM.

---

**Last Updated:** 2026-04-07  
**Test Suite Version:** 1.0
