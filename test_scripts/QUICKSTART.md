# EAMv2 Water Isotope Test Suite - Quick Reference

## Quick Start

```bash
# Navigate to repository
cd /path/to/EAMv2-wiso

# Run tests (simplest method - auto-detects paths)
./test_scripts/run_isotope_tests_docker.sh

# Or with custom paths
./test_scripts/run_isotope_tests_docker.sh \
  --input-data ~/e3sm_data/input \
  --output ~/e3sm_data/output

# Check configuration first (dry run)
./test_scripts/run_isotope_tests_docker.sh --dry-run

# Or start interactive session
./test_scripts/run_isotope_tests_docker.sh --interactive
```

## Configuration Options

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input-data PATH` | Input data directory | `/data/inputdata` |
| `--output PATH` | Output directory | `/data/eamv2` |
| `--repo PATH` | Repository path | Auto-detected |
| `--image IMAGE` | Docker image | `rfiorella/model-containers:e3sm-openmpi-dev-latest` |
| `--dry-run` | Show config without running | - |
| `--interactive` | Start interactive shell | - |

### Environment Variables

```bash
export E3SM_INPUT_DATA=~/e3sm_data/input
export E3SM_OUTPUT=~/e3sm_data/output
export E3SM_REPO_PATH=/custom/path/to/EAMv2-wiso
export E3SM_DOCKER_IMAGE=custom/image:tag

./test_scripts/run_isotope_tests_docker.sh
```

## Files Created

```
test_scripts/
├── README.md                          # Full documentation
├── test_isotopes.sh                   # Main test suite (runs in container)
├── validate_isotope_output.py         # Output validation script
├── run_isotope_tests_docker.sh        # Docker wrapper (run on host)
├── FUTURE_VALIDATION_PHASES.md        # Future testing phases documentation
└── QUICKSTART.md                      # This file
```

## Test Configuration

| Aspect | Setting |
|--------|---------|
| **Compset** | Fi2000climo (2000 climatology with water isotopes) |
| **Resolution** | ne30pg2_oECv3 (standard E3SM low-res) |
| **Duration** | 5 days per test |
| **Packages Tested** | h2o_hdo, h2o_h218o, h2o_hdo_h218o, all_stable_wiso |
| **Total Tests** | 4 |
| **Estimated Time** | 2-6 hours (first run slower due to downloads) |

## Directory Setup

The script auto-detects the repository and allows configurable paths:

| Path Type | Default | Configurable Via |
|-----------|---------|------------------|
| **Repository** | Auto-detected from script location | `--repo` or `E3SM_REPO_PATH` |
| **Input Data** | `/data/inputdata` | `--input-data` or `E3SM_INPUT_DATA` |
| **Output** | `/data/eamv2` | `--output` or `E3SM_OUTPUT` |

**Note:** Directories will be created automatically if they don't exist.

## Output Location

Results are saved to the configured output directory (default: `/data/eamv2`):

```
${OUTPUT_DIR}/isotope_tests_YYYYMMDD_HHMMSS/
├── create_newcase_*.log               # Case creation logs
├── isotope_test_h2o_hdo/
├── isotope_test_h2o_h218o/
├── isotope_test_h2o_hdo_h218o/
├── isotope_test_all_stable_wiso/
└── test_summary.txt                   # Overall results
```

## Checking Results

### Quick Check
```bash
# View summary report
cat /data/eamv2/isotope_tests_*/test_summary.txt
```

### Detailed Validation
```bash
# Validate specific test case
docker run -it \
  -v /data/eamv2:/home/e3smuser/output \
  rfiorella/model-containers:e3sm-openmpi-dev-latest \
  python3 /home/e3smuser/EAMv2-wiso/test_scripts/validate_isotope_output.py \
    /home/e3smuser/output/isotope_tests_YYYYMMDD_HHMMSS/isotope_test_h2o_hdo \
    h2o_hdo
```

## Common Issues

### Docker not running
```bash
# Start Docker Desktop or daemon
# On macOS: open Docker Desktop app
# On Linux: sudo systemctl start docker
```

### Permission errors on /data directories
```bash
sudo chown -R $USER:$USER /data/inputdata /data/eamv2
```

### Input data download timeout
```bash
# The first run downloads 100-500 GB of data
# This can take several hours - be patient!
# Data is cached in /data/inputdata for subsequent runs
```

## Test Success Criteria

- ✓ Model builds with isotope tracers
- ✓ Simulation completes without errors  
- ✓ Output files contain isotope variables
- ✓ Tracers have non-zero values

## What's Tested

### Tracer Packages

1. **h2o_hdo** (14 tracers)
   - H₂O and HDO (deuterium)
   - Tests: Q_H2O, Q_HDO, CLDLIQ_H2O, CLDLIQ_HDO, etc.

2. **h2o_h218o** (14 tracers)
   - H₂O and H₂¹⁸O (oxygen-18)
   - Tests: Q_H2O, Q_H218O, CLDLIQ_H2O, CLDLIQ_H218O, etc.

3. **h2o_hdo_h218o** (21 tracers)
   - H₂O, HDO, and H₂¹⁸O
   - Tests: Q_H2O, Q_HDO, Q_H218O, etc.

4. **all_stable_wiso** (35 tracers)
   - All stable isotopes: H₂O, H₂¹⁶O, HDO, H₂¹⁸O, H₂¹⁷O
   - Tests: Q_H2O, Q_H216O, Q_HDO, Q_H218O, Q_H217O, etc.

## Next Steps

After successful Phase 1 tests, consider implementing:

- **Phase 2:** Physical realism checks (δD, δ¹⁸O patterns)
- **Phase 3:** Mass conservation validation
- **Phase 4:** Performance benchmarking
- **Phase 5:** Scientific validation vs observations

See `FUTURE_VALIDATION_PHASES.md` for details.

## Getting Help

- **Full Documentation:** `test_scripts/README.md`
- **Future Plans:** `test_scripts/FUTURE_VALIDATION_PHASES.md`
- **Issues/Questions:** rfiorella@lanl.gov

## Example Session

```bash
# Complete workflow from start to finish

# 1. Ensure Docker is running
docker info

# 2. Navigate to repository
cd /Users/rfiorella/code/E3SM/EAMv2-wiso

# 3. Run tests
./test_scripts/run_isotope_tests_docker.sh

# 4. Wait 2-6 hours...

# 5. Check results
cat /data/eamv2/isotope_tests_*/test_summary.txt

# 6. Review detailed output
ls -lh /data/eamv2/isotope_tests_*/*/run/*.eam.h0.*.nc
```

## Monitoring Progress

### During Build
```bash
# In another terminal, watch build progress
docker ps  # Get container ID
docker exec -it <container-id> tail -f /home/e3smuser/output/isotope_tests_*/isotope_test_*/build.log
```

### During Run
```bash
# Watch run progress
docker exec -it <container-id> tail -f /home/e3smuser/output/isotope_tests_*/isotope_test_*/run/atm.log.*
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Docker not found | Install Docker Desktop or Docker Engine |
| Permission denied | Run with sudo or fix directory permissions |
| Out of disk space | Clear space or use different mount points |
| Build fails | Check `build.log` in test case directory |
| Run fails | Check `atm.log.*` and `cpl.log.*` in run directory |
| Missing tracers | Verify water_tracer was set correctly |

## Performance Tips

- **First run:** Budget 4-6 hours (includes downloads)
- **Subsequent runs:** ~2-3 hours (cached input data)
- **Parallel execution:** Tests run sequentially by design (resource management)
- **Memory:** Ensure Docker has at least 8 GB RAM allocated
- **CPU:** 4+ cores recommended

## Key Configuration Points

### Build-time
- Water tracer package set via: `CAM_CONFIG_OPTS="-water_tracer <package>"`
- Preprocessor defines: `SET_WTRC_MAX_CNST` (automatic)

### Run-time (namelists)
```fortran
trace_water = .true.        ! Enable tracers
wisotope = .true.           ! Enable fractionation
wtrc_alpha_kinetic = .true. ! Kinetic fractionation
```

### Coupler
```xml
FLDS_WISO = TRUE  ! Pass isotope fields between components
```

---

**Last Updated:** 2026-04-07  
**Version:** 1.0
