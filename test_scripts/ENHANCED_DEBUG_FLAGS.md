# Enhanced Debug Flags Implementation

**Date**: 2026-04-09  
**Purpose**: Comprehensive debugging to diagnose segmentation fault in water isotope model  
**Status**: ✅ Implemented

---

## Problem Statement

The EAMv2 water isotope model was experiencing segmentation faults at runtime startup that were **not being caught** by the existing debug flags:

**Previous flags:**
```cmake
-g -Wall -fbacktrace -fcheck=bounds -ffpe-trap=zero,overflow
```

**Issue**: Silent crashes during model initialization with no diagnostic information.

---

## Solution Implemented

### Enhanced Debug Configuration (Maximum Diagnostics)

Implemented comprehensive debugging using:
1. **AddressSanitizer** - Memory error detection
2. **UndefinedBehaviorSanitizer** - Catches undefined behavior
3. **Complete initialization** - All variables initialized to trap values
4. **Comprehensive runtime checks** - Full array/pointer/memory checking

---

## Files Modified

### 1. `cime_config/machines/cmake_macros/gnu.cmake`

**Location**: Lines 19-48  
**Change**: Replaced minimal DEBUG flags with comprehensive diagnostics

#### Enhanced Flags Added:

**For C Code:**
```cmake
-fsanitize=address              # Detect heap/stack buffer overflow, use-after-free
-fsanitize=undefined            # Detect integer overflow, null pointer, alignment
-fstack-protector-all           # Stack buffer overflow protection
-fno-omit-frame-pointer         # Better stack traces
```

**For C++ Code:**
```cmake
-fsanitize=address              # Memory error detection
-fsanitize=undefined            # Undefined behavior detection
-fstack-protector-all           # Stack protection
-fno-omit-frame-pointer         # Better stack traces
```

**For Fortran Code:**
```cmake
-fcheck=all                     # Check bounds, pointer, memory, recursion, do loops
-ffpe-trap=invalid,zero,overflow,underflow  # Trap ALL floating point exceptions
-finit-real=snan                # ⭐ Initialize reals to signaling NaN (traps on use)
-finit-integer=-2147483647      # Initialize integers to obvious bad value
-finit-logical=true             # Initialize logicals
-finit-character=0              # Initialize characters
-Wuninitialized                 # Compile-time warnings for uninitialized vars
-Wmaybe-uninitialized           # Warnings for possibly uninitialized
-fstack-protector-all           # Stack protection
-fsanitize=address              # Memory error detection
-fno-omit-frame-pointer         # Better stack traces
```

**For Linking:**
```cmake
-fsanitize=address -fsanitize=undefined  # Link with sanitizer libraries
```

#### Key Improvements:

1. **`-fcheck=all`** vs `-fcheck=bounds`: 
   - Old: Only checks array bounds
   - New: Checks array bounds, pointer validity, memory allocation, recursion depth, do loop bounds

2. **`-finit-real=snan`** (Most Important):
   - Initializes ALL real/float variables to Signaling NaN
   - **Any use of uninitialized variable will immediately trap**
   - Produces exact line number and stack trace

3. **`-ffpe-trap=...underflow`**:
   - Old: Only trapped zero, overflow
   - New: Also traps underflow and invalid operations
   - Catches gradual underflow that can lead to NaN propagation

4. **AddressSanitizer (`-fsanitize=address`)**:
   - Detects out-of-bounds heap access
   - Detects stack buffer overflow
   - Detects use-after-free
   - Detects memory leaks
   - **Provides detailed error messages with exact location**

5. **UndefinedBehaviorSanitizer (`-fsanitize=undefined`)**:
   - Detects integer overflow
   - Detects null pointer dereference
   - Detects misaligned memory access
   - Detects division by zero

---

### 2. `test_scripts/test_isotopes.sh`

**Location**: `run_case()` function, lines 447-493  
**Change**: Added sanitizer environment setup and enhanced error reporting

#### Runtime Environment Configuration:

```bash
# AddressSanitizer options
export ASAN_OPTIONS="abort_on_error=1:detect_leaks=1:detect_stack_use_after_return=1:print_stacktrace=1:log_path=${case_dir}/run/asan.log"

# UndefinedBehaviorSanitizer options
export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:log_path=${case_dir}/run/ubsan.log"

# Increase stack size
ulimit -s unlimited

# Enable core dumps
ulimit -c unlimited
```

#### Enhanced Error Reporting:

- Automatically checks for sanitizer log files on failure
- Displays last 20 lines of sanitizer logs in error output
- Provides full paths to detailed log files

---

## Expected Behavior

### During Build:

**Build Time Impact:**
- Previous: ~15-20 minutes
- With sanitizers: ~45-60 minutes (3x slower)
- **Reason**: Compiler adds instrumentation to every memory access

**Compile Warnings:**
- More warnings expected (especially `-Wuninitialized`)
- These are **helpful** - they show potential issues at compile time

**Disk Space:**
- Additional ~2-3 GB for debug symbols and instrumentation

### During Runtime:

**Runtime Impact:**
- Previous: ~10-15 minutes for 5-day simulation
- With sanitizers: ~25-40 minutes (2-3x slower)
- **Reason**: Every memory access is checked at runtime

**Memory Usage:**
- Increased by 50-100% (sanitizer overhead)

### When Segfault Occurs:

**Scenario 1: Uninitialized Variable (Most Likely)**

The `-finit-real=snan` flag will cause the program to trap immediately when an uninitialized variable is used:

```
Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.

Backtrace for this error:
#0  0x7f8a1b2c3d00 in ???
#1  0x7f8a1b2c2e35 in ???
#2  0x7f8a19f5120f in ???
#3  0x559a8c4e7890 in water_tracer_init_
    at /home/e3smuser/EAMv2-wiso/components/eam/src/physics/cam/water_tracers.F90:456
#4  0x559a8c3f2a10 in phys_init_
    at /home/e3smuser/EAMv2-wiso/components/eam/src/physics/cam/physpkg.F90:234
```

**This tells you:**
- Exact file: `water_tracers.F90`
- Exact line: 456
- Exact function: `water_tracer_init_`
- Call chain: phys_init → water_tracer_init

**Scenario 2: Buffer Overflow**

AddressSanitizer will catch and report:

```
==12345==ERROR: AddressSanitizer: stack-buffer-overflow on address 0x7ffca2b4c8d0
WRITE of size 8 at 0x7ffca2b4c8d0 thread T0
    #0 0x559a8c4e7890 in water_tracer_init_ /path/to/water_tracers.F90:456
    #1 0x559a8c3f2a10 in phys_init_ /path/to/physpkg.F90:234

Address 0x7ffca2b4c8d0 is located in stack of thread T0 at offset 336 in frame
    #0 0x559a8c4e7500 in water_tracer_init_

This frame has 3 object(s):
    [32, 144) 'array1'
    [176, 288) 'array2'
    [320, 400) 'array3' <== Memory access at offset 336 is 0 bytes to the left of this variable
```

**This tells you:**
- Exact memory address and offset
- Which array is being overflowed
- How many bytes out of bounds

**Scenario 3: Undefined Behavior**

UndefinedBehaviorSanitizer will report:

```
/home/e3smuser/EAMv2-wiso/components/eam/src/physics/cam/water_tracers.F90:456:12: 
runtime error: load of value 12345678, which is not a valid value for type 'bool'
```

---

## Log Files Generated

After running tests, check these locations:

### 1. Sanitizer Logs
```
/data/eamv2/isotope_tests_YYYYMMDD_HHMMSS/isotope_test_PACKAGE/run/
├── asan.log.*    # AddressSanitizer output (if memory errors detected)
└── ubsan.log.*   # UndefinedBehaviorSanitizer output (if undefined behavior detected)
```

### 2. Standard Logs
```
/data/eamv2/isotope_tests_YYYYMMDD_HHMMSS/isotope_test_PACKAGE/
├── build.log     # Build output with enhanced warnings
├── run.log       # Runtime output
└── run/
    ├── e3sm.log.*   # Main model log
    ├── atm.log.*    # Atmosphere component log
    └── cpl.log.*    # Coupler log
```

---

## How to Use These Flags

### Running Tests:

```bash
cd /Users/rfiorella/code/E3SM/EAMv2-wiso
./test_scripts/run_isotope_tests_docker.sh
```

The enhanced flags are **automatically enabled** when running in DEBUG mode.

### Interpreting Results:

1. **Check build.log for compile-time warnings**
   - Look for `-Wuninitialized` warnings
   - These show variables that might be used before initialization

2. **If test fails, check sanitizer logs first**
   - `asan.log.*` - Memory errors (buffer overflow, use-after-free)
   - `ubsan.log.*` - Undefined behavior (integer overflow, null pointer)

3. **Look for SIGFPE (signal 8)**
   - This indicates use of uninitialized variable (NaN trap)
   - Stack trace shows exact location

4. **Check for specific error patterns**
   - "stack-buffer-overflow" → Array too small for data
   - "heap-buffer-overflow" → Heap allocation too small
   - "use-after-free" → Accessing freed memory
   - "memory leak" → Allocations not freed

---

## After Diagnosis: Optimization Options

Once the bug is found and fixed, you can optimize the flags for faster development:

### Option 1: Keep Full Diagnostics (Recommended during development)
- Keep all flags as-is
- Ensures new bugs are caught immediately

### Option 2: Moderate Diagnostics (Faster builds/runs)
Remove sanitizers but keep initialization:
```cmake
if (DEBUG)
  string(APPEND FFLAGS " -g -Wall -fbacktrace")
  string(APPEND FFLAGS " -fcheck=all")
  string(APPEND FFLAGS " -ffpe-trap=invalid,zero,overflow")
  string(APPEND FFLAGS " -finit-real=snan")  # Keep this!
  # Remove -fsanitize flags
endif()
```

### Option 3: Minimal (Production-ready)
Switch to non-DEBUG build:
```bash
./xmlchange DEBUG=FALSE
```

---

## Troubleshooting

### Issue: Build fails with sanitizer errors

**Solution**: Check if your GCC version supports sanitizers:
```bash
gcc --version  # Should be >= 4.8 for AddressSanitizer
```

If not supported, use Option 2 (Moderate Diagnostics) instead.

### Issue: "Cannot allocate memory" during runtime

**Solution**: Sanitizers increase memory usage. Try:
```bash
# Increase Docker memory limit
docker run --memory=16g ...
```

### Issue: False positives from sanitizers

**Solution**: Suppress specific issues using suppression files (consult AddressSanitizer documentation).

### Issue: Build is too slow

**Solution**: 
1. Use parallel builds: `./case.build -j 8`
2. Switch to Option 2 (Moderate Diagnostics)
3. Build without sanitizers for routine testing

---

## References

### GCC Documentation:
- [AddressSanitizer](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html#index-fsanitize_003daddress)
- [UndefinedBehaviorSanitizer](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html#index-fsanitize_003dundefined)
- [Fortran Debugging Options](https://gcc.gnu.org/onlinedocs/gfortran/Debugging-Options.html)

### Sanitizer Runtime Options:
- [ASAN_OPTIONS](https://github.com/google/sanitizers/wiki/AddressSanitizerFlags)
- [UBSAN_OPTIONS](https://github.com/google/sanitizers/wiki/UndefinedBehaviorSanitizerFlags)

---

## Summary

**What Changed:**
- Enhanced debug flags in `gnu.cmake`
- Sanitizer environment in `test_isotopes.sh`

**What You Get:**
- Exact line number of segfault
- Type of error (uninitialized, overflow, etc.)
- Complete stack trace
- Memory addresses and context

**Next Steps:**
1. Clean previous builds: `rm -rf /data/eamv2/isotope_tests_*`
2. Run tests: `./test_scripts/run_isotope_tests_docker.sh`
3. Wait for build (~60 min) and test run (~30 min)
4. Check sanitizer logs for detailed diagnostics
5. Fix the identified issue
6. Verify fix with another test run

**Success Criteria:**
Model runs without segfault and produces output files with isotope tracers.

---

**Implemented by**: OpenCode AI Assistant  
**Date**: 2026-04-09  
**Version**: 1.0
