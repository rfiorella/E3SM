# C++ and Kokkos Primer for Fortran Developers

This guide explains C++ and Kokkos concepts for someone familiar with CAM's Fortran codebase. It covers the most important differences you'll encounter when working with EAMxx water tracers.

## Table of Contents
1. [Basic Syntax Differences](#basic-syntax-differences)
2. [Modules vs Namespaces](#modules-vs-namespaces)
3. [Arrays and Indexing](#arrays-and-indexing)
4. [Subroutines vs Functions](#subroutines-vs-functions)
5. [Templates (The Big One!)](#templates-the-big-one)
6. [Kokkos Device Functions](#kokkos-device-functions)
7. [Pack Types for Vectorization](#pack-types-for-vectorization)
8. [Memory and References](#memory-and-references)
9. [Constants and Parameters](#constants-and-parameters)
10. [Quick Reference Table](#quick-reference-table)

---

## Basic Syntax Differences

### Comments
```fortran
! Fortran comment
```
```cpp
// C++ comment (single line)
/* C++ comment (multi-line) */
```

### Variable Declaration
```fortran
! Fortran
real(r8) :: temperature
integer :: count
logical :: is_frozen
```
```cpp
// C++
Real temperature;              // Real is typedef for double/float
int count;
bool is_frozen;               // bool instead of logical

// With initialization
Real temperature = 280.0;
int count = 0;
bool is_frozen = false;       // false/true instead of .false./.true.
```

### Logical Operators
```fortran
! Fortran
if (a .and. b .or. .not. c) then
  ...
end if
```
```cpp
// C++
if (a && b || !c) {
  ...
}
```

| Fortran | C++ |
|---------|-----|
| `.and.` | `&&` |
| `.or.` | `\|\|` |
| `.not.` | `!` |
| `.eq.` or `==` | `==` |
| `.ne.` or `/=` | `!=` |
| `.lt.` or `<` | `<` |
| `.le.` or `<=` | `<=` |

### Control Flow
```fortran
! Fortran IF
if (condition) then
  statement1
else if (condition2) then
  statement2
else
  statement3
end if
```
```cpp
// C++ if
if (condition) {
  statement1;
} else if (condition2) {
  statement2;
} else {
  statement3;
}
// Note: semicolons at end of statements!
//       curly braces {} group statements
```

```fortran
! Fortran DO loop
do i = 1, 10
  array(i) = i * 2
end do
```
```cpp
// C++ for loop (NOTE: 0-indexed!)
for (int i = 0; i < 10; ++i) {
  array[i] = i * 2;
}
// i starts at 0, not 1!
// ++i increments i (like i = i + 1)
```

---

## Modules vs Namespaces

Fortran modules and C++ namespaces organize code similarly, but with different syntax.

### Fortran Module
```fortran
module water_isotopes
  use shr_kind_mod, only: r8 => shr_kind_r8
  
  implicit none
  private
  
  ! Public interface
  public :: wiso_alpl
  public :: wiso_alpi
  
  integer, parameter :: pwtspec = 6
  
contains

  function wiso_alpl(isp, tk) result(alpha)
    integer, intent(in) :: isp
    real(r8), intent(in) :: tk
    real(r8) :: alpha
    alpha = ...
  end function wiso_alpl
  
end module water_isotopes
```

### C++ Namespace
```cpp
namespace scream {
namespace water_isotopes {

// No "implicit none" - C++ always requires declarations
// No "private/public" - everything in namespace is accessible

// Constants
constexpr int pwtspec = 6;

// Function declaration (like "public ::")
// Templates replace Fortran interfaces
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk);

// Function implementation
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk) {
  ScalarT alpha = ...;
  return alpha;
}

} // namespace water_isotopes
} // namespace scream
```

### Using Modules/Namespaces
```fortran
! Fortran
use water_isotopes, only: wiso_alpl, pwtspec

call wiso_alpl(isp, tk, alpha)
```
```cpp
// C++
using scream::water_isotopes::wiso_alpl;
using scream::water_isotopes::pwtspec;

alpha = wiso_alpl(isp, tk);  // Function returns value

// OR use namespace for everything
using namespace scream::water_isotopes;
alpha = wiso_alpl(isp, tk);

// OR use explicit qualification
alpha = scream::water_isotopes::wiso_alpl(isp, tk);
```

---

## Arrays and Indexing

### **CRITICAL DIFFERENCE: C++ arrays start at index 0, not 1!**

```fortran
! Fortran: 1-indexed
real(r8), dimension(6) :: array
array(1) = 1.0  ! First element
array(6) = 6.0  ! Last element (sixth element)

do i = 1, 6
  print *, array(i)
end do
```

```cpp
// C++: 0-indexed
std::array<Real,6> array;
array[0] = 1.0;  // First element
array[5] = 6.0;  // Last element (sixth element)

// IMPORTANT: Use square brackets [], not parentheses ()

for (int i = 0; i < 6; ++i) {
  std::cout << array[i] << std::endl;
}
```

### Converting Fortran Array Index to C++
If Fortran uses index `i` (1 to N), C++ uses index `i-1` (0 to N-1):

```fortran
! Fortran
integer, parameter :: isphdo = 3
real(r8) :: value
value = difrm(isphdo)  ! Gets third element
```
```cpp
// C++
constexpr int isphdo_fortran = 3;  // Keep Fortran value for reference
constexpr int isphdo = 2;          // C++ index (3-1 = 2)
Real value = difrm[isphdo];        // Gets third element (index 2)

// OR directly use:
Real value = difrm[2];  // Third element
```

### Array Types in C++
```cpp
// Fixed-size array (like Fortran dimension(N))
std::array<Real, 6> array;

// Dynamic-size array (like Fortran allocatable)
std::vector<Real> array(6);  // 6 elements

// C-style array (avoid! use std::array instead)
Real array[6];  // Old style, less safe
```

---

## Subroutines vs Functions

Fortran has both subroutines (no return value) and functions (return value).  
C++ only has functions, which can return values or return nothing (`void`).

### Fortran Subroutine with Output
```fortran
subroutine wiso_kmol(isp, rbot, zbot, ustar, alpkn)
  integer, intent(in) :: isp
  real(r8), intent(in) :: rbot, zbot, ustar
  real(r8), intent(out) :: alpkn
  
  alpkn = ...  ! Set output value
end subroutine wiso_kmol
```

### C++ Function Returning Value
```cpp
template<typename ScalarT>
ScalarT wiso_kmol(int isp, const ScalarT& rbot,
                  const ScalarT& zbot, const ScalarT& ustar) {
  ScalarT alpkn = ...;
  return alpkn;  // Explicit return
}

// Usage:
Real result = wiso_kmol(isp, rbot, zbot, ustar);
```

### Fortran Function
```fortran
function wiso_alpl(isp, tk) result(alpha)
  integer, intent(in) :: isp
  real(r8), intent(in) :: tk
  real(r8) :: alpha
  alpha = ...
end function wiso_alpl
```

### C++ Equivalent
```cpp
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk) {
  ScalarT alpha = ...;
  return alpha;
}
```

### Parameter Intent Mapping

| Fortran | C++ | Meaning |
|---------|-----|---------|
| `intent(in)` | `const Type&` | Input only, passed by reference |
| `intent(out)` | Return value or `Type&` | Output, return or pass by reference |
| `intent(inout)` | `Type&` | Modified, pass by reference |

---

## Templates (The Big One!)

**Templates are C++'s most powerful feature for generic programming.** They're more powerful than Fortran interfaces.

### The Problem Templates Solve

In Fortran, you'd write:
```fortran
! One version for scalars
function wiso_alpl_scalar(isp, tk) result(alpha)
  real(r8), intent(in) :: tk
  real(r8) :: alpha
  alpha = exp(...)
end function

! Another version for arrays (element-wise)
function wiso_alpl_array(isp, tk) result(alpha)
  real(r8), intent(in) :: tk(:)
  real(r8) :: alpha(size(tk))
  integer :: i
  do i = 1, size(tk)
    alpha(i) = exp(...)
  end do
end function

! Interface to combine them
interface wiso_alpl
  module procedure wiso_alpl_scalar
  module procedure wiso_alpl_array
end interface
```

In C++ with templates:
```cpp
// ONE function works for scalars AND packs!
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk) {
  ScalarT alpha = ekat::impl::exp(...);
  return alpha;
}

// Compiler generates versions automatically:
Real scalar_result = wiso_alpl(isp, tk_scalar);
Pack<Real,16> pack_result = wiso_alpl(isp, tk_pack);
```

### Template Syntax Explained

```cpp
template<typename ScalarT>
//       ^^^^^^^^ "typename" declares a type parameter
//               ^^^^^^^ ScalarT is our chosen name (could be T, Type, etc.)

ScalarT wiso_alpl(int isp, const ScalarT& tk) {
// ^^^^^^ Return type is ScalarT (whatever type was passed in)
//                             ^^^^^^^ Parameter type is also ScalarT
  
  ScalarT alpha = 1.0;  // Variable type is ScalarT
  return alpha;
}
```

**Think of templates like this:**
- Template parameter is like a "type variable"
- Compiler substitutes actual types when you use the function
- One definition → many specialized versions generated automatically

### Why Templates Are Better Than Fortran Interfaces

1. **Less code duplication**: Write once, works for many types
2. **Type safety**: Compiler checks everything
3. **Works with user-defined types**: Not just built-in types
4. **Optimized per type**: Each type gets its own optimized version

---

## Kokkos Device Functions

**Goal**: Write code once that runs on both CPU and GPU.

### The Problem

```fortran
! Fortran with OpenACC
!$acc routine seq
function wiso_alpl(isp, tk) result(alpha)
  !$acc declare device_resident(...)
  real(r8), intent(in) :: tk
  real(r8) :: alpha
  alpha = exp(...)
end function
```

### Kokkos Solution

```cpp
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION  // <-- Magic happens here!
ScalarT wiso_alpl(int isp, const ScalarT& tk) {
  ScalarT alpha = ekat::impl::exp(...);  // GPU-compatible exp
  return alpha;
}
```

### What `KOKKOS_INLINE_FUNCTION` Does

- **On CPU**: Expands to `inline` (optimization hint)
- **On CUDA GPU**: Expands to `__device__ __host__ inline`
- **On HIP GPU**: Expands to `__device__ __host__ inline`
- **Result**: Function works on both CPU and GPU!

### Rules for Device Functions

1. **Mark with** `KOKKOS_INLINE_FUNCTION`
2. **Use Kokkos/EKat math functions**:
   - ❌ `std::exp(x)` - CPU only
   - ✅ `ekat::impl::exp(x)` - CPU and GPU
   - ✅ `Kokkos::exp(x)` - CPU and GPU

3. **Can't use**:
   - `printf`, `std::cout` (use `Kokkos::printf` on device)
   - `std::vector`, `std::string` (use Kokkos::View instead)
   - File I/O

4. **Can use**:
   - Math operations (`+`, `-`, `*`, `/`)
   - Conditionals (`if`, `for`)
   - `constexpr` constants
   - Other `KOKKOS_INLINE_FUNCTION` functions

---

## Pack Types for Vectorization

**Packs enable SIMD (Single Instruction Multiple Data) operations.**

### Fortran Approach
You operate on arrays element-by-element, and compiler tries to vectorize:
```fortran
real(r8) :: tk(100), alpha(100)
do i = 1, 100
  alpha(i) = exp(alpal/tk(i)**3 + ...)
end do
! Compiler may vectorize this loop (or may not)
```

### C++ Pack Approach
Explicitly operate on "packs" of values:
```cpp
// Pack of 16 doubles processed together
using Pack16 = ekat::Pack<Real, 16>;

Pack16 tk = ...;  // 16 temperature values
Pack16 alpha = ekat::impl::exp(alpal / (tk*tk*tk) + ...);
// This processes 16 values in parallel with SIMD instructions!
```

### What Is a Pack?

```cpp
ekat::Pack<Real, N>  // N scalars processed together
// N is typically 1, 4, 8, or 16 depending on architecture
```

**Think of it like this:**
- Fortran array: `real(r8), dimension(N)`
- C++ Pack: `Pack<Real,N>`

**But**:
- Fortran array: N elements processed sequentially (maybe vectorized by compiler)
- C++ Pack: N elements processed in PARALLEL using SIMD instructions

### Why Templates + Packs = Awesome

Same function works for both!
```cpp
template<typename ScalarT>
ScalarT wiso_alpl(int isp, const ScalarT& tk) {
  ScalarT alpha = ekat::impl::exp(...);
  return alpha;
}

// With scalars
Real tk_scalar = 280.0;
Real alpha_scalar = wiso_alpl(2, tk_scalar);

// With packs (16 values at once!)
Pack<Real,16> tk_pack = ...;  // 16 temperatures
Pack<Real,16> alpha_pack = wiso_alpl(2, tk_pack);
// Processes 16 temperatures in one call!
```

### Pack Math Operations

Packs support all standard operations:
```cpp
Pack<Real,16> a, b, c;
c = a + b;          // Add 16 pairs
c = a * b;          // Multiply 16 pairs
c = a / b;          // Divide 16 pairs
c = exp(a);         // Exp of 16 values
c = min(a, b);      // Min of 16 pairs
```

---

## Memory and References

### Pass by Value vs Pass by Reference

```fortran
! Fortran: intent controls copying
subroutine foo(x, y)
  real(r8), intent(in) :: x     ! Can't modify (but might be copied)
  real(r8), intent(out) :: y    ! Will modify
end subroutine
```

```cpp
// C++: & controls references

// Pass by value (copies the variable)
void foo(Real x) {
  x = x + 1;  // Modifies COPY, not original
}

// Pass by reference (no copy, can modify original)
void foo(Real& x) {
  x = x + 1;  // Modifies ORIGINAL
}

// Pass by const reference (no copy, can't modify)
void foo(const Real& x) {
  // x = x + 1;  // ERROR! Can't modify const reference
  Real y = x + 1;  // OK, read x to make new variable
}
```

### Best Practice: Use const&

For input parameters, always use `const &`:
```cpp
// GOOD: Efficient, safe
Real wiso_alpl(int isp, const Real& tk) {
  ...
}

// BAD: Copies large objects (inefficient)
Real wiso_alpl(int isp, Real tk) {
  ...
}
```

**Rule of thumb:**
- **Small types** (int, bool): pass by value
- **Large types** (Real, Pack, arrays): pass by `const &` for input

---

## Constants and Parameters

### Compile-Time Constants

```fortran
! Fortran
integer, parameter :: pwtspec = 6
real(r8), parameter :: dkfac = 0.58_r8
real(r8), dimension(pwtspec), parameter :: &
    difrm = (/ 1.0_r8, 1.0_r8, 0.9757_r8, ... /)
```

```cpp
// C++
constexpr int pwtspec = 6;
constexpr Real dkfac = 0.58;
constexpr std::array<Real, pwtspec> difrm = {
  1.0, 1.0, 0.9757, ...
};
```

**`constexpr`** means:
- Computed at compile time (like Fortran PARAMETER)
- Can be used in template parameters and array sizes
- Very efficient (no runtime cost)

### When to Use const vs constexpr

```cpp
constexpr int SIZE = 100;     // Compile-time constant
const Real PI = 3.14159;      // Runtime constant (but known immediately)
const Real x = some_function();  // Runtime constant (computed once)

// In functions:
void foo() {
  constexpr int N = 10;       // Compile-time constant
  const Real PI = 3.14159;    // Runtime constant
  Real x = get_value();       // Variable (can change)
}
```

---

## Quick Reference Table

| Feature | Fortran | C++ |
|---------|---------|-----|
| **Module/Namespace** | `module name` | `namespace name {` |
| **Use module** | `use name` | `using namespace name;` |
| **Comment** | `! comment` | `// comment` |
| **Variable decl** | `real(r8) :: x` | `Real x;` |
| **Array** | `dimension(N)` or `(N)` | `std::array<Type,N>` |
| **Array index** | `array(i)` [1-indexed] | `array[i]` [0-indexed] |
| **Logical** | `logical` | `bool` |
| **True/False** | `.true.` / `.false.` | `true` / `false` |
| **And/Or/Not** | `.and.` / `.or.` / `.not.` | `&&` / `\|\|` / `!` |
| **Parameter** | `parameter ::` | `constexpr` |
| **Function** | `function name(...) result(x)` | `Type name(...) { return x; }` |
| **Subroutine** | `subroutine name(...)` | `void name(...)` or return value |
| **Intent(in)** | `intent(in) :: x` | `const Type& x` |
| **Intent(out)** | `intent(out) :: x` | Return value or `Type& x` |
| **Intent(inout)** | `intent(inout) :: x` | `Type& x` |
| **Do loop** | `do i = 1, N` | `for (int i = 0; i < N; ++i)` |
| **If-then** | `if (...) then` | `if (...) {` |
| **End if** | `end if` | `}` |
| **Exp function** | `exp(x)` | `ekat::impl::exp(x)` (device) |
| **Print** | `print *, x` | `std::cout << x;` (host) |
| **Module variable** | `module :: var` | `namespace::var` |
| **Interface** | `interface name ... end interface` | `template<typename T> ...` |
| **Generic prog** | Module procedures + interface | Templates |
| **GPU compatible** | `!$acc routine seq` | `KOKKOS_INLINE_FUNCTION` |

---

## Common Pitfalls for Fortran Developers

### 1. Forgetting Zero-Indexing
```fortran
! Fortran: First element is 1
do i = 1, N
  array(i) = ...
```
```cpp
// C++: First element is 0
for (int i = 0; i < N; ++i) {  // Note: < not <=
  array[i] = ...;
}
```

### 2. Using () Instead of []
```cpp
// WRONG
array(i) = ...;  // Error! This tries to call array as a function!

// RIGHT
array[i] = ...;  // Access element i
```

### 3. Forgetting Semicolons
```cpp
// WRONG
int x = 5
int y = 10

// RIGHT
int x = 5;
int y = 10;
```

### 4. Using std:: Functions on GPU
```cpp
// WRONG - std:: functions don't work on GPU
KOKKOS_INLINE_FUNCTION
Real foo(Real x) {
  return std::exp(x);  // ERROR on GPU!
}

// RIGHT - Use ekat::impl:: or Kokkos:: versions
KOKKOS_INLINE_FUNCTION
Real foo(Real x) {
  return ekat::impl::exp(x);  // Works on GPU and CPU
}
```

### 5. Not Understanding References
```cpp
// Inefficient (copies large object)
Real foo(Pack<Real,16> x) {
  return x[0];
}

// Efficient (no copy)
Real foo(const Pack<Real,16>& x) {
  return x[0];
}
```

---

## Learning Path

1. ✅ **Start with syntax**: Get comfortable with `{}`, `;`, `[]`, `&&`, etc.
2. ✅ **Understand zero-indexing**: Always remember to subtract 1!
3. ✅ **Learn const references**: Use `const &` for read-only parameters
4. ✅ **Master templates basics**: Understand `template<typename T>`
5. ✅ **Get comfortable with KOKKOS_INLINE_FUNCTION**: Mark device functions
6. ✅ **Use ekat::impl:: math functions**: For device compatibility
7. ⏭️ **Advanced**: Packs, Kokkos Views, parallel patterns

---

## Additional Resources

- **Kokkos Documentation**: https://kokkos.github.io/kokkos-core-wiki/
- **EKat Documentation**: In EAMxx source tree at `externals/ekat/`
- **C++ Reference**: https://en.cppreference.com/
- **Ask the Team**: EAMxx developers are friendly and helpful!

---

## Summary

**Key Takeaways:**
1. **C++ is stricter than Fortran** - explicit types, semicolons, etc.
2. **Zero-indexing is everywhere** - arrays, loops, everything starts at 0
3. **Templates are powerful** - one function works for many types
4. **Kokkos enables portability** - same code runs on CPU and GPU
5. **Use const& liberally** - efficient and safe for input parameters
6. **Pack types enable vectorization** - explicit SIMD programming

**Most Important Rule:**
> When in doubt, look at existing EAMxx code and follow the patterns!

The water_isotopes.hpp and water_types.hpp files have extensive comments showing these patterns in action.
