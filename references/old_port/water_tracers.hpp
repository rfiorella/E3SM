#ifndef WATER_TRACERS
#define WATER_TRACERS

#include <string>
#include <array>
#include <Kokkos_Core.hpp>

namespace scream {
namespace WaterTracers {

// Constants (should be set from shared config headers)
constexpr int pcnst   = 7; /* value set elsewhere */;
constexpr int pwtspec = 7; /* value set elsewhere */;
constexpr int pwtype  = 7; /* value set elsewhere */;

// WTRC_MAX_CNST is provided by CMake as a preprocessor define
// When SCREAM_TRACE_WATER=OFF: WTRC_MAX_CNST=1 (enforced at configure time)
// When SCREAM_TRACE_WATER=ON:  WTRC_MAX_CNST=SCREAM_NUM_WATER_TRACERS
#ifndef WTRC_MAX_CNST
#error "WTRC_MAX_CNST must be defined by CMake"
#endif

constexpr int WTRC_WSET_STD = 1;

// Namelist Variables (runtime-configurable)
inline bool trace_water              = false;
inline bool wisotope                 = false;
inline bool wtrc_lh2oadj             = false;
inline bool wtrc_lzmlin              = true;
inline bool wtrc_warn_only           = true;
inline bool wtrc_add_cvprecip        = false;
inline bool wtrc_add_stprecip        = false;
inline bool wtrc_alpha_kinetic       = false;
inline bool wtrc_check_total_h2o     = false;
inline bool wtrc_check_show_types    = false;
inline bool wtrc_detrain_in_macrop   = false;
inline bool wtrc_use_ice_supsat      = false;

// Numerical Configuration
inline int wtrc_niter = 10;
inline int wtrc_citer = 20;

inline double wtrc_qchkmin = 1.e-14;
inline double wtrc_qmin    = 1.e-18;

// These could be `Kokkos::View` if they need device access
inline std::array<double, pwtspec> wtrc_fixed_alpha = {1.0};
inline std::array<double, pwtspec> wtrc_fixed_rstd  = {1.0};

inline std::string water_tracer_model = "none";

// DEPRECATED: The following arrays are legacy scaffolding that is replaced
// by the compile-time water tracer registry (water_tracer_registry.hpp).
// Use the registry query API instead:
//   - tracer_name(i)          replaces wtrc_names[i]
//   - tracer_isotopologue(i)  replaces wtrc_species[i]
//   - tracer_is_tag(i)        replaces wtrc_is_tag[i]
// These arrays are retained temporarily for backwards compatibility but are
// not populated and should not be used in new code.

// Tracer name arrays (fixed-length strings like Fortran)
using NameArray = std::array<std::string, WTRC_MAX_CNST>;

inline NameArray wtrc_names           = {};
inline NameArray wtrc_species_names  = {};
inline NameArray wtrc_type_names     = {};
inline NameArray wtrc_srfvap_names   = {};
inline NameArray wtrc_srfpcp_names   = {};
inline NameArray wtrc_tag_names      = {};
inline NameArray wtrc_out_names      = {};

// Bulk water names and indices
inline constexpr std::array<const char*, pwtype> wtrc_bulk_names = {
  "Q", "CLDLIQ", "CLDICE", "RAINQM", "SNOWQM", "RAINQC", "SNOWQC"
};

inline std::array<int, pwtype> wtrc_bulk_indices = {};

// Derived water tracer configuration
inline int wtrc_ncnst = 0;
inline std::array<int, WTRC_MAX_CNST> wtrc_indices = {};
inline std::array<int, WTRC_MAX_CNST> wtrc_species = {};
inline std::array<int, WTRC_MAX_CNST> wtrc_types   = {};
inline std::array<bool, WTRC_MAX_CNST> wtrc_is_tag = {};

// Index arrays (grouped by traits)
inline std::array<std::array<int, WTRC_MAX_CNST / pwtype>, pwtype> wtrc_iatype = {};
inline std::array<int, pwtspec> wtrc_nspec = {};
inline std::array<std::array<int, WTRC_MAX_CNST>, pwtspec> wtrc_iaspec = {};
inline int wtrc_nwset = 0;
inline std::array<std::array<int, WTRC_MAX_CNST / pwtype>, pwtype> wtrc_iawset = {};
inline std::array<std::array<int, WTRC_MAX_CNST / pwtype>, pwtype> wtrc_srfpcp_indices = {};

// Surface tracer indices
inline int wtrc_nsrfvap = 0;
inline std::array<int, WTRC_MAX_CNST> wtrc_iasrfvap = {};
inline int wtrc_nsrfpcp = 0;
inline std::array<int, WTRC_MAX_CNST> wtrc_iasrfpcp = {};

// Constituent flags (based on pcnst)
inline std::array<int, pcnst> iwater = {};
inline std::array<int, pcnst> iwspec = {};
inline std::array<bool, pcnst> iwistag = {};

// Subview accessor for bulk water tracer (CMP index 0)
// Returns a (COL, LEV) view from a rank-3 (COL, CMP, LEV) water field.
// Compiled in all builds - with WTRC_MAX_CNST=1 (default), this extracts
// the single CMP slice; with WTRC_MAX_CNST>1, it extracts the bulk tracer.
//
// Usage: auto qv_bulk = get_bulk_water_subview(qv_rank3);
//
// Template parameter ViewType should be a Kokkos::View<Real***, ...>
template<typename ViewType>
KOKKOS_INLINE_FUNCTION
auto get_bulk_water_subview(const ViewType& water_field_rank3) {
  return Kokkos::subview(water_field_rank3, Kokkos::ALL(), 0, Kokkos::ALL());
}

} // namespace WaterTracers
} // namespace scream

#endif // WATER_TRACERS