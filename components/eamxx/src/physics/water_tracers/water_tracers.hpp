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
//constexpr int SET_WTRC_MAX_CNST = /* define as in Fortran environment */;
constexpr int WTRC_MAX_CNST = 1;
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

} // namespace water_tracers
} // namespace scream

#endif // WATER_TRACERS