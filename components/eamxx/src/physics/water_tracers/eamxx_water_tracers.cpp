#include "water_tracers.hpp"
#include "water_isotopes.hpp"

#include <vector>
#include <string>


#include "share/eamxx_types.hpp"
#include "physics_constants.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/logging/ekat_logger.hpp"

namespace scream { 
namespace WaterTracers { 
// =========================================================================================
// WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
//   : AtmosphereProcess(comm, params)
// {
//   // Nothing to do here
// }

// void WtypeGetItype(name) {
//   // Retrieve WaterType based on type name
// };


void WtypeGetAlpha() {
  // retrieve fractionation factor for process that converts
  // source water type (isrctype) to desination water type (idsttype)

};

// // Set of functions to determine the type of water tracer
// bool trc_is_wtrc(tracer_num) {
//   // returns true if tracer is a water tracer
// };




void OceanTracerFlux() {
  // compute water tracer exchange from ocean
  // see L. 531-635 of wiso_flxoce in CAM
};

}
}
