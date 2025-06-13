#ifndef WATER_TYPES
#define WATER_TYPES

#include "share/eamxx_types.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/logging/ekat_logger.hpp"

#include <string>
#include <vector>
#include <array>

namespace scream {
namespace WaterTypes {

/* 
* Isotopic constants used to determine fractionation factors, etc.
*/

// Define all water types that can be specified
struct WaterTypes {
  // this defines the various types water types that will be queried by other functions
  static constexpr int pwtype = 7;      // number of water types

  enum class WaterType : int {
    Undefined = 0,
    Vapor = 1,
    CloudLiquid = 2,
    CloudIce = 3,
    StratiformRain = 4, 
    StratiformSnow = 5,
    ConvectiveRain = 6,
    ConvectiveSnow = 7
  };

  // Do these exist elsewhere? RPF
  static const std::vector<std::string> wtype_names; // names of water types in model
  static const std::vector<std::string> wtype_suffix; // suffixes of water types

};

static const std::vector<std::string> wtype_names = {"VAPOR", "LIQUID", "ICE", "STRAT_RAIN", "STRAT_SNOW", "CONV_RAIN", "CONV_SNOW"};
static const std::vector<std::string> wtype_suffix = {"_v", "_l", "_i", "_R", "_S", "_r", "_s"};

}
}

#endif // WATER_TYPES