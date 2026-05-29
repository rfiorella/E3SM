#ifndef WATER_TRACER_METADATA_HPP
#define WATER_TRACER_METADATA_HPP

#include <string>

namespace scream {
namespace water_tracers {

enum class WaterTracerKind {
  bulk,                  // Slot 0, canonical total water
  evaporation,          // Tracer with modified surface evaporation flux (Group 1 testing)
  // Future campaigns will add:
  // isotope, tagged_partition, tagged_diagnostic, auxiliary
};

struct WaterTracerMetadata {
  std::string name;
  std::string longname;
  WaterTracerKind kind;
  // Units implicitly kg/kg for all water tracers in Group 1

  // Fields below reserved for future campaigns
  bool conserved_independently = true;
  std::string ratio_standard;  // For isotope tracers
  std::string partition_group_id;  // For tagged_partition tracers

  WaterTracerMetadata() : kind(WaterTracerKind::bulk), conserved_independently(true) {}

  WaterTracerMetadata(const std::string& n, const std::string& ln, WaterTracerKind k)
    : name(n), longname(ln), kind(k), conserved_independently(true) {}
};

} // namespace water_tracers
} // namespace scream

#endif // WATER_TRACER_METADATA_HPP
