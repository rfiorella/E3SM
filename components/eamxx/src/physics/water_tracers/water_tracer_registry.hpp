#ifndef WATER_TRACER_REGISTRY_HPP
#define WATER_TRACER_REGISTRY_HPP

#include "water_tracer_metadata.hpp"
#include <vector>
#include <map>
#include <stdexcept>

namespace scream {
namespace water_tracers {

class WaterTracerRegistry {
public:
  // Singleton access
  static WaterTracerRegistry& instance() {
    static WaterTracerRegistry registry;
    return registry;
  }

  // Register a new tracer
  void register_tracer(const WaterTracerMetadata& meta) {
    if (tracers_by_name_.find(meta.name) != tracers_by_name_.end()) {
      throw std::runtime_error("Tracer already registered: " + meta.name);
    }
    tracers_.push_back(meta);
    tracers_by_name_[meta.name] = tracers_.size() - 1;
  }

  // Get tracer metadata by name
  const WaterTracerMetadata& get(const std::string& name) const {
    auto it = tracers_by_name_.find(name);
    if (it == tracers_by_name_.end()) {
      throw std::runtime_error("Tracer not found: " + name);
    }
    return tracers_[it->second];
  }

  // Get tracer metadata by index
  const WaterTracerMetadata& get(int index) const {
    if (index < 0 || index >= static_cast<int>(tracers_.size())) {
      throw std::runtime_error("Tracer index out of range: " + std::to_string(index));
    }
    return tracers_[index];
  }

  // Get total number of registered tracers
  int count() const {
    return static_cast<int>(tracers_.size());
  }

  // Get all tracer names
  std::vector<std::string> get_names() const {
    std::vector<std::string> names;
    names.reserve(tracers_.size());
    for (const auto& tracer : tracers_) {
      names.push_back(tracer.name);
    }
    return names;
  }

private:
  WaterTracerRegistry() {
    // Register bulk water (slot 0) by default
    WaterTracerMetadata bulk_meta("bulk_H2O", "Bulk Water (Total)", WaterTracerKind::bulk);
    tracers_.push_back(bulk_meta);
    tracers_by_name_["bulk_H2O"] = 0;
  }

  ~WaterTracerRegistry() = default;

  // Prevent copying
  WaterTracerRegistry(const WaterTracerRegistry&) = delete;
  WaterTracerRegistry& operator=(const WaterTracerRegistry&) = delete;

  std::vector<WaterTracerMetadata> tracers_;
  std::map<std::string, int> tracers_by_name_;
};

} // namespace water_tracers
} // namespace scream

#endif // WATER_TRACER_REGISTRY_HPP
