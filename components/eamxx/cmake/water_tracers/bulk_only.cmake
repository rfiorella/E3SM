# Water tracer configuration: bulk H2O only
# Usage: -DSCREAM_WATER_TRACERS_FILE=<path-to-eamxx>/cmake/water_tracers/bulk_only.cmake

add_water_tracer(
  NAME bulk
  ISOTOPOLOGUE H2O
)
