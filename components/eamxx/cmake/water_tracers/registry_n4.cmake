# Water tracer configuration: 4-tracer test case
# Usage: -DSCREAM_WATER_TRACERS_FILE=<path-to-eamxx>/cmake/water_tracers/registry_n4.cmake
#
# Tracers:
#  0: bulk - bulk H2O (physically fractionates, represents total water)
#  1: passive - passive H2O tag (no fractionation, source attribution)
#  2: hdo - HDO isotopologue
#  3: h218o - H218O isotopologue

add_water_tracer(
  NAME bulk
  ISOTOPOLOGUE H2O
)

add_water_tracer(
  NAME passive
  ISOTOPOLOGUE H2O
  TAG
)

add_water_tracer(
  NAME hdo
  ISOTOPOLOGUE HDO
)

add_water_tracer(
  NAME h218o
  ISOTOPOLOGUE H218O
)
