# Water tracer infrastructure tests

include(EkatCreateUnitTest)

# Test 1: Tracer field layout infrastructure
EkatCreateUnitTest(test_tracer_field_infrastructure
  SOURCES test_tracer_field_infrastructure.cpp
  LIBS scream_share
  THREADS 1
)

# Test 2: Tracer access patterns (BFB validation)
EkatCreateUnitTest(test_tracer_access_patterns
  SOURCES test_tracer_access_patterns.cpp
  LIBS scream_share
  THREADS 1
)
