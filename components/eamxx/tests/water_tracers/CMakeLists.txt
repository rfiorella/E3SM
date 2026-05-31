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

# Test 3: SHOC qv field access with tracer dimension
EkatCreateUnitTest(test_qv_shoc
  SOURCES test_qv_shoc.cpp
  LIBS scream_share
  THREADS 1
)

# Test 4: P3 qv field access with tracer dimension
EkatCreateUnitTest(test_qv_p3
  SOURCES test_qv_p3.cpp
  LIBS scream_share
  THREADS 1
)

# Test 5: P3 precipitation field access with tracer dimension
EkatCreateUnitTest(test_precip_tracer_access
  SOURCES test_precip_tracer_access.cpp
  LIBS scream_share
  THREADS 1
)
