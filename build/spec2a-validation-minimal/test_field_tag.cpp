#include "../../src/share/field/field_tag.hpp"
#include <iostream>

int main() {
    using namespace scream;

    // Test that TRACER FieldTag exists and compiles
    FieldTag tag = FieldTag::Tracer;

    std::cout << "PASS: TRACER FieldTag compiles" << std::endl;
    std::cout << "FieldTag value: " << static_cast<int>(tag) << std::endl;

    return 0;
}
