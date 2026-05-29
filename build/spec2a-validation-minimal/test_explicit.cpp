#include <iostream>

// Simulate 3D tracer field access with explicit indexing
template<typename T>
struct Field3D {
    T data[1][10][20];  // [tracer][col][lev]

    T& operator()(int t, int c, int l) { return data[t][c][l]; }
    const T& operator()(int t, int c, int l) const { return data[t][c][l]; }
};

int main() {
    Field3D<double> qv;

    // Initialize with test pattern
    for (int c = 0; c < 10; ++c) {
        for (int l = 0; l < 20; ++l) {
            qv(0, c, l) = c * 100.0 + l;  // Explicit slot-0 indexing
        }
    }

    // Verify
    bool pass = true;
    for (int c = 0; c < 10; ++c) {
        for (int l = 0; l < 20; ++l) {
            double expected = c * 100.0 + l;
            if (qv(0, c, l) != expected) {
                pass = false;
                break;
            }
        }
    }

    if (pass) {
        std::cout << "PASS: Explicit indexing qv(0, icol, ilev) works correctly" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL: Explicit indexing produced wrong values" << std::endl;
        return 1;
    }
}
