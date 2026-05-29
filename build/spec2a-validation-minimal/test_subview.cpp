#include <iostream>

// Simulate subview accessor pattern
template<typename T>
struct Field3D {
    T data[1][10][20];

    T& operator()(int t, int c, int l) { return data[t][c][l]; }
    const T& operator()(int t, int c, int l) const { return data[t][c][l]; }
};

template<typename T>
struct SubView2D {
    T (*data)[20];  // Pointer to [col][lev] slice

    SubView2D(T (*ptr)[20]) : data(ptr) {}

    T& operator()(int c, int l) { return data[c][l]; }
    const T& operator()(int c, int l) const { return data[c][l]; }
};

template<typename T>
SubView2D<T> get_tracer_bulk_subview(Field3D<T>& field) {
    return SubView2D<T>(field.data[0]);  // Return slice at tracer=0
}

int main() {
    Field3D<double> qv;

    // Get subview of slot-0
    auto qv_bulk = get_tracer_bulk_subview(qv);

    // Initialize via subview (2D indexing)
    for (int c = 0; c < 10; ++c) {
        for (int l = 0; l < 20; ++l) {
            qv_bulk(c, l) = c * 100.0 + l;
        }
    }

    // Verify via original 3D field
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
        std::cout << "PASS: Subview accessor pattern works correctly" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL: Subview pattern produced wrong values" << std::endl;
        return 1;
    }
}
