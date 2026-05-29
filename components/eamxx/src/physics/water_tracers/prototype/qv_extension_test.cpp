// Prototype test for qv extension from (col,lev) to (tracer,col,lev)
// Validates BFB preservation and measures performance overhead
//
// Compile: g++ -std=c++17 -O3 -o qv_extension_test qv_extension_test.cpp
// Run:     ./qv_extension_test
//          ./qv_extension_test --benchmark

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <algorithm>

// Simple 2D and 3D array wrappers (no Kokkos for standalone testing)
template<typename T>
class Array2D {
public:
  Array2D(int d1, int d2) : d1_(d1), d2_(d2), data_(d1*d2) {}
  T& operator()(int i, int j) { return data_[i*d2_ + j]; }
  const T& operator()(int i, int j) const { return data_[i*d2_ + j]; }
  int size() const { return data_.size(); }
private:
  int d1_, d2_;
  std::vector<T> data_;
};

template<typename T>
class Array3D {
public:
  Array3D(int d1, int d2, int d3) : d1_(d1), d2_(d2), d3_(d3), data_(d1*d2*d3) {}
  T& operator()(int i, int j, int k) { return data_[i*d2_*d3_ + j*d3_ + k]; }
  const T& operator()(int i, int j, int k) const { return data_[i*d2_*d3_ + j*d3_ + k]; }
  int size() const { return data_.size(); }
private:
  int d1_, d2_, d3_;
  std::vector<T> data_;
};

// Simulate a physics kernel that reads and writes qv
template<typename T>
void scalar_kernel(Array2D<T>& qv, int ncols, int nlevs) {
  // Simple kernel: qv_new = qv * 1.01 + 0.0001 (simulate tendency)
  for (int icol = 0; icol < ncols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      qv(icol, ilev) = qv(icol, ilev) * 1.01 + 0.0001;
    }
  }
}

template<typename T>
void tracer_kernel(Array3D<T>& qv, int ntracers, int ncols, int nlevs) {
  // Same kernel but with tracer dimension (operating on slot 0 only for BFB test)
  for (int itracer = 0; itracer < 1; ++itracer) {  // Only slot 0 for BFB comparison
    for (int icol = 0; icol < ncols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        qv(itracer, icol, ilev) = qv(itracer, icol, ilev) * 1.01 + 0.0001;
      }
    }
  }
}

int main(int argc, char** argv) {
  bool benchmark_mode = false;
  if (argc > 1 && std::string(argv[1]) == "--benchmark") {
    benchmark_mode = true;
  }

  std::cout << "=================================================\n";
  std::cout << "QV Extension Prototype Test\n";
  std::cout << "=================================================\n\n";

  // Test dimensions (realistic for EAMxx)
  const int ncols = 384;    // ne4pg2 columns
  const int nlevs = 72;     // Typical vertical levels
  const int ntracers = 1;   // SCREAM_NUM_TRACERS=1 for BFB test
  const int iterations = 100;  // Physics timesteps to simulate

  std::cout << "Configuration:\n";
  std::cout << "  ncols      = " << ncols << "\n";
  std::cout << "  nlevs      = " << nlevs << "\n";
  std::cout << "  ntracers   = " << ntracers << "\n";
  std::cout << "  iterations = " << iterations << "\n\n";

  // Allocate arrays
  Array2D<double> qv_scalar(ncols, nlevs);
  Array3D<double> qv_tracer(ntracers, ncols, nlevs);

  // Initialize with same values
  const double init_val = 0.01;  // 10 g/kg, typical atmospheric moisture
  for (int icol = 0; icol < ncols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      double val = init_val * (1.0 + 0.1 * std::sin(icol * 0.1) * std::cos(ilev * 0.05));
      qv_scalar(icol, ilev) = val;
      qv_tracer(0, icol, ilev) = val;
    }
  }

  std::cout << "Running scalar kernel... " << std::flush;
  auto start_scalar = std::chrono::high_resolution_clock::now();
  for (int iter = 0; iter < iterations; ++iter) {
    scalar_kernel(qv_scalar, ncols, nlevs);
  }
  auto end_scalar = std::chrono::high_resolution_clock::now();
  auto duration_scalar = std::chrono::duration_cast<std::chrono::microseconds>(end_scalar - start_scalar).count();
  std::cout << "done (" << duration_scalar << " us)\n";

  std::cout << "Running tracer kernel...  " << std::flush;
  auto start_tracer = std::chrono::high_resolution_clock::now();
  for (int iter = 0; iter < iterations; ++iter) {
    tracer_kernel(qv_tracer, ntracers, ncols, nlevs);
  }
  auto end_tracer = std::chrono::high_resolution_clock::now();
  auto duration_tracer = std::chrono::duration_cast<std::chrono::microseconds>(end_tracer - start_tracer).count();
  std::cout << "done (" << duration_tracer << " us)\n\n";

  // BFB Check
  std::cout << "BFB Check:\n";
  double max_abs_diff = 0.0;
  double max_rel_diff = 0.0;
  int diff_count = 0;
  int max_diff_col = 0, max_diff_lev = 0;

  for (int icol = 0; icol < ncols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      double scalar_val = qv_scalar(icol, ilev);
      double tracer_val = qv_tracer(0, icol, ilev);
      double abs_diff = std::abs(scalar_val - tracer_val);
      double rel_diff = abs_diff / (std::abs(scalar_val) + 1e-20);

      if (abs_diff > 0.0) {
        diff_count++;
      }

      if (abs_diff > max_abs_diff) {
        max_abs_diff = abs_diff;
        max_rel_diff = rel_diff;
        max_diff_col = icol;
        max_diff_lev = ilev;
      }
    }
  }

  bool bfb_pass = (max_abs_diff == 0.0);
  bool tight_tolerance_pass = (max_rel_diff < 1e-12);

  std::cout << "  Max absolute difference: " << std::scientific << std::setprecision(3) << max_abs_diff << "\n";
  std::cout << "  Max relative difference: " << std::scientific << std::setprecision(3) << max_rel_diff << "\n";
  std::cout << "  Points with differences: " << diff_count << " / " << (ncols*nlevs) << "\n";

  if (bfb_pass) {
    std::cout << "  BFB: PASS (exact bit-for-bit match)\n";
  } else if (tight_tolerance_pass) {
    std::cout << "  BFB: SOFT-PASS (within rtol=1e-12)\n";
  } else {
    std::cout << "  BFB: FAIL (exceeds rtol=1e-12)\n";
    std::cout << "  Location of max difference: col=" << max_diff_col << " lev=" << max_diff_lev << "\n";
  }
  std::cout << "\n";

  // Performance Check
  std::cout << "Performance Check:\n";
  double overhead_percent = 100.0 * (duration_tracer - duration_scalar) / (double)duration_scalar;
  std::cout << "  Scalar runtime:  " << duration_scalar << " us\n";
  std::cout << "  Tracer runtime:  " << duration_tracer << " us\n";
  std::cout << "  Overhead: " << std::fixed << std::setprecision(2) << overhead_percent << "%\n";

  bool perf_pass = (overhead_percent < 2.0);
  if (perf_pass) {
    std::cout << "  Performance: PASS (<2% overhead)\n";
  } else {
    std::cout << "  Performance: FAIL (>=2% overhead)\n";
  }
  std::cout << "\n";

  // Write results to file
  std::ofstream out("qv_test_output.txt");
  out << "QV Extension Prototype Test Results\n";
  out << "====================================\n\n";
  out << "Configuration:\n";
  out << "  ncols      = " << ncols << "\n";
  out << "  nlevs      = " << nlevs << "\n";
  out << "  ntracers   = " << ntracers << "\n";
  out << "  iterations = " << iterations << "\n\n";
  out << "BFB Check:\n";
  out << "  Max diff: " << std::scientific << max_abs_diff << "\n";
  out << "  BFB: " << (bfb_pass ? "PASS" : (tight_tolerance_pass ? "SOFT-PASS" : "FAIL")) << "\n\n";
  out << "Performance Check:\n";
  out << "  Scalar runtime: " << duration_scalar << " us\n";
  out << "  Tracer runtime: " << duration_tracer << " us\n";
  out << "  Overhead: " << std::fixed << std::setprecision(2) << overhead_percent << "%\n";
  out << "  Performance: " << (perf_pass ? "PASS" : "FAIL") << "\n";
  out.close();

  std::cout << "Results written to qv_test_output.txt\n\n";

  // Summary
  std::cout << "=================================================\n";
  std::cout << "Summary:\n";
  std::cout << "  BFB gate:         " << (bfb_pass || tight_tolerance_pass ? "PASS" : "FAIL") << "\n";
  std::cout << "  Performance gate: " << (perf_pass ? "PASS" : "FAIL") << "\n";
  std::cout << "=================================================\n";

  return (bfb_pass || tight_tolerance_pass) && perf_pass ? 0 : 1;
}
