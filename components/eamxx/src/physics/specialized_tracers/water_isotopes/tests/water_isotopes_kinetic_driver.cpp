// Prints the kinetic fractionation factor alpkn of alpha_kinetic for every
// case of a CSV file, for comparison against the Fortran wiso_kmol.
//
//   water_isotopes_kinetic_driver <cases.csv>
//
// Input rows are "species,rbot,zbot,ustar" (species HDO or H218O; a header
// line and lines starting with '#' are skipped). Output rows are
// "species,rbot,zbot,ustar,alpkn" with %.17e, using the default
// WaterIsotopeParameters (Merlivat 1978 diffusivities). Not a ctest test.

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_kinetic.hpp"
#include "share/core/eamxx_session.hpp"

#include <ekat_pack.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main (int argc, char** argv)
{
  using namespace scream;
  using water_isotopes::WaterIsotopologues;
  using Pack1 = ekat::Pack<Real, 1>;

  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <cases.csv>\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "error: cannot open " << argv[1] << "\n";
    return 2;
  }

  initialize_eamxx_session(false);
  int status = 0, nrows = 0;
  {
    const water_isotopes::WaterIsotopeParameters<Real> params;
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line.rfind("species", 0) == 0) continue;
      std::stringstream ss(line);
      std::string sp_name, field;
      double v[3];
      std::getline(ss, sp_name, ',');
      bool ok = true;
      for (double& x : v) {
        ok = ok && static_cast<bool>(std::getline(ss, field, ','));
        if (ok) x = std::stod(field);
      }
      WaterIsotopologues sp;
      if (sp_name == "HDO") {
        sp = WaterIsotopologues::HDO;
      } else if (sp_name == "H218O") {
        sp = WaterIsotopologues::H218O;
      } else {
        ok = false;
      }
      if (!ok) {
        std::cerr << "error: bad row: " << line << "\n";
        status = 1;
        break;
      }
      // Host evaluation; alpha_kinetic is the same inline function on device.
      const Real a = water_isotopes::WaterIsotopeKinetic::alpha_kinetic(
          Pack1(v[0]), Pack1(v[1]), Pack1(v[2]), sp, params)[0];
      std::printf("%s,%.17e,%.17e,%.17e,%.17e\n", sp_name.c_str(), v[0], v[1], v[2], double(a));
      ++nrows;
    }
  }
  finalize_eamxx_session();

  if (status == 0 && nrows == 0) {
    std::cerr << "error: no cases in " << argv[1] << "\n";
    status = 1;
  }
  return status;
}
