#ifndef C6C12_H
#define C6C12_H
#include <string>
#include <utility>
#include <vector>

using namespace std;
struct para_base {
  double c6_sigma, c12_epsilon, mass, radius, charge;
  string atom_type, code_name;
  para_base()
      : c6_sigma{}, c12_epsilon{}, mass{1.0}, radius{}, charge{}, code_name{} {}
  template <typename T1, typename T11, typename T2, typename T3>
  para_base(T1 &&_c6_sigma, T1 &&_c12_epsilon, T11 &&_mass, T3 &&_radius,
            T2 &&_atom_type)
      : c6_sigma{forward<T1>(_c6_sigma)},
        c12_epsilon{forward<T1>(_c12_epsilon)}, mass{forward<T11>(_mass)},
        radius{forward<T1>(_radius)}, atom_type{forward<T2>(_atom_type)} {}

  template <typename T1, typename T11, typename T2, typename T3>
  para_base(T1 &&_c6_sigma, T1 &&_c12_epsilon, T11 &&_mass, T3 &&_radius,
            T1 &&_charge, T2 &&_atom_type, T2 &&_code_name)
      : c6_sigma{forward<T1>(_c6_sigma)},
        c12_epsilon{forward<T1>(_c12_epsilon)}, mass{forward<T11>(_mass)},
        radius{forward<T1>(_radius)}, charge{_charge},
        atom_type{forward<T2>(_atom_type)}, code_name{_code_name} {}
};
class para {
public:
  vector<para_base> base;
  para() { base.reserve(200); }
};
#endif