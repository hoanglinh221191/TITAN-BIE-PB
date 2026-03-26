#ifndef MOL_PARA_H
#define MOL_PARA_H
#include <cmath>
#include <utility>

using namespace std;
class Multi_MOL_parameters {
  double charge{}, radius{}, sigma{}, epsilon{}, sqrt_epsilon{}, mass{1.};

public:
  Multi_MOL_parameters()
      : charge{}, radius{}, sigma{}, epsilon{}, sqrt_epsilon{}, mass{1.} {};
  Multi_MOL_parameters(const Multi_MOL_parameters &rhs)
      : charge{rhs.charge}, radius{rhs.radius}, sigma{rhs.sigma},
        epsilon{rhs.epsilon}, sqrt_epsilon{rhs.sqrt_epsilon}, mass{rhs.mass} {}
  Multi_MOL_parameters(Multi_MOL_parameters &&rhs) noexcept
      : charge{move(rhs.charge)}, radius{move(rhs.radius)},
        sigma{move(rhs.sigma)}, epsilon{move(rhs.epsilon)},
        sqrt_epsilon{move(rhs.sqrt_epsilon)}, mass{move(rhs.mass)} {}
  Multi_MOL_parameters &operator=(const Multi_MOL_parameters &rhs) noexcept {
    charge = rhs.charge;
    radius = rhs.radius;
    sigma = rhs.sigma;
    epsilon = rhs.epsilon;
    sqrt_epsilon = rhs.sqrt_epsilon;
    mass = rhs.mass;
    return *this;
  }
  Multi_MOL_parameters &operator=(Multi_MOL_parameters &&rhs) noexcept {
    charge = move(rhs.charge);
    radius = move(rhs.radius);
    sigma = move(rhs.sigma);
    epsilon = move(rhs.epsilon);
    sqrt_epsilon = move(rhs.sqrt_epsilon);
    mass = move(rhs.mass);
    return *this;
  }
  ~Multi_MOL_parameters() = default;
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  Multi_MOL_parameters(T1 &&_charge, T2 &&_radius, T3 &&_sigma, T4 &&_epsilon,
                       T5 &&_mass)
      : charge{forward<T1>(_charge)}, radius{forward<T2>(_radius)},
        sigma{forward<T3>(_sigma)}, epsilon{forward<T4>(_epsilon)},
        mass{forward<T5>(_mass)} {
    sqrt_epsilon = std::sqrt(epsilon);
  }
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void set_data(T1 &&_charge, T2 &&_radius, T3 &&_sigma, T4 &&_epsilon,
                T5 &&_mass) {
    charge = forward<T1>(_charge);
    radius = forward<T2>(_radius);
    sigma = forward<T3>(_sigma);
    epsilon = forward<T4>(_epsilon);
    mass = forward<T5>(_mass);
    sqrt_epsilon = std::sqrt(epsilon);
  }
  template <typename T> void set_charge(T &&_charge) {
    charge = forward<T>(_charge);
  }
  template <typename T> void set_mass(T &&_mass) { mass = forward<T>(_mass); }
  template <typename T, typename T2, typename T3>
  void set_epsilon_sigma_radius(T &&_epsilon, T2 &&_sigma, T3 &&_radius) {
    epsilon = _epsilon;
    sigma = _sigma;
    radius = _radius;
    sqrt_epsilon = std::sqrt(epsilon);
  }
  const double get_mass() const {
    if (mass > 0.01)
      return mass;
    else
      return 0.01;
  }
  const double get_charge() const { return charge; }
  const double get_radius() const { return radius; }
  const double get_epsilon() const { return epsilon; }
  const double get_sqrt_epsilon() const { return sqrt_epsilon; }
  const double get_sigma() const { return sigma; }
};
#endif