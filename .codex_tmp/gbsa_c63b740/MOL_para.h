#ifndef MOL_PARA_H
#define MOL_PARA_H
#include <cmath>
#include <cstdint> // Added for gcc 13 strictness
#include <cstdlib> // For aligned_alloc on Linux
#ifdef _WIN32
#include <malloc.h>
#endif
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

// Custom C++20 standard-compliant aligned allocator for AVX2/AVX-512
template <typename T, std::size_t Alignment = 32> struct AlignedAllocator {
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using propagate_on_container_move_assignment = std::true_type;

  template <class U> struct rebind {
    using other = AlignedAllocator<U, Alignment>;
  };

  AlignedAllocator() noexcept = default;
  template <typename U>
  constexpr AlignedAllocator(const AlignedAllocator<U, Alignment> &) noexcept {}

  [[nodiscard]] T *allocate(std::size_t n) {
    if (n == 0)
      return nullptr;
    if (n > std::size_t(-1) / sizeof(T))
      throw std::bad_alloc();
    void *ptr = nullptr;
#ifdef _WIN32
    ptr = _aligned_malloc(n * sizeof(T), Alignment);
#else
    std::size_t size = (n * sizeof(T) + Alignment - 1) & ~(Alignment - 1);
    ptr = aligned_alloc(Alignment, size);
#endif
    if (!ptr)
      throw std::bad_alloc();
    return static_cast<T *>(ptr);
  }

  void deallocate(T *p, std::size_t n) noexcept {
#ifdef _WIN32
    _aligned_free(p);
#else
    free(p);
#endif
  }

  template <typename U>
  bool operator==(const AlignedAllocator<U, Alignment> &) const noexcept {
    return true;
  }
  template <typename U>
  bool operator!=(const AlignedAllocator<U, Alignment> &) const noexcept {
    return false;
  }
};

template <typename T, typename U, std::size_t A>
bool operator==(const AlignedAllocator<T, A> &,
                const AlignedAllocator<U, A> &) {
  return true;
}
template <typename T, typename U, std::size_t A>
bool operator!=(const AlignedAllocator<T, A> &,
                const AlignedAllocator<U, A> &) {
  return false;
}

using AlignedDoubleVector = std::vector<double, AlignedAllocator<double, 32>>;

class Multi_MOL_parameters {
public:
  AlignedDoubleVector charge;
  AlignedDoubleVector radius;
  AlignedDoubleVector sigma;
  AlignedDoubleVector epsilon;
  AlignedDoubleVector sqrt_epsilon;
  AlignedDoubleVector mass;

  Multi_MOL_parameters() = default;
  explicit Multi_MOL_parameters(size_t size) { resize(size); }

  struct Proxy {
    Multi_MOL_parameters *ptr;
    size_t idx;

    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    void set_data(T1 &&_charge, T2 &&_radius, T3 &&_sigma, T4 &&_epsilon,
                  T5 &&_mass) {
      ptr->charge[idx] = std::forward<T1>(_charge);
      ptr->radius[idx] = std::forward<T2>(_radius);
      ptr->sigma[idx] = std::forward<T3>(_sigma);
      ptr->epsilon[idx] = std::forward<T4>(_epsilon);
      ptr->mass[idx] = std::forward<T5>(_mass);
      ptr->sqrt_epsilon[idx] = std::sqrt(ptr->epsilon[idx]);
    }
    void set_charge(double c) { ptr->charge[idx] = c; }
    void set_mass(double m) { ptr->mass[idx] = m; }
    void set_epsilon_sigma_radius(double e, double s, double r) {
      ptr->epsilon[idx] = e;
      ptr->sigma[idx] = s;
      ptr->radius[idx] = r;
      ptr->sqrt_epsilon[idx] = std::sqrt(e);
    }

    double get_mass() const {
      return (ptr->mass[idx] > 0.01) ? ptr->mass[idx] : 0.01;
    }
    double get_charge() const { return ptr->charge[idx]; }
    double get_radius() const { return ptr->radius[idx]; }
    double get_epsilon() const { return ptr->epsilon[idx]; }
    double get_sqrt_epsilon() const { return ptr->sqrt_epsilon[idx]; }
    double get_sigma() const { return ptr->sigma[idx]; }
  };

  struct ConstProxy {
    const Multi_MOL_parameters *ptr;
    size_t idx;

    double get_mass() const {
      return (ptr->mass[idx] > 0.01) ? ptr->mass[idx] : 0.01;
    }
    double get_charge() const { return ptr->charge[idx]; }
    double get_radius() const { return ptr->radius[idx]; }
    double get_epsilon() const { return ptr->epsilon[idx]; }
    double get_sqrt_epsilon() const { return ptr->sqrt_epsilon[idx]; }
    double get_sigma() const { return ptr->sigma[idx]; }
  };

  Proxy operator[](size_t i) { return Proxy{this, i}; }
  ConstProxy operator[](size_t i) const { return ConstProxy{this, i}; }

  size_t size() const { return charge.size(); }
  bool empty() const { return charge.empty(); }
  void clear() {
    charge.clear();
    radius.clear();
    sigma.clear();
    epsilon.clear();
    sqrt_epsilon.clear();
    mass.clear();
  }
  void resize(size_t n) {
    charge.resize(n, 0.0);
    radius.resize(n, 0.0);
    sigma.resize(n, 0.0);
    epsilon.resize(n, 0.0);
    sqrt_epsilon.resize(n, 0.0);
    mass.resize(n, 1.0);
  }

  // Support range-based for loops
  struct Iterator {
    Multi_MOL_parameters *ptr;
    size_t idx;
    Proxy operator*() { return Proxy{ptr, idx}; }
    Iterator &operator++() {
      ++idx;
      return *this;
    }
    bool operator!=(const Iterator &other) const { return idx != other.idx; }
  };

  struct ConstIterator {
    const Multi_MOL_parameters *ptr;
    size_t idx;
    ConstProxy operator*() const { return ConstProxy{ptr, idx}; }
    ConstIterator &operator++() {
      ++idx;
      return *this;
    }
    bool operator!=(const ConstIterator &other) const {
      return idx != other.idx;
    }
  };

  Iterator begin() { return Iterator{this, 0}; }
  Iterator end() { return Iterator{this, charge.size()}; }
  ConstIterator begin() const { return ConstIterator{this, 0}; }
  ConstIterator end() const { return ConstIterator{this, charge.size()}; }
};
#endif