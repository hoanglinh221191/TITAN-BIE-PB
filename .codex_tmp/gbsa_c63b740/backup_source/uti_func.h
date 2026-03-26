#ifndef UTI_FUNCTION_H
#define UTI_FUNCTION_H
#include<utility>
using namespace std;
template<typename T>
inline double square(T&& base) { long double t{ base*base }; return t; }
template<typename T>
inline double cube(T&& base) { return square(forward<T>(base))*base; }
template<typename T>
inline double c6(T&& base) { long double temp{ cube(base*base) }; return temp; }
template<typename T>
inline double c12(T&& base) { long double temp{ c6(std::forward<T>(base)) }; return temp*temp; }
template<typename T1, typename T2, typename T3>
inline double norm3(T1&& x, T2&& y, T3&& z) { return sqrt(square(forward<T1>(x)) + square(forward<T2>(y)) + square(forward<T3>(z))); }
#endif