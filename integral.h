#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <functional>

[[maybe_unused]] double integrate_trap(const std::function<double(double)>& f, double a, double b, size_t n);
double integrate_Simpson(const std::function<double(double)>& f, double a, double b, size_t n);

#endif // INTEGRAL_H
