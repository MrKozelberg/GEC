#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <functional>

double integrate_trap(std::function<double(double)> f, double a, double b, size_t n);
double integrate_Simpson(std::function<double(double)> f, double a, double b, size_t n);

#endif // INTEGRAL_H
