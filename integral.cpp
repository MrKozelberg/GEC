#include "integral.h"
#include <cmath>
#include <functional>
#include <iostream>

double integrate_trap(const std::function<double(double)>& f, double a, double b, size_t n)
{
    double func[n], help[n];
    const double step = (b - a) / double(n);
    // filling vectors with values
    for (size_t i = 0; i < n; ++i) {
        func[i] = f(a + i * step);
        help[i] = 2;
    }
    double ans = 0.0;
    for (size_t i = 0; i < n; ++i) {
        help[i] *= (step / 2);
        ans += func[i] * help[i];
    }
    return ans;
}

double integrate_Simpson( const std::function<double(double)>& f, double a, double b, size_t n )
{
    try {
        if (n % 2 != 1) {
            throw n;
        }
    } catch (int i) {
        std::cout << "integrate_Simpson has a wrong argument" << std::endl;
        exit(-1);
    }
    const double step = (b - a) / double(n);
    double func[n], help[n];
    double A[2];
    A[0] = 2;
    A[1] = 4;
    for (size_t i = 0; i < n; ++i) {
        func[i] = f(a + i * step);
        if (i != 0 or i != n) {
            help[i] = A[i % 2];
        } else
            help[i] = 1;
    }
    double ans = 0.0;
    for (size_t i = 0; i < n; ++i) {
        help[i] *= ( step / 3 );
        ans += func[i] * help[i];
    }
    return ans;
}
