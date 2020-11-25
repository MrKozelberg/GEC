#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <cmath>
#include <functional>
#include "integral.h"

double integrate_trap(std::function<double(double)> f, double a, double b, size_t n) {
    using namespace boost::numeric::ublas;
    vector<double> func(n), help(n);
    double step = (b - a) / double (n);
    //filling vectors with values
    for (size_t i = 0; i < n; ++ i)
    {
        func(i) = f(i * step);
        if (i != 0 or i != n){
            help(i) = 2;
        } else {
            help(i) = 2;
        }
    }
    help *= (step / 2);
    return inner_prod(func, help);
}

double integrate_Simpson(std::function<double(double)> f, double a, double b, size_t n) {
    using namespace boost::numeric::ublas;
    try {
        if (n % 2 != 1){
            throw n;
        }
    } catch (int i) {
        std::cout << "integrate_Simpson has a wrong argument" << std::endl;
        exit(-1);
    }
    double step = (b - a) / double (n);
    vector<double> func(n), help(n);
    vector<double> A(2);
    A(0) = 2; A(1) = 4;
    for (size_t i = 0; i < n; ++ i)
    {
        func[i] = f(i * step);
        if (i != 0 or i != n){
            help(i) = A(i % 2);
        } else help(i) = 1;
    }
    help *= (step / 3);
    return inner_prod(func, help);
}

