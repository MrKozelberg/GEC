#ifndef SIGMA_H
#define SIGMA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <iomanip>

#include "integral.h"
#include "StdAtm.h"

/// COMPUTATION OF THE ION-PAIR PRODUCTION RATE ASSOCIATED WITH GALACTIC COSMIC RAYS
class SMZ15: private StdAtm{
private:
    double L(double lymbda, double lymbda_0){
        if(std::abs(lymbda) < lymbda_0){
            return std::abs(lymbda);
        } else {
            return lymbda_0;
        }
    }
    /// constants
    static constexpr double T_q = 297.15;            ///[K]
    static constexpr double p_q = 98658.55232;       ///[Pa]
    std::vector<double> A = {49.0 * M_PI / 180.0,
                             49.0 * M_PI / 180.0,
                             57.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0};
    std::vector<double> a = {8.0 * M_PI / 180.0,
                             8.0 * M_PI / 180.0,
                             8.0 * M_PI / 180.0,
                             10.0 * M_PI / 180.0,
                             10.0 * M_PI / 180.0,
                             9.0 * M_PI / 180.0};
    std::vector<double> H = {0.0, 5.0, 10.0, 16.0, 23.0, 37.0};                     ///[km]
    std::vector<double> B = {0.0, 0.0, 5.0, 21.0, 14.0, 0.0};                       ///[km]
    std::vector<double> b = {0.0, 0.0, 0.0, 14.0, 9.0, 0.0};                        ///[km]
    std::vector<double> C = {1.4E6, 15.0E6, 64.0E6, 93.0E6, 74.0E6, 45.0E6};        ///[m^(-3)*s^(-1)]
    std::vector<double> c = {0.1E6, 0.5E6, 2.0E6, 5.0E6, 3.0E6, 2.0E6};             ///[m^(-3)*s^(-1)]
    std::vector<double> D = {0.8E6, 10.0E6, 236.0E6, 402.0E6, 421.0E6, 450.0E6};    ///[m^(-3)*s^(-1)]
    std::vector<double> d = {0.1E6, 2.5E6, 83.0E6, 225.0E6, 236.0E6, 243.0E6};      ///[m^(-3)*s^(-1)]
    std::vector<double> g = {0.0, 0.0, 4.0, 6.0, 5.0, 0.0};
    std::vector<double> h = {1.7, 1.7, 3.9, 3.2, 3.4, 4.0};
    std::vector<double> K(double xi){
        std::vector<double> vec(6);
        for(size_t i = 0; i < 6; ++i){
            vec[i] = A[i] - a[i] * xi;
        }
        return vec;
    }
    std::vector<double> deltaH(double xi){
        std::vector<double> vec(6);
        vec[0] = 0.0;
        for(size_t i = 1; i < 6; ++i){
            vec[i] = B[i] - b[i] * xi;
        }
        return vec;
    }
    std::vector<double> U(double xi){
        std::vector<double> vec(6);
        for(size_t i = 0; i < 6; ++i){
            vec[i] = C[i] - c[i] * xi;
        }
        return vec;
    }
    std::vector<double> deltaU(double xi){
        std::vector<double> vec(6);
        for(size_t i = 0; i < 6; ++i){
            vec[i] = D[i] - d[i] * xi;
        }
        return vec;
    }
    std::vector<double> Z(double lymbda, double xi){
        std::vector<double> vec(6);
        vec[0] = 0.0;
        for(size_t i = 1; i < 6; ++i){
            vec[i] = H[i] + deltaH(xi)[i] * pow(sin(L(lymbda, K(xi)[i])) / sin(K(xi)[i]), g[i]);
        }
        return vec;
    }
    std::vector<double> Q(double lymbda, double xi){
        std::vector<double> vec(6);
        for(size_t i = 0; i < 6; ++i){
            vec[i] = U(xi)[i] + deltaU(xi)[i] * pow(sin(L(lymbda, K(xi)[i])) / sin(K(xi)[i]), h[i]);
        }
        return vec;
    }
    std::vector<double> P(double lymbda, double xi){
        std::vector<double> vec(5);
        vec[0] = 0.0; vec[3] = 0.0;
        vec[1] = Q(lymbda, xi)[1] * log(Q(lymbda, xi)[1] / Q(lymbda, xi)[0]) / Z(lymbda, xi)[1];
        vec[2] =  2.0 * (Q(lymbda, xi)[2] - Q(lymbda, xi)[1]) / (Z(lymbda, xi)[2] - Z(lymbda, xi)[1])
               - vec[1];
        if (Z(lymbda, xi)[4] != Z(lymbda, xi)[5]){
            vec[4] =  3.0 * (Q(lymbda, xi)[4] - Q(lymbda, xi)[5]) / (Z(lymbda, xi)[5] - Z(lymbda, xi)[4]);
        } else vec[4] = 0.0;
        return vec;
    }
    double A_Q(double lymbda, double xi){
        return ((Q(lymbda, xi)[2] - Q(lymbda, xi)[1]) - P(lymbda, xi)[1] * (Z(lymbda, xi)[2] -
                Z(lymbda, xi)[1])) / pow((Z(lymbda, xi)[2] - Z(lymbda, xi)[1]), 2.0);
    }
    double B_Q(double lymbda, double xi){
        return Q(lymbda, xi)[1] - pow(P(lymbda, xi)[1] * (Z(lymbda, xi)[2] - Z(lymbda, xi)[1]), 2.0)
               / (4.0 * ((Q(lymbda, xi)[2] - Q(lymbda, xi)[1]) - P(lymbda, xi)[1] * (Z(lymbda, xi)[2] -
                Z(lymbda, xi)[1])));
    }
    double C_Q(double lymbda, double xi){
        return (2.0 * (Q(lymbda, xi)[2] - Q(lymbda, xi)[1]) * Z(lymbda, xi)[1] - P(lymbda, xi)[1]
                * (pow(Z(lymbda, xi)[2], 2.0) - pow(Z(lymbda, xi)[1], 2.0))) / (2.0 *
                ((Q(lymbda, xi)[2] - Q(lymbda, xi)[1]) - P(lymbda, xi)[1] * (Z(lymbda, xi)[2] -
                Z(lymbda, xi)[1])));
    }
    double Q_STP(double z, double lymbda, double xi){
        std::vector<double> Z_ = Z(lymbda, xi), P_ = P(lymbda, xi), Q_ = Q(lymbda, xi);
        double A_Q_ = A_Q(lymbda, xi), B_Q_ = B_Q(lymbda, xi), C_Q_ = C_Q(lymbda, xi);
        if(z < Z_[1]){
            return Q_[0] * pow(Q_[1] / Q_[0], z / Z_[1]);
        } else if(z < Z_[2]){
            return B_Q_ + A_Q_ * pow(z - C_Q_, 2.0);
        } else if(z < Z_[3]){
            return Q_[3] - (Q_[3] - Q_[2]) * pow((Z_[3] - z) / (Z_[3] - Z_[2]), P_[2] *
                   (Z_[3] - Z_[2]) / (Q_[3] - Q_[2]));
        } else if(z < Z_[4]){
            return Q_[3] - (Q_[3] - Q_[4]) * pow((z - Z_[3]) / (Z_[4] - Z_[3]), P_[4] *
                   (Z_[4] - Z_[3]) / (Q_[3] - Q_[4]));
        } else if(z < Z_[5]){
            return Q_[5] + (Q_[4] - Q_[5]) * pow((Z_[5] - z) / (Z_[5] - Z_[4]), 3.0);
        } else{
            return Q_[5];
        }
    }

public:
    SMZ15(){};
    ///ion-pair production rate
    ///z is a geometrical altitude [km] (LOOK DEF ABOVE)
    ///p [Pa] | T [K] | lymbda [deg] | q, Q [m^(-3)*s^(-1)] | xi = (sin(pi*t))^2
    /// takes the altitude z in km, the latitude lat in degrees and the parameter xi representing the solar cycle phase
    /// xi is confined between 0 and 1; xi = 0 corresponds to solar minima, xi = 1 corresponds to solar maxima
    /// returns the ion-pair production rate q in m^(-3) s^(-1)
    /// uses the parameterisation of Slyunyaev et al. (2015)
    double q(double z, double lymbda, double xi){
        return Q_STP(z, lymbda, xi) * StdAtm::pressure(z) * T_q / (p_q * StdAtm::temperature(z));
    }
};

/// COMPUTATION OF THE ION-ION RECOMBINATION COEFFICIENT [TZ 06]
class TZ06: private StdAtm{
private:
    /// constants
    static constexpr double A = 6.0E-14;                                                  ///[m^3 / s]
    static constexpr double a = 0.5;
    std::vector<double> B = {0.0, 1.702E-12, 1.035E-12, 6.471E-12};                      ///[m^3 / s]
    static constexpr double T_0 = 300.0;                                                 ///[K]
    static constexpr double N_0 = 2.69E25;                                               ///[m^(-3)]
    std::vector<double> b = {0.0, -1.984, 4.374, -0.191};
    std::vector<double> Z = {0.0, 10.0, 20.0};                                           ///[km]
    std::vector<double> c = {0.0, -0.451, 0.769, 0.901};
    /// Boltzmann constant
    static constexpr double k = 1.380649E-23;                                            ///[J/K]
    /// concentration of air molecules
    double N(double z){
        return StdAtm::pressure(z) / (k * StdAtm::temperature(z));
    }
public:
    TZ06(){};
    /// ION-ION RECOMBINATION COEFFICIENT alpha
    /// z [km] | p [Pa] | T [K]
    double alpha(double z){
        double var = A * pow(T_0 / StdAtm::temperature(z), a);
        if(z < Z[1]){
            var += B[1] * pow(T_0 / StdAtm::temperature(z), b[1]) * pow(N(z) / N_0, c[1]);
        } else if(z < Z[2]){
            var += B[2] * pow(T_0 / StdAtm::temperature(z), b[2]) * pow(N(z) / N_0, c[2]);
        } else{
            var += B[3] * pow(T_0 / StdAtm::temperature(z), b[3]) * pow(N(z) / N_0, c[3]);
        }
        return var;
    }
};

/// COMPUTATION OF THE ION MOBILITIES
class ZT07: private StdAtm{
private:
    /// constants
    static constexpr double mu_0_plus = 1.4E-4;              /// [m^2/(V*s)]
    static constexpr double mu_0_minus = 1.9E-4;             /// [m^2/(V*s)]
    static constexpr double T_mu = 120;                      /// [K]
    static constexpr double T_0 = 288.15;                    /// [K]
    static constexpr double p_0 = 101'325;                   /// [Pa]
public:
    ZT07(){};
    /// ion mobilities [m^2/(V*s)]
    /// p [Pa] | T [K]
    double mu_plus(double p, double T){
        return mu_0_plus * p_0 / p * pow(T / T_0, 1.5) * (T_0 + T_mu) / (T + T_mu);
    }
    double mu_minus(double p, double T){
        return mu_0_minus * p_0 / p * pow(T / T_0, 1.5) * (T_0 + T_mu) / (T + T_mu);
    }
    double mu_sum(double z){
        return mu_plus(StdAtm::pressure(z), StdAtm::temperature(z)) +
                mu_minus(StdAtm::pressure(z), StdAtm::temperature(z));
    }
};


#endif // SIGMA_H