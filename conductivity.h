//
// Created by mrkozelberg on 12.04.2021.
//

#ifndef GLOBAL_ELECTRIC_CIRCUIT_CONDUCTIVITY_H
#define GLOBAL_ELECTRIC_CIRCUIT_CONDUCTIVITY_H


#include <cmath>
#include <cassert>
#include "std_atm.h"

class IonMobility {
private:
    constexpr static double mu0_p = 1.4E-4; // m^2/(V*s)
    constexpr static double mu0_m = 1.9E-4; // m^2/(V*s)
    constexpr static double T_mu = 120; // K
    constexpr static double T_0 = 288.15; // K
    constexpr static double p_0 = 101325; // Pa
    static double common_factor(double T, double p) {
        return p_0 / p * std::pow(T / T_0, 1.5) * (T_0 + T_mu) / (T + T_mu);
    }

public:
    IonMobility() = default;

    static double mu_p_get(double T, double p) {
        return mu0_p * common_factor(T, p);
    }

    static double mu_m_get(double T, double p) {
        return mu0_m * common_factor(T, p);
    }

};

class IonIonRecombCoeff {
private:
    constexpr static double A = 6E-14; // m^3/s
    constexpr static double B_1 = 1.702E-12; // m^3/s
    constexpr static double B_2 = 1.035E-12; // m^3/s
    constexpr static double B_3 = 6.471E-12; // m^3/s
    constexpr static double T_0 = 300.0; // K
    constexpr static double N_0 = 2.69E25; // m^(-3)
    constexpr static double b_1 = -1.984;
    constexpr static double b_2 = 4.374;
    constexpr static double b_3 = -0.191;
    constexpr static double z_1 = 10.0; // km
    constexpr static double z_2 = 20.0; // km
    constexpr static double c_1 = -0.451;
    constexpr static double c_2 = 0.769;
    constexpr static double c_3 = 0.901;
    constexpr static double k = 1.380649E-23; // J/K
public:
    IonIonRecombCoeff() = default;

    static double alpha(double z, double p, double T) {
        assert(z >= 0);
        if (z < z_1) {
            return A * std::sqrt(T_0 / T) + B_1 * std::pow(T_0 / T, b_1) * std::pow(p / (k * T * N_0), c_1);
        }
        if (z >= z_1 && z < z_2) {
            return A * std::sqrt(T_0 / T) + B_2 * std::pow(T_0 / T, b_2) * std::pow(p / (k * T * N_0), c_2);
        }
        if (z >= z_2) {
            return A * std::sqrt(T_0 / T) + B_3 * std::pow(T_0 / T, b_3) * std::pow(p / (k * T * N_0), c_3);
        }
        return EXIT_FAILURE;
    }
};

class IonPairProdRate {
private:
    constexpr static double A_0 = 49 * M_PI / 180;
    constexpr static double A_1 = 49 * M_PI / 180;
    constexpr static double A_2 = 57 * M_PI / 180;
    constexpr static double A_3 = 63 * M_PI / 180;
    constexpr static double A_4 = 63 * M_PI / 180;
    constexpr static double A_5 = 63 * M_PI / 180;
    constexpr static double a_0 = 8 * M_PI / 180;
    constexpr static double a_1 = 8 * M_PI / 180;
    constexpr static double a_2 = 8 * M_PI / 180;
    constexpr static double a_3 = 10 * M_PI / 180;
    constexpr static double a_4 = 10 * M_PI / 180;
    constexpr static double a_5 = 9 * M_PI / 180;
    constexpr static double H_1 = 5.0; // km
    constexpr static double H_2 = 10.0; // km
    constexpr static double H_3 = 16.0; // km
    constexpr static double H_4 = 23.0; // km
    constexpr static double H_5 = 37.0; // km
    constexpr static double B_1 = 0.0; // km
    constexpr static double B_2 = 5.0; // km
    constexpr static double B_3 = 21.0; // km
    constexpr static double B_4 = 14.0; // km
    constexpr static double B_5 = 0.0; // km
    constexpr static double b_1 = 0.0; // km
    constexpr static double b_2 = 0.0; // km
    constexpr static double b_3 = 14.0; // km
    constexpr static double b_4 = 9.0; // km
    constexpr static double b_5 = 0.0; // km
    constexpr static double C_0 = 1.4E6; // s^(-1)
    constexpr static double C_1 = 15E6; // s^(-1)
    constexpr static double C_2 = 64E6; // s^(-1)
    constexpr static double C_3 = 93E6; // s^(-1)
    constexpr static double C_4 = 74E6; // s^(-1)
    constexpr static double C_5 = 45E6; // s^(-1)
    constexpr static double c_0 = 0.1E6; // m^(-3) * s^(-1)
    constexpr static double c_1 = 0.5E6; // m^(-3) * s^(-1)
    constexpr static double c_2 = 2E6; // m^(-3) * s^(-1)
    constexpr static double c_3 = 5E6; // m^(-3) * s^(-1)
    constexpr static double c_4 = 3E6; // m^(-3) * s^(-1)
    constexpr static double c_5 = 2E6; // m^(-3) * s^(-1)
    constexpr static double D_0 = 0.8E6; // m^(-3) * s^(-1)
    constexpr static double D_1 = 10E6; // m^(-3) * s^(-1)
    constexpr static double D_2 = 236E6; // m^(-3) * s^(-1)
    constexpr static double D_3 = 402E6; // m^(-3) * s^(-1)
    constexpr static double D_4 = 421E6; // m^(-3) * s^(-1)
    constexpr static double D_5 = 450E6; // m^(-3) * s^(-1)
    constexpr static double d_0 = 0.1E6; // m^(-3) * s^(-1)
    constexpr static double d_1 = 2.5E6; // m^(-3) * s^(-1)
    constexpr static double d_2 = 83E6; // m^(-3) * s^(-1)
    constexpr static double d_3 = 225E6; // m^(-3) * s^(-1)
    constexpr static double d_4 = 236E6; // m^(-3) * s^(-1)
    constexpr static double d_5 = 243E6; // m^(-3) * s^(-1)
    constexpr static double g_1 = 0;
    constexpr static double g_2 = 4;
    constexpr static double g_3 = 6;
    constexpr static double g_4 = 5;
    constexpr static double g_5 = 0;
    constexpr static double h_0 = 1.7;
    constexpr static double h_1 = 1.7;
    constexpr static double h_2 = 3.9;
    constexpr static double h_3 = 3.2;
    constexpr static double h_4 = 3.4;
    constexpr static double h_5 = 4;

    static double Lambda(double lambda, double lambda_0) {
        return (std::fabs(lambda) < lambda_0) ? std::fabs(lambda) : lambda_0;
    }

    static double Q(double z, double lambda, double xi) {
        assert(z >= 0);
        // mb optimize this place
        const double K_0 = A_0 - a_0 * xi;
        const double K_1 = A_1 - a_1 * xi;
        const double K_2 = A_2 - a_2 * xi;
        const double K_3 = A_3 - a_3 * xi;
        const double K_4 = A_4 - a_4 * xi;
        const double K_5 = A_5 - a_5 * xi;
        const double delta_H_1 = B_1 - b_1 * xi;
        const double delta_H_2 = B_2 - b_2 * xi;
        const double delta_H_3 = B_3 - b_3 * xi;
        const double delta_H_4 = B_4 - b_4 * xi;
        const double delta_H_5 = B_5 - b_5 * xi;
        const double U_0 = C_0 - c_0 * xi;
        const double U_1 = C_1 - c_1 * xi;
        const double U_2 = C_2 - c_2 * xi;
        const double U_3 = C_3 - c_3 * xi;
        const double U_4 = C_4 - c_4 * xi;
        const double U_5 = C_5 - c_5 * xi;
        const double delta_U_0 = D_0 - d_0 * xi;
        const double delta_U_1 = D_1 - d_1 * xi;
        const double delta_U_2 = D_2 - d_2 * xi;
        const double delta_U_3 = D_3 - d_3 * xi;
        const double delta_U_4 = D_4 - d_4 * xi;
        const double delta_U_5 = D_5 - d_5 * xi;
        const double Z_1 = H_1 + delta_H_1 * std::pow(std::sin(Lambda(lambda, K_1)) / std::sin(K_1), g_1);
        const double Z_2 = H_2 + delta_H_2 * std::pow(std::sin(Lambda(lambda, K_2)) / std::sin(K_2), g_2);
        const double Z_3 = H_3 + delta_H_3 * std::pow(std::sin(Lambda(lambda, K_3)) / std::sin(K_3), g_3);
        const double Z_4 = H_4 + delta_H_4 * std::pow(std::sin(Lambda(lambda, K_4)) / std::sin(K_4), g_4);
        const double Z_5 = H_5 + delta_H_5 * std::pow(std::sin(Lambda(lambda, K_5)) / std::sin(K_5), g_5);
        const double Q_0 = U_0 + delta_U_0 * std::pow(std::sin(Lambda(lambda, K_0)) / std::sin(K_0), h_0);
        const double Q_1 = U_1 + delta_U_1 * std::pow(std::sin(Lambda(lambda, K_1)) / std::sin(K_1), h_1);
        const double Q_2 = U_2 + delta_U_2 * std::pow(std::sin(Lambda(lambda, K_2)) / std::sin(K_2), h_2);
        const double Q_3 = U_3 + delta_U_3 * std::pow(std::sin(Lambda(lambda, K_3)) / std::sin(K_3), h_3);
        const double Q_4 = U_4 + delta_U_4 * std::pow(std::sin(Lambda(lambda, K_4)) / std::sin(K_4), h_4);
        const double Q_5 = U_5 + delta_U_5 * std::pow(std::sin(Lambda(lambda, K_5)) / std::sin(K_5), h_5);
        const double P_1 = Q_1 * std::log(Q_1 / Q_0) / Z_1;
        const double P_2 = 2 * (Q_2 - Q_1) / (Z_2 - Z_1) - P_1;
        const double P_4 = 3 * (Q_4 - Q_5) / (Z_5 - Z_4);
        const double A_Q = ((Q_2 - Q_1) - P_1 * (Z_2 - Z_1)) / (Z_2 - Z_1) / (Z_2 - Z_1);
        const double B_Q = Q_1 - P_1 * P_1 * (Z_2 - Z_1) * (Z_2 - Z_1) / 4 / (Q_2 - Q_1 - P_1 * (Z_2 - Z_1));
        const double C_Q =
                (2 * (Q_2 - Q_1) * Z_1 - P_1 * (Z_2 * Z_2 - Z_1 * Z_1)) / 2 / (Q_2 - Q_1 - P_1 * (Z_2 - Z_1));
        if (z < Z_1) {
            return Q_0 * std::pow(Q_1 / Q_0, z / Z_1);
        }
        if (z >= Z_1 && z < Z_2) {
            return B_Q + A_Q * (z - C_Q) * (z - C_Q);
        }
        if (z >= Z_2 && z < Z_3) {
            return Q_3 - (Q_3 - Q_2) * std::pow((Z_3 - z) / (Z_3 - Z_2), P_2 * (Z_3 - Z_2) / (Q_3 - Q_2));
        }
        if (z >= Z_3 && z < Z_4) {
            return Q_3 - (Q_3 - Q_4) * std::pow((z - Z_3) / (Z_4 - Z_3), P_4 * (Z_4 - Z_3) / (Q_3 - Q_4));
        }
        if (z >= Z_4 && z < Z_5) {
            return Q_5 + (Q_4 - Q_5) * std::pow((Z_5 - z) / (Z_5 - Z_4), 3.0);
        }
        if (z >= Z_5) {
            return Q_5;
        }
        return EXIT_FAILURE;
    }

    constexpr static double T_Q = 297.15; // K
    constexpr static double p_Q = 98658.55232; // Pa
public:
    IonPairProdRate() = default;

    static double q(double z, double lambda, double xi, double p, double T) {
        assert(xi >= 0.0 && xi <= 1.0);
        return Q(z, lambda, xi) * p * T_Q / p_Q / T;
    }
};

#endif //GLOBAL_ELECTRIC_CIRCUIT_CONDUCTIVITY_H
