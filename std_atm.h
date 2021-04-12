#ifndef STDATM_H
#define STDATM_H

#include <cmath>

///It calculates a geopotential altitude from a geometric altitude [km]
/// |DEF| Geometric height is an altitude above mean sea level.
/// |DEF| Geopotential is an adjustment to geometric height that
///       accounts for the variation of gravity with latitude and altitude.
double gp_from_geom(double H)
{
    double earth_radius = 6356.766; ///[km] according to the U.S. Standard Atmosphere (1976)
    return H * earth_radius / (H + earth_radius);
}

///altitude is from 0 km to 70 km
///temperature and pressure can be calculated
class StdAtm {
private:
    constexpr static double T_0 = 288.15;
    constexpr static double p_0 = 1.01325e5;
    constexpr static double R = 8.31432;
    constexpr static double M = 28.9644;
    constexpr static double g = 9.80665;
    constexpr static double z[7] = { 0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 70.0 };
    constexpr static double T[7] = { T_0, 216.65, 216.65, 228.65, 270.65, 270.65, 217.45 };
    constexpr static double p[7] = { p_0, 22632.1, 5474.89, 868.019, 110.906, 66.9389, 4.63422 };
    constexpr static double gamma[7] = { 0.0, -6.5, 0.0, 1.0, 2.8, 0.0, -2.8 };


public:
    StdAtm() = default;;

    static double temperature(double H)
    {
        double zeta = gp_from_geom(H);
        size_t n = 0;
        for (size_t k = 1; k < 7; ++k) {
            if (zeta >= z[k - 1] and z[k] > zeta) {
                n = k;
                break;
            } else {
                n = 6;
            }
        }
        return T[n - 1] + gamma[n] * (zeta - z[n - 1]);
    }
    static double pressure(double H)
    {
        double alt = gp_from_geom(H);
        size_t n = 6;
        for (size_t k = 1; k < 7; ++k) {
            if (alt >= z[k - 1] and z[k] > alt) {
                n = k;
                break;
            } else {
                n = 6;
            }
        }
        if (gamma[n] == 0.0) {
            return p[n - 1] * exp(-g * M * (alt - z[n - 1]) / R / T[n - 1]);
        } else {
            return p[n - 1] * pow((1.0 + gamma[n] * (alt - z[n - 1]) / T[n - 1]), -g * M / gamma[n] / R);
        }
    }
};

#endif // STDATM_H
