#include <cmath>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "integral.h"
#include "sigma.h"
#include "geomagnetic.h"

/// This is a struct of columns
/// You enter them in main() and take an ionosphere potential
struct Column {
    double area;
    std::function<double(double)> conductivity;
    std::function<double(double)> current;
    double bound_val;
};

/// explicit conductivity functions
class ParentCond : protected SMZ15, protected TZ06, protected ZT07 {
protected:
    static constexpr double sigma_0 = 5.0e-14;
    static constexpr double H_0 = 6.0; /// [km]
    static constexpr double e_0 = 1.602176634e-19; /// [C]  
};

///This class based upon the 'sigma.h'
///z in kilometres
///lambda is geomagnetic longitude [rad]
class [[maybe_unused]] Conductivity : protected ParentCond {
public:
    static double sigma(double z, double lambda, double xi)
    {
        ZT07 m1; SMZ15 m2; TZ06 m3;
        return e_0 * m1.mu_sum(z) * sqrt(m2.q(z, lambda, xi) / m3.alpha(z));
    }
};

///Exponential conductivity function is for fast tests
///z in kilometres
class [[maybe_unused]] ExpSigma: protected ParentCond {
public:
    static double sigma(double z, ...)
    {
        return sigma_0 * exp(z / H_0);
    };
};

class [[maybe_unused]] ConstSigma: protected ParentCond {
public:
    static double sigma(...)
    {
        return sigma_0;
    }
};


/// explicit current functions

class ParentCurr {
protected:
    static constexpr double j_0 = 1.2e-10;
};

class StepJ : protected ParentCurr {
public:
    static double j(double z, ...)
    {
        return (z >= 6.0 and z < 11.0) ? j_0 : 0.0;
    }
};

class ZeroJ : protected ParentCurr {
public:
    static double j(...)
    {
        return 0.0;
    }
};

class [[maybe_unused]] SimpleGeoJ : protected ParentCurr, private ZeroJ, private StepJ {
public:
    static double j(double z, double lat)
    {
        return (std::abs(lat) <= 5) ? StepJ::j(z) : ZeroJ::j(z);
    }
};

/// explicit boundary values function
///it is an empty parent class for a while
class ParentBoundVal {
};

class [[maybe_unused]] ZeroPhiS : private ParentBoundVal {
public:
    static double phi_s(...)
    {
        return 0.0;
    }
};

/// It is the central class, here you can find any parameters of the model
class GECModel {
protected:
    static constexpr double H = 70.0; ///in km
    static constexpr double pot_step = 1.0; ///in km
    static constexpr unsigned int_points = 701;
    static constexpr unsigned int_pot_points = 11;
    static constexpr size_t N = H / pot_step + 1;
    double V = 0.0;
    bool isVCalculated = false;
    std::vector<Column> model;

    /// This is a function that calculates an ionosphere potential
    double calc_ip()
    {
        std::vector<double> area_by_int(model.size());
        std::vector<double> int_of_cur_by_cond(model.size());
        for (unsigned i = 0; i < model.size(); ++i) {
            area_by_int[i] = model[i].area / integrate_Simpson([this, i](double z)
                        { return 1.0 / model[i].conductivity(z); },
                    0.0,
                    H,
                    int_points);
            int_of_cur_by_cond[i] = integrate_Simpson([this, i](double z)
                        { return model[i].current(z) / model[i].conductivity(z); },
                    0.0,
                    H,
                    int_points);
        }
        double up = 0.0, down = 0.0;
        for (unsigned i = 0; i < model.size(); ++i) {
            up += area_by_int[i] * (int_of_cur_by_cond[i] - model[i].bound_val);
            down += area_by_int[i];
        }
        return up / down;
    }

    /// This is a function that calculates the potential on the i * pot_step
    std::array<double, N> calc_pot(unsigned column_num)
    {
        std::array<double, N> vec{};
        double I1 = 0.0, I2 = 0.0, C1, C2;
        C1 = integrate_Simpson([this, column_num](double z)
                { return 1.0 / model[column_num].conductivity(z); },
                0.0, H, int_points);
        C2 = integrate_Simpson([this, column_num](double z)
                { return model[column_num].current(z) / model[column_num].conductivity(z); },
                0.0, H, int_points);
        std::array<double, N> h{};
        for (unsigned n = 0; n < N; ++n) {
            h[n] = n * pot_step;
        }
        for (unsigned n = 1; n < N; ++n) {
            I1 += integrate_Simpson([this, column_num](double z)
                    { return model[column_num].current(z) / model[column_num].conductivity(z); },
                    h[n - 1], h[n], int_pot_points);
            I2 += integrate_Simpson([this, column_num](double z)
                    { return 1.0 / model[column_num].conductivity(z); },
                    h[n - 1], h[n], int_pot_points);
            vec[n] = I1 - I2 * (C2 - model[column_num].bound_val - getIP()) / C1;
        }
        return vec;
    }

public:
    explicit GECModel() = default;
    double getIP()
    {
        if (not isVCalculated) {
            V = calc_ip();
            isVCalculated = true;
        }
        return V;
    }

    [[maybe_unused]] void getPot(const std::string& filename, unsigned column_num)
    {
        std::array<double, N> vec = calc_pot(column_num);
        std::ofstream fout(filename);
        if (!fout.is_open()) {
            std::cout << "Impossible to find a file" << std::endl;
            exit(-1);
        }
        for (unsigned n = 1; n < N; ++n) {
            fout << n * pot_step << "\t" << vec[n] << std::endl;
        }
        fout.close();
    }
};

/// It is the latitude and longitude grid with the parametrization you want
template <class Cond, class Curr, class BoundVal>
class [[maybe_unused]] GeoModel : public GECModel {
protected:
    static constexpr double earth_radius2 = 6370.0*6370.0; ///< km^2
    unsigned N, M;
    double delta_lat, delta_lon;
    double cell_area(unsigned n, unsigned m, double d_lat, double d_lon)
    {
        double lat_n = -90.0 + d_lat * n;
        if (n != N-1 and m != M-1) {
            return fabs(earth_radius2 * M_PI / 180.0 * d_lon * (sin(M_PI / 180.0 * (lat_n + d_lat)) - sin(M_PI / 180.0 * lat_n)));
        } else {
            if (m == M-1) {
                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - m * d_lon) * (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            } else {
                return fabs(earth_radius2 * M_PI / 180.0 * d_lon * (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            }
        }
    }

    double lat_arg(unsigned n, double d_lat)
    {
        double lat_n = -90.0 + d_lat * n;
        if (n == N-1) {
            return 0.5 * (lat_n + 90.0);
        } else {
            return lat_n + 0.5 * d_lat;
        }
    }

    double lon_arg(unsigned m, double d_lon)
    {
        double lon_m = d_lon * m;
        if (m == M-1) {
            return 0.5 * (lon_m + 360.0);
        } else {
            return lon_m + 0.5 * d_lon;
        }
    }

public:
    /// \param K
    /// K == true => the arguments are coordinate steps
    /// K != false => the arguments define the number of steps
    [[maybe_unused]] GeoModel(double arg1, double arg2, bool K)
    {
        if (K) {
            delta_lat = arg1;
            delta_lon = arg2;
            N = std::ceil(180.0 / delta_lat);
            M = std::ceil(360.0 / delta_lon);

        } else {
            N = (unsigned) arg1;
            M = (unsigned) arg2;
            delta_lat = 180.0 / N;
            delta_lon = 360.0 / M;
        }
        model.reserve(N * M);
        for (unsigned n = 0; n < N; ++n) {
            for (unsigned m = 0; m < M; ++m) {
                model.push_back({ cell_area(n, m, delta_lat, delta_lon),
                                  [this, n](double z){return Cond::sigma(z, lat_arg(n, delta_lat), 0.5);},
                                  [n, this](double z){return Curr::j(z, lat_arg(n, delta_lat));},
                                  BoundVal::phi_s(lat_arg(n, delta_lat), lon_arg(m, delta_lon))
                                });
            }
        }
    }
};


/// \brief main
/// \return IP
int main()
{
    GeoModel<Conductivity, SimpleGeoJ, ZeroPhiS> m(180.0, 360.0, false);
    //m.getPot("plots/potential_2_columns.txt", 180*180);
    std::cout << "Ionosphere potential is " << m.getIP() << " [kV]" << std::endl;
    /*double latm, longm, altm;
    gdz_to_mag(2020.2, 50.0, 50.0, 10.0, latm, longm, altm);
    std::cout << "latm = " << latm << "\n"
              << "longm = " << longm << "\n"
              << "altm = " << altm << "\n";*/
    return EXIT_SUCCESS;
}
