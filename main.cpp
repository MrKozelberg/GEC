#include <cmath>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <fstream>
#include "cnpy/cnpy.h"
#include <iostream>
#include <string>
#include <vector>

#include "integral.h"
#include "sigma.h"
#include "geomagnetic.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-comma-subscript"
/// This is a struct of columns
/// You enter them in main() and take an ionosphere potential
struct Column {
    double area;
    std::function<double(double)> conductivity;
    std::function<double(double)> current;
    double bound_val;
};

/// It is the central class, here you can find any parameters of the model
class GECModel {
protected:
    static constexpr double H = 70.0; ///in km
    static constexpr double pot_step = 1.0; ///in km
    static constexpr unsigned int_points = 271;
    static constexpr unsigned int_pot_points = 5;
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
    //input DATA
// py_array[n,m,k] is cpp_array[n*180*360 + m*360 + k]
// You should check file type
// '<f4' is equavalent of float
// '<f8' --- double
    cnpy::npz_t data = cnpy::npz_load("data/DATA-2015-12-31-00.npz");
    cnpy::NpyArray cape_arr = data["cape"];
    cnpy::NpyArray cbot_arr = data["cbot"];
    cnpy::NpyArray ctop_arr = data["ctop"];
    cnpy::NpyArray alpha_arr = cnpy::npy_load("data/alpha.npy");
    float_t* cape = cape_arr.data<float>();
    double_t* cbot = cbot_arr.data<double>();
    double_t* ctop = ctop_arr.data<double>();
    double_t* alpha_ = alpha_arr.data<double>();
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
class Conductivity : protected ParentCond {
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
    static constexpr double j_0 = 1.72e-10;
};

class StepJ : protected ParentCurr {
public:
    static double j(double z, ...)
    {
        return (z >= 6.0 and z < 11.0) ? j_0 : 0.0;
    }
};

class [[maybe_unused]] ZeroJ : protected ParentCurr {
public:
    static double j(...)
    {
        return 0.0;
    }
};

class [[maybe_unused]] SimpleGeoJ : protected ZeroJ, private StepJ {
public:
    static double j(double z, double lat)
    {
        return (std::abs(lat) <= 10) ? StepJ::j(z) : ZeroJ::j(z);
    }
};

class GeoJ : public ZeroJ, public GECModel  {
private:
    static constexpr double cape_0 = 1'000; // J/kg
public:
    double j(double z, int t, unsigned lat, unsigned lon) {
        return GECModel::cape[t*180*360 + lat*360 + lon] < cape_0
               || GECModel::cbot[t*180*360 + lat*360 + lon] >= z * 1'000
               || GECModel::ctop[t*180*360 + lat*360 + lon] <= z * 1'000 ? 0 : j_0;
    }
};

/// explicit boundary values function
///it is an empty parent class for a while
class ParentBoundVal {
protected:
    static constexpr double phi_0 = 15.0; //[kV]
    static constexpr double theta_0 = 20.0; //[deg]
};

class [[maybe_unused]] ZeroPhiS : private ParentBoundVal {
public:
    static double phi_s(...)
    {
        return 0.0;
    }
};

class [[maybe_unused]] VollandPhiS : protected ParentBoundVal {
public:
    //theta is a longitude [deg]
    static double k(double theta)
    {
        try {
            if (theta < theta_0) {
                throw theta; // NOLINT(hicpp-exception-baseclass)
            }
        } catch (double i) {
            std::cout << "k(theta) has a wrong argument" << std::endl;
            exit(-2);
        }
        if (theta < 30.0) {
            return 1.0;
        } else if (theta >= 30.0 && theta < 90.0) {
            return (1 + 2 * sin(3.0 * theta)) / 2;
        }
        return 0.0;
    }
    //theta = geomagnetic longitude
    //lambda = latitude [deg]
    static double phi_s(double theta, double lambda)
    {
        return phi_0 * sin(lambda) * ((theta < theta_0) ?
                                       sin(theta) / sin(theta_0) :
                                       k(theta) * pow(sin(theta_0) / sin(theta), 4) );
    }
};

class [[maybe_unused]] ParentGrid {
protected:
    static constexpr double earth_radius2 = 6370.0 * 6370.0; ///< km^2
    unsigned N, M;
    double delta_lat, delta_lon;

public:
    void set_N(unsigned n){
        N = n;
    }
    void set_M(unsigned m){
        M = m;
    }
    void set_delta_lat(double lat){
        delta_lat = lat;
    }
    void set_delta_lon(double lon){
        delta_lon = lon;
    }

    unsigned get_N() const{
        return N;
    }
    unsigned get_M() const{
        return M;
    }
    double get_delta_lat() const{
        return delta_lat;
    }
    double get_delta_lon() const{
        return delta_lon;
    }

    double cell_area(unsigned n, unsigned m, double d_lat, double d_lon) {
        double lat_n = -90.0 + d_lat * n;
        if (n != N - 1 and m != M - 1) {
            return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
                        (sin(M_PI / 180.0 * (lat_n + d_lat)) - sin(M_PI / 180.0 * lat_n)));
        } else {
            if (m == M - 1) {
                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - m * d_lon) *
                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            } else {
                return fabs(
                        earth_radius2 * M_PI / 180.0 * d_lon * (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            }
        }
    }

    double lat_arg(unsigned n, double d_lat) {
        double lat_n = -90.0 + d_lat * n;
        if (n == N - 1) {
            return 0.5 * (lat_n + 90.0);
        } else {
            return lat_n + 0.5 * d_lat;
        }
    }

    double lon_arg(unsigned m, double d_lon) {
        double lon_m = d_lon * m;
        if (m == M - 1) {
            return 0.5 * (lon_m + 360.0);
        } else {
            return lon_m + 0.5 * d_lon;
        }
    }

    double lon_arg_m(double lon_arg_, double lat_arg_) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(2015.9, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lon_m;
    }

    double lat_arg_m(double lon_arg_, double lat_arg_) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(2015.9, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lat_m;
    }
};


/// It is the latitude and longitude grid with the parametrization you want
template <class Cond, class Curr, class BoundVal>
class [[maybe_unused]] GeoGrid : public ParentGrid, public GECModel {
public:
    GeoGrid(double arg1, double arg2, bool K)
    {
            if (K) {
                set_delta_lat(arg1);
                set_delta_lon(arg2);
                set_N(std::ceil(180.0 / get_delta_lat()));
                set_M(std::ceil(360.0 / get_delta_lon()));
                model.reserve(get_N() * get_M());
            } else {
                set_N(arg1);
                set_M(arg2);
                model.reserve(get_N() * get_M());
                set_delta_lat(180.0 / get_N());
                set_delta_lon(360.0 / get_M());
            }
            double lat_n = -90.0;
            for (unsigned n = 0; n < get_N(); ++n) {
                lat_n =+ delta_lat * n;
                for (unsigned m = 0; m < get_M(); ++m) {
                    model.push_back({ cell_area(n, m, get_delta_lat(), get_delta_lon()),
                                      [this, n, m](double z)
                                        {return Cond::sigma(z,
                                                            lat_arg_m(lon_arg(m, delta_lon),
                                                                      lat_arg(n, delta_lat)),
                                                            0.5);
                                        },
                                      [n, this](double z){return Curr::j(z, lat_arg(n, delta_lat));},
                                      BoundVal::phi_s(lat_arg(n, delta_lat), lon_arg(m, delta_lon))
                                    });

                }
            }
    }
};

class GeoGridNickolay : public ParentGrid, public Conductivity, public GeoJ, public ZeroPhiS {
private:
    static constexpr double k = 1.0;
public:
    GeoGridNickolay(int t)
    {
        set_delta_lat(1.0);
        set_delta_lon(1.0);
        set_N(180);
        set_M(360);
        model.reserve(2 * get_N() * get_M());
        double lat_n = -90.0;
        for (unsigned n = 0; n < get_N(); ++n) {
            lat_n = +delta_lat * n;
            for (unsigned m = 0; m < get_M(); ++m) {
                model.push_back({cell_area(n, m, get_delta_lat(), get_delta_lon()) * (1 - k * alpha_[t*180*360 + n*360 + m]),
                                 [this, n, m](double z) {
                                     return Conductivity::sigma(z,
                                                                lat_arg_m(lon_arg(m, delta_lon),
                                                                          lat_arg(n, delta_lat)),
                                                                0.5);
                                 },
                                 [n, this](double z) { return ZeroJ::j(); },
                                 ZeroPhiS::phi_s(lat_arg(n, delta_lat), lon_arg(m, delta_lon))
                                });
                model.push_back({cell_area(n, m, get_delta_lat(), get_delta_lon()) * k * alpha_[t*180*360 + n*360 + m],
                                 [this, n, m](double z) {
                                     return Conductivity::sigma(z,
                                                                lat_arg_m(lon_arg(m, delta_lon),
                                                                          lat_arg(n, delta_lat)),
                                                                0.5);
                                 },
                                 [n, m, t, this](double z) { return GeoJ::j(z, t, n, m); },
                                 ZeroPhiS::phi_s(lat_arg(n, delta_lat), lon_arg(m, delta_lon))
                                });
            }
        }
    }
};

int main()
{
    GeoGrid<Conductivity, SimpleGeoJ, ZeroPhiS> m(1.0, 1.0, true);
    //m.getPot("plots/potential_2_columns_sourse.txt", 90*360);
    std::cout << "Ionosphere potential is " << m.getIP() << " [kV]" << std::endl;

/// Nickolay's data
    /*std::ofstream fout("plots/diurnal_variation.txt");
    if (!fout.is_open()) {
        std::cout << "Impossible to find a file" << std::endl;
        exit(-1);
    }
    int i = 1;
    GeoGridNickolay m(i);
    for (int i = 1; i <= 48; i++) {
        GeoGridNickolay m(i);
        fout << i << '\t' << m.getIP() << '\n';
    }
    fout.close();*/
    return EXIT_SUCCESS;
}
    //GeoGrid<Conductivity, SimpleGeoJ, ZeroPhiS> m(1.0, 1.0, true);
    //m.getPot("plots/potential_2_columns.txt", 180*180);
    //std::cout << "Ionosphere potential is " << m.getIP() << " [kV]" << std::endl;
    /*double latm, longm, altm;
    gdz_to_mag(2020.2, 50.0, 50.0, 10.0, latm, longm, altm);
    std::cout << "latm = " << latm << "\n"
              << "longm = " << longm << "\n"
              << "altm = " << altm << "\n";*/
