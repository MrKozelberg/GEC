#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "integral.h"
#include "sigma.h"
#include "std_atm.h"

/// This is srtuct of columns;
/// You enter them in main() and take an ionosphere potential;
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

///z in kilometres
class Conductivity : protected ParentCond {
public:
    double sigma(double z, double lambda, double xi)
    {
        return e_0 * mu_sum(z) * sqrt(q(z, lambda, xi) / alpha(z));
    }
};

class ExpSigma: protected ParentCond {
public:
    static double sigma(double z, ...)
    {
        return sigma_0 * exp(z / H_0);
    };
};

class ConstSigma: protected ParentCond {
public:
    static double sigma(double z, ...)
    {
        return sigma_0;
    }
};

/// explicit current functions
class ParentCurr {
protected:
    static constexpr double j_0 = 1e-9;
};

class StepJ : protected ParentCurr {
public:
    static double j(double z, ...)
    {
        return (z >= 5.0 and z < 10.0) ? j_0 : 0.0;
    }
};

class ZeroJ : protected ParentCurr {
public:
    static double j(double z, ...)
    {
        return 0.0;
    }
};

class SimpleGeoJ : protected ParentCurr, private ZeroJ, private StepJ {
public:
    static double j(double z, double lat)
    {
        return (std::abs(lat) <= 10) ? StepJ::j(z) : ZeroJ::j(z);
    }
};

class ParentBoundVal {
//it is an empty class for a while
};

class ZeroPhiS : private ParentBoundVal {
public:
    static double phi_s(double lat, double lot)
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

    /// This is a function that calculates an ionosphere potention
    double calc_ip()
    {
        std::vector<double> area_by_int(model.size());
        std::vector<double> int_of_cur_by_cond(model.size());
        for (unsigned i = 0; i < model.size(); ++i) {
            area_by_int[i] = model[i].area / integrate_Simpson([this, i](double z) { return 1.0 / model[i].conductivity(z); }, 0.0, H, int_points);
            int_of_cur_by_cond[i] = integrate_Simpson([this, i](double z) { return model[i].current(z) / model[i].conductivity(z); },
                                                      0.0, H, int_points);
        }
        double up = 0.0, down = 0.0;
        for (unsigned i = 0; i < model.size(); ++i) {
            up += area_by_int[i] * (int_of_cur_by_cond[i] - model[i].bound_val);
            down += area_by_int[i];
        }
        return up / down;
    }

    /// This is a function that calculates the potention on the i * pot_step
    std::array<double, N> calc_pot(unsigned column_num)
    {
        std::array<double, N> vec;
        double I1 = 0.0, I2 = 0.0, C1, C2;
        C1 = integrate_Simpson([this, column_num](double z) { return 1.0 / model[column_num].conductivity(z); }, 0.0, H, int_points);
        C2 = integrate_Simpson([this, column_num](double z) { return model[column_num].current(z) / model[column_num].conductivity(z); }, 0.0, H, int_points);
        std::array<double, N> h;
        for (unsigned n = 0; n < N; ++n) {
            h[n] = n * pot_step;
        }
        for (unsigned n = 1; n < N; ++n) {
            I1 += integrate_Simpson([this, column_num](double z) { return model[column_num].current(z) / model[column_num].conductivity(z); },
                h[n - 1], h[n], int_pot_points);
            I2 += integrate_Simpson([this, column_num](double z) { return 1.0 / model[column_num].conductivity(z); },
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
    void getPot(std::string filename, unsigned column_num)
    {
        std::array<double, N> vec = calc_pot(column_num);
        std::ofstream fout(filename);
        if (fout.is_open() == false) {
            std::cout << "Impossible to find a file" << std::endl;
            exit(-1);
        }
        for (unsigned n = 1; n < N; ++n) {
            fout << n * pot_step << "\t" << vec[n] << std::endl;
        }
        fout.close();
    }
};

///This is the simplest test parametrization
template <class Cond, class Curr1, class Curr2, class BoundVal>
class SimpliestModel : public GECModel, private Cond, private ParentCurr, private ParentBoundVal {
public:
    SimpliestModel() : GECModel()
    {
        model.reserve(2);
        model.push_back({ 1.0, [this](double z) { return Cond::sigma(z, 0.0, 0.0); },
                          [this](double z){return Curr1::j(z);},
                          BoundVal::phi_s(0.5,0.5) });

        model.push_back({ 1.0, [this](double z) { return Cond::sigma(z, 0.0, 0.0); },
                          [this](double z){return Curr2::j(z);},
                          BoundVal::phi_s(1.5,0.5) });
    }
};


template <class Cond, class Curr, class BoundVal>
class GeoModel : public GECModel, private Cond, private ParentCurr {
protected:
    static constexpr double earth_radius2 = 6371.0*6371.0; //km^2
    unsigned N, M;
    double delta_lat, delta_lon;
    double cell_area(unsigned n, unsigned m, double delta_lat, double delta_lon)
    {
        double lat_n = -90.0 + delta_lat * n;
        if (n != N-1 and m != M-1) {
            return fabs(earth_radius2 * M_PI / 180.0 * delta_lon * (sin(M_PI / 180.0 * (lat_n + delta_lat)) - sin(M_PI / 180.0 * lat_n)));
        } else {
            if (m == M-1) {
                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - m * delta_lon) * (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            } else {
                return fabs(earth_radius2 * M_PI / 180.0 * delta_lon * (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            }
        }
    }

    double lat_arg(unsigned n, double delta_lat)
    {
        double lat_n = -90.0 + delta_lat * n;
        if (n == N-1) {
            return 0.5 * (lat_n + 90.0);
        } else {
            return lat_n + 0.5 * delta_lat;
        }
    }

    double lon_arg(unsigned m, double delta_lon)
    {
        double lon_m = delta_lon * m;
        if (m == M-1) {
            return 0.5 * (lon_m + 360.0);
        } else {
            return lon_m + 0.5 * delta_lon;
        }
    }

public:
    ///K == true => the arguments are coorginate steps
    ///K != false => the arguments define the number of steps
    GeoModel(double arg1, double arg2, bool K)
    {
        if (K==true) {
            delta_lat = arg1;
            delta_lon = arg2;
            N = std::ceil(180.0 / delta_lat);
            M = std::ceil(360.0 / delta_lon);
            model.reserve(N * M);
        } else {
            N = arg1;
            M = arg2;
            model.reserve(N * M);
            delta_lat = 180.0 / N;
            delta_lon = 360.0 / M;
        }                 //mb reserve( (N-1) * (M-1) )
        double lat_n = -90.0;
        for (unsigned n = 0; n < N; ++n) {
            lat_n =+ delta_lat * n;
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

int main()
{
    ///Degrees
    GeoModel<ExpSigma, SimpleGeoJ, ZeroPhiS> m(180.0, 360.0, false);
    //SimpliestModel<ExpSigma, StepJ, ZeroJ, ZeroPhiS> m;
    m.getPot("plots/potential_2_columns.txt", 180*180);
    //std::cout << "Ionosphere potential is " << m.getIP() << "\t[kV]" << std::endl;
    return 0;
}
