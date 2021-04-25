#include <cmath>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>

#include "conductivity.h"
#include "geomagnetic.h"
#include "cnpy/cnpy.h"

template<typename T>
T sqr(T x) {
    return x * x;
}

/// some constants
constexpr size_t number_of_uniform_points_per_km = 19; // points per kilometre
constexpr size_t max_high_of_uniform_grid = 20; // maximum high of linear grid
constexpr size_t number_of_points_per_upper_atm = 7; // points per upper atmosphere
constexpr size_t steps = number_of_uniform_points_per_km * max_high_of_uniform_grid + number_of_points_per_upper_atm;

/// This is a struct of columns
struct Column {
    double _area{};
    double _altitude[steps]{}; // last point should be in 70 km
    double _conductivity[steps]{}; //должно быть обязательно нечетным числом
    double _current[steps]{};
    double _bound_value{};
};

class GECModelCalculation {
protected:
    size_t N, M;
    double delta_lat{}, delta_lon{};
    std::vector<Column> model;
private:
    double IP = 0.0;
    bool isIPCalculated = false;

    double calc_IP() {

        double area_by_int[2*N*M];
        double int_curr_by_cond[2*N*M];

        /// initialization of help
        double help[steps]{};
        try {
            if (steps % 2 != 1) throw steps;
        } catch (size_t k) {
            std::cout << "'steps' has to be uneven" << std::endl;
            exit(-1);
        }
        double A[2];
        A[0] = 2;
        A[1] = 4;
        for (size_t i = 0; i < steps; ++i) {
            help[i] = ( (i != 0) or (i != steps - 1) ? A[i % 2] : 1);
        }

        /// calculation of integrals
        for (size_t i = 0; i < 2*N*M; ++i) {
            double int_reverse_cond = 0.0;
            double int_curr_by_cond_ = 0.0;
            for (size_t q = 0; q < steps - 1; ++q) {
                double h = model[i]._altitude[q + 1] - model[i]._altitude[q];
                int_reverse_cond += 1 / model[i]._conductivity[q] * h / 3;
                int_curr_by_cond_ += model[i]._current[q] / model[i]._conductivity[q] * h / 3;
            }
            area_by_int[i] = model[i]._area / int_reverse_cond;
            int_curr_by_cond[i] = int_curr_by_cond_;
        }

        /// calculation of sums
        double up = 0.0, down = 0.0;
        for (size_t i = 0; i < 2*N*M; ++i) {
            up += area_by_int[i] * (int_curr_by_cond[i] - model[i]._bound_value);
            down += area_by_int[i];
        }
        return up / down;
    }
public:
    GECModelCalculation(size_t N_, size_t M_) {
        N = N_;
        M = M_;
        delta_lat = 180.0 / N;
        delta_lon = 360.0 / M;
    }

    double get_IP() {
        if (not isIPCalculated) {
            IP = calc_IP();
            isIPCalculated = true;
        }
        return IP;
    }
};

/// This set of classes provides some functions for calculation a column area
class ParentArea {
protected:
    static constexpr double earth_radius2 = 6370.0 * 6370.0; ///< km^2
public:
    ParentArea() = default;
};

class StupidArea : public ParentArea {
public:
    StupidArea() = default;

    static double area(...) {
        return 1;
    }
};

class GeoArea : public ParentArea {
public:
    GeoArea() = default;

    static double area(size_t n, size_t N, size_t m, size_t M, double d_lat, double d_lon) {
        double lat_n = -90.0 + d_lat * n;
        if (n != N - 1 and m != M - 1) {
            return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
                        (sin(M_PI / 180.0 * (lat_n + d_lat)) - sin(M_PI / 180.0 * lat_n)));
        } else {
            if (m == M - 1) {
                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - m * d_lon) *
                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            } else {
                return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            }
        }
    }
};

class ParentAltitude {
protected:
    static constexpr double max_altitude = 70.0;
public:
    double altitude[steps]{};

    ParentAltitude() = default;

    ~ParentAltitude() = default;
};

class UniformAltitude : public ParentAltitude {
public:
    UniformAltitude() {
        for (size_t i = 0; i < steps; ++i) {
            altitude[i] = double(i) / (steps - 1) * max_altitude;
        }
    }

    ~UniformAltitude() = default;

};

/* IT DOES NOT WORK! */
/*class ExpAltitude : public ParentAltitude {
private:
    static constexpr double A = (max_altitude - max_high_of_uniform_grid) /
               (steps - max_high_of_uniform_grid * number_of_uniform_points_per_km) *
               (steps - max_high_of_uniform_grid * number_of_uniform_points_per_km);
    static constexpr double x_0 = max_high_of_uniform_grid * number_of_uniform_points_per_km - 1;
public:
    ExpAltitude() = default;
    ~ExpAltitude() = default;

    void altitude_calc() {
        for (size_t i = 0; i < steps; ++i) {
            if (i <= max_high_of_uniform_grid * number_of_uniform_points_per_km) {
                altitude[i] =
                        i / (max_high_of_uniform_grid * number_of_uniform_points_per_km - 1) * max_high_of_uniform_grid;
            } else {
                altitude[i] = A * sqr(i - x_0) + max_high_of_uniform_grid;
            }
        }
    }
};*/

class ParentConductivity : public IonMobility, public IonPairProdRate, public IonIonRecombCoeff, public StdAtm {
protected:
    static constexpr double sigma_0 = 5.0e-14;
    static constexpr double H_0 = 6.0; // km
    constexpr static double e_0 = 1.602176634e-19; // C
public:
    double sigma[steps];

    ParentConductivity() = default;
};

template<class Alt>
class Conductivity : public ParentConductivity {
public:
    static double sigma_func(double z, double lambda, double xi) {
        const double T = StdAtm::temperature(z);
        const double p = StdAtm::pressure(z);
        const double n = std::sqrt(IonPairProdRate::q(z, lambda, xi, p, T) / IonIonRecombCoeff::alpha(z, p, T));
        return e_0 * (IonMobility::mu_p_get(T, p) + IonMobility::mu_m_get(T, p)) * n;
    }

    Conductivity(double lambda, double xi) : ParentConductivity() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            sigma[i] = sigma_func(a.altitude[i], lambda, xi);
        }
    }
};

template<class Alt>
class [[maybe_unused]] ExpCond : public ParentConductivity {
public:
    static double sigma_func(double z, ...) {
        return sigma_0 * std::exp(z / H_0);
    };

    ExpCond(double lambda, double xi) : ParentConductivity() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            sigma[i] = sigma_func(a.altitude[i], lambda, xi);
        }
    }
};


template<class Alt>
class [[maybe_unused]] ConstSigma : protected ParentConductivity {
public:
    static double sigma_func(...) {
        return sigma_0;
    }

    ConstSigma(double lambda, double xi) : ParentConductivity() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            sigma[i] = sigma_func(a.altitude[i], lambda, xi);
        }
    }
};

class ParentCurrent {
protected:
    static constexpr double j_0 = 6.4e-9;
public:
    ParentCurrent() = default;

    double j[steps];
};

template<class Alt>
class StepCurrent : public ParentCurrent {
public:
    static double current_func(double z, ...) {
        return (z >= 6.0 and z < 11.0) ? j_0 : 0.0;
    }

    StepCurrent(unsigned lat, unsigned lon, double cbot, double ctop, double cape) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat, lon, cbot, ctop, cape);
        }
    }
};

template<class Alt>
class [[maybe_unused]] ZeroCurrent : public ParentCurrent {
public:
    static double current_func(...) {
        return 0.0;
    }

    ZeroCurrent(...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func();
        }
    }
};

template<class Alt>
class [[maybe_unused]] SimpleGeoCurrent : public ParentCurrent {
public:
    static double current_func(double z, double lat, ...) {
        return (std::abs(lat) <= 10) ? StepCurrent<Alt>::current_func(z) : ZeroCurrent<Alt>::current_func(z);
    }

    SimpleGeoCurrent(unsigned lat, unsigned lon, double cbot, double ctop, double cape) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat, lon, cbot, ctop, cape);
        }
    }
};

template<class Alt>
class GeoCurrent : public ParentCurrent {
private:
    static constexpr double cape_0 = 1'000; // J/kg
public:
    double current_func(double z, unsigned lat, unsigned lon, double cbot, double ctop, double cape) {
        return cape < cape_0
               || cbot >= z * 1'000
               || ctop <= z * 1'000 ? 0 : j_0;
    }

    GeoCurrent(unsigned lat, unsigned lon, double cbot, double ctop, double cape) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat, lon, cbot, ctop, cape);
        }
    }
};

/// explicit boundary values function
///it is an empty parent class for a while
class ParentBoundValue {
protected:
    static constexpr double phi_0 = 15.0; //[kV]
    static constexpr double theta_0 = 20.0; //[deg]
public:
    ParentBoundValue() = default;
};

class [[maybe_unused]] ZeroPhiS : public ParentBoundValue {
public:
    static double phi_s(...) {
        return 0.0;
    }
};

class [[maybe_unused]] VollandPhiS : public ParentBoundValue {
public:
    //theta is a longitude [deg]
    static double k(double theta) {
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
    static double phi_s(double theta, double lambda) {
        return phi_0 * sin(lambda) * ((theta < theta_0) ?
                                      sin(theta) / sin(theta_0) :
                                      k(theta) * pow(sin(theta_0) / sin(theta), 4));
    }
};

class ParentLatitudeLongitudeGrid : public GECModelCalculation {
protected:
    static constexpr double earth_radius2 = 6371.0 * 6371.0; ///< km^2
    //input DATA
    // py_array[n,m,k] is cpp_array[n*180*360 + m*360 + k]
    // CHECK datatypes of npy-files!
    // '<f8' is an equivalent of double
    // '<f4' is an equivalent of float
    cnpy::npz_t data = cnpy::npz_load("data/DATA-2015-12-31-00.npz");
    cnpy::NpyArray cape_arr = data["cape"];
    cnpy::NpyArray cbot_arr = data["cbot"];
    cnpy::NpyArray ctop_arr = data["ctop"];
    cnpy::NpyArray alpha_arr = cnpy::npy_load("data/alpha.npy");
    float_t *cape = cape_arr.data<float>();
    double_t *cbot = cbot_arr.data<double>();
    double_t *ctop = ctop_arr.data<double>();
    double_t *alpha_ = alpha_arr.data<double>();
public:
    ParentLatitudeLongitudeGrid() : GECModelCalculation(90, 180) {};
    double lat_arg(unsigned n, double d_lat) const {
        double lat_n = -90.0 + d_lat * n;
        if (n == N - 1) {
            return 0.5 * (lat_n + 90.0);
        } else {
            return lat_n + 0.5 * d_lat;
        }
    }

    double lon_arg(unsigned m, double d_lon) const {
        double lon_m = d_lon * m;
        if (m == M - 1) {
            return 0.5 * (lon_m + 360.0);
        } else {
            return lon_m + 0.5 * d_lon;
        }
    }

    static double lon_arg_m(double lon_arg_, double lat_arg_) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(2015.9, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lon_m;
    }

    static double lat_arg_m(double lon_arg_, double lat_arg_) {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(2015.9, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lat_m;
    }
};

template<class Alt, class PhiS>
class DataLatitudeLongitudeGrid : public ParentLatitudeLongitudeGrid {
private:
    static constexpr double k = 1.0; // coefficient defines how many
public:
    explicit DataLatitudeLongitudeGrid(size_t hour) : ParentLatitudeLongitudeGrid() {
        model.reserve(2 * N * M);
        double lat_n = -90.0;
        for (size_t n = 0; n < N; ++n) {
            lat_n += delta_lat;
            for (size_t m = 0; m < M; ++m) {
                Alt alt;
                Conductivity<Alt> cond(lat_arg_m(lon_arg(m, delta_lon), lat_arg(n, delta_lat)),0.5);
                ZeroCurrent<Alt> zero_j;
                GeoCurrent<Alt> geo_j(n, m, cbot[hour * 180 * 360 + n * 360 + m], ctop[hour * 180 * 360 + n * 360 + m],
                                      cape[hour * 180 * 360 + n * 360 + m]);

                /// Part without sources
                size_t n_1 = n * M + m;
                model[n_1]._area =
                        GeoArea::area(n, N, m, M, delta_lat, delta_lon) *
                        (1 - k * alpha_[hour * 360 * 180 + n * 360 + m]);
                std::copy(std::begin(alt.altitude), std::end(alt.altitude),
                          std::begin(model[n_1]._altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma),
                          std::begin(model[n_1]._conductivity));
                std::copy(std::begin(zero_j.j), std::end(zero_j.j),
                          std::begin(model[n_1]._current));
                model[n_1]._bound_value = PhiS::phi_s(lat_arg(n, delta_lat),
                                                                      lon_arg(m, delta_lon));

                /// Part with sources
                size_t n_2 = 2 * N * M - 1 - n * M - m;
                model[n_2]._area =
                        GeoArea::area(n, N, m, M, delta_lat, delta_lon) *
                        k * alpha_[hour * 360 * 180 + n * 360 + m];
                std::copy(std::begin(alt.altitude), std::end(alt.altitude),
                           std::begin(model[n_2]._altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma),
                          std::begin(model[n_2]._conductivity));
                std::copy(std::begin(geo_j.j), std::end(geo_j.j),
                          std::begin(model[n_2]._current));
                model[n_2]._bound_value = PhiS::phi_s(lat_arg(n, delta_lat),
                                                                      lon_arg(m, delta_lon));
            }
        }
    }
};

int main() {
    std::ofstream fout("plots/diurnal_variation.txt");
    if (!fout.is_open()) {
        std::cout << "Impossible to find a file" << std::endl;
        exit(-1);
    }
    for (size_t hour = 0; hour < 49; ++hour) {
        DataLatitudeLongitudeGrid<UniformAltitude, ZeroPhiS> m(hour);
        std::cout << hour << '\t' << m.get_IP() << std::endl;
    }
    fout.close();
    return EXIT_SUCCESS;
}
