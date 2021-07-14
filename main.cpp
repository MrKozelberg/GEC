#include <cmath>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>

#include "conductivity.h"
#include "geomagnetic.h"
#include "cnpy/cnpy.h"

template<typename T>
T sqr(T x) {
    return x * x;
}

template<typename T>
T cub(T x) {
    return x * x * x;
}

constexpr double z_0 = 0.0;
constexpr double z_1 = 20.0; ///< the break of height grid (integer)
constexpr double z_max = 70.0;
constexpr size_t n_1 = 11; ///< points per kilometer lower than z_1
constexpr size_t n_2 = 19; ///< points per upper atmosphere
constexpr size_t steps = n_1 * size_t(z_1 - z_0) + n_2; ///< total number of points

/*!
 \brief This is a container of column data

 Atmosphere is divorced into columns (add the scheme!)
 */
struct Column {
    double area{};
    double altitude[steps]{}; ///< array of height points
    double sigma[steps]{}; ///< conductivity
    double j_s[steps]{}; ///< source current
    double phi_s{}; ///< additive IP from different ionospheric sources
};

/*!
 \brief (Parent Height Grid) Parent class for classes of height grids

 It is possible to create several parameterization of height grid
 */
class ParentHG {
public:
    double altitude[steps]{};
    ParentHG() = default;
    ~ParentHG() = default;
};

/*!
 \brief Linear height grid with a break of a derivative

 (Add a picture!)
 */
class LinHG: public ParentHG {
public:
    LinHG(): ParentHG() {
        for (size_t i = 0; i < steps; i++) {
            if (i <= n_1 * size_t(z_1)) {
                altitude[i] = double(i) / double(n_1);
            } else {
                altitude[i] = (double(i) - double(n_1) * z_1) * (z_max - z_1) / n_2 + z_1;
            };
        }
    }
};

/*!
 \brief Parent class for classes of area parameterizations

 It is possible to work with different parameterizations
 */
class ParentArea {
protected:
    static constexpr double earth_radius2 = 6370.0 * 6370.0; ///< km^2
public:
    ParentArea() = default;
    ~ParentArea() = default;
};

/*!
 \brief This class provides a function that calculates area of geographical cell
 */
class GeoArea : public ParentArea {
public:
    GeoArea() : ParentArea() {};

    /*!
    \brief This compute areas of geographical cells (like on a globe)
    \param n Cell number by latitude
    \param m Cell number by longitude
    \param N Dimension by latitude
    \param M Dimension by longitude
    \param d_lat Cell dimension by latitude
    \param d_lon Cell dimension by longitude
    */
    static double area(size_t n, size_t N, size_t m, size_t M, double d_lat, double d_lon) {
        double lat_n = -90.0 + d_lat * double(n);
        if (n != N - 1 and m != M - 1) {
            return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
                        (sin(M_PI / 180.0 * (lat_n + d_lat)) - sin(M_PI / 180.0 * lat_n)));
        } else {
            if (m == M - 1) {
                return fabs(earth_radius2 * M_PI / 180.0 * (360.0 - double(m) * d_lon) *
                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            } else {
                return fabs(earth_radius2 * M_PI / 180.0 * d_lon *
                            (sin(M_PI / 180.0 * 90.0) - sin(M_PI / 180.0 * lat_n)));
            }
        }
    }
};

/*!
 \brief Parent class for classes of conductivity parameterizations
 */
class ParentConductivity : public IonMobility, public IonPairProdRate, public IonIonRecombCoeff, public StdAtm {
protected:
    static constexpr double sigma_0 = 6.0e-14;
    static constexpr double H_0 = 6.0; // km
    constexpr static double e_0 = 1.602176634e-19; // C
public:
    double sigma[steps];///<

    ParentConductivity() = default;
};

/*!
 \brief parameterization of conductivity which we mainly work with

 This class provides function that creates array of conductivity value оn the considered grid
 \tparam Alt Height grid
 */
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

/*!
 \brief Simple parameterization of conductivity

 This class provides function that creates array of conductivity value оn the considered grid
 \tparam Alt Height grid
 */
template<class Alt>
class [[maybe_unused]] ExpCond : public ParentConductivity {
public:
    static double sigma_func(double z, ...) {
        return sigma_0 * std::exp(z / H_0);
    };

    ExpCond(...) : ParentConductivity() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) sigma[i] = sigma_func(a.altitude[i]);
    }
};

/*!
 \brief Constant parameterization of conductivity, it is needed for testing

 This class provides function that creates array of conductivity value оn the considered grid
 \tparam Alt Height grid
 */
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

/*!
 \brief Parent class for classes of parameterizations of source currents
 */
class ParentCurrent {
protected:
    static constexpr double j_0 = 6.4e-9;
public:
    ParentCurrent() = default;

    double j[steps];
};

/*!
 \brief Simplest parameterization of current

 This class provides function that creates array of current value оn the considered grid
 \tparam Alt Height grid
 */
template<class Alt>
class StepCurrent : public ParentCurrent {
public:
    static double current_func(double z, ...) {
        return (z >= 5.0 and z <= 10.0) ? j_0 : 0.0;
    }

    StepCurrent(...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i]);
        }
    }
};

/*!
 \brief Zero parameterization of current, it is needed for testing

 This class provides function that creates array of current value оn the considered grid
 \tparam Alt Height grid
 */
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

/*!
 \brief Simple parameterization of current

 This class provides function that creates array of current value оn the considered grid
 \tparam Alt Height grid
 */
template<class Alt>
class [[maybe_unused]] SimpleGeoCurrent : public ParentCurrent {
public:
    static double current_func(double z, double lat, ...) {
        return (std::abs(lat) <= 5 and std::abs(z-7.5) <= 2.5) ? j_0 : 0.0;
    }

    SimpleGeoCurrent(unsigned lat, ...) : ParentCurrent() {
        Alt a{};
        for (size_t i = 0; i < steps; ++i) {
            j[i] = current_func(a.altitude[i], lat);
        }
    }
};

/*!
 \brief Parameterization of current which we mainly work with

 This class provides function that creates array of current value оn the considered grid
 \tparam Alt Height grid
 */
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

/*!
 \brief Parent class for classes that provide parameterization of ionospheric sources

 This is not a boundary value, it is something like a boundary value
 */
class ParentBoundValue {
protected:
    static constexpr double phi_0 = 15.0; //[kV]
    static constexpr double theta_0 = 20.0; //[deg]
public:
    ParentBoundValue() = default;
};

/*!
 \brief Class for zero parameterization of phi_s

 We mainly use this parameterization
 */
class [[maybe_unused]] ZeroPhiS : public ParentBoundValue {
public:
    static double phi_s(...) {
        return 0.0;
    }
};

/*
 \brief Volland's parameterization of ionospheric sources

 \warning Mb it doesn't work, it was writen off some article
 */
class [[maybe_unused]] VollandPhiS : public ParentBoundValue {
public:

    /*!
     \param theta it is a longitude in deg
     \return Some k from the article
     */
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

    /*!
     \param theta It is a geomagnetic longitude
     \param lambda it is a latitude in deg
     \return The value of ionospheric sources potential
     */
    static double phi_s(double theta, double lambda) {
        return phi_0 * sin(lambda) * ((theta < theta_0) ?
                                      sin(theta) / sin(theta_0) :
                                      k(theta) * pow(sin(theta_0) / sin(theta), 4));
    }
};

/*!
 \brief This class provides possibility to calculate IP and the electric potential vs altitude

 \warning The electric potential vs altitude is temporarily unavailable
 */
class GECModel {
protected:
    size_t N, M; ///< dimensions of latitude-longitude grid
    double delta_lat{}, delta_lon{}; ///< dimensions of latitude-longitude cell
    std::vector<Column> model;
private:
    double IP = 0.0;
    bool isIPCalculated = false;
    double help[steps]{}; ///< Help-array is used when calculating the integrals

    /*!
     * \brief IP calculator
     * \todo Add formula!
     */    

    double calc_IP_Simp() {
        double up = 0.0;
        double down = 0.0;
        double h = 0.0;
        for (size_t i = 0; i < model.capacity(); ++i) {
            double int_curr_by_cond = 0.0;
            double int_rev_cond = 0.0;

            /*!
             * \brief Calculation of integrals with the help of Simpson's rule
             *
             * On each step integrands is approximated by parapola with coefficients A, B, C
             */
            double x1{}, x2{}, x3{}, y1{}, y2{}, y3{}, denom{};
            double A{}, B{}, C{};
            for (size_t k = 0; k + 2 < n_1 * size_t(z_1); k+=2){
                x1 = model[i].altitude[k], x2 = model[i].altitude[k+1], x3 = model[i].altitude[k+2];
                y1 = 1 / model[i].sigma[k], y2 = 1 / model[i].sigma[k+1], y3 = 1 / model[i].sigma[k+2];
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                B = (sqr(x3) * (y1 - y2) + sqr(x2) * (y3 - y1) + sqr(x1) * (y2 - y3)) / denom;
                C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                int_rev_cond += A * (cub(x3) - cub(x1)) / 3 + B * (sqr(x3) - sqr(x1)) / 2 + C * (x3 - x1);

                y1 = model[i].j_s[k] / model[i].sigma[k], y2 = model[i].j_s[k+1] / model[i].sigma[k+1], y3 = model[i].j_s[k] / model[i].sigma[k];
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                B = (sqr(x3) * (y1 - y2) + sqr(x2) * (y3 - y1) + sqr(x1) * (y2 - y3)) / denom;
                C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                int_curr_by_cond += A * (cub(x3) - cub(x1)) / 3 + B * (sqr(x3) - sqr(x1)) / 2 + C * (x3 - x1);
            }

            up += model[i].area * (int_curr_by_cond - model[i].phi_s) / int_rev_cond;
            down += model[i].area / int_rev_cond;
        }
        return up / down;
    }

    double calc_IP_trap() {
        double up = 0.0;
        double down = 0.0;
        double h = 0.0;

        /*!
         \warning .capacity() shows how much memory has been reserved for this vector and .size() shows how much memory is being used
         */
        for (size_t i = 0; i < model.capacity(); ++i) {
            double int_curr_by_cond = 0.0;
            double int_rev_cond = 0.0;
            // int calc
            if (model[i].altitude[2] - model[i].altitude[1] == 0){
                std::cout << i << "\n";
            }
            for (size_t q = 1; q < steps; ++q) {
                h = model[i].altitude[q] - model[i].altitude[q - 1];
                int_curr_by_cond += (model[i].j_s[q] / model[i].sigma[q] + model[i].j_s[q - 1] / model[i].sigma[q - 1]) * h / 2.0;
                int_rev_cond += (1 / model[i].sigma[q] + 1 / model[i].sigma[q - 1]) * h / 2;
            }
            up += model[i].area * (int_curr_by_cond - model[i].phi_s) / int_rev_cond;
            down += model[i].area / int_rev_cond;
        }
        return up / down;
    }

public:

    /*!
     It takes N and M, initializes dimensions of cell and help-array
     */
    GECModel(size_t N_, size_t M_) {
        assert(steps % 2 != 0);
        N = N_;
        M = M_;
        delta_lat = 180.0 / double(N);
        delta_lon = 360.0 / double(M);
        /// initialization of help array
        double A[2];
        A[0] = 2;
        A[1] = 4;
        for (size_t i = 0; i < steps; ++i) {
            help[i] = ((i != 0) or (i != steps - 1) ? A[i % 2] : 1);
        }
    }

    double get_IP() {
        if (not isIPCalculated) {
            IP = calc_IP_Simp();
            isIPCalculated = true;
        }
        return IP;
    }
};

/*!
 \brief Parent class for latitude-longitude partitioning classes, using 1x1 deg grid

 Works with NPZ-files taken from colleagues

 \warning py_array[n,m,k] is cpp_array[n*180*360 + m*360 + k]
          CHECK datatypes of npy-files!
          '<f8' is an equivalent of double
          '<f4' is an equivalent of float
          PW is cumulative data of precipitation!
 */
class ParentLatLonGrid : public GECModel {
protected:
    static constexpr double earth_radius2 = 6371.0 * 6371.0; ///< Earth radius in km^2
    /*static constexpr*/ double year = 2015.9; ///< {year}.{the day number from the year beginning by the number of days per this year}
    // input DATA
    cnpy::npz_t data = cnpy::npz_load("data/DATA-2015-12-31-00.npz");
    cnpy::NpyArray cape_arr = data["cape"];
    cnpy::NpyArray cbot_arr = data["cbot"];
    cnpy::NpyArray ctop_arr = data["ctop"];

    /*!
     \brief The portion of the area of the grid column occupied by GEC generators

     It equals the ratio of "rainc" to non-cumulative "pw"
     */
    cnpy::NpyArray alpha_arr = cnpy::npy_load("data/alpha.npy");
    float_t *cape = cape_arr.data<float>();
    double_t *cbot = cbot_arr.data<double>();
    double_t *ctop = ctop_arr.data<double>();
    double_t *alpha = alpha_arr.data<double>();

public:
    /*!
     1x1 deg grid, if you want to change the dimension of the cell, change arguments of GECModel constructor
     */
    ParentLatLonGrid() : GECModel(180, 360) {};

    [[nodiscard]] double lat_arg(unsigned n, double d_lat) const {
        return -90.0 + 0.5 + d_lat * n;
        /*
        if (n == N - 1) {
            return 0.5 * (lat_n + 90.0);
        } else {
            return lat_n + 0.5 * d_lat;
        }
        */
    }

    [[nodiscard]] double lon_arg(unsigned m, double d_lon) const {
        double lon_m = d_lon * m;
        if (m == M - 1) {
            return 0.5 * (lon_m + 360.0);
        } else {
            return lon_m + 0.5 * d_lon;
        }
    }

    [[nodiscard]] double lon_arg_m(double lon_arg_, double lat_arg_) const {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(year, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lon_m;
    }

    [[nodiscard]] double lat_arg_m(double lon_arg_, double lat_arg_) const {
        double lat_m, lon_m, alt_m;
        gdz_to_mag(year, lat_arg_, lon_arg_, 10.0, lat_m, lon_m, alt_m);
        return lat_m;
    }

    void set_year(double inputYear) {
        if (year != inputYear) {
            year = inputYear;
        }
    }
};

/*!
 Child class for latitude-longitude partitioning considering date and time
 \tparam Alt Height grid
 \tparam PhiS Ionospheric sources
 */
template<class Alt, class PhiS>
class DateLatLonGrid : public ParentLatLonGrid {
private:
    static constexpr double coef = 0.34; ///< some proportionality coefficient that defines what area occupied by sources
public:
    explicit DateLatLonGrid(double thisYear, size_t hour) : ParentLatLonGrid() {
        set_year(thisYear);
        model.reserve(2 * N * M);

        /*!
         This variable serves as the lower boundary of the cells when calculating their areas;
         when finding the values of functions from coordinates, use a special function for the argument
         */
        double lat_n = -90.0;
        for (size_t n = 0; n < N; ++n) {
            for (size_t m = 0; m < M; ++m) {
                Alt alt;
                Conductivity<Alt> cond(lat_arg_m(lon_arg(m, delta_lon), lat_arg(n, delta_lat)),0.5);
                ZeroCurrent<Alt> zero_j;
                GeoCurrent<Alt> geo_j(n, m, cbot[hour * 180 * 360 + n * 360 + m], ctop[hour * 180 * 360 + n * 360 + m],
                                      cape[hour * 180 * 360 + n * 360 + m]);

                /// Part without sources
                size_t n1 = n * 2 * M + 2 * m;
                model[n1].area =
                        GeoArea::area(n, N, m, M, delta_lat, delta_lon) *
                        (1 - coef * alpha[hour * 360 * 180 + n * 360 + m]);
                std::copy(std::begin(alt.altitude), std::end(alt.altitude),
                          std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma),
                          std::begin(model[n1].sigma));
                std::copy(std::begin(zero_j.j), std::end(zero_j.j),
                          std::begin(model[n1].j_s));
                model[n1].phi_s = PhiS::phi_s(lat_arg(n, delta_lat),
                                              lon_arg(m, delta_lon));

                /// Part with sources
                size_t n2 = n * 2 * M + 2 * m + 1;
                model[n2].area =
                        GeoArea::area(n, N, m, M, delta_lat, delta_lon) *
                        coef * alpha[hour * 360 * 180 + n * 360 + m];
                std::copy(std::begin(alt.altitude), std::end(alt.altitude),
                          std::begin(model[n2].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma),
                          std::begin(model[n2].sigma));
                std::copy(std::begin(geo_j.j), std::end(geo_j.j),
                          std::begin(model[n2].j_s));
                model[n2].phi_s = PhiS::phi_s(lat_arg(n, delta_lat),
                                              lon_arg(m, delta_lon));
            }
            lat_n += delta_lat;
        }
    }
};

/*!
 \brief This is very simplistic GEC model, using just to test the programme

 1x1, exp-sigma, simple-geo-j, zero-phi_s
 */
template<class Alt>
class Test1 : public ParentLatLonGrid {
public:
    Test1() : ParentLatLonGrid() {
        model.reserve(N * M);
        double lat_n = -90.0;
        for (size_t n = 0; n < N; ++n) {
            //std::cout << n << "\t" << lat_n << "\t" << lat_arg(n, delta_lat) << "\n";
            for (size_t m = 0; m < M; ++m) {
                Alt alt{};
                ExpCond<Alt> cond{};
                StepCurrent<Alt> j_test{};
                ZeroCurrent<Alt> j_zero{};

                size_t n1 = n * M + m;
                model[n1].area = 1.0;
                std::copy(std::begin(alt.altitude), std::end(alt.altitude), std::begin(model[n1].altitude));
                std::copy(std::begin(cond.sigma), std::end(cond.sigma), std::begin(model[n1].sigma));
                model[n1].phi_s = ZeroPhiS::phi_s();

                if (std::fabs(lat_arg(n,delta_lat)) < 5.0) {
                    std::copy(std::begin(j_test.j), std::end(j_test.j), std::begin(model[n1].j_s));
                } else {
                    std::copy(std::begin(j_zero.j), std::end(j_zero.j), std::begin(model[n1].j_s));
                }
            }
            lat_n += delta_lat;
        }

    }
};

int main() {
    std::ofstream basicOfstream("plots/diurnal_variation.txt");
    if (!basicOfstream.is_open()) {
        std::cout << "Impossible to find a file" << std::endl;
        exit(-1);
    }

    for (size_t hour = 18; hour <= 42; ++hour) {
        DateLatLonGrid<LinHG, ZeroPhiS> m(2015.9, hour);
        basicOfstream << hour << '\t' << m.get_IP() << std::endl;
    }

    basicOfstream.close();

    /*
    DateLatLonGrid<LinHG, ZeroPhiS> m(2015.9, 28);
    std::cout << m.get_IP() << std::endl;
    */

    /*!
     Test 1, the analytic answer is 8736.80375713836, while the programme gives 8927.4
     */
    /*
    Test1<LinHG> m{};
    std::cout << m.get_IP() << std::endl;
    */
    return EXIT_SUCCESS;
}
