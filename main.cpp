#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <iomanip>

#include "integral.h"

/// This is srtuct of columns;
/// You enter them in main() and take an ionosphere potential;
struct Column {
    double area;
    std::function<double(double)> conductivity;
    std::function<double(double)> current;
};

///It calculates a geopotential altitude from a geometric altitude [km]
/// |DEF| Geometric height is an altitude above mean sea level.
/// |DEF| Geopotential is an adjustment to geometric height that
///       accounts for the variation of gravity with latitude and altitude.
double gp_from_geom(double H){
    double earth_radius = 6356.766;                 ///[km] according to the U.S. Standard Atmosphere (1976)
    return H * earth_radius / (H + earth_radius);
}

///altitude is from 0 km to 70 km
///temperature and pressure can be calculated
class StdAtm{
private:
    double T_0 = 288.15;
    double p_0 = 1.01325e5;
    double R = 8.31432;
    double M = 28.9644;
    double g = 9.80665;
    constexpr static double z[7] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 70.0};
    double T[7], p[7];
    constexpr static double gamma[7] = {0.0, -6.5, 0.0, 1.0, 2.8, 0.0, -2.8};
public:
    StdAtm(){
        T[0] = T_0; p[0] = p_0;
        for(size_t n = 1; n < 7; ++n){
            T[n] = T[n-1] + gamma[n]*(z[n] - z[n-1]);
            if (gamma[n] == 0.0){
                p[n] = p[n-1]*exp(-g*M*(z[n] - z[n-1])/(R*T[n-1]));
            } else {
                p[n] = p[n-1]*pow((1.0 + gamma[n]*(z[n] - z[n-1])/T[n-1]), -g*M/(gamma[n]*R));
            }
            //std::cout << n << " " << std::setprecision(17) << p[n] << std::endl;
        }
    };
    /*
1 22632.063973462926+
2 5474.8886696777781+
3 868.01868475522883+
4 110.90630555496605+
5 66.938873118687368+
6 4.634221541689616+

    1 22632.063973462926
    2 5474.888669677778
    3 868.0186847552288
    4 110.90630555496605
    5 66.93887311868737
    6 4.634221541689616
    */
    double temperature(double H){
        double zeta = H;//gp_from_geom(H);
        size_t n = 1;
        for(size_t k = 1; k < 7; ++k) {
            if (zeta >= z[k-1] and z[k] >= zeta){n = k;}
        }
        return T[n-1] + gamma[n]*(zeta - z[n-1]);
    }
    double pressure(double H){
        double alt = H;//gp_from_geom(H);
        size_t n = 6;
        for(size_t k = 1; k < 7; ++k) {
            if (alt >= z[k-1] and z[k] > alt){n = k;}
        }
        if (gamma[n] == 0.0){
            return p[n-1]*exp(-g*M*(alt - z[n-1])/R/T[n-1]);
        } else {
            return p[n-1]*pow((1.0 + gamma[n]*(alt - z[n-1])/T[n-1]), -g*M/gamma[n]/R);
        }
    }
};

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

/// explicit conductivity functions
class Conductivities: private SMZ15, private TZ06, private ZT07{
protected:
    static constexpr double sigma_0 = 5.0e-14;
    static constexpr double H_0 = 6.0;                      /// [km]
    static constexpr double e_0 = 1.602176634e-19;          /// [C]
public:
    ///z in kilometres
    static double exp_conductivity(double z){
             return sigma_0*exp(z/H_0);
    };
    static double const_conductivity(double z){
        return sigma_0;
    }
    double conductivity(double z, double lymbda, double xi){
        return e_0 * mu_sum(z) * sqrt(q(z, lymbda, xi) / alpha(z));
    }
};

/// explicit current functions
class Currents{
protected:
    static constexpr double j_0 = 1e-9;
public:
    static double step_current(double z){

        if(z>=5.0 and z<10.0){
            return j_0;
        } else{
            return 0.0;
        }
    }
    static double zero_current(double z){
        return 0.0;
    }
};

/// It is the central class, here you can find any parameters of the model
class GECModel{
protected:
    static constexpr double H = 70.0;
    double V = 0.0;
    std::vector<std::vector<double>> potential;
    std::vector<double> altitude;
    unsigned pot_points = 70;
    bool isVCalculated = false;
    bool isPhiCalculated = false;
    std::vector<Column> model;
    unsigned int_points = 10'001;

    /// This is a function that calculates an ionosphere potention
    double calc_ip() {
        std::vector<double> area_by_int(model.size());
        std::vector<double> int_of_cur_by_cond(model.size());
        for(unsigned i = 0; i < model.size(); ++i){
            area_by_int[i] = model[i].area / integrate_Simpson(
                        [this,i](double z) {return 1 / model[i].conductivity(z);},
                        0.0, H, int_points);
            int_of_cur_by_cond[i] = integrate_Simpson(
                        [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                        0.0, H, int_points);
        }
        double up = 0.0, down = 0.0;
        for(unsigned i = 0; i < model.size(); ++i){
            up += area_by_int[i] * int_of_cur_by_cond[i];
            down += area_by_int[i];
        }
        return up / down;
    }

    /// This is a function that calculates dependence of the potention on the varying altitude
    void calc_pot(){
        //unsigned new_int_points = 2 * unsigned(int_points / pot_points) + 1;
        potential.reserve(model.size());
        double I1 = 0.0, I2 = 0.0, C1, C2;
        altitude.reserve(pot_points);
        altitude[0] = 0.0;
        for(size_t j = 1; j < pot_points; ++j){
             altitude[j] = j * H / pot_points;
        }
        potential.reserve(model.size());
        for(unsigned i = 0; i < model.size(); ++i){
            potential[i].reserve(pot_points);
        }
        for(unsigned i = 0; i < model.size(); ++i){
            C1 = integrate_Simpson(
                        [this,i](double z) {return 1 / model[i].conductivity(z);},
                        0.0, H, int_points);
            C2 = integrate_Simpson(
                        [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                        0.0, H, int_points);
            potential[i][0] = 0.0;
            for(size_t j = 1; j < pot_points; ++j){
                /// WARNING! Integrate from 0 only!! (?)
                I1 = integrate_Simpson(
                            [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                            0.0, altitude[j], int_points);
                I2 = integrate_Simpson(
                            [this,i](double z) {return 1 / model[i].conductivity(z);},
                            0.0, altitude[j], int_points);
                potential[i][j] = I1 - (I2 / C1) * (C2 - getIP());
            }
        }
    }

public:
    GECModel(std::vector<Column> aModel): model(std::move(aModel)) {
    }
    GECModel() {
    }
    double getIP() {
        if (not isVCalculated) {
            V = calc_ip();
            isVCalculated = true;
        }
        return V;
    }
    double getPot(size_t col, size_t j) {
        if (not isPhiCalculated) {
            calc_pot();
            isPhiCalculated = true;
        }
        return potential[col][j];
    }
    double getAlt(size_t j) {
        if (not isPhiCalculated) {
            calc_pot();
            isPhiCalculated = true;
        }
        return altitude[j];
    }
    size_t getPhiPoints(){return pot_points;}
    double getH(){return H;}
};

class ConcreteGECModel: public GECModel, private Conductivities, private Currents, private StdAtm {
public:
    ConcreteGECModel(): GECModel() {
        ///This is the simplest parametrization
        model.reserve(2);
        model.push_back({1.0, [this](double z){return Conductivities::conductivity(z, 1.2, 1.0);}, Currents::step_current});
        model.push_back({1.0, [this](double z){return Conductivities::conductivity(z, 0.9, 1.0);}, Currents::zero_current});
    }
};

int main(){

    /// The simplest parametrization (2 columns, considered conductivity, step and zero currents)
    /*ConcreteGECModel m;
    std::ofstream fout("plots/potential_2_columns.txt");
    if(fout.is_open() == false){
        std::cout << "Impossible to find a file" << std::endl;
        return 1;
    }
    for(size_t i = 0; i < m.getPhiPoints(); ++i){
        fout << m.getAlt(i) << "\t" << m.getPot(0,i) << "\t" <<  m.getPot(1,i) << std::endl;
    }
    fout << 70 << "\t" << m.getIP() << "\t" << m.getIP() << std::endl;
    fout.close();*/

    /// US Standard Atmosphere 1976
    /// This part of programme works correctly
    /*std::ofstream fout("plots/my_atm.txt");
    if(fout.is_open() == false){
        std::cout << "Impossible to find a file" << std::endl;
        return 1;
    }
    StdAtm atm;
    for(double z = 0.0; z <= 70.1; z += 0.1){
        fout << z << "\t" << atm.temperature(z) << "\t" << atm.pressure(z) << "\n";
    }
    fout.close();*/

    /// Conductivity parametrization
    /*std::ofstream fout("plots/cond.txt");
    if(fout.is_open() == false){
        std::cout << "Impossible to find a file" << std::endl;
        return 1;
    }
    Conductivities sigma;
    for(double z = 0.0; z <= 70.1; z += 0.1){
        /// z is a geometric altitude
        fout << z << "\t" << sigma.conductivity(z, 1.0, 0.0) << "\t"
             << sigma.conductivity(z, 1.0, 0.5) << "\t"
             << sigma.conductivity(z, 1.0, 1.0) << "\n";
    }
    fout.close();*/

    std::ofstream fout("plots/test2.txt");
    if(fout.is_open() == false){
        std::cout << "Impossible to find a file" << std::endl;
        return 1;
    }
    Conductivities sigma;
    StdAtm atm;
    for(double z = 0.0; z <= 70.1; z = z+0.1){
        fout << std::fixed  << std::setprecision(1) << z << "\t" << std::setprecision(16)<< atm.pressure(z) << "\n";
    }
    fout.close();

    std::cout << std::setprecision(1) << 1.9 << "\t" << std::setprecision(15)<< atm.pressure(1.9) << "\t" << std::setprecision(15)<< atm.temperature(1.9) << "\n";

    return 0;
}
